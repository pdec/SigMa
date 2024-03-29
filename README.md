# SigMa
SigMa was designed for iterative prophage identification in bacterial genomes.

It is meant to allow incorporation of validated (manually or automatically) prophage predictions of analyzed sets of genomes prior next iteration(s).

## Installation
Before you start using SigMa, you should have programs and packages installed on your system or conda env. 
We provide a conda configuration file `sigma_env.yml` to make the process easier, yet you will still need to install MeShClust yourself as well as configure CheckV database.

```bash
git clone git@github.com:pdec/SigMa.git
cd SigMa
# use conda or mamba to create env
conda env create --file sigma_env.yml

# Install MeShClust
# go to location you want to install MeShClust
git clone git@github.com:BioinformaticsToolsmith/Identity.git
cd Identity
mkdir bin
cd bin
cmake ..
make
# either add current directory to PATH 
PATH=`pwd`"/${PATH:+:${PATH}}"; export PATH;
# or make symbolic link to meshclust into directory from you PATH
ln -s `pwd` <directory_in_your_PATH>

# configure CheckV
conda activate sigma
checkv download_database <path_to_checkv_db>
# you might also want to add this path to your system env or you might as well skip that and provide path while running SigMa
export CHECKVDB=<path_to_checkv_db>
```
### Requirements
  - mmseqs==14.7e284
  - pandas>=1.4.3
  - numpy>=1.22.3
  - biopython>=1.78
  - blast>=2.13.0
  - checkv>=1.0.1 (requires database initialization)
  - MeShClust (requires a separate installation)

## Reference datasets
By default SigMa relies on the following datasets:
- [PHROGs](https://phrogs.lmge.uca.fr/) MMSEQs HMM profile database
- prophages we have manually identified
- [INPHARED](https://github.com/RyanCook94/inphared) database

We cluster prophages and INPHARED database separatey to track the advantage of prophage discovery and reuse.

During the run, SigMa considers all reference datasets separately but can combine that if `--combine` or `--combine --sig_sources combined` is used.
After candidate regions are picked from each signal groups (all reference datasets and potentially 'combined' group) they are all overlayed and merged to pick unique ranges within each sequence to verify further with CheckV and/or manually. Coordinates of all candidates are still written to summary files.

### Set up reference databases
The repo comes with a set of prophages we have manually identified - `./data/pps.gb.gz` and a subset of PHROGS database `./data/PHROGS_SigMa.tar.gz` with profiles for tail, head and packaging, lysis, connector, and integration and excision functional categories. This way we're reducing a lot of false positives. 
We'll put all reference databases to `./dbs/` directory within the repo.

Let's create the db directory.

`if [ ! -z ./dbs ]; then mkdir ./dbs; fi`

Now let's prepare pps database running provided script. It takes GenBank and FASTA files as input and generates protein and nucleotide clustered files produced with MMSeqs and [MeShClust](https://github.com/BioinformaticsToolsmith/Identity).

```bash
conda activate sigma
python scripts/prepdb.py --gb ./data/pps.gb.gz --dbdir ./dbs --threads 4 --tmp ./tmp`
```

Extract provided PHROGS MMseqs2 profiles:
```bash
tar xzvf ./data/PHROGS_SigMa.tar.gz -C ./dbs
```

To download all PHROGS MMSeqs2 profiles you can use these commands:

```bash
curl https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz
tar zxvf phrogs_mmseqs_db.tar.gz
mv phrogs_mmseqs_db/phrogs_profile_db* ./dbs/
rm -r phrogs_mmseqs_db*
```

However, we have noticed that to use the newest MMSeqs v.14 profile database needs to be recreated.
You can use the following commands to achieve that.

```bash
%%bash
# download MSA FASTA files
wget https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz
mmseqs tar2db MSA_phrogs.tar.gz phrogsMSA --output-dbtype 11 --tar-include '.+\.fma$'
mmseqs msa2profile phrogsMSA ./dbs/phrogsProfileDB
rm MSA_phrogs.tar.gz phrogsMSA*
```

To download INPHARED latest database you cas use these commands but remember to modify the date.
We'll download all phage genbanks as we'll cluster them anyway and it's easier to download the whole file but you may just download all genbanks based on <date>_data_excluding_refseq.tsv file.


```bash
conda activate sigma
curl http://inphared.s3.climb.ac.uk/1Oct2022_phages_downloaded_from_genbank.gb -o inphared.gb
# consider more threads as this one might take a while
python scripts/prepdb.py --gb inphared.gb --dbdir ./dbs --threads 8 --tmp ./tmp
rm inphared.gb
```

You can do that with your own GenBank or FASTA files as well. Input file name without extention will be further used as reference dataset name, therefore avoid using whitespace characters.

# Workflow
SigMa is supposed to run in a single on multiple iterations mode with and without breaks for manual validation of predictions from certain iterations.

## Overiew
1. Search input GenBank files against a set of reference datasets --> determine candidate regions
2. Validate cadidate regions --> verified regions are written in `regions.verified.tsv` and `regions.verified/` directory where GenBanks of those are stored
    1. Automatically using CheckV - regions will be deduplicated and written to the final file
    2. Manually - regions need to be written by user
    3. Semi-automatic - run CheckV and let User further verify that
3. Perform another iteration
    1. Automatically using CheckV regions
    2. Manually - after calling SigMa with `--resume` argument
4. Check for change in verified regions
    1. Update databases
    2. Update CheckV reference datadabase
## Single run steps
Each iteration consists of the following steps:
1. Set up reference datasets --> `reference/`
    - set up datasets
2. Prepare query sequences --> `query/`
    - extract protein and nucleotide sequences
3. Search --> `search/`
    - PhiSpy --> `phispy/`
    - nucleotide-based search
    - protein-based search
    - mmseqs profile-based search
4. Evaluate results
    - read results
    - map results/signal presence on sequences
    - merge adjacent/overlapping signals 
    - determine candidate regions
    - merge overlapping regions and consider as candidates
    - write candidate regions --> `regions/candidate.fasta`
    - write Artemis plot files --> `plots/`
5. Validate predictions --> `regions/verified.fasta` and `regions/verified.gb`
    - manually - go throught predictions and select phage and non-phage regions
    - automatically - run CheckV on picked regions and select 
       - Complete
       - High-quality
       - Medium-quality but of length >= 20kb and fraction >= 0.85

## Citation
If you use SigMa remember to cite the tools it incorporates:
- Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T.L. BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421, (2009). https://doi.org/10.1186/1471-2105-10-421.
- Steinegger M. and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35, 1026–1028 (2017). https://doi.org/10.1038/nbt.3988.
- Cook R, Brown N, Redgwell T, Rihtman B, Barnes M, Clokie M, Stekel DJ, Hobman JL, Jones MA, Millard A. INfrastructure for a PHAge REference Database: Identification of Large-Scale Biases in the Current Collection of Cultured Phage Genomes. PHAGE, 214-223 (2021). http://doi.org/10.1089/phage.2021.0007.
- Terzian P., Olo Ndela E., Galiez C., Lossouarn J., Pérez Bucio R.E., Mom R., Toussaint A., Petit M.A., Enault F. PHROG : families of prokaryotic virus proteins clustered using remote homology", NAR Genomics and Bioinformatics, Volume 3, Issue 3, September 2021, lqab067. https://doi.org/10.1093/nargab/lqab067
- Hani Z. Girgis. MeShClust v3.0: High-quality clustering of DNA sequences using the mean shift algorithm and alignment-free identity scores BMC Genomics 23, 423 (2022). https://doi.org/10.1186/s12864-022-08619-0
## Funding
- The National Science Centre PRELUDIUM 15 grant no. 2018/29/N/NZ8/00228.
- The Polish National Agency for Academic Exchange Bekker Programme Fellowship no. BPN/BEK/2021/1/00416.