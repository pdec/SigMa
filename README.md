# SigMa
SigMa was designed for iterative prophage identification in bacterial genomes.

It is meant to allow incorporation of validated (manually or autoamtically) prophage predictions of analyzed set of genomes prior next iteration.

## TODO
### Major things to do
- [x] - dereplication of candidate regions in automatic run mode
- [x] - write GenBank files
- [x] - write summary TSV
- [ ] - managing config files or input mapping files
- [x] - default reference datasets configuration
- [ ] - managing iterations
- [ ] - determine SigMa requirements and prepare envs

### Optional
- [ ] - implement PhiSpy
- [ ] - rewrite diamond to MMSEQs for protein searches
- [ ] - implement HMMER (pVOGs, VOGDB)
## Installation
### Requirements
- python3
- ncbi-blast
- mmseqs2
- diamond
- MeShClust
- CheckV

## Reference datasets
By default SigMa relies on the following datasets:
- [PHROGs](https://phrogs.lmge.uca.fr/) MMSEQs HMM profile database
- prophages we have manually identified
- [INPHARED](https://github.com/RyanCook94/inphared) database

We cluster prophages and INPHARED database separatey to track the advantage of prophage discovery and reuse.

During the run, SigMa considers all reference datasets separately but can combine that if `--combine` or `--combine --sig_sources combined` is used.
After candidate regions are picked from each signal groups (all reference datasets and potentially 'combined' group) they are all overlayed and merged to pick unique ranges within each sequence to verify further with CheckV and/or manually. Coordinates of all candidates are still written to summary files.

### Set up reference databases
The repo comes with a set of prophages we have manually identified - `./data/pps.gb.gz`.
We'll put all reference databases to `./dbs/` directory within the repo.

Let's create the db directory.

`if [ ! -z ./dbs ]; then mkdir ./dbs; fi`

Now let's prepare pps database running provided script. It takes GenBank and FASTA files as input and generates protein and nucleotide clustered files produced with MMSeqs and [MeShClust](https://github.com/BioinformaticsToolsmith/Identity).

`python3 scripts/prepdb.py --gb ./data/pps.gb.gz --dbdir ./dbs --threads 4 --tmp ./tmp`

To download PHROGS MMSeqs HMM profiles you can use these commands:

```bash
curl https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz
tar zxvf phrogs_mmseqs_db.tar.gz
mv phrogs_mmseqs_db/phrogs_profile_db* ./dbs/
rm -r phrogs_mmseqs_db*
```

To download INPHARED latest database you cas use these commands but remember to modify the date.
We'll download all phage genbanks as we'll cluster them anyway and it's easier to download the whole file but you may just download all genbanks based on <date>_data_excluding_refseq.tsv file.


```bash
curl http://inphared.s3.climb.ac.uk/1Oct2022_phages_downloaded_from_genbank.gb -o inphared.gb
# consider more threads as this one might take a while
python3 scripts/prepdb.py --gb inphared.gb --dbdir ./dbs --threads 4 --tmp ./tmp
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
    - [x] set up datasets
2. Prepare query sequences --> `query/`
    - [x] extract protein and nucleotide sequences
3. Search --> `search/`
    - [ ] PhiSpy
    - [x] nucleotide-based search
    - [x] protein-based search
    - [x] mmseqs profile-based search
    - [ ] hmmer profile-based search
4. Evaluate results
    - [x] read results
    - [x] map results/signal presence on sequences
    - [x] merge adjacent/overlapping signals 
    - [x] determine candidate regions
    - [x] merge overlapping regions and consider as candidates
    - [x] write candidate regions --> `regions/candidate.fasta`
    - [x] write Artemis plot files --> `artemis_plots/`
5. Validate predictions --> `regions/verified.fasta` and `regions/verified.gb`
    - [ ] manually - go throught predictions and select phage and non-phage regions
    - [x] automatically - run CheckV on picked regions and select 
      - [x] Complete
      - [x] High-quality
      - [x] Medium-quality but of length >= 20kb and fraction >= 0.85

