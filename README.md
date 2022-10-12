# SigMa
SigMa was designed for iterative prophage identification in bacterial genomes.

It is meant to allow incorporation of validated (manually or autoamtically) prophage predictions of analyzed set of genomes prior next iteration.

## TODO
### Major things to do
- [x] - dereplication of candidate regions in automatic run mode
- [ ] - managing config files or input mapping files
- [ ] - default reference datasets configuration
- [ ] - managing iterations
- [ ] - determine SigMa requirements and prepare envs

### Optional
- [ ] - implement PhiSpy
- [ ] - rewrite diamond to MMSEQs for protein searches
- [ ] - implement HMMER (pVOGs, VOGDB)
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
    - [x] write candidate regions --> `regions.candidate/`
    - [x] write Artemis plot files --> `artemis_plots/`
5. Validate predictions --> `regions.verified`
    - [ ] manually - go throught predictions and select phage and non-phage regions
    - [x] automatically - run CheckV on picked regions and select Complete or High-quality ones only

