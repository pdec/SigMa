# SigMa
SigMa was designed for iterative prophage identification in bacterial genomes.

It is meant to allow incorporation of validated (manually or autoamtically) prophage predictions of analyzed set of genomes prior next iteration.

# Workflow
Each iteraction consists of the following steps:
1. Set up reference database
    - 1st - set up default databases
    - 2nd and later - make databases from predictions
2. Prepare query sequences
    - 1st - pick a subset of query sequences and extract protein and nucleotide sequences
    - 2nd and up to N - pick another subset of query sequences
    - N - re-search last N subsets
3. Search
    - PhiSpy
    - nucleotide-based search
    - protein-based search
    - profile-based search
4. Evaluate results
    - read results
    - map results/signal presence on sequences
    - merge adjacent/overlapping signals 
    - prioritise signals
5. Validate predictions
    - manually - go throught predictions and select phage and non-phage regions
    - automatically - run CheckV on picked regions
    - manually and automatically - run CheckV only on manually picked regions

