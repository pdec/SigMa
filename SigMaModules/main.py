
import sys
import os
import argparse
import numpy as np
from typing import List

from .models import SigMa
from .utils import create_logger, log_progress, CustomHelpFormatter
from .version import __version__

def main():

    parser = argparse.ArgumentParser(description='SigMa - a tool for iterative mapping of phage signatures and prophage recovery.', epilog='Example: SigMa.py run -r reference.fasta -q query.fasta -o output_dir\n', formatter_class=CustomHelpFormatter)

    parser_setup = parser.add_argument_group('Setup')
    parser_search = parser.add_argument_group('Search')
    parser_evaluate = parser.add_argument_group('Evaluate')
    parser_validate = parser.add_argument_group('Validate')

    REFERENCE_TYPES = [
        'fasta_nt',
        'fasta_aa',
        'genbank',
        'hmm',
        'mmseqs_db',
    ]

    QUERY_TYPES = [
        'genbank'
    ]

    SIG_SOURCES = [
        'all',
        'combined',
    ]

    LOGGING = [
        'DEBUG',
        'INFO',
        'WARNING',
        'ERROR',
        'CRITICAL',
    ]

    ### Run
    parser.add_argument('-o', '--outdir', help='Output directory', metavar = '<path>', required=True)
    parser.add_argument('-b', '--batches', help='Number of batches to run on input files [%(default)i]', default = 1, metavar = '<num>', type = int)
    parser.add_argument('-t', '--threads', help='Number of threads to use for data preparation [%(default)i]', default = 4, metavar = '<num>', type = int)
    parser.add_argument('-v', '--version', help='Show program\'s version and exit.', action='version', version='%(prog)s 0.1')
    parser.add_argument('-l', '--logging', help='Output logging level [%(default)s]. Allowed options: %(choices)s', default = 'INFO', choices = LOGGING, metavar = ' ', type = str)
    
    ### Setup
    parser_setup.add_argument('-r', '--reference', nargs = '+', help='Reference dataset(s)', metavar = '<path>', type = str, action = 'extend', required = True)
    parser_setup.add_argument('-R', '--reference_type', nargs = '+', choices = REFERENCE_TYPES, help='Reference dataset type(s). Allowed types: %(choices)s', metavar = '<type>', type = str, action = 'extend', required = True)
    parser_setup.add_argument('-q', '--query', nargs = '+', help='Query dataset(s)', metavar = '<path>', type = str, action = 'extend', required = True)
    parser_setup.add_argument('-Q', '--query_type', nargs = '+', choices = QUERY_TYPES, help='Query dataset type(s). Allowed types: %(choices)s', metavar = '<type>', type = str, action = 'extend', required = True)

    ### Search
    parser_search.add_argument('--reuse', help='Reuse existing search results', action = 'store_true')
    # nt-based
    parser_search.add_argument('--nt_pident', help='Minimum nucleotide identity [%(default)i]', default = 60, metavar = ' ', type = float)
    parser_search.add_argument('--nt_evalue', help='Maximum nucleotide e-value [%(default)1.0e]', default = 1e-30, metavar = ' ', type = float)
    parser_search.add_argument('--nt_length', help='Minimum nucleotide length [%(default)i]', default = 1000, metavar = ' ', type = int)
    # aa-based
    parser_search.add_argument('--aa_pident', help='Minimum amino acid identity [%(default)i]', default = 75, metavar = ' ', type = int)
    parser_search.add_argument('--aa_evalue', help='Maximum amino acid e-value [%(default)1.0e]', default = 1e-5, metavar = ' ', type = float)
    parser_search.add_argument('--aa_qscovs', help='Minimum query and subject HSP coverages [%(default)i]', default = 90, metavar = ' ', type = float)
    # hmm-based
    parser_search.add_argument('--hmm_evalue', help='Maximum hmm e-value [%(default)1.0e]', default = 1e-10, metavar = ' ', type = float)
    # mmseqs2-based
    parser_search.add_argument('--mmseqs_sens', help='MMseqs2 search sensitivity [%(default).1f]', default = 5.7, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_evalue', help='MMseqs2 maximum e-value [%(default)1.0e]', default = 1e-5, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_pident', help='MMseqs2 minimum amino acid identity [%(default).2f]', default = 0.7, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_cov', help='MMseqs2 minimum amino acid coverage for clustering [%(default).2f]', default = 0.7, metavar = ' ', type = float)

    ### Evaluate
    parser_evaluate.add_argument('--combine', help='Combine all signals as a separate category', action='store_true')
    parser_evaluate.add_argument('--max_nt_gap', help='Maximum nucleotide distance between signal [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_nt_sig', help='Minimum nucleotide signal within region [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--max_aa_gap', help='Maximum distance between proteins with signal [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_aa_sig', help='Minimum number of proteins with signal within region [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_sig_frac', help='Minimum fraction of signal within region [%(default).2f]', default = 0.5, metavar = ' ', type = float)

    ### Validate
    parser_validate.add_argument('--checkv_env', help='Name of the conda env with CheckV installed if not installed system-wide.', metavar = ' ', type = str)
    parser_validate.add_argument('--checkv_db', help='Path to CheckV database or if no CHECKVDB was set.', metavar = ' ', type = str)
    parser_validate.add_argument('--checkv_len', help='Minimum length of all CheckV candidates to consider [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_validate.add_argument('--checkv_sig_frac', help='Minimum fraction of signal of all CheckV candidates to consider [%(default).2f]', default = 0.5, metavar = ' ', type = float)
    parser_validate.add_argument('--checkv_mq_len', help='Minimum length of CheckV *Medium-quality* candidates [%(default)i]', default = 20000, metavar = ' ', type = int)
    parser_validate.add_argument('--checkv_mq_sig_frac', help='Minimum fraction of signal of CheckV *Medium-quality* candidates [%(default).2f]', default = 0.85, metavar = ' ', type = float)
    parser_validate.add_argument('--artemis_plots', help='Generate Artemis plots for query records', action = 'store_true')
    parser_validate.add_argument('--sig_sources', help='Signal sources to consider for validation. Allowed types: %(choices)s [%(default)s] ', choices = SIG_SOURCES, default = 'all', metavar = ' ', type = str)
    
    args = parser.parse_args()
    if len(sys.argv) < 2:
        args.print_usage()
        sys.exit(1)

    # make working directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    args.logger = create_logger(os.path.join(args.outdir, 'log.txt'), log_level=args.logging)
    
    log_progress(f"Starting SigMa v.{__version__} analysis", loglevel = "INFO")
    log_progress(f"Command: " + ' '.join(sys.argv), loglevel = "INFO")
    log_progress(f"Output directory: {args.outdir}", loglevel = "INFO")
    log_progress(f"# of threads: {args.threads}", loglevel = "INFO")

    # check if reference and query datasets are the same length
    if len(args.reference) != len(args.reference_type):
        if len(args.reference_type) == 1:
            args.reference_type = args.reference_type * len(args.reference)

    if len(args.query) != len(args.query_type):
        if len(args.query_type) == 1:
            args.query_type = args.query_type * len(args.query)
    
    log_progress(f"# of reference datasets: {len(args.reference)}", loglevel = "INFO")
    log_progress(f"# of query datasets: {len(args.query)}", loglevel = "INFO")
    
    # declare SigMa object
    sigma = SigMa(args)

    # add reference datasets
    sigma.prepare_targets()

    # add query datasets
    sigma.prepare_queries()

    # search query datasets
    sigma.search_queries()

    # evaluate query datasets
    sigma.evaluate_signals()

    # filter merged regions
    candidate_regions = sigma.filter_regions(sig_group = 'merged')
    if len(candidate_regions) > 0:
        # list candidate regions
        log_progress(f"Selected {len(candidate_regions)} unique candidate regions", msglevel = 0, loglevel = "INFO")

        # write cadidate regions
        sigma.write_regions(candidate_regions, 'candidate')

        # write Artemis plot files
        if args.artemis_plots: sigma.write_artemis_plots()

        # run CheckV
        sigma.run_checkv()

        # filter only high-quality regions
        sigma.filter_hq_regions()

        # list high-quality regions
        log_progress(f"Selected {len(sigma.hq_regions)} high-quality regions", msglevel = 0, loglevel = "INFO")
        sigma.list_regions(sigma.hq_regions)

        # write summary
        sigma.write_summary()

        # write verified regions
        sigma.write_regions(sigma.hq_regions, 'verified', format = ['fasta', 'genbank'])
    else:
        log_progress('No candidate regions found.')

def run():
    """Run SigMa analysis."""
    main()