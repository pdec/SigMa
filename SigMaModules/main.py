
import sys
import os
import argparse
import glob
from typing import List

from .models import SigMa
from .read import read_batch_file
from .utils import call_process, create_logger, log_progress, make_batches, post_batch, read_batches_done, CustomHelpFormatter
from .version import __version__

def run_sigma(args):
    """
    Runs a single SigMa instance
    """
    # declare SigMa object
    sigma = SigMa(args)

    # add reference datasets
    sigma.prepare_targets()

    # analyze query datasets one by one
    sigma.handle_queries()

    # filter merged regions
    candidate_regions = sigma.filter_regions(sig_groups = 'merged', sig_sources = 'merged')
    if len(candidate_regions) > 0:
        # list candidate regions
        log_progress(f"Considering {len(candidate_regions)} candidate regions for validation", msglevel = 0, loglevel = "INFO")

        # write cadidate regions
        sigma.write_regions(candidate_regions, 'candidate')

        # write Artemis plot files
        if args.artemis_plots: sigma.write_artemis_plots()

        if args.skip_checkv:
            log_progress('CheckV validation skipped.', msglevel = 0, loglevel = "INFO")
            
            # write summary
            sigma.write_summary()

            # write candidate regions
            sigma.write_regions(candidate_regions, 'candidate', format = ['fasta', 'genbank'])
        else:
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

    EXT_TOOLS = [
        'phispy',
    ]

    SIG_GROUPS = [
        'all',
        'combined',
        'nt_based',
        'aa_based'
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
    parser.add_argument('-c', '--batch_size', help='Number of query files to process in each batch [%(default)i]', default = 0, metavar = '<num>', type = int)
    parser.add_argument('-B', '--batches_file', help='File with batches to run on input files. Each batch starts with a "batch_" followed by file query file names one in each line', metavar = '<path>', type = str)
    parser.add_argument('-t', '--threads', help='Number of threads to use for data preparation [%(default)i]', default = 4, metavar = '<num>', type = int)
    parser.add_argument('-v', '--version', help='Show program\'s version and exit.', action='version', version='%(prog)s 0.1')
    parser.add_argument('-l', '--logging', help='Output logging level [%(default)s]. Allowed options: %(choices)s', default = 'INFO', choices = LOGGING, metavar = ' ', type = str)
    parser.add_argument('--resume', help='Resume works on the last processed batch. Based on batch_status.tsv entries.', action = 'store_true')
    
    ### Setup
    parser_setup.add_argument('-r', '--reference', nargs = '+', help='Reference dataset(s)', metavar = '<path>', type = str, action = 'extend', required = True)
    parser_setup.add_argument('-R', '--reference_type', nargs = '+', choices = REFERENCE_TYPES, help='Reference dataset type(s). Allowed types: %(choices)s', metavar = '<type>', type = str, action = 'extend', required = True)
    parser_setup.add_argument('-i', '--indir', help='Input directory with query datasets', metavar = '<path>', type = str)
    parser_setup.add_argument('-q', '--query', nargs = '+', help='Query dataset(s)', metavar = '<path>', type = str, action = 'extend')
    parser_setup.add_argument('-Q', '--query_type', nargs = '+', choices = QUERY_TYPES, help='Query dataset type(s). Allowed types: %(choices)s', metavar = '<type>', type = str, action = 'extend', default = ['genbank'])

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
    # parser_search.add_argument('--hmm_evalue', help='Maximum hmm e-value [%(default)1.0e]', default = 1e-10, metavar = ' ', type = float)
    # mmseqs2-based
    parser_search.add_argument('--mmseqs_sens', help='MMseqs2 search sensitivity [%(default).1f]', default = 5.7, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_evalue', help='MMseqs2 maximum e-value [%(default)1.0e]', default = 1e-5, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_pident', help='MMseqs2 minimum amino acid identity [%(default).2f]', default = 0.7, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_cov', help='MMseqs2 minimum amino acid coverage for clustering [%(default).2f]', default = 0.7, metavar = ' ', type = float)
    # other tools
    parser_search.add_argument('--ext_tools', help='External prophage identification programs to use. Choices: %(choices)s', choices = EXT_TOOLS, nargs = '+', action = 'extend', metavar = ' ', type = str)

    ### Evaluate
    parser_evaluate.add_argument('--combine', help='Combine all signals as a separate category', action='store_true')
    parser_evaluate.add_argument('--max_nt_gap', help='Maximum nucleotide distance between signal [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_nt_sig', help='Minimum nucleotide signal within region [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--max_aa_gap', help='Maximum distance between proteins with signal [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_aa_sig', help='Minimum number of proteins with signal within region [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_sig_frac', help='Minimum fraction of signal within region [%(default).2f]', default = 0.5, metavar = ' ', type = float)

    ### Validate
    parser_validate.add_argument('--skip_checkv', help='Do not verify candidate regions with CheckV', action = 'store_true')
    parser_validate.add_argument('--checkv_env', help='Name of the conda env with CheckV installed if not installed system-wide.', metavar = ' ', type = str)
    parser_validate.add_argument('--checkv_db', help='Path to CheckV database or if no CHECKVDB was set.', metavar = ' ', type = str)
    parser_validate.add_argument('--checkv_len', help='Minimum length of all CheckV candidates to consider [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_validate.add_argument('--checkv_sig_frac', help='Minimum fraction of signal of all CheckV candidates to consider [%(default).2f]', default = 0.5, metavar = ' ', type = float)
    parser_validate.add_argument('--checkv_mq_len', help='Minimum length of CheckV *Medium-quality* candidates [%(default)i]', default = 20000, metavar = ' ', type = int)
    parser_validate.add_argument('--checkv_mq_sig_frac', help='Minimum fraction of signal of CheckV *Medium-quality* candidates [%(default).2f]', default = 0.85, metavar = ' ', type = float)
    parser_validate.add_argument('--artemis_plots', help='Generate Artemis plots for query records', action = 'store_true')
    parser_validate.add_argument('--sig_groups', help='Signal groups to consider for validation. Allowed types: %(choices)s [%(default)s] ', choices = SIG_GROUPS, default = 'all', metavar = ' ', type = str, nargs = '+')
    
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
    if args.indir: log_progress(f"Input directory: {args.indir}", loglevel = "INFO")
    log_progress(f"# of threads: {args.threads}", loglevel = "INFO")

    # check input
    if not args.query and not args.batches_file and not args.indir:
        log_progress(f"You must provide at least one query dataset with either --query, --batches_file or --indir", loglevel = "ERROR")
        exit(1)
    elif args.query and args.batches_file:
        log_progress(f"You must provide either --query or --batches_file, not both", loglevel = "ERROR")
        exit(1)
    elif args.query and args.indir:
        log_progress(f"You must provide either --query or --indir, not both", loglevel = "ERROR")
        exit(1)
    elif args.batches_file and args.indir:
        log_progress(f"You must provide either --batches_file or --indir, not both", loglevel = "ERROR")
        exit(1)
    elif args.indir:
        args.query = glob.glob(os.path.join(args.indir, '*'))
        log_progress(f"Found {len(args.query)} query datasets in {args.indir}", loglevel = "DEBUG")

    # batches of query datasets
    if (args.batches > 1) or (args.batches_file) or (args.batch_size > 0):
        if args.batches > 1 or args.batch_size > 0:
            log_progress(f"Splitting query datasets based on input order", loglevel = "INFO")
            args.query_batches = make_batches(args.query, args.batches, args.batch_size)
            if args.batch_size > 0:
                args.batches = len(args.query_batches)
        elif args.batches_file:
            log_progress("Reading batches", loglevel = "INFO")
            args.query, args.query_batches = read_batch_file(args.batches_file)
            args.batches = len(args.query_batches)
        log_progress(f"# of query batches: {args.batches}", loglevel = "INFO")
      
    # check reference input
    if len(args.reference) != len(args.reference_type):
        if len(args.reference_type) == 1:
            args.reference_type = args.reference_type * len(args.reference)

    # check query input
    if len(args.query) != len(args.query_type):
        if len(args.query_type) == 1:
            args.query_type = args.query_type * len(args.query)
    
    log_progress(f"# of reference datasets: {len(args.reference)}", loglevel = "INFO")
    log_progress(f"# of query datasets: {len(args.query)}", loglevel = "INFO")
        
    # check if all files exist and are not empty
    log_progress("Checking input query files", loglevel = "INFO")
    for infile in args.query:
        if not os.path.exists(infile):
            log_progress(f"Input file {infile} does not exist. Quiting.", loglevel = "ERROR")
            sys.exit(1)
        if os.path.getsize(infile) == 0:
            log_progress(f"Input file {infile} is empty. Quiting.", loglevel = "ERROR")
            sys.exit(1)
    
    # run analyses
    if args.batches > 1:
        log_progress("Running in batch mode.", loglevel = "INFO")

        bi = 0
        args.main_outdir = args.outdir
        args.main_query = args.query
        args.main_checkv_db = args.checkv_db
        args.batch_refs = args.reference
        args.batch_ref_types = args.reference_type
        args.prepdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts', 'prepdb.py')

        # track which batches have been completed
        args.batches_done = read_batches_done(args.main_outdir)

        for batch, batch_query in args.query_batches.items():
            bi += 1
            log_progress(f">>> Batch [{bi}/{args.batches}]: {batch}", loglevel = "INFO")
            # check if batch has already been completed
            if (args.resume) and (batch in args.batches_done):
                log_progress(f"Batch {batch} has already been completed. Skipping.", loglevel = "INFO")
                continue
            # update output directory
            args.outdir = os.path.join(args.main_outdir, batch)
            log_progress(f"Output directory: {args.outdir}", loglevel = "DEBUG")
            # update query datasets
            args.query = batch_query
            # run sigma
            run_sigma(args)
            # write batch status
            verified, candidates = post_batch(args.main_outdir, args.outdir, batch, bi)
            verified_gbs = os.path.join(args.main_outdir, f"all_verified.gb")
            verified_db = os.path.join(args.main_outdir, f"batch_dbs")
            verified_nt_derep = os.path.join(verified_db, f"all_verified.nt.fasta")
            verified_aa_derep = os.path.join(verified_db, f"all_verified.aa.fasta")
            verified_checkv_db = os.path.join(verified_db, f"checkv_db")

            log_progress(f">>> Batch [{bi}/{args.batches}]: {candidates} candidates and {verified} verified regions", loglevel = "INFO")

            # whenever there are verified regions, update the reference databases
            if verified > 0:
                log_progress(f"Using verified regions as new references", loglevel = "INFO")
                # call prepdb
                cmd = f"python {args.prepdb} --gb {verified_gbs} --dbdir {verified_db} --threads {args.threads} --tmp {os.path.join(verified_db, 'tmp')} --logging {args.logging} --simple_log"
                call_process(cmd, "prepdb")
                
                # update checkv database
                checkv_env = '' if not args.checkv_env else f"conda run -n {args.checkv_env}"
                cmd = f"{checkv_env} checkv update_database {args.main_checkv_db} {verified_checkv_db} {verified_nt_derep} --threads {args.threads} --quiet --restart"
                call_process(cmd, program="checkv")
                
                # update reference databases
                if verified_nt_derep not in args.reference and verified_aa_derep not in args.reference:
                    # nt
                    args.reference.append(verified_nt_derep)
                    args.reference_type.append('fasta_nt')
                    # aa
                    args.reference.append(verified_aa_derep)
                    args.reference_type.append('fasta_aa')
                    # checkv - change the database if updated at least once
                    args.checkv_db = verified_checkv_db

            # update batch status
            args.batches_done.append(batch)
            
    else:
        run_sigma(args)

def run():
    """Run SigMa analysis."""
    main()