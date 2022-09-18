
import sys
import os
import argparse

from .query import SigMaQuery
from .reference import SigMaRefNT, SigMaRefAA, SigMaRefHMM, SigMaRefMMSEQS
from .utils import create_logger, log_progress, list_databases
from .write import write_fasta
from .version import __version__


def main():

    class CustomHelpFormatter(argparse.HelpFormatter):
        def __init__(self, prog):
            super().__init__(prog, max_help_position=60, width=120)

        def _format_action_invocation(self, action):
            if not action.option_strings or action.nargs == 0:
                return super()._format_action_invocation(action)
            default = self._get_default_metavar_for_optional(action)
            args_string = self._format_args(action, default)
            return ', '.join(action.option_strings) + ' ' + args_string

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

    ### Run
    parser.add_argument('-o', '--outdir', help='Output directory', metavar = '<path>', required=True)
    parser.add_argument('-t', '--threads', help='Number of threads to use for data preparation [%(default)i]', default = 4, metavar = '<num>', type = int)
    parser.add_argument('-v', '--version', help='Show program\'s version and exit.', action='version', version='%(prog)s 0.1')
    
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
    parser_evaluate.add_argument('--max_nt_gap', help='Maximum nucleotide distance between signal [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_nt_sig', help='Minimum nucleotide signal within region [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--max_aa_gap', help='Maximum distance between proteins with signal [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_aa_sig', help='Minimum number of proteins with signal within region [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_sig_frac', help='Minimum fraction of signal within region [%(default).2f]', default = 0.5, metavar = ' ', type = float)

    ### Validate
    parser_validate.add_argument('--checkv_db_path', help='Path to CheckV database or if no CHECKVDB was set.', metavar = ' ', type = str)

    
    args = parser.parse_args()
    if len(sys.argv) < 2:
        args.print_usage()
        sys.exit(1)

    # make working directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    args.logger = create_logger(os.path.join(args.outdir, 'log.txt'))
    
    log_progress(f"Starting SigMa analysis", loglevel = "INFO")
    log_progress(f"Command: " + ' '.join(sys.argv), loglevel = "ERROR")
    log_progress(f"Output directory: {args.outdir}", loglevel = "WARNING")
    log_progress(f"Number of threads: {args.threads}", loglevel = "CRITICAL")
    log_progress(f"Version: {__version__}", loglevel = "DEBUG")

    # check if reference and query datasets are the same length
    if len(args.reference) != len(args.reference_type) or len(args.query) != len(args.query_type):
        log_progress('Reference and query datasets must be the same length.')
        sys.exit(1)

    # make reference directory
    ref_dir = os.path.join(args.outdir, 'reference')
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
    
    # prepare reference databases
    log_progress(f"Preparing reference datasets...")
    refs = []
    for ref_dataset_path, ref_type in zip(args.reference, args.reference_type):
        if ref_type == 'fasta_nt': refs.append(SigMaRefNT(ref_dataset_path, ref_type, ref_dir))
        elif ref_type == 'fasta_aa': refs.append(SigMaRefAA(ref_dataset_path, ref_type, ref_dir))
        elif ref_type == 'hmm': refs.append(SigMaRefHMM(ref_dataset_path, ref_type, ref_dir))
        elif ref_type == 'mmseqs_db': refs.append(SigMaRefMMSEQS(ref_dataset_path, ref_type, ref_dir))

    # prepare query datasets
    log_progress(f"Preparing query datasets...")
    queries = []
    for query_type, query_dataset_path in zip(args.query_type, args.query):
        queries.append(SigMaQuery(query_dataset_path, query_type))
    
    # make query directory
    query_dir = os.path.join(args.outdir, 'query')
    if not os.path.exists(query_dir):
        os.makedirs(query_dir)

    # parse query datasets
    log_progress(f"Parsing query datasets...")
    query_nt_path = os.path.join(query_dir, 'query_nt.fasta')
    query_aa_path = os.path.join(query_dir, 'query_aa.fasta')

    # write FASTA files
    log_progress(f"Writing query FASTA files")
    write_fasta([q.records_to_fasta() for q in queries], query_nt_path)
    write_fasta([q.cdss_to_fasta() for q in queries], query_aa_path)

    # search query datasets
    log_progress(f"Searching query datasets...")
    search_dir = os.path.join(args.outdir, 'search')
    if not os.path.exists(search_dir):
        os.makedirs(search_dir)

    for ref in sorted(refs, key = lambda x: x.type):
        search_out = ref.get_output_path(search_dir)
        if args.reuse and os.path.exists(search_out):
            log_progress(f"reusing {search_out}", msglevel=1)
            continue
        if ref.type == 'fasta_nt':
            log_progress(f"searching query nucleotide sequences against {ref.name}...", msglevel=1)
            ref.search(query_nt_path, search_out, args.nt_evalue, args.nt_pident, args.threads)
        elif ref.type == 'fasta_aa':
            log_progress(f"searching query amino acid sequences against {ref.name}...", msglevel=1)
            ref.search(query_aa_path, search_out, args.aa_evalue, args.aa_pident, args.aa_qscovs, args.threads)
        elif ref.type == 'hmm':
            log_progress(f"searching query amino acid sequences against {ref.name}...", msglevel=1)
            ref.search(query_aa_path, search_out, args.hmm_evalue, args.threads)
        elif ref.type == 'mmseqs_db':
            log_progress(f"searching query nucleotide sequences against {ref.name}...", msglevel=1)
            ref.search(query_aa_path, search_out, args.mmseqs_sens, args.mmseqs_evalue, args.mmseqs_pident, args.mmseqs_cov, args.threads)

    # parse search results
    log_progress(f"Parsing search results...")
    log_progress(f"NOT IMPLEMENTED YET")
    

def run():
    """Run SigMa analysis."""
    main()