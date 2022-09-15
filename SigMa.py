__author__ = "Przemysław Decewicz"
__email__ = "p.decewicz@uw.edu.pl"
__status__ = "Development"


import sys
import os
import argparse
import logging

from version import __version__


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
    parser.add_argument('-o', '--output', help='Output directory', metavar = '<path>', required=True)
    parser.add_argument('-t', '--threads', help='Number of threads to use for data preparation [%(default)i]', default = 4, metavar = '<num>', type = int)
    parser.add_argument('-v', '--version', help='Show program\'s version and exit.', action='version', version='%(prog)s 0.1')
    
    ### Setup
    parser_setup.add_argument('-r', '--reference', nargs = '+', help='Reference dataset(s)', metavar = '<path>', type = str, action = 'append', required = True)
    parser_setup.add_argument('-R', '--reference_type', nargs = '+', choices = REFERENCE_TYPES, help='Reference dataset type(s). Allowed types: %(choices)s', metavar = '<type>', type = str, action = 'append', required = True)
    parser_setup.add_argument('-q', '--query', nargs = '+', help='Query dataset(s)', metavar = '<path>', type = str, action = 'append', required = True)
    parser_setup.add_argument('-Q', '--query_type', nargs = '+', choices = QUERY_TYPES, help='Query dataset type(s). Allowed types: %(choices)s', metavar = '<type>', type = str, action = 'append', required = True)

    ### Search
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
    parser_search.add_argument('--mmseqs_sensivitiy', help='MMseqs2 search sensitivity [%(default).1f]', default = 5.7, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_evalue', help='MMseqs2 maximum e-value [%(default)1.0e]', default = 1e-5, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_pident', help='MMseqs2 minimum amino acid identity [%(default).2f]', default = 0.7, metavar = ' ', type = float)
    parser_search.add_argument('--mmseqs_coverage', help='MMseqs2 minimum amino acid coverage for clustering [%(default).2f]', default = 0.7, metavar = ' ', type = float)

    ### Evaluate
    parser_evaluate.add_argument('--max_nt_gap', help='Maximum nucleotide distance between signal [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_nt_sig', help='Minimum nucleotide signal within region [%(default)i]', default = 5000, metavar = ' ', type = int)
    parser_evaluate.add_argument('--max_aa_gap', help='Maximum amino acid distance between signal [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_aa_sig', help='Minimum amino acid signal within region [%(default)i]', default = 5, metavar = ' ', type = int)
    parser_evaluate.add_argument('--min_sig_frac', help='Minimum fraction of signal within region [%(default).2f]', default = 0.5, metavar = ' ', type = float)

    ### Validate
    parser_validate.add_argument('--checkv_db_path', help='Path to CheckV database or if no CHECKVDB was set.', metavar = ' ', type = str)

    
    args = parser.parse_args()
    if len(sys.argv) < 2:
        args.print_usage()
        sys.exit(1)

    # make working directory
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # setup logger
    logFormatter = logging.Formatter(fmt='%(asctime)s [%(levelname)-8.8s]: %(message)s', datefmt='%Y-%m-%d %I:%M:%S')#, filename = os.path.join(args.output, 'log.txt'))
    logger = logging.getLogger()
    fileHandler = logging.FileHandler(os.path.join(args.output, 'log.txt'))
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    logger.setLevel(logging.DEBUG)

    logger.info(f"Starting SigMa analysis")
    logger.debug(f"Command: " + ' '.join(sys.argv))
    logger.error(f"Output directory: {args.output}")
    logger.warning(f"Number of threads: {args.threads}")
    logger.critical(f"Version: {__version__}")

    # check if reference and query datasets are the same length
    if len(args.reference) != len(args.reference_type) or len(args.query) != len(args.query_type):
        logger.error('Reference and query datasets must be the same length.')
        sys.exit(1)
    
    # prepare reference datasets
    logger.info(f"Preparing reference datasets...")
    ref_datasets = []


    



if __name__ == "__main__":
    main()