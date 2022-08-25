__author__ = "Przemysław Decewicz"
__email__ = "p.decewicz@uw.edu.pl"
__status__ = "Development"


import sys
import os
import argparse
import subprocess
import re

def setup(args):
    """
    Setup SigMa reference databases and query sequences
    :param args:
    """
    pass

def search(args):
    """
    Search for signatures in a reference sequence
    :param args:
    """
    pass

def evaluate(args):
    """
    Evaluate the significance of a phage signal
    :param args:
    """
    pass

def validate(args):
    """
    Validate a signature
    :param args:
    """
    pass

def run(args):
    """
    Run SigMa
    :param args:
    """
    pass

def main():
    parser = argparse.ArgumentParser(description='SigMa v- a tool for iterative mapping of phage signatures and prophage recovery.', epilog='Example: SigMa.py run -r reference.fasta -q query.fasta -o output_dir')
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')

    parser_run = subparsers.add_parser('run', help='Complete SigMa run')
    parser_setup = subparsers.add_parser('setup', help='Setup SigMa reference databases and query sequences')
    parser_search = subparsers.add_parser('search', help='Search for signatures in query sequence')
    parser_evaluate = subparsers.add_parser('evaluate', help='Evaluate the observed phage signatures')
    parser_validate = subparsers.add_parser('validate', help='Validate signatures')
    
    args = parser.parse_args()
    if args.subparser_name == None:
        parser.print_help()
        sys.exit(1)
    elif args.subparser_name == 'setup':
        setup(args)
    elif args.subparser_name == 'search':
        search(args)
    elif args.subparser_name == 'evaluate':
        evaluate(args)
    elif args.subparser_name == 'validate':
        validate(args)
    elif args.subparser_name == 'run':
        run(args)  


if __name__ == "__main__":
    main()