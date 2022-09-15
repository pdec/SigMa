"""
A module to read various files
"""

import binascii
import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import List


def is_gzip_file(filename: str) -> bool:
    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param filename: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(filename, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'


def parse_genbank(filename: str, verbose=False) -> List[SeqRecord]:
    """
    Parse a genbank file and return a Bio::Seq object
    """
    try:
        if is_gzip_file(filename):
            handle = gzip.open(filename, 'rt')
        else:
            handle = open(filename, 'r')
    except IOError as e:
        print(f"There was an error opening {filename}", file=sys.stderr)
        sys.exit(20)

    return list(SeqIO.parse(handle, "genbank"))


def parse_fasta(filename: str, verbose=False) -> List[SeqRecord]:
    """
    Parse a fasta file and return a list of SeqRecord objects
    :param filename: the fasta file to parse
    :param verbose: print verbose output
    """
    try:
        if is_gzip_file(filename):
            handle = gzip.open(filename, 'rt')
        else:
            handle = open(filename, 'r')
    except IOError as e:
        print(f"There was an error opening {filename}", file=sys.stderr)
        sys.exit(20)

    return list(SeqIO.parse(handle, "fasta"))
