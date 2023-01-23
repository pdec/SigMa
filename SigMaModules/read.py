"""
A module to read various files
"""

import binascii
import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import List
from .utils import log_progress


def is_gzip_file(file_path: str) -> bool:
    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param file_path: the path to file to test
    :return: True if the file is gzip compressed else false
    """
    with open(file_path, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'


def parse_genbank(file_path: str, verbose=False) -> List[SeqRecord]:
    """
    Parse a genBank file into a list of SeqRecord objects
    :param file_path: the path to GenBank file
    :return: list of SeqRecords
    """
    try:
        if is_gzip_file(file_path):
            handle = gzip.open(file_path, 'rt')
        else:
            handle = open(file_path, 'r')
    except IOError as e:
        log_progress(f"There was an error opening {file_path}", msglevel = 0, loglevel = "ERROR")
        sys.exit(20)

    return list(SeqIO.parse(handle, "genbank"))


def parse_fasta(file_path: str, verbose=False) -> List[SeqRecord]:
    """
    Parse a fasta file and return a list of SeqRecord objects
    :param file_path: the fasta file to parse
    :param verbose: print verbose output
    """
    try:
        if is_gzip_file(file_path):
            handle = gzip.open(file_path, 'rt')
        else:
            handle = open(file_path, 'r')
    except IOError as e:
        log_progress(f"There was an error opening {file_path}", msglevel = 0, loglevel = "ERROR")
        sys.exit(20)

    return list(SeqIO.parse(handle, "fasta"))
