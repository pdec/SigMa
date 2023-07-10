"""
A module to read various files
"""

from .utils import log_progress
from typing import List, Dict
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import binascii
import sys
import gzip
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


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
        log_progress(
            f"There was an error opening {file_path}", msglevel=0, loglevel="ERROR")
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
        log_progress(
            f"There was an error opening {file_path}", msglevel=0, loglevel="ERROR")
        sys.exit(20)

    return list(SeqIO.parse(handle, "fasta"))


def read_batch_file(batch_file: str, verbose=False) -> Dict[str, List[str]]:
    """
    Read a batch file and return a list of file paths
    Format of batch file is:
    'batch_name1
    file1
    file2
    batch_name2
    file3'
    ...
    :param batch_file: the batch file to read
    :param verbose: print verbose output
    :return: list of file paths

    """
    try:
        if is_gzip_file(batch_file):
            handle = gzip.open(batch_file, 'rt')
        else:
            handle = open(batch_file, 'r')
    except IOError as e:
        log_progress(
            f"There was an error opening {batch_file}", msglevel=0, loglevel="ERROR")
        sys.exit(20)

    batch_dict = {}
    query = []
    for line in handle:
        line = line.strip()
        if line.startswith('#'):
            continue
        if line.startswith('batch'):
            batch_name = line
            batch_dict[batch_name] = []
        else:
            batch_dict[batch_name].append(line)
            query.append(line)

    for batch, files in batch_dict.items():
        log_progress(
            f"Batch {batch} contains {len(files)} files", msglevel=1, loglevel="INFO")

    return query, batch_dict
