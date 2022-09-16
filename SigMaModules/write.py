from typing import Union, List, Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re

from .utils import colours

__author__ = 'Przemyslaw Decewicz'


def format_seq(seq : Union[str, Seq, SeqRecord], width : int = 60) -> str:
    """
    Splits seq across multiple lines based on width.
    :param seq: sequence to format
    :param width: width of the line
    """

    return re.sub(r'(.{' + re.escape(str(width)) + '})', '\\1\n', str(seq).upper(), 0, re.DOTALL).strip()

def write_fasta(seqs : Union[List, Dict], file_path : str, width : int = 60):
    """
    Write sequences to FASTA file
    :param seqs: list of SeqRecord objects or dictionary of SeqRecord objects
    :param file_path: path to FASTA file
    :return:
    """

    with open(file_path, 'w') as f:
        if isinstance(seqs, list):
            for seq in seqs:
                if isinstance(seq, SeqRecord):
                    f.write(f">{seq.id}\n{format_seq(seq.seq)}\n")
                else:
                    f.write(f">{seq[0]}\n{format_seq(seq[1])}\n")
        elif isinstance(seqs, dict):
            for header, seq in seqs.items():
                if isinstance(seq, SeqRecord):
                    seq = seq.seq
                f.write(f">{header}\n{format_seq(seq)}\n")
        else:
            raise TypeError('seqs must be a list tuples or dictionary of sequence header and sequence/SeqRecord objects')

def write_df_to_artemis(df : pd.DataFrame,  filename : str,  colours : List[str] = colours()):
    """
    Write and clean dataframe to Artemis graph file
    :param df: Pandas dataframe with base position in first column and values in other columns
    :param filename: name of the file to write to
    :return:
    """
    df.to_csv(filename, sep=' ', index=False)
    with open(filename, 'r') as f:
        lines = f.readlines()

    with open(filename, 'w') as f:
        header = lines[0].split()[1:]
        f.write('# BASE {}\n'.format(' '.join(header)))

        colours_rgb = []
        colours_names = []
        for rgb in colours[:len(header)]:
            colours_rgb.append(':'.join(rgb[0].split()))
            colours_names.append(rgb[1])
        f.write('# colour {}\n'.format (' '.join(colours_rgb)))
        f.write('# label {}\n'.format(' '.join(header))) # the label must go after colour
        f.write('# name {}\n'.format (' '.join(colours_names)))
        for line in lines[1:]:
            if len(line.split()) == 1:
                continue
            f.write(line)

