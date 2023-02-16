from typing import Union, List, Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import gzip
import json
import re

from .utils import colours, log_progress

__author__ = 'Przemyslaw Decewicz'


def format_seq(seq : Union[str, Seq, SeqRecord], width : int = 60) -> str:
    """
    Splits seq across multiple lines based on width.
    :param seq: sequence to format
    :param width: width of the line
    """

    return re.sub(r'(.{' + re.escape(str(width)) + '})', '\\1\n', str(seq).upper(), 0, re.DOTALL).strip()

def write_fasta(
    seqs : Union[List[str], List[SeqRecord], Dict], file_path : str, width : int = 60):
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
                    f.write(f"{seq}")
        elif isinstance(seqs, dict):
            for header, seq in seqs.items():
                if isinstance(seq, SeqRecord):
                    seq = seq.seq
                f.write(f">{header}\n{format_seq(seq)}\n")
        else:
            raise TypeError('seqs must be a list tuples or dictionary of sequence header and sequence/SeqRecord objects')

def write_df_to_artemis(df : pd.DataFrame,  file_path : str,  colours : List[str] = colours()):
    """
    Write and clean dataframe to Artemis graph file
    :param df: Pandas dataframe with base position in first column and values in other columns
    :param file_path: name of the file to write to
    :return:
    """

    with gzip.open(file_path, 'wt', compresslevel = 6) as f:

        header = df.columns
        f.write('# BASE {}\n'.format(' '.join(header)))

        colours_rgb = []
        colours_names = []
        for rgb in colours[:len(header)]:
            colours_rgb.append(':'.join(rgb[0].split()))
            colours_names.append(rgb[1])
        f.write('# colour {}\n'.format (' '.join(colours_rgb)))
        f.write('# label {}\n'.format(' '.join(header))) # the label must go after colour
        f.write('# name {}\n'.format (' '.join(colours_names)))

    # add dataframe content to the file
    df.dropna(axis=0, how='all').to_csv(file_path, sep=' ', index=True, compression={'method': 'gzip', 'compresslevel': 6, 'mtime': 1}, mode='a')

def write_df_to_plotly(df : pd.DataFrame, file_path : str):
    """
    Write and clean dataframe to Plotly graph file
    :param df: Pandas dataframe with base position in first column and values in other columns
    :param file_path: name of the file to write to
    :return:
    """

    # NOTE
    # single `null` value needs to be inserted between each consecutive pair of non-null 
    # values to make breaks in the line in both x and y axes in plotly.js
    
    j = {"data": [], "layout": {}} # initialise the plotly json object
    for col in df.columns:
        tdf = df[[col]]
        tdf = tdf.loc[~(tdf==0).all(axis=1)]
        j["data"].append({"x": tdf.index.to_list(), "y": tdf[col].to_list(), "name": col, "type": "line"})

    with open(file_path, 'w') as f: 
        json.dump(j, f)


