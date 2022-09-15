
"""
Module to work of features in input GenBank or FASTA files
"""
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from typing import List, Union
import re


def format_seq(seq : Union[str, Seq, SeqRecord], width : int = 60) -> str:
    """
    Splits seq across multiple lines based on width.
    :param seq: sequence to format
    :param width: width of the line
    """

    return re.sub(r'(.{' + re.escape(str(width)) + '})', '\\1\n', str(seq).upper(), 0, re.DOTALL).strip()

def get_features_of_type(seqiorec: SeqRecord, ftype: str) -> List[SeqFeature]:
    """
    Get features of a given type from SeqRecord
    :param seqiorec: a SeqRecord object
    :param ftype: type of a feature
    :return:
    """

    flist = []
    for fcnt, feature in enumerate(seqiorec.features, 1):
        if feature.type == ftype:
            if ftype == 'CDS':
                if 'translation' not in feature.qualifiers:
                    feature.qualifiers['translation'] = [feature.extract(seqiorec.seq).translate(table = 11, to_stop = True)]
                if 'protein_id' not in feature.qualifiers:
                    feature.qualifiers['protein_id'] = [f'{seqiorec.id}_ft_{fcnt:06d}']
            flist.append(feature)
    
    return flist

def cds_to_fasta(cds: SeqFeature) -> str:
    """
    Get amino acid sequence of CDS feature
    :param cds: a SeqFeature object
    :return: amino acid sequence in FASTA format
    """
    
    return f">{cds.qualifiers['protein_id'][0]}\n{format_seq(cds.qualifiers['translation'][0])}\n"