
"""
Module to work of features in input GenBank or FASTA files
"""
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from typing import List

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

def get_cds_aa(cds: SeqFeature) -> List[str, str]:
    """
    Get amino acid sequence of CDS feature
    :param cds: a SeqFeature object
    :return: a list of header and amino acid sequence
    """
    
    return [cds.qualifiers['protein_id'][0], cds.qualifiers['translation'][0]]