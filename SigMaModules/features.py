
"""
Module to work of features in input GenBank or FASTA files
"""
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from typing import List,Union

def get_features_of_type(seqiorec: Union[SeqRecord, List[SeqRecord]], ftype: str) -> List[SeqFeature]:
    """
    Get features of a given type from SeqRecord
    :param seqiorec: a SeqRecord object
    :param ftype: type of a feature
    :return:
    """

    if isinstance(seqiorec, SeqRecord):
        seqiorec = [seqiorec]
    flist = []
    for seqrec in seqiorec:
        for fcnt, feature in enumerate(seqrec.features, 1):
            if feature.type == ftype:
                if ftype == 'CDS':
                    if 'translation' not in feature.qualifiers:
                        feature.qualifiers['translation'] = [feature.extract(seqrec.seq).translate(table = 11, to_stop = True)]
                    if 'protein_id' not in feature.qualifiers:
                        feature.qualifiers['protein_id'] = [f'{seqrec.id}_ft_{fcnt:06d}']
                    if 'record_id' not in feature.qualifiers:
                        feature.qualifiers['record_id'] = [seqrec.id]
                flist.append(feature)
    
    return flist

def get_cds_header_and_aa(cds: SeqFeature) -> List[str]:
    """
    Get amino acid sequence of CDS feature
    :param cds: a SeqFeature object
    :return: a list of header and amino acid sequence
    """
    
    return [f"{cds.qualifiers['record_id'][0]}|{cds.qualifiers['protein_id'][0]}|{int(cds.location.nofuzzy_start)}..{cds.location.nofuzzy_end}..{cds.strand}", cds.qualifiers['translation'][0]]

def get_feature_coords(feature: SeqFeature) -> List[int]:
    """
    Get coordinates of a feature
    :param feature: a SeqFeature object
    :return: a list of coordinates
    """
    return [feature.location.nofuzzy_start, feature.location.nofuzzy_end]