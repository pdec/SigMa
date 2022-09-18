"""
Module defining SigMa query class
"""

from .read import parse_fasta, parse_genbank
from .features import get_features_of_type
from .utils import log_progress
from .write import format_seq

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from typing import List, Dict, OrderedDict, Tuple

class SigMaQuery:
    """This is class for storing information about the query input"""
    
    def __init__(self, file_path : str, input_type : str) -> None:
        """
        SigMaQuery object constructor
        """
        self.file_path = file_path
        self.type = input_type
        self.records = []
        self.cdss = []
        self.singal = {}
    
        if self.type == 'fasta':
            self.records = parse_fasta(self.file_path)
        elif self.type == 'genbank':
            for record in parse_genbank(self.file_path):
                self.records.append(record)
                rec_cdss = get_features_of_type(record, 'CDS')
                self.cdss.extend(rec_cdss)
                log_progress(f"{record.id}: {len(rec_cdss)} CDSs", msglevel = 1)

    def __str__(self):
        """
        SigMaQuery string representation
        """

        return f"Query: {self.file_path} [{self.type}]: {len(self.records)} record{'s' if len(self.records) > 1 else ''} and {len(self.cdss)} CDS{'s' if len(self.cdss) > 1 else ''}"

        
    def records_to_fasta(self) -> str:
        """
        Return a fasta string of the query records
        :return: fasta string
        """
        fasta = ""
        for record in self.records:
            fasta += f">{record.id}|{len(record.seq)}\n{format_seq(record.seq)}"

        return fasta

    def cdss_to_fasta(self) -> str:
        """
        Return a fasta string of the query CDSs
        :return: fasta string
        """
        fasta = ""
        for cds in self.cdss:
            fasta += f">{cds.qualifiers['record_id'][0]}|{cds.qualifiers['protein_id'][0]}|{int(cds.location.nofuzzy_start)}..{cds.location.nofuzzy_end}..{cds.strand}\n{format_seq(cds.qualifiers['translation'][0])}\n"

        return fasta

    def has_record(self, record_id : str) -> bool:
        """
        Returns True if record_id is in records.
        :return: bool
        """

        return record_id in [record.id for record in self.records]

    def get_record(self, record_id : str) -> SeqRecord:
        """
        Returns a SeqRecord object
        :param record_id: SeqRecord id
        :return: SeqRecord object
        """

        for record in self.records:
            if record.id == record_id:
                return record

    def get_record_lengt(self, record_id) -> int:
        """
        Returns the length of a record
        :param record_id: SeqRecord id
        :return: int with length
        """

        for record in self.records:
            if record.id == record_id:
                return len(record.seq)

    def get_records_ids_and_size(self) -> List[Tuple[str, int]]:
        """
        Returns a list of record ids and their length
        :return: list of tuples with id and length"""

        return [(record.id, len(record.seq)) for record in self.records]

    def get_records_ids(self) -> List[str]:
        """
        Returns a list of available records.
        :return: list with records
        """
        
        return [record.id for record in self.records]

    def get_cdss_num_per_record(self, record_id : str) -> int:
        """
        Returns the numbero fo CDSs features in record.
        :return: num of CDS features
        """

        return sum([1 if record_id == cds.qualifiers['record_id'][0] else 0 for cds in self.cdss])

    def get_cdss_per_record(self) -> OrderedDict:
        """
        Returns a dict of CDSs per record.
        :return: dict with lists of CDS protein_id per record.id
        """

        cdss_per_record = OrderedDict()
        for cds in self.cdss:
            try:
                cdss_per_record[cds.qualifiers['record_id'][0]].append(cds.qualifiers['protein_id'][0])
            except KeyError:
                cdss_per_record[cds.qualifiers['record_id'][0]] = [cds.qualifiers['protein_id'][0]]

        return cdss_per_record
