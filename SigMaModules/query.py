"""
Module defining SigMa query class
"""

from .read import parse_fasta, parse_genbank
from .features import get_features_of_type
from .utils import log_progress
from .write import format_seq

from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from typing import List, Dict, OrderedDict, Tuple
import numpy as np

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
        self.signal = {'nt_based': {}, 'aa_based': {}}
    
        if self.type == 'fasta':
            self.records = parse_fasta(self.file_path)
        elif self.type == 'genbank':
            for record in parse_genbank(self.file_path):
                self.records.append(record)
                rec_cdss = get_features_of_type(record, 'CDS')
                self.cdss.extend(rec_cdss)
                if len(rec_cdss) == 0:
                    log_progress(f"{record.id}: {len(rec_cdss)} CDSs", msglevel = 1, loglevel = "WARNING")
                else:
                    log_progress(f"{record.id}: {len(rec_cdss)} CDSs", msglevel = 1)

    ### overwriting built in methods ###
    def __str__(self):
        """
        SigMaQuery string representation
        """

        return f"Query: {self.file_path} [{self.type}]: {len(self.records)} record{'s' if len(self.records) > 1 else ''} and {len(self.cdss)} CDS{'s' if len(self.cdss) > 1 else ''}"

    def __repr__(self):
        """
        SigMaQuery representation
        """

        return f"SigMaQuery('{self.file_path}', '{self.type}')"

    ### to FASTA ###        
    def records_to_fasta(self) -> str:
        """
        Return a fasta string of the query records
        :return: fasta string
        """
        fasta = ""
        for record in self.records:
            fasta += f">{record.id}|{len(record.seq)}\n{format_seq(record.seq)}\n"

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

    ### lookup attributes methods ###
    def has_record(self, record_id : str) -> bool:
        """
        Returns True if record_id is in records.
        :return: bool
        """

        return record_id in [record.id for record in self.records]

    def has_signal(self, signal_type : str) -> bool:
        """
        Returns True if signal_type is in signal.
        :return: bool
        """

        return signal_type in self.signal.keys()

    ### get methods ###
    def get_record(self, record_id : str) -> SeqRecord:
        """
        Returns a SeqRecord object
        :param record_id: SeqRecord id
        :return: SeqRecord object
        """

        for record in self.records:
            if record.id == record_id:
                return record

    def get_record_length(self, record_id) -> int:
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

    def get_cds_order_num(self, protein_id) -> int:
        """
        Returns the order number of a CDS feature in the query.
        :param protein_id: CDS protein_id
        :return: int with order number
        """

        for i, cds in enumerate(self.cdss):
            if cds.qualifiers['protein_id'][0] == protein_id:
                return i

    def print_signal_summary(self) -> None:
        """
        Returns a summary of the signals in the query.
        """
        
        for signal_group in self.signal.keys():
            if len(self.signal[signal_group]) == 0: continue
            log_progress(f"{signal_group} signal group:", msglevel = 1)
            for record_id, refs in self.signal[signal_group].items():
                log_progress(f"{record_id}:", msglevel = 2)
                for ref, signal in refs.items():
                    log_progress(f"{ref}: {sum([1 if x else 0 for x in signal])}/{signal.size} of total {sum(signal)}", msglevel = 3)

    ### modify attributes methods ###
    def add_signal(self, signal_group : str, signal_name : str, signal_arrays : Dict[str, np.ndarray]) -> None:
        """
        Add signal array to the query.
        :param signal_group: signal group
        :param signal_name: name of the signal
        :param signal_array: NumPy array with signal values
        """
        
        for record_id, signal_array in signal_arrays.items():
            if self.has_record(record_id):
                try:
                    self.signal[signal_group][record_id][signal_name] = signal_array
                except KeyError:
                    try:
                        self.signal[signal_group][record_id] = {signal_name: signal_array}
                    except KeyError:
                        self.signal[signal_group] = {record_id: {signal_name: signal_array}}

