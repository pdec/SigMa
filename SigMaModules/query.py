"""
Module defining SigMa query class
"""

from .read import parse_fasta, parse_genbank
from .features import get_features_of_type, get_feature_coords
from .utils import log_progress
from .write import format_seq, write_df_to_artemis

from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from typing import List, Dict, OrderedDict, Tuple
import os
import numpy as np
import pandas as pd

class SigMaRegion:
    """Class for handling predicted signal regions"""

    def __init__(self, record : SeqRecord, start : int, end : int, reference : str, signal : np.ndarray, signal_group : str):
        self.record = record
        self.start = start
        self.end = end
        self.reference = reference
        self.signal = signal
        self.signal_group = signal_group
        self.category = 'prophage'
        self.status = 'candidate'

    def __repr__(self):
        return f"{self.record.id} [{self.start}..{self.end}] ({self.get_len()}; {self.get_sig_frac()})"

    def __str__(self):
        return self.__repr__()

    def __len__(self):
        return len(self.record.seq)
    
    def get_len(self):
        return self.__len__()

    def get_sig_frac(self):
        return np.count_nonzero(self.signal) / len(self.signal)

    def record_to_fasta(self) -> str:
        """
        Write a fasta file of the query records
        :return: FASTA string
        """
        header = f"{self.record.id}|{self.start}..{self.end}|{self.reference}|{self.signal_group}|{self.category}|{self.status}"
        return f">{header}\n{format_seq(self.record.seq)}"


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
        self.regions = {'nt_based': [], 'aa_based': []}
    
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

    def get_regions_nt_fasta(self) -> List[str]:
        """
        Returns a list of FASTA string of the query regions
        :return: list of FASTA string
        """

        fastas = []
        for signal_group, regions in self.regions.items():
            for region in regions:
                fastas.extend(region.record_to_fasta())

        return fastas
    
    ### print or write methods ###
    def print_signal_summary(self) -> None:
        """
        Returns a summary of the signals in the query.
        """
        
        for signal_group in self.signal.keys():
            if len(self.signal[signal_group]) == 0: continue
            log_progress(f"{signal_group} signal group for {self.file_path}:", msglevel = 1)
            for record_id, refs in self.signal[signal_group].items():
                log_progress(f"{record_id}:", msglevel = 2)
                for ref, signal in refs.items():
                    log_progress(f"{ref}: {sum([1 if x else 0 for x in signal])}/{signal.size} of total {sum(signal)}", msglevel = 3)

    def _get_record_signal_df(self, record_id : str) -> pd.DataFrame:
        """
        Internal function to make a pandas DataFrame with the signal data.
        :param record_id: SeqRecord id
        :return: pandas DataFrame
        """

        signal_df = pd.DataFrame()
        signal_group = 'nt_based'
        for ref, signal in self.signal[signal_group][record_id].items():
            colname = f"{signal_group}_{ref}"
            signal_df[colname] = signal

        return signal_df
        
    def regions_to_artemis_plot(self, output_dir : str) -> None:
        """
        Writes Artemis plot files of the query regions.
        :param output_dir: output directory
        """

        for record_id in self.get_records_ids():
            output_path = os.path.join(output_dir, f"{record_id}.artemis.plot")
            write_df_to_artemis(self._get_record_signal_df(record_id), output_path)

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

    def evaluate(self, max_nt_gap : int, min_nt_signals : int, max_aa_gap : int, min_aa_signals : int, min_sig_frac : float) -> None:
        """
        Evaluate signal for each approach and database #TODO allow for combining signal from different approaches
        :param max_nt_gap: max nt gap between signal
        :param min_nt_sig: min nt signal
        :param max_aa_gap: max aa gap between signal
        :param min_aa_sig: min aa signal
        :param min_sig_frac: min signal fraction
        """
        for signal_group in self.signal.keys():
            # setup thresholds for signal considerations and regions merging
            if signal_group == 'nt_based':
                max_gap_size = max_nt_gap
                min_sig_signals = min_nt_signals
            elif signal_group == 'aa_based':
                max_gap_size = max_aa_gap
                min_sig_signals = min_aa_signals
            min_signal_frac = min_sig_frac
        
            log_progress(f"{signal_group} evaluation of {self.file_path}...", msglevel = 1)
            for record_id, refs in self.signal[signal_group].items():
                for ref, signal_array in refs.items():
                    cand_cnt = 0
                    log_progress(f"{record_id} {np.count_nonzero(signal_array)} {'positions' if signal_group == 'nt_based' else 'proteins'} based on {ref}", msglevel = 2)
                    i_pos = -1
                    i_gap = -1
                    pos_len = 0
                    gap_size = 0
                    region_length = 0
                    region_values = []

                    # start searching
                    for i, v in enumerate(signal_array):
                        # count positive signal
                        if v > 0: # >= min_signal_value:
                            if i_pos == -1: # if that's the first value after negative
                                i_pos = i
                                pos_len = 0
                            pos_len += 1
                            gap_size = 0
                        else:
                            i_gap = i
                            gap_size += 1
                            # if the gap is too big check if region can be considered
                            if (gap_size > max_gap_size):
                                # if other thresholds are met, consider region
                                if (pos_len >= min_sig_signals) and (pos_len / (len(region_values) - gap_size + 1) >= min_signal_frac):
                                    if signal_group == 'nt_based':
                                        region_start = i_pos
                                        region_end = i_gap - max_gap_size + 1
                                    elif signal_group == 'aa_based':
                                        region_start = get_feature_coords(self.cdss[i_pos])[0]
                                        region_end = get_feature_coords(self.cdss[i_gap - max_gap_size + 1])[1]
                                    candidate_region = SigMaRegion(
                                        record = self.get_record(record_id)[region_start : region_end],
                                        start = region_start,
                                        end = region_end,
                                        reference = ref,
                                        signal = signal_array[i_pos : i_gap - max_gap_size + 1],
                                        signal_group = signal_group
                                        )
                                    
                                    self.regions[signal_group].append(candidate_region)
                                    cand_cnt += 1
                                    log_progress(f"{cand_cnt}: {candidate_region}", msglevel = 3)
                                else:
                                    # thresholds unmet
                                    pass
                                # reset, as maximum gap size reached
                                i_pos = -1
                                pos_len = 0
                                region_values = []

                        region_values.append(v)

                    # when the gap was not big enough at the end, consider the last region
                    if (gap_size < max_gap_size) and (pos_len >= min_sig_signals) and (pos_len / (len(region_values) - gap_size + 1) >= min_signal_frac):
                            if signal_group == 'nt_based':
                                region_start = i_pos
                                region_end = i_gap - max_gap_size + 1
                            elif signal_group == 'aa_based':
                                region_start = get_feature_coords(self.cdss[i_pos])[0]
                                region_end = get_feature_coords(self.cdss[i_gap - max_gap_size + 1])[1]
                            candidate_region = SigMaRegion(
                                record = self.get_record(record_id)[region_start : region_end],
                                start = region_start,
                                end = region_end,
                                reference = ref,
                                signal = signal_array[i_pos : i_gap - max_gap_size + 1],
                                signal_group = signal_group
                                )
                                    
                            self.regions[signal_group].append(candidate_region)
                            cand_cnt += 1
                            log_progress(f"{cand_cnt}: {candidate_region}", msglevel = 3)

