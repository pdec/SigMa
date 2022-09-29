"""
Module defining SigMa models
"""

from .read import parse_fasta, parse_genbank
from .features import get_features_of_type, get_feature_coords
from .utils import log_progress, call_process
from .write import format_seq, write_df_to_artemis

from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Tuple
import os
import numpy as np
import pandas as pd
import shutil


class SigMa():
    """
    Main SigMa model
    """

    def __init__(self, args) -> None:
        self.targets: List[Target] = []
        self.queries: List[Query] = []
        self.record_queries: List[RecordQuery] = []
        self.regions: List[Region] = []
        self.args = args
    


class Input():
    """
    General input data class
    """

    def __init__(self, file_path : str, type : str, name : str = None):
        self.file_path = file_path
        self.type = type
        self.name = name if name else os.path.basename(file_path)
    
class Target(Input):
    """
    Target data class
    """

    def __init__(self, file_path : str, type : str, params : Dict, name : str = None):
        super().__init__(file_path, type, name)
        self.params = params
        self.db_path = ''
        self.output_path = ''

    ### default methods ###
    def __str__(self):
        """
        String representation of SigMaRef
        """

        return f"Reference: {self.db} [{self.type}]"

    ### custom methods ###
    def prepare(self, outdir_path : str):
        """
        Prepare target data for searching and set paths
        :param outdir_path: path to output directory
        """

        log_progress(f"Creating reference nucleotide database from {self.name}...", 1)
        # set db_path and results_path
        self.db_path = os.path.join(outdir_path, f"{self.name}.db")
        self.outfile_path = os.path.join(outdir_path,f"{self.name}.{self.type}.tsv")

        if self.type == 'fasta_nt':
            cmd = 'makeblastdb -in {} -dbtype nucl -out {}'.format(self.file_path, self.db_path)
            call_process(cmd)
        elif self.type == 'fasta_aa':
            cmd = 'diamond makedb --in {} --db {} --quiet'.format(self.file_path, self.db_path)
            call_process(cmd)
        elif self.type == 'mmsqes_db':
            self.db_path = self.file_path        

    def search(self, query_path : str):
        """
        Run search
        :param query_path: path to FASTA file
        """
        
        log_progress(f"Searching {self.name} with {query_path}...", 1)

        if self.type == 'fasta_nt':
            cmd = 'blastn -query {} -db {} -out {} -outfmt "6 qaccver saccver pident length qstart qend evalue" -evalue {} -perc_identity {} -num_threads {} -max_target_seqs 10000'.format(query_path, self.db_path, self.output_path, self.params.nt_evalue, self.params.nt_pident, self.params.threads)
            call_process(cmd)

        elif self.type == 'fasta_aa':
            cmd = 'diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid sseqid evalue pident qcovhsp scovhsp --evalue {} --id {} --query-cover {} --subject-cover {} --very-sensitive -c1 --threads {} --max-target-seqs 0 --quiet'.format(query_path, self.db_path, self.output_path, self.params.aa_evalue, self.params.aa_pident, self.params.aa_qscovs, self.params.aa_qscovs, self.params.threads)
            call_process(cmd)

        elif self.type == 'mmseqs_db':
            outdir = os.path.dirname(self.output_path)
            tmp_path = os.path.join(outdir, f"mmseqs_tmp")
            protein_db_path = os.path.join(tmp_path, 'proteindb')
            results_path = os.path.join(tmp_path, 'mmseqs.out')

            # create tmp directory
            if not os.path.exists(tmp_path):
                os.makedirs(tmp_path)

            # create proteindb
            cmd = f"mmseqs createdb {query_path} {protein_db_path}"
            call_process(cmd)

            # make a search
            cmd = f"mmseqs search -s {self.params.mmseqs_sens} --threads {self.params.threads} -e {self.params.mmseqs_evalue} {self.db_path} {protein_db_path} {results_path} {tmp_path}"
            call_process(cmd)

            # create a results tsv file
            cmd = f"mmseqs createtsv {self.db} {protein_db_path} {results_path} {self.output_path}"
            call_process(cmd)

            # delete tmp directory
            shutil.rmtree(tmp_path)

        return

    def read_output(self,queries : List) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """
        Read diamond output and return information about regions with signal
        :param file_path: path to diamond output file
        :param queries: list of SigMaQuery objectss
        :return: two dictionaries for nucleotide-based signal and proteins with list of similar reference proteins
        """

        nt_signal_arrays = {}
        aa_signal_arrays = {}

        if self.type == 'fasta_nt':
            for line in open(outfile_path, 'r'):
                if line.startswith('#'):
                    continue
                else:
                    qaccver, saccver, pident, length, qstart, qend, evalue = line.strip().split('\t')
                    signal = 1 # TODO: add option to use pident or evalue as signal
                    record_id, qlength = qaccver.split("|")
                    qlength = int(qlength)
                    
                    if record_id not in nt_signal_arrays:
                        nt_signal_arrays[record_id] = np.zeros(int(qlength))

                    # record nt signal
                    if int(length) >= self.params.nt_length:
                        nt_signal_arrays[record_id][int(qstart) - 1 : int(qend)] += signal

        elif self.type == 'fasta_aa':
            for line in open(outfile_path, 'r'):
                if line.startswith('#'):
                    continue
                else:
                    qseqid, sseqid, evalue, pident, qcovhsp, scovhsp = line.strip().split('\t')
                    signal = 1 # TODO: add option to use pident, qcovhsp, scovhsp or evalue as signal
                    record_id, protein_id, protein_coords = qseqid.split('|')
                    start, end, strand = map(int, protein_coords.split('..'))

                    # get q
                    for q in queries:
                        if q.has_record(record_id):
                            break
                    else:
                        log_progress(f"Record {record_id} not found in queries", loglevel = "WARNING", msglevel = 1)

                    # record nt equivalent aa signal
                    if record_id not in nt_signal_arrays:
                        nt_signal_arrays[record_id] = np.zeros(q.get_record_length(record_id))
                        aa_signal_arrays[record_id] = np.zeros(q.get_cdss_num_per_record(record_id))

                    # record nt signal
                    nt_signal_arrays[record_id][int(start) - 1 : int(end)] += signal
                    # record aa signal
                    aa_signal_arrays[record_id][q.get_cds_order_num(protein_id)] += signal

        elif self.type == 'mmseqs_db':
            for line in open(outfile_path, 'r'):
                if line.startswith('#'):
                    continue
                else:
                    qseqid, sseqid, alnScore, seqIdentity, eVal, qStart, qEnd, qLen, tStart, tEnd, tLen = line.strip().split('\t')
                    signal = 1 # TODO: add option to use pident, qcovhsp, scovhsp or evalue as signal
                    record_id, protein_id, protein_coords = sseqid.split('|')
                    start, end, strand = map(int, protein_coords.split('..'))

                    # get q
                    for q in queries:
                        if q.has_record(record_id):
                            break
                    else:
                        log_progress(f"Record {record_id} not found in queries", loglevel = "WARNING", msglevel = 1)

                    # record nt equivalent aa signal
                    if record_id not in nt_signal_arrays:
                        nt_signal_arrays[record_id] = np.zeros(q.get_record_len(record_id))
                        aa_signal_arrays[record_id] = np.zeros(q.get_record_cds_num(record_id))

                    # record nt signal
                    nt_signal_arrays[record_id][int(start) - 1 : int(end)] += signal
                    # record aa signal
                    aa_signal_arrays[record_id][q.get_record_cds_order_num(protein_id)] += signal
                
        return nt_signal_arrays, aa_signal_arrays

class Query(Input):
    """
    Query data class
    """

    def __init__(self, file_path : str, type : str, name : str):
        super().__init__(file_path, type, name)
        self.records = []
        self.cdss = []

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

    ### default methods ###
    def __str__(self):
        """
        Query string representation
        """

        return f"Query: {self.file_path} [{self.type}]: {len(self.records)} record{'s' if len(self.records) > 1 else ''} and {len(self.cdss)} CDS{'s' if len(self.cdss) > 1 else ''}"

    def __repr__(self):
        """
        Query representation
        """

        return f"SigMaQuery('{self.file_path}', '{self.type}')"

    ### lookup methods ###
    def has_record(self, record_id : str) -> bool:
        """
        Returns True if record_id is in records.
        :return: bool
        """

        return record_id in [record.id for record in self.records]

    ### get methods ###
    def get_record(self, record_id : str) -> SeqRecord:
        """
        Returns record with record_id.
        :return: SeqRecord
        """

        for record in self.records:
            if record.id == record_id:
                return record
        else:
            raise ValueError(f"Record {record_id} not found in query")

    def get_record_len(self, record_id : str) -> int:
        """
        Returns record length with record_id.
        :return: int
        """

        return len(self.get_record(record_id))

    def get_record_cds_num(self, record_id : str) -> int:
        """
        Returns number of CDSs in record with record_id.
        :return: int
        """

        return len(self.get_record(record_id))

    def get_record_cds_order_num(self, record_id : str, cds_id : str) -> int:
        """
        Returns CDS order number in record with record_id.
        :return: int
        """

        return self.get_record(record_id).cdss.index(cds_id)

    def get_records_ids(self) -> List[str]:
        """
        Returns a list of available records.
        :return: list with records
        """
        
        return [record.record.id for record in self.records]

    def get_records_ids_lens(self) -> List[Tuple[str, int]]:
        """
        Returns a list of record ids and their length
        :return: list of tuples with id and length"""

        return [(record.record.id, record.len()) for record in self.records]

class Record():
    """
    General record class
    """

    def __init__(self, record : SeqRecord):
        self.record = record
        self.name = f"{self.record.id}|{self.len()}"

    ### default methods ###
    def __len__(self):
        return len(self.record.seq)
    
    def len(self):
        return self.__len__()

    ### lookup methods ###
    def has_signal(self, signal_type : str) -> bool:
        """
        Returns True if signal_type is in signal.
        :return: bool
        """

        return signal_type in self.signal.keys()

    ### get methods ###
    def get(self) -> SeqRecord:
        """
        Returns a SeqRecord object
        :return: SeqRecord object
        """

        return self.record
    
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
    
    ### custom methods ###
    def to_fasta_nt(self) -> str:
        """
        Returns a fasta string of the Record
        :return: fasta string
        """

        return f">{self.name}\n{format_seq(self.record.seq)}\n"

    def cdss_to_fasta(self) -> str:
        """
        Return a fasta string of the query CDSs
        :return: fasta string
        """
        fasta = ""
        for cds in self.cdss:
            fasta += f">{cds.qualifiers['record_id'][0]}|{cds.qualifiers['protein_id'][0]}|{int(cds.location.nofuzzy_start)}..{cds.location.nofuzzy_end}..{cds.strand}\n{format_seq(cds.qualifiers['translation'][0])}\n"

        return fasta
    
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


class RecordQuery(Record):
    """
    Single record query class
    """

    def __init__(self, record : SeqRecord):
        super().__init__(record)
        self.regions = {'nt_based': {}, 'aa_based': {}}
        self.signals = {'nt_based': {}, 'aa_based': {}}

    ### get methods ###
    def _get_signal_df(self) -> pd.DataFrame:
        """
        Internal function to make a pandas DataFrame with the signal data.
        :param record_id: SeqRecord id
        :return: pandas DataFrame
        """

        signal_df = pd.DataFrame()
        signal_group = 'nt_based'
        for ref, signal in self.signal[signal_group].items():
            colname = f"{signal_group}_{ref}"
            signal_df[colname] = signal

        return signal_df
        
    ### custom methods ###
    def artemis_plot(self, output_dir : str) -> None:
        """
        Writes Artemis plot files of the query regions.
        :param output_dir: output directory
        """

        for record_id in self.get_records_ids():
            output_path = os.path.join(output_dir, f"{record_id}.artemis.plot")
            write_df_to_artemis(self._get_record_signal_df(record_id), output_path)

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
                                    candidate_region = Region(
                                        record = self.get_record(record_id)[region_start : region_end],
                                        name = f"{record_id}|{ref}|{cand_cnt}",
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
                            candidate_region = Region(
                                record = self.get_record(record_id)[region_start : region_end],
                                name = f"{record_id}|{ref}|{cand_cnt}",
                                start = region_start,
                                end = region_end,
                                reference = ref,
                                signal = signal_array[i_pos : i_gap - max_gap_size + 1],
                                signal_group = signal_group
                                )
                                    
                            self.regions[signal_group].append(candidate_region)
                            cand_cnt += 1
                            log_progress(f"{cand_cnt}: {candidate_region}", msglevel = 3)


class Region(Record):
    """
    Region class
    """

    def __init__(self, record : SeqRecord, name : str, start : int, end : int, signal_source : str, category : str = 'prophage', status : str = 'candidate'):
        self.record = record
        self.name = name
        self.start = start
        self.end = end
        self.signal_source = signal_source
        self.category = category
        self.status = status

    ### default methods ###
    def __repr__(self):
        return f"{self.record.id}[{self.start}..{self.end}] {self.name} ({self.len()}; {self.get_sig_frac()})"

    def __str__(self):
        return self.__repr__()

    ### get methods ###
    def get_sig_frac(self):
        return np.count_nonzero(self.signal) / len(self.signal)

    ### convert methods ###
    def to_fasta_nt(self) -> str:
        """
        Write a fasta file of the query records
        :return: FASTA string
        """
        header = f"{self.record.id}|{self.start}..{self.end}|{self.reference}|{self.signal_group}|{self.category}|{self.status}"
        return f">{header}\n{format_seq(self.record.seq)}"


