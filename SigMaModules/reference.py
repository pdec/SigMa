"""
Module defining SigMaRef class for handling reference datasets and their use
"""

from .utils import log_progress, call_process

import os
import shutil
import numpy as np

from typing import Dict, List, Tuple

class SigMaRef():
    """This is class for storing and working on reference datasets"""

    def __init__(self, file_path : str, input_type : str, ref_dir) -> None:
        """
        SigMaRef object constructor
        """
        self.file_path = file_path
        self.type = input_type
        self.name = os.path.basename(file_path)
        self.db = os.path.join(ref_dir, self.name)

    def __str__(self):
        """
        String representation of SigMaRef
        """

        return f"Reference: {self.db} [{self.type}]"

    def get_output_path(self, out_path : str) -> str:
        """
        Get the path for the output file
        :param out_path: output directory path
        :return: path to the output file
        """

        return os.path.join(out_path,f"{self.name}.{self.type}.tsv")
class SigMaRefNT(SigMaRef):
    """This is class for storing and working on nucleotide reference datasets"""

    def __init__(self, file_path : str, input_type : str, ref_dir) -> None:
        """
        SigMaRefNT object constructor
        """
        super().__init__(file_path, input_type, ref_dir)

        log_progress(f"Creating reference nucleotide database from {self.name}...", 1)
        cmd = 'makeblastdb -in {} -dbtype nucl -out {}'.format(self.file_path, self.db)
        call_process(cmd)

    def search(self, query_path : str, outfile_path : str, evalue : float, pident : int, threads : int):
        """
        Run BLASTN search
        :param query_path: path to FASTA file
        :param db_path: path to BLAST database
        :param outfile_path: path to output file
        :param evalue: e-value threshold
        :param pident: percent identity threshold
        :param threads: number of threads
        """
        
        cmd = 'blastn -query {} -db {} -out {} -outfmt "6 qaccver saccver pident length sstart send evalue" -evalue {} -perc_identity {} -num_threads {} -max_target_seqs 10000'.format(query_path, self.db, outfile_path, evalue, pident, threads)
        call_process(cmd)

    
    def read_output(self, outfile_path : str, min_nt_length : int) -> Dict[str, np.ndarray]:
        """
        Read blastn output and return information about regions with signal
        :param file_path: path to blastn output file
        :param min_nt_length: minimum length of the alignment
        :return: dictionary of nucleotide records ids and numpy arrays with signal covered by the alignment
        """
        
        nt_signal_arrays = {}
        for line in open(outfile_path, 'r'):
            if line.startswith('#'):
                continue
            else:
                qaccver, saccver, pident, length, sstart, send, evalue = line.strip().split('\t')
                signal = 1 # TODO: add option to use pident or evalue as signal
                record_id, qlength = qaccver.split("|")
                qlength = int(qlength)
                
                if record_id not in nt_signal_arrays:
                    nt_signal_arrays[record_id] = np.zeros(int(qlength))

                # record nt signal
                if int(length) >= min_nt_length:
                    nt_signal_arrays[record_id][int(sstart) - 1 : int(send)] += signal

                        

        return nt_signal_arrays


class SigMaRefAA(SigMaRef):
    """This is class for storing and working on amino acid reference datasets"""

    def __init__(self, file_path : str, input_type : str, ref_dir) -> None:
        """
        SigMaRefAA object constructor
        """
        super().__init__(file_path, input_type, ref_dir)

        log_progress(f"Creating reference amino acid database from {self.name}...", 1)
        cmd = 'diamond makedb --in {} --db {} --quiet'.format(self.file_path, self.db)
        call_process(cmd)

    def search(self, query_path : str, outfile_path : str, evalue : float, pident : int, qscovs : int, threads : int):
        """
        Run DIAMOND search
        :param query_path: path to FASTA file
        :param outfile_path: path to output file
        :param evalue: e-value threshold
        :param pident: percent identity threshold
        :param qscovs: query and subject coverage threshold
        :param threads: number of threads
        """
        
        cmd = 'diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid sseqid evalue pident qcovhsp scovhsp --evalue {} --id {} --query-cover {} --subject-cover {} --very-sensitive -c1 --threads {} --max-target-seqs 0 --quiet'.format(query_path, self.db, outfile_path, evalue, pident, qscovs, qscovs, threads)
        call_process(cmd)

        return

    def read_output(self, outfile_path : str, queries : List) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """
        Read diamond output and return information about regions with signal
        :param file_path: path to diamond output file
        :param queries: list of SigMaQuery objectss
        :return: two dictionaries for nucleotide-based signal and proteins with list of similar reference proteins
        """
        
        nt_signal_arrays = {}
        aa_signal_arrays = {}
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
                

        return nt_signal_arrays, aa_signal_arrays


class SigMaRefHMM(SigMaRef):
    """This is class for storing and working on HMM reference datasets"""

    def __init__(self, file_path, input_type, ref_dir):
        """
        SigMaRefHMM object constructor
        """
        super().__init__(file_path, input_type, ref_dir)
        self.db = self.file_path

        log_progress(f"Using provided HMM database from {self.name}...", 1)

    ### HMMER
    def search(self, query_path : str, outfile_path : str, evalue : float, threads : int):
        """
        Run HMMSCAN search
        :param file_path: path to FASTA file
        :param db_path: path to HMM database
        :param outfile_path: path to output file
        :param evalue: e-value threshold
        :param threads: number of threads
        """
        
        cmd = 'hmmsearch --cpu {} --domtblout {} --noali --notextw -E {} --domE {} {} {}'.format(threads, outfile_path, evalue, evalue, self.db, query_path)
        call_process(cmd)

        return

class SigMaRefMMSEQS(SigMaRef):
    """This is class for storing and working on MMSEQS reference datasets"""

    def __init__(self, file_path, input_type, ref_dir):
        """
        SigMaRefMMSEQS object constructor
        """
        super().__init__(file_path, input_type, ref_dir)
        self.db = self.file_path

        log_progress(f"Using provided MMSEQS HMM database from {self.name}...", 1)

    ### MMSEQS
    def search(self, query_path : str, out_path : str, sensitivity : float, evalue : float, pident : float, coverage : float, threads : int):
        """
        Run MMSEQS search
        :param query_path: path to FASTA file with protein sequences
        :param db_path: path to MMSEQS database
        :param out_path: path to output directory
        :param sensitivity: sensitivity of the search
        :param evalue: e-value threshold
        :param pident: percent identity threshold
        :param coverage: query and subject coverage threshold
        :param threads: number of threads
        """
        outdir = os.path.dirname(out_path)
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
        cmd = f"mmseqs search -s {sensitivity} --threads {threads} -e {evalue} {self.db} {protein_db_path} {results_path} {tmp_path}"
        call_process(cmd)

        # create a results tsv file
        cmd = f"mmseqs createtsv {self.db} {protein_db_path} {results_path} {out_path}"
        call_process(cmd)

        # delete tmp directory
        shutil.rmtree(tmp_path)

        return

    def read_output(self, outfile_path : str, queries : List) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """
        Read MMSEQS output and return information about regions with signal
        :param file_path: path to MMSEQS output file
        :param queries: list of SigMaQuery objectss
        :return: two dictionaries for nucleotide-based signal and proteins with list of similar reference proteins
        """
        
        nt_signal_arrays = {}
        aa_signal_arrays = {}
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
                    nt_signal_arrays[record_id] = np.zeros(q.get_record_length(record_id))
                    aa_signal_arrays[record_id] = np.zeros(q.get_cdss_num_per_record(record_id))

                # record nt signal
                nt_signal_arrays[record_id][int(start) - 1 : int(end)] += signal
                # record aa signal
                aa_signal_arrays[record_id][q.get_cds_order_num(protein_id)] += signal
                

        return nt_signal_arrays, aa_signal_arrays

