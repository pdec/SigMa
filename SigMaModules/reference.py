"""
Module defining SigMaRef class for handling reference datasets and their use
"""

from .utils import log_progress, call_process

import os
import shutil

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