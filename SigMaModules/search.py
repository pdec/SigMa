"""
Module to handling search stage. Includes making reference databases and running searches.
"""
import os
import shutil
from typing import List

from .utils import call_process

### BLASTN
def make_db_from_fasta(file_path : str, db_path : str, db_type : str):
    """
    Create a BLAST database from a FASTA file
    :param file_path: path to FASTA file
    :param db_path: path to BLAST database
    :param db_type: type of BLAST database
    """
    
    cmd = 'makeblastdb -in {} -dbtype {} -out {}'.format(file_path, db_type, db_path)
    call_process(cmd)

    return

def call_blastn(file_path : str, db_path : str, outfile_path : str, evalue : float, pident : int, threads : int):
    """
    Run BLASTN search
    :param file_path: path to FASTA file
    :param db_path: path to BLAST database
    :param outfile_path: path to output file
    :param evalue: e-value threshold
    :param pident: percent identity threshold
    :param threads: number of threads
    """
    
    cmd = 'blastn -query {} -db {} -out {} -outfmt "6 qaccver saccver pident length sstart send evalue" -evalue {} -perc_identity {} -num_threads {} -max_target_seqs 10000'.format(file_path, db_path, outfile_path, evalue, pident, threads)
    call_process(cmd)

    return 

### DIAMOND
def make_diamond_db_from_fasta(file_path : str, db_path : str):
    """
    Create a DIAMOND database from a FASTA file
    :param file_path: path to FASTA file
    :param db_path: path to DIAMOND database
    :param db_type: type of DIAMOND database
    """
    
    cmd = 'diamond makedb --in {} --db {} --quiet'.format(file_path, db_path)
    call_process(cmd)

    return

def call_diamond(file_path : str, db_path : str, outfile_path : str, evalue : float, pident : int, qscovs : int, threads : int):
    """
    Run DIAMOND search
    :param file_path: path to FASTA file
    :param db_path: path to DIAMOND database
    :param outfile_path: path to output file
    :param evalue: e-value threshold
    :param pident: percent identity threshold
    :param qscovs: query and subject coverage threshold
    :param threads: number of threads
    """
    
    cmd = 'diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid sseqid evalue pident qcovhsp scovhsp --evalue {} --id {} --query-cover {} --subject-cover {} --very-sensitive -c1 --threads {} --max-target-seqs 0 --quiet'.format(file_path, db_path, outfile_path, evalue, pident, qscovs, qscovs, threads)
    call_process(cmd)

    return

### HMMER
def call_hmmsearch(file_path : str, db_path : str, out_path : str, evalue : float, threads : int):
    """
    Run HMMSCAN search
    :param file_path: path to FASTA file
    :param db_path: path to HMM database
    :param out_path: path to output file
    :param evalue: e-value threshold
    :param threads: number of threads
    """
    
    cmd = 'hmmsearch --cpu {} --domtblout {} --noali --notextw -E {} --domE {} {} {}'.format(threads, out_path, evalue, evalue, db_path, file_path)
    call_process(cmd)

    return

### MMSEQS
def call_mmseqs(file_path : str, db_path : str, out_path : str, sensitivity : float, evalue : float, pident : float, coverage : float, threads : int):
    """
    Run MMSEQS search
    :param file_path: path to FASTA file with protein sequences
    :param db_path: path to MMSEQS database
    :param out_path: path to output directory
    :param evalue: e-value threshold
    :param threads: number of threads
    """
    tmp_path = os.path.join(f"{out_path}.tmp")
    results_tsv_path = os.path.join(f"{out_path}.tsv")
    protein_db_path = os.path.join(tmp_path, 'proteindb')
    results_path = os.path.join(tmp_path, 'mmseqs.out')

    # create tmp directory
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

    # create proteindb
    cmd = f"mmseqs createdb {file_path} {protein_db_path}"
    call_process(cmd)

    # make a search
    cmd = f"mmseqs search -s {sensitivity} --threads {threads} -e {evalue} {db_path} {protein_db_path} {results_path} {tmp_path}"
    call_process(cmd)

    # create a results tsv file
    cmd = f"mmseqs createtsv {db_path} {protein_db_path} {results_path} {results_tsv_path}"
    call_process(cmd)

    # delete tmp directory
    shutil.rmtree(tmp_path)

    return results_tsv_path
