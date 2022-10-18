"""
Utility script to prepare reference databases for SigMa using GenBank files.
"""

import argparse
import os
import sys
import shutil
from typing import List
from Bio.Seq import UndefinedSequenceError

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from SigMaModules.read import parse_fasta, parse_genbank
from SigMaModules.write import write_fasta
from SigMaModules.utils import call_process, create_logger, CustomHelpFormatter, log_progress

def get_file_name(path : str) -> str:
    """
    Simply gets file name without extenstions.
    :param path: path to file
    :return: file name
    """
    if path.endswith('.gz'):
        path = path[:-3]

    return os.path.basename(path).rsplit('.', 1)[0]

def mmseqs_cluster(dbname : str, fasta_path : str, tmp_dir : str, sens : float, cov : float, cov_mode : int, id : float, evalue : float, threads : int) -> None:
    """
    Calls MMseqs protein clustering command
    :param dbname: name of the database
    :param fasta_path: path to fasta file
    :param tmp_dir: path to temporary directory
    :param sens: sensitivity
    :param cov: coverage threshold
    :param cov_mode: coverage mode
    :param id: identity threshold
    :param evalue: e-value threshold
    :param threads: number of threads
    :return: None
    """

    # call MMseqs command
    cmd = f"mmseqs easy-cluster -v 2 -s {sens} --cov-mode {cov_mode} -c {cov} --min-seq-id {id} -e {evalue} --threads {threads} {fasta_path} {os.path.join(tmp_dir, dbname)} {os.path.join(tmp_dir, 'protcltml')}"
    call_process(cmd)

    return os.path.join(tmp_dir, f"{dbname}_rep_seq.fasta")

def meshclust_cluster(fasta_path : str, clusters_path :str, nt_id : float, threads : int) -> List[str]:
    """
    Runs MeShClust to cluster nuclotide sequences.
    :param fasta_path: path to fasta file
    :param clusters_path: path to output file
    :param nt_id: sequence identity threshold
    :param threads: number of threads
    """

    # call MeShClust command
    cmd = f"meshclust -d {fasta_path} -o {clusters_path} -t {nt_id} -a n -c {threads}"
    call_process(cmd)

    clusters_ids = []
    C_cnt = 0
    M_cnt = 0
    E_cnt = 0
    with open(clusters_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.endswith('C'):
                clusters_ids.append(line.split()[1][1:])
                C_cnt += 1
            elif line.endswith('M'):
                M_cnt += 1
            elif line.endswith('E'):
                E_cnt += 1

    log_progress(f"MeShClust: {C_cnt} clusters created, {M_cnt} cluster members, {E_cnt} cluster extended members", loglevel="INFO", msglevel=1)

    return clusters_ids

def main():

    parser = argparse.ArgumentParser(
        description='prepdb.py - a script for handling ref datasets preparation for SigMa.', 
        epilog='Example: prepdb.py --gb ./data/pps.gb.gz --dbdir ./dbs --cov 0.8 --cov_mode 0 --evalue 1e-3 --threads 4 --tmp ./tmp\n', 
        formatter_class=CustomHelpFormatter)

    parser.add_argument('--gb', help='GenBank file', nargs = '+', metavar = '<path>')
    parser.add_argument('--fasta_aa', help='Protein FASTA file', nargs = '+', metavar = '<path>')
    parser.add_argument('--fasta_nt', help='Nucleotide FASTA file', nargs = '+', metavar = '<path>')
    parser.add_argument('--dbdir', help='Output directory', metavar = '<path>', required = True)
    parser.add_argument('--sens', help='Sensitivity (0.0-7.5) [%(default)0.1f]', default = 7.5, metavar = ' ', type = float)
    parser.add_argument('--id', help='Identity threshold (0.0-1.0) [%(default)0.2f]', default = 0.9, metavar = ' ', type = float)
    parser.add_argument('--nt_id', help='Nucleotide identity threshold (0.0-1.0) [%(default)0.2f]', default = 0.8, metavar = ' ', type = float)
    parser.add_argument('--cov', help='Minimum coverage (0.0-1.0) [%(default)0.2f]', default = 0.8, metavar = ' ', type = float)
    parser.add_argument('--cov_mode', help='Coverage mode [%(default)i]', default = 0, metavar = ' ', type = int)
    parser.add_argument('--evalue', help='Maximum e-value [%(default)1.0e]', default = 1e-3, metavar = ' ', type = float)
    parser.add_argument('--threads', help='Number of threads [%(default)i]', default = 4, metavar = ' ', type = int)
    parser.add_argument('--tmp', help='Temporary directory', metavar = '<path>', required = True)
    parser.add_argument('--log', help='Log file', metavar = '<path>')

    args = parser.parse_args()
    if len(sys.argv) < 2:
            args.print_usage()
            sys.exit(1)

    # check for required arguments
    if not any([args.gb, args.fasta_aa, args.fasta_nt]):
        print("Any input file was provided. Use one of --gb, --fasta_aa or --fasta_nt.")
        sys.exit(1)

    # create dbdir
    if not os.path.exists(args.dbdir):
        os.makedirs(args.dbdir)
    
    # create tmpdir
    if not os.path.exists(args.tmp):
        os.makedirs(args.tmp)

    # create logger
    if args.log:
        args.logger = create_logger(args.log)
    else:
        args.logger = create_logger(os.path.join(args.dbdir, "prepdb.log"))
    
    log_progress("Preparing reference datasets for SigMa.", loglevel="INFO", msglevel=0)

    # input files for clustering
    aa_fastas = []
    nt_fastas = []
    
    if args.fasta_aa:
        aa_fastas = args.fasta_aa

    if args.fasta_nt:
        nt_fastas = args.fasta_nt
    
    # parse GenBank file
    if args.gb:
        for gb in args.gb:
            dbname = get_file_name(gb)
            aa_file = os.path.join(args.tmp, dbname + '.aa.fasta')
            nt_file = os.path.join(args.tmp, dbname + '.nt.fasta')
            aa_cnt = 0
            nt_cnt = 0
            log_progress(f"Parsing GenBank file: {gb}", loglevel = "INFO")
            log_progress(f"Using name '{dbname}' for the database", loglevel = "INFO", msglevel=1)
            records = parse_genbank(gb)

            # write protein sequences
            log_progress(f"Writing output data", loglevel = "INFO", msglevel=1)
            aa_out = open(aa_file, 'w')
            nt_out = open(nt_file, 'w')
            for record in records:
                try:
                    if len(record.seq) == 0:
                        log_progress(f"Missing sequence for record {record.id}. Skipping.", loglevel="WARNING", msglevel=1)
                        continue
                except UndefinedSequenceError:
                    log_progress(f"Problem with nucleotide sequence for record {record.id}. Skipping.", loglevel="WARNING", msglevel=1)
                    continue
                # write nt
                nt_cnt += 1
                nt_head = f">{dbname}|{record.id}|{nt_cnt}"
                nt_out.write(f"{nt_head}\n{str(record.seq)}\n")
                # write aa
                r_aa_cnt = 0
                for f in record.features:
                    if f.type == 'CDS':
                        aa_cnt += 1
                        r_aa_cnt += 1
                        aa_head = f"{nt_head}|{aa_cnt}"
                        aa_seq = f.qualifiers['translation'][0] if 'translation' in f.qualifiers else f.extract(record.seq).translate(table = 11, to_stop = True)
                        aa_out.write(f"{aa_head}\n{aa_seq}\n")
                log_progress(f"Record {record.id} had {r_aa_cnt} CDS features", loglevel="INFO", msglevel=1)
            # close files                
            aa_out.close()
            nt_out.close()

            if aa_cnt > 0:
                log_progress(f"{aa_cnt} protein sequences written to {aa_file}", loglevel="INFO", msglevel=1)
                aa_fastas.append(aa_file)
            else:
                log_progress(f"There were no CDS features in provided GenBank file.", loglevel="WARNING", msglevel=1)
            if nt_cnt > 0:
                log_progress(f"{nt_cnt} nucleotide sequences written to {nt_file}", loglevel = "INFO", msglevel=1)
                nt_fastas.append(nt_file)

    # cluster protein sequences
    if aa_fastas:
        log_progress(f"Clustering protein sequences", loglevel = "INFO")
        for aa_fasta in aa_fastas:
            dbname = get_file_name(aa_fasta)
            aa_cl_fasta = os.path.join(args.dbdir, f"{dbname}.fasta")
            log_progress(f"Clustering {aa_fasta}", loglevel = "INFO", msglevel=1)
            rep_seq_fasta = mmseqs_cluster(dbname, aa_fasta, args.tmp, args.sens, args.cov, args.cov_mode, args.id, args.evalue, args.threads)
            records = parse_fasta(rep_seq_fasta)
            log_progress(f"MMSeqs: {len(records)} clusters created", loglevel = "INFO", msglevel=1)
            write_fasta(records, aa_cl_fasta)

    # cluster nucleotide sequences
    if nt_fastas:
        log_progress(f"Clustering nucleotide sequences", loglevel = "INFO")
        for nt_fasta in nt_fastas:
            dbname = get_file_name(nt_fasta)
            nt_cl_fasta = os.path.join(args.dbdir, f"{dbname}.fasta")
            log_progress(f"Clustering {nt_fasta}", loglevel = "INFO", msglevel=1)
            clusters_file = os.path.join(args.tmp, f"{dbname}.meshclusters.tsv")
            clusters_ids = meshclust_cluster(nt_fasta, clusters_file, args.nt_id, args.threads)
            write_fasta([record for record in parse_fasta(nt_fasta) if record.id in clusters_ids], nt_cl_fasta)

    log_progress(f"Removing temporary directory {args.tmp}", loglevel = "INFO")
    shutil.rmtree(args.tmp)

    log_progress("Done! You can now run use output files for SigMa!", loglevel = "INFO", msglevel=0)

if __name__ == '__main__':
    main()