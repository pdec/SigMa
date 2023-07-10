"""
A helper module with small utility commands.
"""

import argparse
import datetime
import logging
import os
import subprocess
import sys
from collections import OrderedDict
from typing import Dict, List, Tuple


class CustomHelpFormatter(argparse.HelpFormatter):
    def __init__(self, prog):
        super().__init__(prog, max_help_position=60, width=120)

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


def create_logger(log_path: str, log_level: str = "INFO", simple: bool = False):
    """
    Create a logger.
    :param log_path: path to log file
    :param log_level: the log level to use
    :param simple: use simple log format
    """
    # logger
    if simple:
        logFormatter = logging.Formatter(fmt='%(message)s')
    else:
        logFormatter = logging.Formatter(
            fmt='%(asctime)s [%(levelname)-1.1s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger('SigMa')
    logger.setLevel(log_level)
    # file logger
    fileHandler = logging.FileHandler(log_path)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    # console logger
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    return logger


def log_progress(msg: str, stderr: bool = False, stdout: bool = False, quiet: bool = False, loglevel: str = "INFO", msglevel: int = 0, msglevelstr: str = "--"):
    """
    Log program execution progress to log file and stdout.
    :param msg: the message to log
    :param stderr: log to stderr
    :param stdout: log to stdout
    :param quiet: don't log anything
    :param loglevel: the log level to use
    :param msglevel: additional indentation level
    """

    if msglevel > 0:
        indent = msglevelstr * msglevel
        msg = f"{indent} {msg}"

    if loglevel == "INFO":
        logging.getLogger('SigMa').info(msg)
    elif loglevel == "DEBUG":
        logging.getLogger('SigMa').debug(msg)
    elif loglevel == "WARNING":
        logging.getLogger('SigMa').warning(msg)
    elif loglevel == "ERROR":
        logging.getLogger('SigMa').error(msg)
    elif loglevel == "CRITICAL":
        logging.getLogger('SigMa').critical(msg)

    if not quiet:
        if stderr:
            print(msg, file=sys.stderr)
        elif stdout:
            print(msg, file=sys.stdout)


def colours() -> List[Tuple[str, str]]:
    """
    Get list of colours for Artemis graph file
    :return:
    """

    colours = [
        ('0 0 0', 'black'),
        ('230 25 75', 'red'),
        ('60 180 75', 'green'),
        ('0 130 200', 'blue'),
        ('255 225 25', 'yellow'),
        ('145 30 180', 'purple'),
        ('70 240 240', 'cyan'),
        ('245 130 48', 'orange'),
        ('240 50 230', 'pink'),
        ('210 245 60', 'lime'),
        ('250 190 212', 'peach'),
        ('0 128 128', 'teal'),
        ('220 190 255', 'lavender'),
        ('170 110 40', 'brown'),
        ('255 250 200', 'beige'),
        ('128 0 0', 'maroon'),
        ('170 255 195', 'mint'),
        ('128 128 0', 'olive'),
        ('255 215 180', 'coral'),
        ('0 0 128', 'navy'),
        ('128 128 128', 'gray'),
        ('255 255 255', 'white'),
    ]
    return colours


def call_process(cmd, program: str = ''):
    """
    Call a subprocess.
    :param cmd: a command string to run
    """

    log_progress(f"running: {cmd}", msglevel=2, loglevel="DEBUG")
    if not program:
        program = cmd.split()[0]
    start_time = datetime.datetime.now()
    try:
        with subprocess.Popen(
            cmd,
            shell=True,
            universal_newlines=True,
            bufsize=1,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        ) as process:

            with process.stdout:
                for line in process.stdout:
                    if not line.strip():
                        continue
                    log_progress(f"{program}: {line.strip()}",
                                 msglevel=2, loglevel="DEBUG")

    except (OSError, subprocess.CalledProcessError) as exception:
        with process.stderr:
            for line in process.stderr:
                if not line.strip():
                    continue
                log_progress(f"{program}: {line.strip()}",
                             msglevel=2, loglevel="ERROR")
        log_progress(f"{program}: Exception occured: {exception}",
                     msglevel=2, loglevel="ERROR")
        exit(1)
    else:
        # no exception was raised
        pass

    end_time = datetime.datetime.now()
    time_taken = end_time - start_time
    log_progress(f"{program}: Finished. Time taken: {time_taken}",
                 msglevel=2, loglevel='DEBUG')


def list_databases(dbs: Dict[str, str]) -> None:
    """
    List available databases.
    :param dbs: a dictionary of databases
    :return: a string of databases
    """

    for db_type, db_paths in dbs.items():
        log_progress(db_type, msglevel=1, loglevel="INFO")
        for db_path in db_paths:
            log_progress(db_path, msglevel=2, loglevel="INFO")


def make_batches(items: List, batches: int = 1, batch_size: int = 0) -> Dict[str, List[str]]:
    """
    Make batches of items based on either the number of batches or batch size.
    :param items: a list of items
    :param batches: the number of batches to create
    :param batche_size: the number of items per batch
    :return: a list of batches
    """

    batch_dict = OrderedDict()
    if batch_size > 0:
        batches = int(len(items) / batch_size)
        if batch_size * batches < len(items):
            batches += 1
    elif batches > 1:
        batch_size = int(len(items) / batches)
        if batch_size * batches < len(items):
            batch_size += 1

    log_progress(f"Batch size: {batch_size}", msglevel=1, loglevel="INFO")
    log_progress(f"Number of batches: {batches}", msglevel=1, loglevel="INFO")

    for bi, si in enumerate(range(0, len(items), batch_size), 1):
        batch_name = f"batch_{bi:04d}"
        batch_dict[batch_name] = items[si:si + batch_size]

    log_progress(
        f"Batch {list(batch_dict.keys())[0]} contains {len(list(batch_dict.values())[0])} files", msglevel=1, loglevel="INFO")
    log_progress(
        f"Batch {list(batch_dict.keys())[-1]} contains {len(list(batch_dict.values())[-1])} files", msglevel=1, loglevel="INFO")

    return batch_dict


def read_batches_done(main_outdir: str) -> List[str]:
    """
    Get a list of batches.
    :param main_outdir: the main output directory
    :return: a list of batches
    """

    status_tsv = os.path.join(main_outdir, 'batch_status.tsv')
    batches = []
    if os.path.exists(status_tsv):
        with open(status_tsv) as inf:
            inf.readline()
            for line in inf:
                if not line.strip():
                    continue
                batch = line.strip().split()[1]
                batches.append(batch)
    log_progress(f"Found {len(batches)} completed batches",
                 msglevel=1, loglevel="INFO")

    return batches


def post_batch(main_outdir: str, outdir: str, batch: str, bi: int) -> Tuple[int, int]:
    """
    If a batch is complete, combine output batch files of verified and candidate regions into combined files.
    :param main_outdir: the main output directory
    :param outdir: the batch output directory
    :param batch: the batch name
    :param bi: batch number
    """

    # combined files
    status_tsv = os.path.join(main_outdir, 'batch_status.tsv')
    candidates_fasta = os.path.join(main_outdir, 'all_candidates.fasta')
    verified_gb = os.path.join(main_outdir, 'all_verified.gb')
    candidates_tsv = os.path.join(main_outdir, 'all_candidates.tsv')
    verified_tsv = os.path.join(main_outdir, 'all_verified.tsv')
    verified_regions = 0
    candidate_regions = 0

    # initialize combined files
    if bi == 1:
        with open(status_tsv, 'w') as outf:
            outf.write('batch_num\tbatch\tstatus\tverified\tcandidate\n')
        with open(candidates_fasta, 'w') as outf:
            pass
        with open(verified_gb, 'w') as outf:
            pass
        with open(candidates_tsv, 'w') as outf:
            pass
        with open(verified_tsv, 'w') as outf:
            pass

    # batch outputs
    verified_regions_gb = os.path.join(outdir, 'regions', 'verified.gb')
    verified_regions_tsv = os.path.join(outdir, 'summary.verified.tsv')
    candidate_regions_fasta = os.path.join(
        outdir, 'regions', 'candidate.fasta')
    candidate_regions_tsv = os.path.join(outdir, 'summary.candidates.tsv')

    # combine output files if they exist
    if os.path.exists(verified_regions_gb):
        with open(verified_regions_gb, 'r') as inf, open(verified_gb, 'a') as outf:
            outf.write(inf.read())
    if os.path.exists(verified_regions_tsv):
        with open(verified_regions_tsv, 'r') as inf, open(verified_tsv, 'a') as outf:
            for line in inf:
                verified_regions += 1
                outf.write(line)

    if os.path.exists(candidate_regions_fasta):
        with open(candidate_regions_fasta, 'r') as inf, open(candidates_fasta, 'a') as outf:
            outf.write(inf.read())
    if os.path.exists(candidate_regions_tsv):
        with open(candidate_regions_tsv, 'r') as inf, open(candidates_tsv, 'a') as outf:
            for line in inf:
                candidate_regions += 1
                outf.write(line)

    # write batch status
    with open(os.path.join(main_outdir, 'batch_status.tsv'), 'a') as outf:
        outf.write(
            f"{bi}\t{batch}\tdone\t{verified_regions}\t{candidate_regions}\n")

    return verified_regions, candidate_regions
