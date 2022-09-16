"""
A helper module with small utility commands.
"""

import datetime
import logging
import subprocess
import sys
from typing import Dict, List, Tuple

def create_logger(log_path):
    """
    Create a logger.
    :param log_path: path to log file
    """
    # logger
    logFormatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d [%(levelname)-8.8s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger('SigMa')
    logger.setLevel(logging.DEBUG)
    # file logger
    fileHandler = logging.FileHandler(log_path)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    # console logger
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    return logger

def log_progress(msg : str, stderr : bool = False, stdout : bool = False, quiet : bool = False, loglevel = "INFO"):
    """
    Log program execution progress to log file and stdout.
    :param msg: the message to log
    :param stderr: log to stderr
    :param stdout: log to stdout
    :param quiet: don't log anything
    :param loglevel: the log level to use
    """

    
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
        ('0 0 0','black'),
        ('230 25 75','red'),
        ('60 180 75','green'),
        ('0 130 200','blue'),
        ('255 225 25','yellow'),
        ('145 30 180','purple'),
        ('70 240 240','cyan'),
        ('245 130 48','orange'),
        ('240 50 230','pink'),
        ('210 245 60','lime'),
        ('250 190 212','peach'),
        ('0 128 128','teal'),
        ('220 190 255','lavender'),
        ('170 110 40','brown'),
        ('255 250 200','beige'),
        ('128 0 0','maroon'),
        ('170 255 195','mint'),
        ('128 128 0','olive'),
        ('255 215 180','coral'),
        ('0 0 128','navy'),
        ('128 128 128','gray'),
        ('255 255 255','white'),
    ]
    return colours

def call_process(cmd):
    """
    Call a subprocess.
    :param cmd: a command string to run
    """

    log_progress(f"Running: {cmd}")
    start_time = datetime.datetime.now()
    subprocess.call(cmd, shell=True, universal_newlines=True, bufsize=1, stdout=subprocess.PIPE)
    end_time = datetime.datetime.now()
    time_taken = end_time - start_time
    log_progress(f"Time taken: {time_taken}")
    return

def list_databases(dbs : Dict[str, str]) -> str:
    """
    List available databases.
    :param dbs: a dictionary of databases
    :return: a string of databases
    """

    for db_type, db_paths in dbs.items():
        log_progress(f"-- {db_type}")
        for db_path in db_paths:
            log_progress(f"   -- {db_path}")
    