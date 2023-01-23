"""
A helper module with small utility commands.
"""

import argparse
import datetime
import logging
import os
import subprocess
import sys
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

def create_logger(log_path : str, log_level : str = logging.INFO, simple : bool = False):
    """
    Create a logger.
    :param log_path: path to log file
    """
    # logger
    if simple:
        logFormatter = logging.Formatter(fmt='%(message)s')
    else:
        logFormatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d [%(levelname)-8.8s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
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

def log_progress(msg : str, stderr : bool = False, stdout : bool = False, quiet : bool = False, loglevel : str = "INFO", msglevel : int = 0, msglevelstr : str = "--"):
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

def call_process(cmd, program : str = ''):
    """
    Call a subprocess.
    :param cmd: a command string to run
    """

    log_progress(f"running: {cmd}", msglevel=2, loglevel="DEBUG")
    if not program: program = cmd.split()[0]
    start_time = datetime.datetime.now()
    with subprocess.Popen(cmd, shell=True, universal_newlines=True, bufsize=1, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
        with process.stdout:
            for line in process.stdout:
                if not line.strip(): continue
                log_progress(f"{program}: {line.strip()}", msglevel = 2, loglevel = "DEBUG")
    end_time = datetime.datetime.now()
    time_taken = end_time - start_time
    log_progress(f"time taken: {time_taken}", msglevel=2, loglevel='DEBUG')

def list_databases(dbs : Dict[str, str]) -> str:
    """
    List available databases.
    :param dbs: a dictionary of databases
    :return: a string of databases
    """

    for db_type, db_paths in dbs.items():
        log_progress(db_type, msglevel = 1, loglevel = "INFO")
        for db_path in db_paths:
            log_progress(db_path, msglevel = 2, loglevel = "INFO")
    