# create a module named SigMaModules and put the following code in it:

from .main import run
from .query import SigMaQuery
from .utils import create_logger, log_progress, call_process, colours
from .read import parse_genbank, parse_fasta
from .features import get_features_of_type, get_cds_header_and_aa
from .write import format_seq, write_fasta, write_df_to_artemis
from .search import make_db_from_fasta, call_blastn, make_diamond_db_from_fasta, call_diamond, call_hmmsearch, call_mmseqs
from .version import __version__

__all__ = [
    'run',
    'SigMaQuery',
    'create_logger', 'log_progress', 'call_process', 'colours',
    'parse_genbank', 'parse_fasta',
    'get_features_of_type', 'get_cds_header_and_aa',
    'format_seq', 'write_fasta', 'write_df_to_artemis', 
    'make_db_from_fasta', 'call_blastn', 'make_diamond_db_from_fasta', 'call_diamond', 'call_hmmsearch', 'call_mmseqs',
    '__version__'
    ]
