"""
Module defining SigMa models
"""

from .read import parse_fasta, parse_genbank
from .utils import log_progress, call_process
from .write import format_seq, write_df_to_artemis

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from collections import OrderedDict
from typing import List, Dict, Tuple, Union
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
        self.hq_regions: List[Region] = []
        self.args = args
        self.dirs = {}

        # setup output directories
        directories = ['reference', 'query', 'search', 'artemis_plots', 'regions', 'checkv', 'phispy']
        for dir_name in directories:
            if dir_name == 'artemis_plots' and not args.artemis_plots:
                continue
            elif dir_name == 'phispy' and args.ext_tools and ('phispy' not in args.ext_tools):
                continue
            dir_path = os.path.join(self.args.outdir, dir_name)
            self.dirs[dir_name] = dir_path
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)


    ### custom methods ###
    def prepare_targets(self):
        log_progress("Preparing references", msglevel = 0, loglevel = "INFO")
    
        for ref_dataset_path, ref_type in zip(self.args.reference, self.args.reference_type):
            # create target object
            self.targets.append(
                Target(
                    file_path = ref_dataset_path, 
                    type = ref_type, 
                    params = self.args
                    )
                )
            # prepare database if needed
            self.targets[-1].prepare(self.dirs['reference'])

        log_progress("Reference prepared", msglevel=0, loglevel="INFO")
        for target in self.targets:
            log_progress(f"{target}", msglevel = 1, loglevel = "INFO")

    def handle_queries(self):
        """
        Runs analyses for each input query from reading to signal evaluation
        :return: None
        """
        log_progress("Handling queries", msglevel = 0, loglevel = "INFO")

        # prepare counters
        q_num = len(self.args.query)
        r_num = 0 # regions
        rm_num = 0 # regions merged
        qi = 0
        for query_type, query_dataset_path in zip(self.args.query_type, self.args.query):
            # PARSE QUERY
            qi += 1
            done_q = f"{qi}/{q_num}"
            log_progress(f"[{done_q}] Query: {query_dataset_path}", msglevel = 0, loglevel = "INFO")
            self.queries.append(
                # create Query object
                Query(
                    file_path = query_dataset_path,
                    type = query_type
                    )
                )

            # PARSE QUERY RECORDS
            rqi = 0
            rq_num = len(self.queries[-1].get_record_queries())
            
            # prepare query files
            query_nt_path = os.path.join(self.dirs['query'], f'q{qi}_nt.fasta')
            query_aa_path = os.path.join(self.dirs['query'], f'q{qi}_aa.fasta')
            output_prefix = os.path.join(self.dirs['search'], f'q{qi}')

            # combine sequnces from a single Query
            query_nts = []
            query_aas = []

            # combine sequnces from a single Query
            for record_query in self.queries[-1].get_record_queries():
                rqi += 1
                # just link RecordQueries to SigMa
                self.record_queries.append(record_query)
                # counter for progress
                done_rq = f"{rqi}/{rq_num}"
                if not record_query.has_nt() and not record_query.has_aa():
                    log_progress(f"[{done_q} - {done_rq}] {record_query} skipped - no sequences", msglevel = 1, loglevel = "DEBUG")
                    continue
                elif record_query.len() < self.args.min_nt_sig:
                    log_progress(f"[{done_q} - {done_rq}] {record_query} skipped - too short)", msglevel = 1, loglevel = "DEBUG")
                    continue
                else:
                    log_progress(f"[{done_q} - {done_rq}] {record_query}", msglevel = 1, loglevel = "DEBUG")

                # write query files
                if record_query.has_nt():
                    query_nts.extend(record_query.get_fasta_nt())
                else:
                    log_progress(f"[{done_q} - {done_rq}] {record_query} has no nucleotide sequences", msglevel = 2, loglevel = "WARNING")
                if record_query.has_aa():
                    query_aas.extend(record_query.get_fasta_aa())
                else:
                    log_progress(f"[{done_q} - {done_rq}] {record_query} has no proteins", msglevel = 2, loglevel = "WARNING")
            
            # write query files
            if len(query_nts) > 0:
                self.write_fastas(query_nts, query_nt_path)
            else:
                log_progress(f"0 nucleotide sequences in {self.queries[-1]}. Not writing.", msglevel = 2, loglevel = "WARNING")
            if len(query_aas) > 0:
                self.write_fastas(query_aas, query_aa_path)
            else:
                log_progress(f"0 proteins in {self.queries[-1]}. Not writing.", msglevel = 2, loglevel = "WARNING")
            
            # Search combined sequenced against Targets
            for target in sorted(self.targets, key = lambda x: {'fasta_nt': 0, 'fasta_aa': 1, 'hmm': 2, 'mmseqs_db': 3}[x.type]):
                # search query against target
                output_path = ''
                if target.type == 'fasta_nt':
                    if len(query_nts) > 0:
                        output_path = target.search(query_nt_path, output_prefix)
                else:
                    if len(query_aas) > 0:
                        output_path = target.search(query_aa_path, output_prefix)

                # add signal only if output file exists and is not empty
                if not output_path or (not os.path.exists(output_path) and os.path.getsize(output_path) > 0):
                    continue

                # parse search results
                nt_signal_array, aa_signal_array = target.read_output(self.queries[-1], output_path) # returns a dict of arrays per record_id
                for record_query in self.queries[-1].get_record_queries():
                    rqid = record_query.get_id()
                    signal_group = 'nt_based'
                    if rqid in nt_signal_array and nt_signal_array[rqid] is not None:
                        record_query.add_signal(signal_group, target.name, nt_signal_array[rqid])
                    if rqid in aa_signal_array and aa_signal_array[rqid] is not None and target.type != 'fasta_nt':
                        signal_group = 'aa_based'
                        record_query.add_signal(signal_group, target.name, aa_signal_array[rqid])
            
            # search against specific tools
            if self.args.ext_tools and 'phispy' in self.args.ext_tools:
                nt_signal_array = self.run_phispy(self.queries[-1], qi)
                for record_query in self.queries[-1].get_record_queries():
                    rqid = record_query.get_id()
                    if rqid in nt_signal_array and nt_signal_array[rqid] is not None:
                        record_query.add_signal('nt_based', 'phispy', nt_signal_array[rqid])

            # combine all signals together as a separate category
            if self.args.combine: 
                for record_query in self.queries[-1].get_record_queries():
                    record_query.combine_signals()

            # print signal summary
            for record_query in self.queries[-1].get_record_queries():
                record_query.print_signal_summary()

                # EVALUATE SIGNALS
                log_progress(f"[{done_q} - {done_rq}] Evaluating {record_query}", msglevel = 1, loglevel = "DEBUG")
                record_query.evaluate(self.args.max_nt_gap, self.args.min_nt_sig, self.args.max_aa_gap, self.args.min_aa_sig, self.args.min_sig_frac, self.args.sig_groups)
                record_query.merge_regions()
                for regions in record_query.get_regions().values():
                    self.regions.extend(regions)


            r_num = len(self.regions)
            rm_num = len(self.filter_regions(sig_groups = 'merged'))
            log_progress(f"{rm_num} merged non-overlapping regions were identified based on {r_num - rm_num} regions", msglevel = 1, loglevel = "INFO")

    # not used - as of now
    # def prepare_queries(self):
    #     log_progress("Preparing queries", msglevel = 0, loglevel = "INFO")
    #     for query_type, query_dataset_path in zip(self.args.query_type, self.args.query):
    #         self.queries.append(
    #             # create Query object
    #             Query(
    #                 file_path = query_dataset_path,
    #                 type = query_type
    #                 )
    #             )
    #         # link RecordQueries
    #         self.record_queries.extend(self.queries[-1].get_record_queries())

    # not used - as of now
    # def search_queries(self) -> None:
    #     """
    #     Search queries against references
    #     :return: None
    #     """
    #     log_progress("Searching queries", msglevel = 0, loglevel = "INFO")

    #     # search queries against reference databases
    #     q_num = len(self.queries)
    #     for qi, query in enumerate(self.queries, 1):
    #         recq_num = len(query.get_record_queries())
    #         if recq_num == 0:
    #             log_progress(f"Query {record_query} has no records. Skipping.", msglevel = 1, loglevel = "WARNING")
    #             continue
            
    #         for rqi, record_query in enumerate(query.get_record_queries(), 1):
    #             # sanity checks to avoid unnecessary searches
    #             done = f"{qi}/{q_num} - {rqi}/{recq_num}"
    #             if not record_query.has_nt() and not record_query.has_aa():
    #                 log_progress(f"[{done}] Searching {record_query} skipped - no sequences", msglevel = 0, loglevel = "DEBUG")
    #                 continue
    #             elif record_query.len() < self.args.min_nt_sig:
    #                 log_progress(f"[{done}] Searching {record_query} skipped - too short)", msglevel = 0, loglevel = "DEBUG")
    #                 continue
    #             else:
    #                 log_progress(f"[{done}] Searching {record_query}", msglevel = 0, loglevel = "INFO")

    #             # prepare query files
    #             query_nt_path = os.path.join(self.dirs['query'], f'q{qi}_r{rqi}_nt.fasta')
    #             query_aa_path = os.path.join(self.dirs['query'], f'q{qi}_r{rqi}_aa.fasta')
    #             output_prefix = os.path.join(self.dirs['search'], f'q{qi}_r{rqi}')
                
    #             # write query files
    #             if record_query.has_nt():
    #                 self.write_fastas(record_query.get_fasta_nt(), query_nt_path)
    #             else:
    #                 log_progress(f"Query {record_query} has no nucleotide sequences. Not writing.", msglevel = 1, loglevel = "WARNING")
    #             if record_query.has_aa():
    #                 self.write_fastas(record_query.get_fasta_aa(), query_aa_path)
    #             else:
    #                 log_progress(f"Query {record_query} has no proteins. Not writing.", msglevel = 1, loglevel = "WARNING")

    #             # search query against targets
    #             for target in sorted(self.targets, key = lambda x: {'fasta_nt': 0, 'fasta_aa': 1, 'hmm': 2, 'mmseqs_db': 3}[x.type]):
    #                 output_path = ''
    #                 if target.type == 'fasta_nt':
    #                     if record_query.has_nt():
    #                         output_path = target.search(query_nt_path, output_prefix)
    #                 else:
    #                     if record_query.has_aa():
    #                         output_path = target.search(query_aa_path, output_prefix)

    #                 # add signal only if output file exists and is not empty
    #                 if not output_path or (not os.path.exists(output_path) and os.path.getsize(output_path) > 0):
    #                     continue

    #                 # parse search results
    #                 nt_signal_array, aa_signal_array = target.read_output(record_query, output_path)
    #                 signal_group = 'nt_based'
    #                 if nt_signal_array is not None:
    #                     record_query.add_signal(signal_group, target.name, nt_signal_array)
    #                 if aa_signal_array is not None and target.type != 'fasta_nt':
    #                     signal_group = 'aa_based'
    #                     record_query.add_signal(signal_group, target.name, aa_signal_array)
                
    #             # combine signals
    #             if self.args.combine: record_query.combine_signals()

    #             # print signal summary
    #             record_query.print_signal_summary()

    # # not used - as of now
    # def evaluate_signals(self) -> None:
    #     """
    #     Evaluate signals
    #     :return: None
    #     """
    #     log_progress("Evaluating signals", msglevel = 0, loglevel = "INFO")
    #     for rqi, record_query in enumerate(self.record_queries, 1):
    #         log_progress(f"[{rqi}/{len(self.record_queries)}] Evaluating {record_query}", msglevel = 1, loglevel = "DEBUG")
    #         record_query.evaluate(self.args.max_nt_gap, self.args.min_nt_sig, self.args.max_aa_gap, self.args.min_aa_sig, self.args.min_sig_frac)
    #         record_query.merge_regions()
    #         for regions in record_query.get_regions().values():
    #             self.regions.extend(regions)

    #     log_progress(f"In total {len(self.regions)} regions were identified", msglevel = 0, loglevel = "INFO")

    def filter_regions(self, sig_groups : Union[str, List[str]] = None, sig_sources : Union[str, List[str]] = None, status : Union[str, List[str]] = None) -> List:
        """
        Filter regions based on provided parameters
        :param sig_groups: signal groups to choose
        :param sig_sources: signal sources to choose
        :param status: status to choose
        :return: list of filtered regions                        
        """

        # make a list if a string is provided
        if isinstance(sig_groups, str):
            sig_groups = [sig_groups]
        if isinstance(sig_sources, str):
            sig_sources = [sig_sources]
        if isinstance(status, str):
            status = [status]

        # set sig_sources
        # at the moment only 'all' or 'combined' are handled
        if (not sig_sources) or (sig_sources == ['all']):
            sig_sources = ['all', 'combined', 'merged']

        # set sig_groups
        if (not sig_groups) or (sig_groups == ['all']):
            sig_groups = ['nt_based', 'aa_based', 'merged']

        regions = []
        for region in self.regions:
            if 'all' in sig_sources:
                pass
            elif region.signal_source not in sig_sources:
                continue
            if region.signal_group not in sig_groups:
                continue
            if status and region.status not in status:
                continue

            regions.append(region)

        return regions

    def filter_hq_regions(self):
        """
        Filter high quality regions only
        :return: list of high quality regions
        """
        
        for region in self.filter_regions(sig_groups = 'merged', sig_sources = 'merged', status = ['CheckV Complete', 'CheckV High-quality', 'CheckV Medium-quality']):
            if region.status in ['CheckV Complete', 'CheckV High-quality']:
                if region.len() >= self.args.checkv_len and region.get_sig_frac() >= self.args.checkv_sig_frac:
                    self.hq_regions.append(region)
            elif region.status in ['CheckV Medium-quality']:
                if region.len() >= self.args.checkv_mq_len and region.get_sig_frac() >= self.args.checkv_mq_sig_frac:
                    self.hq_regions.append(region)
    
    def run_checkv(self) -> None:
        """
        Runs CheckV to automatically evaluate the completeness and contamination of candidate phage regions.
        After that the output file is read non-overlapping high-quality and complete phage regions are picked.
        """

        log_progress("Running CheckV", msglevel = 0, loglevel = "INFO")
        # prepare output file path
        checkv_dir = self.dirs['checkv']
        checkv_env = '' if not self.args.checkv_env else f"conda run -n {self.args.checkv_env} "
        checkv_db = '' if not self.args.checkv_db else f" -d {self.args.checkv_db} "
        cmd = f"{checkv_env}checkv end_to_end {checkv_db} {self.args.outdir}/regions/candidate.fasta {checkv_dir} -t {self.args.threads} --remove_tmp"
        call_process(cmd, program="checkv")

        log_progress("Processing CheckV output data", msglevel = 0, loglevel = "INFO")
        qual_sum_path = os.path.join(checkv_dir, 'quality_summary.tsv')
        cont_path = os.path.join(checkv_dir, 'contamination.tsv')
        # read quality summary file and filter high-quality and complete records
        df = pd.read_csv(qual_sum_path, sep = '\t', keep_default_na=False)
        cdf = pd.read_csv(cont_path, sep = '\t', keep_default_na=False)
        ckeys = list(set(df.columns).intersection(set(cdf.columns)))
        checkv = df.merge(cdf, on = ckeys, how='left').to_dict('records', into=OrderedDict)

        # update checkv information for each region
        for vr in checkv:
            self.get_region(vr['contig_id']).update_checkv(vr)

    def run_phispy(self, query, prefix) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Runs PhiSpy to identify prophages in the input genome.
        :param query: query object
        :param prefix: prefix for the output file
        :return: dictionary with phispy results
        """

        log_progress("Running PhiSpy", msglevel = 1, loglevel = "DEBUG")
        # prepare output file path
        cmd = f"PhiSpy.py --output_dir {self.dirs['phispy']} --threads {self.args.threads} --output_choice 1 --color {query.get_file_path()} -p {prefix}"
        call_process(cmd, program="phispy")

        log_progress("Processing PhiSpy output data", msglevel = 1, loglevel = "DEBUG")
        phispy_output = os.path.join(self.dirs['phispy'], f"{prefix}_prophage_coordinates.tsv")
        if not os.path.exists(phispy_output):
            log_progress(f"PhiSpy output file {phispy_output} not found!", loglevel="WARNING")
            return {}
        else:
            df = pd.read_csv(phispy_output, sep = '\t', names = ['pp', 'contig', 'start', 'end', 'start_att1', 'start_att2', 'end_att1', 'end_att2', 'start_att', 'end_att', 'comment'], keep_default_na=False)          

            nt_signal_array = {}
            for record_id in df['contig'].unique():
                nt_signal_array[record_id] = np.zeros(query.get_record_len(record_id))
                for i, row in df.loc[df['contig'] == record_id].iterrows():
                    nt_signal_array[record_id][int(row['start']) - 1 : int(row['end'])] += 1
        
        return nt_signal_array

    def list_regions(self, regions : List) -> None:
        """
        Logs the list of selected regions.
        :param regions: a list of Region objects
        """

        for region in regions:
            log_progress(f"{region} - {region.status}", msglevel=1, loglevel="INFO")
        
    ### get methods ###
    def get_region(self, header : str):
        """
        Returns region based on a unique header sequence.
        :param header: Region header string
        :return: Region object
        """

        for region in self.regions:
            if region.header == header:
                return region
        else:
            log_progress(f"Missing {header} in regions!", loglevel="ERROR")
            exit(1)

    ### write methods ###
    def write_regions(self, regions : List, group : str, format : Union[str, List[str]] = 'fasta') -> None:
        """
        Write regions to file
        :param regions: list of Region objects
        :param group: group name, either 'candidate' or 'verified'
        :param format: format of the sequence to write: fasta or genbank
        :return: None
        """
        log_progress(f"Writing regions", msglevel = 0, loglevel = "INFO")
        # prepare output file paths
        regions_dir = self.dirs['regions']
        if isinstance(format, str):
            format = [format]
        
        for f in format:
            if f == 'fasta':
                regions_file_path = os.path.join(regions_dir, f"{group}.fasta")
                self.write_fastas([region.to_fasta_nt() for region in regions], regions_file_path)
            elif f == 'genbank':
                regions_file_path = os.path.join(regions_dir, f"{group}.gb")
                self.write_genbanks([region.to_genbank() for region in regions], regions_file_path)

    def write_fastas(self, fastas : Union[str, List[str]], out_path : str) -> None:
        """
        Write FASTA sequences to file
        :param fastas: list of FASTA sequences
        :param out_path: path to output file
        :return: None
        """

        if isinstance(fastas, str):
            fastas = [fastas]

        cnt = 0
        with open(out_path, 'w') as out:
            for fasta in fastas:
                out.write(fasta)
                cnt += fasta.count('>')
        log_progress(f"wrote {cnt} sequence{'s' if cnt > 1 else ''} to {out_path}", msglevel = 1, loglevel = "DEBUG")

    def write_genbanks(self, seqrecs : Union[SeqRecord, List[SeqRecord]], out_path : str) -> None:
        """
        Write SeqRecords of regions to GenBank file
        :param seqrecs: list of SeqRecord objects
        :param out_path: path to output file
        :return: None
        """

        if isinstance(seqrecs, str):
            seqrecs = [seqrecs]

        with open(out_path, 'w') as out:
            SeqIO.write(seqrecs, out, 'genbank')
        
        log_progress(f"wrote {len(seqrecs)} sequence{'s' if len(seqrecs) > 1 else ''} to {out_path}", msglevel = 1, loglevel = "DEBUG")

    def write_summary(self, sep : str = '\t') -> None:
        """
        Write summary of the results
        :param sep: table separator
        :return: None
        """

        log_progress("Writing summary", msglevel = 0, loglevel = "INFO")
        # prepare output file paths
        gen_sum_path = os.path.join(self.args.outdir, 'summary.all_regions.tsv')
        cand_sum_path = os.path.join(self.args.outdir, 'summary.candidates.tsv')
        ver_sum_path = os.path.join(self.args.outdir, 'summary.verified.tsv')

        # write general summary
        with open(gen_sum_path, 'w') as out:
            for region in self.filter_regions(sig_groups=['all'], sig_sources=['all']):
                out.write(sep.join(map(str, region.to_list())) + "\n")
        log_progress(f"wrote summary to {gen_sum_path}", msglevel = 1, loglevel="INFO")    
        
        with open(cand_sum_path, 'w') as out:
            for region in self.filter_regions(sig_groups=['merged'], sig_sources=['merged']):
                out.write(sep.join(map(str, region.to_list())) + "\n")
        log_progress(f"wrote summary to {cand_sum_path}", msglevel = 1, loglevel="INFO")
        
        with open(ver_sum_path, 'w') as out:
            for region in self.hq_regions:
                out.write(sep.join(map(str, region.to_list())) + "\n")
        log_progress(f"wrote summary to {ver_sum_path}", msglevel = 1, loglevel="INFO")
            
    def write_artemis_plots(self) -> None:
        """
        Writes Artemis plot files for all query records.
        :return: None
        """

        log_progress("Writing Artemis plot files", msglevel = 0, loglevel = "INFO")
        # prepare output file path
        artemis_dir = self.dirs['artemis_plots']
        for record_query in self.record_queries:
            record_query.artemis_plot(artemis_dir)

class Input():
    """
    General input data class
    """

    def __init__(self, file_path : str, type : str, name : str = None):
        self.file_path = file_path
        self.type = type
        self.name = name if name else os.path.basename(file_path)

    def get_file_path(self) -> str:
        """
        Get file path
        :return: file path
        """

        return self.file_path
class Target(Input):
    """
    Target data class
    """

    def __init__(self, file_path : str, type : str, params : Dict, name : str = None):
        super().__init__(file_path, type, name)
        self.params = params
        self.db_path = ''
        self.output_paths = []

    ### default methods ###
    def __str__(self):
        """
        String representation of SigMaRef
        """

        return f"Reference: {self.name} [{self.type}]"

    ### custom methods ###
    def prepare(self, outdir_path : str):
        """
        Prepare target data for searching and set paths
        :param outdir_path: path to output directory
        """

        log_progress(f"Creating reference {self.type} database from {self.name}", msglevel=1, loglevel="INFO")
        # set db_path and results_path
        self.db_path = os.path.join(outdir_path, f"{self.name}.db")

        if self.type == 'fasta_nt':
            cmd = 'makeblastdb -in {} -dbtype nucl -out {}'.format(self.file_path, self.db_path)
            call_process(cmd)
        elif self.type == 'fasta_aa':
            cmd = 'diamond makedb --in {} --db {} --quiet'.format(self.file_path, self.db_path)
            call_process(cmd)
        elif self.type == 'mmseqs_db':
            self.db_path = self.file_path        

    def search(self, query_path : str, output_prefix : str) -> str:
        """
        Run search
        :param query_path: path to FASTA file
        :param output_prefix: prefix for output files
        :return: path to output file
        """
        
        self.output_paths.append(f"{output_prefix}.{self.name}.{self.type}.tsv")
        output_path = self.output_paths[-1]

        if self.params.reuse and os.path.exists(output_path):
            log_progress(f"reusing {output_path}", msglevel=1, loglevel="DEBUG")
            return output_path
        else:
            log_progress(f"searching against {self.name}", msglevel = 1, loglevel = "DEBUG")

        if self.type == 'fasta_nt':
            cmd = 'blastn -query {} -db {} -out {} -outfmt "6 qaccver saccver pident length qstart qend evalue" -evalue {} -perc_identity {} -num_threads {} -max_target_seqs 10000'.format(query_path, self.db_path, output_path, self.params.nt_evalue, self.params.nt_pident, self.params.threads)
            call_process(cmd)

        elif self.type == 'fasta_aa':
            cmd = 'diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid sseqid evalue pident qcovhsp scovhsp --evalue {} --id {} --query-cover {} --subject-cover {} --very-sensitive -c1 --threads {} --max-target-seqs 0 --quiet'.format(query_path, self.db_path, output_path, self.params.aa_evalue, self.params.aa_pident, self.params.aa_qscovs, self.params.aa_qscovs, self.params.threads)
            call_process(cmd)

        elif self.type == 'mmseqs_db':
            outdir = os.path.dirname(output_path)
            tmp_path = os.path.join(outdir, f"mmseqs_tmp")
            protein_db_path = os.path.join(tmp_path, 'proteindb')
            results_path = os.path.join(tmp_path, 'mmseqs.out')

            # create tmp directory
            if not os.path.exists(tmp_path):
                os.makedirs(tmp_path)

            # create proteindb
            cmd = f"mmseqs createdb {query_path} {protein_db_path} -v 2"
            call_process(cmd)

            # make a search
            cmd = f"mmseqs search -s {self.params.mmseqs_sens} --threads {self.params.threads} -e {self.params.mmseqs_evalue} {self.db_path} {protein_db_path} {results_path} {tmp_path} -v 2"
            call_process(cmd)

            # create a results tsv file
            cmd = f"mmseqs createtsv {self.db_path} {protein_db_path} {results_path} {output_path} -v 2"
            call_process(cmd)

            # delete tmp directory
            shutil.rmtree(tmp_path)

        return output_path

    def read_output(self, query, output_path : str) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, Dict[str, np.ndarray]]]:
        """
        Read diamond output and return information about regions with signal
        :param file_path: path to diamond output file
        :param query: Query object
        :return: dictionary with information about regions with signal
        :return: a touple with two dictionaries of dictionaries for Query and their nucleotide-based signal and proteins with list of similar reference proteins
        """

        nt_signal_array = {}
        aa_signal_array = {}

        # don't read output if it doesn't exist or is empty
        if not os.path.exists(output_path):
            return nt_signal_array, aa_signal_array
        elif os.path.getsize(output_path) == 0:
            return nt_signal_array, aa_signal_array

        if self.type == 'fasta_nt':
            for line in open(output_path, 'r'):
                if line.startswith('#'):
                    continue
                else:
                    qaccver, saccver, pident, length, qstart, qend, evalue = line.strip().split('\t')
                    record_id, qlength = qaccver.split("|")
                    signal = 1 # TODO: add option to use pident or evalue as signal

                    if record_id not in nt_signal_array:
                        qlength = int(qlength)
                        nt_signal_array[record_id] = np.zeros(int(qlength))

                    # record nt signal
                    if int(length) >= self.params.nt_length:
                        nt_signal_array[record_id][int(qstart) - 1 : int(qend)] += signal

        elif self.type == 'fasta_aa':
            for line in open(output_path, 'r'):
                if line.startswith('#'):
                    continue
                else:
                    qseqid, sseqid, evalue, pident, qcovhsp, scovhsp = line.strip().split('\t')
                    signal = 1 # TODO: add option to use pident, qcovhsp, scovhsp or evalue as signal
                    record_id, protein_id, protein_coords = qseqid.split('|')
                    start, end, strand = map(int, protein_coords.split('..'))

                    # record nt equivalent aa signal
                    if record_id not in nt_signal_array:
                        nt_signal_array[record_id] = np.zeros(query.get_record_len(record_id))
                        aa_signal_array[record_id] = np.zeros(query.get_record_cds_num(record_id))

                    # record nt signal
                    nt_signal_array[record_id][int(start) - 1 : int(end)] += signal
                    # record aa signal
                    aa_signal_array[record_id][query.get_record_cds_order_num(record_id, protein_id)] += signal

        elif self.type == 'mmseqs_db':
            for line in open(output_path, 'r'):
                if line.startswith('#'):
                    continue
                else:
                    qseqid, sseqid, alnScore, seqIdentity, eVal, qStart, qEnd, qLen, tStart, tEnd, tLen = line.strip().split('\t')
                    signal = 1 # TODO: add option to use pident, qcovhsp, scovhsp or evalue as signal
                    record_id, protein_id, protein_coords = sseqid.split('|')
                    start, end, strand = map(int, protein_coords.split('..'))

                    # record nt equivalent aa signal
                    if record_id not in nt_signal_array:
                        nt_signal_array[record_id] = np.zeros(query.get_record_len(record_id))
                        aa_signal_array[record_id] = np.zeros(query.get_record_cds_num(record_id))

                    # record nt signal
                    nt_signal_array[record_id][int(start) - 1 : int(end)] += signal
                    # record aa signal
                    aa_signal_array[record_id][query.get_record_cds_order_num(record_id, protein_id)] += signal
                
        return nt_signal_array, aa_signal_array

class Query(Input):
    """
    Query data class
    """

    def __init__(self, file_path : str, type : str, name : str = None):
        super().__init__(file_path, type, name)
        self.records = []
        self.cdss = []

        if self.type == 'fasta':
            self.records = [RecordQuery(record) for record in parse_fasta(self.file_path)]
        elif self.type == 'genbank':
            cnt = 0
            for record in sorted(parse_genbank(self.file_path), key = lambda x: len(x.seq), reverse=True):
                self.records.append(RecordQuery(record, file_path))
                self.cdss.extend(self.records[-1].get_cdss())
                log_progress(f"{record.id}: {self.records[-1].len()} bps and {len(self.cdss) - cnt} CDSs", msglevel = 1, loglevel = "DEBUG")
                cnt = len(self.cdss)

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
    def get_record(self, record_id : str):
        """
        Returns record with record_id.
        :return: RecordQuery
        """

        for record in self.records:
            if record.get_id() == record_id:
                return record

    def get_record_queries(self) -> List:
        """
        Returns query records.
        :return: List(RecordQuery)
        """

        return self.records

    def get_record_len(self, record_id : str) -> int:
        """
        Returns record length with record_id.
        :return: int
        """

        return self.get_record(record_id).len()

    def get_record_cds_num(self, record_id : str) -> int:
        """
        Returns number of CDSs in record with record_id.
        :return: int
        """

        return self.get_record(record_id).get_cdss_num()

    def get_record_cds_order_num(self, record_id : str, cds_id : str) -> int:
        """
        Returns CDS order number in record with record_id.
        :return: int
        """

        return self.get_record(record_id).get_cds_order_num(cds_id)

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
    
class Record():
    """
    General record class
    """

    def __init__(self, record : SeqRecord, file_path : str):
        self.record : SeqRecord = record
        self.name : str = f"{self.record.id}|{self.len()}"
        self.cdss : List[SeqFeature] = self.get_features_of_type('CDS')
        self.cdss_idx : Dict[str : int] = {cds.qualifiers['protein_id'][0] : i for i, cds in enumerate(self.cdss)}
        self.nt_path : str = ""
        self.aa_path : str = ""
        self.file_path : str = file_path

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

    def has_nt(self) -> bool:
        """
        Returns True if nt_path is not empty.
        :return: bool
        """

        return len(self.get_record().seq) > 0

    def has_aa(self) -> bool:
        """
        Returns True if aa_path is not empty.
        :return: bool
        """

        return len(self.get_cdss()) > 0

    ### get methods ###
    def get_id(self) -> str:
        """
        Returns record id
        :return: str
        """

        return self.record.id

    def get_record(self) -> SeqRecord:
        """
        Returns a SeqRecord object
        :return: SeqRecord object
        """

        return self.record
    
    def get_cdss(self) -> List:
        """
        Returns a list of CDSs
        :return: List(CDS)
        """

        return self.cdss

    def get_cdss_num(self) -> int:
        """
        Returns number of CDSs
        :return: int
        """

        return len(self.cdss)

    def get_cds_order_num(self, cds_id : str) -> int:
        """
        Returns CDS order number in record with record_id.
        :return: int
        """

        try:
            return self.cdss_idx[cds_id]
        except KeyError:
            log_progress(f"Could not find CDS {cds_id} in record {self.record.id}", msglevel = 1, loglevel = "ERROR")
            log_progress(f"Available CDSs: {self.cdss_idx.items()}", msglevel = 1, loglevel = "ERROR")
            log_progress(f"Available CDSs: {[cds.qualifiers['protein_id'][0] for cds in self.cdss]}", msglevel = 1, loglevel = "ERROR")
            raise KeyError

    def get_features_of_type(self, ftype: str, min_cds_len : int = 14) -> List[SeqFeature]:
        """
        Get features of a given type from SeqRecord
        :param ftype: type of a feature
        :param min_cds_len: minimum length of the feature to consider. Set to 14 as this is minimum for MMseqs2 "Maybe the sequences length is less than 14 residues."
        :return:
        """

        flist = []
        # filter out features of NoneType that may mess up the sorting
        pref_feautures = [f for f in self.record.features if f is not None]
        for fcnt, feature in enumerate(sorted(pref_feautures, key = lambda x: (x.location.start, x.location.end - x.location.start), reverse = False), 1):
            if feature.type == ftype:
                if ftype == 'CDS':
                    if 'translation' not in feature.qualifiers:
                        feature.qualifiers['translation'] = [feature.extract(self.record.seq).translate(table = 11, to_stop = True)]
                    if 'protein_id' not in feature.qualifiers:
                        feature.qualifiers['protein_id'] = [f'{self.record.id}_ft_{fcnt:06d}']
                    if 'record_id' not in feature.qualifiers:
                        feature.qualifiers['record_id'] = [self.record.id]
                    # check the CDS length
                    if len(feature.qualifiers['translation'][0]) < min_cds_len:
                        log_progress(f"CDS {feature.qualifiers['protein_id'][-1]} is shorter than {min_cds_len} aa which breaks MMseqs2 run. Skipping.", msglevel = 1, loglevel = "WARNING")
                        continue
                flist.append(feature)
        
        return flist

    def get_cds_header_and_aa(cds: SeqFeature) -> List[str]:
        """
        Get amino acid sequence of CDS feature
        :param cds: a SeqFeature object
        :return: a list of header and amino acid sequence
        """
        
        return [f"{cds.qualifiers['record_id'][0]}|{cds.qualifiers['protein_id'][0]}|{int(cds.location.nofuzzy_start)}..{cds.location.nofuzzy_end}..{cds.strand}", cds.qualifiers['translation'][0]]
        
    def get_fasta_nt(self) -> str:
        """
        Returns a fasta string of the Record
        :return: fasta string
        """

        return f">{self.name}\n{format_seq(self.record.seq)}\n"

    def get_fasta_aa(self) -> str:
        
        """
        Return a fasta string of the query CDSs
        :return: fasta string
        """
        fasta = ""
        for cds in self.cdss:
            fasta += f">{cds.qualifiers['record_id'][0]}|{cds.qualifiers['protein_id'][0]}|{int(cds.location.nofuzzy_start) + 1}..{cds.location.nofuzzy_end}..{cds.strand}\n{format_seq(cds.qualifiers['translation'][0])}\n"

        return fasta
    
    ### write methods ###
    def write_fasta(self, output_path : str, fasta : str) -> None:
        """
        Writes a fasta file of the Record
        :param output_path: path to save fasta file
        :param fasta: fasta string
        """

        with open(output_path, 'w') as f:
            f.write(fasta)

class RecordQuery(Record):
    """
    Single record query class
    """

    def __init__(self, record : SeqRecord, file_path : str):
        super().__init__(record, file_path)
        self.regions = {'nt_based': [], 'aa_based': [], 'merged': []}
        self.signals = {'nt_based': {}, 'aa_based': {}}

    ### default methods ###
    def __str__(self):
        """
        Returns a string representation of the RecordQuery
        :return: str
        """

        return f"{self.name}"

    ### get methods ###
    def get_regions(self):
        """
        Return list of regions
        """
        return self.regions

    def _get_signal_df(self) -> pd.DataFrame:
        """
        Internal function to make a pandas DataFrame with the signal data.
        :param record_id: SeqRecord id
        :return: pandas DataFrame
        """

        signal_df = pd.DataFrame()
        signal_group = 'nt_based'
        for ref, signal in self.signals[signal_group].items():
            colname = f"{signal_group}_{ref}"
            signal_df[colname] = signal

        return signal_df
        
    ### custom methods ###
    def print_signal_summary(self) -> None:
        """
        Returns a summary of the signals in the query.
        """
        
        for signal_group, signal_names in self.signals.items():
            if len(self.signals[signal_group]) == 0: continue
            log_progress(f"{signal_group} signal group for {self.name}:", msglevel = 1, loglevel='DEBUG')
            for signal_name, signal in signal_names.items():
                log_progress(f"{signal_name}: {sum([1 if x else 0 for x in signal])}/{signal.size} of total {sum(signal)}", msglevel = 2, loglevel='DEBUG')

    def artemis_plot(self, output_dir : str) -> None:
        """
        Writes Artemis plot files of the query regions.
        :param output_dir: output directory
        """

        output_path = os.path.join(output_dir, f"{self.record.id}.artemis.plot")
        write_df_to_artemis(self._get_signal_df(), output_path)

    def add_signal(self, signal_group : str, signal_name : str, signal_array : np.ndarray) -> None:
        """
        Add signal array to the query.
        :param signal_group: signal group
        :param signal_name: name of the signal
        :param signal_array: NumPy array with signal values
        """
        
        self.signals[signal_group][signal_name] = signal_array

    def combine_signals(self) -> None:
        """
        Combines all signals to consider separately
        """

        for signal_group in self.signals.keys():
            combined_array = np.zeros(self.len() if signal_group == 'nt_based' else self.get_cdss_num())
            for signal_array in self.signals[signal_group].values():
                combined_array += signal_array
            self.add_signal(signal_group, 'combined', combined_array)

    def evaluate(self, max_nt_gap : int, min_nt_signals : int, max_aa_gap : int, min_aa_signals : int, min_sig_frac : float, sig_groups : List[str]) -> None:
        """
        Evaluate signal for each approach and database #TODO allow for combining signal from different approaches
        :param max_nt_gap: max nt gap between signal
        :param min_nt_sig: min nt signal
        :param max_aa_gap: max aa gap between signal
        :param min_aa_sig: min aa signal
        :param min_sig_frac: min signal fraction
        :param sig_sources: list of signal sources to consider
        """
        log_progress(f"Using " + ", ".join(sig_groups) + " signal groups.", msglevel = 1, loglevel='DEBUG')
        if 'all' in sig_groups or 'combined' in sig_groups:
            sig_groups = ['nt_based', 'aa_based']
        for signal_group in self.signals.keys():
            # skip the signal if not in the list of signal groups to consider
            if signal_group not in sig_groups:
                continue
            # setup thresholds for signal considerations and regions merging
            if signal_group == 'nt_based':
                max_gap_size = max_nt_gap
                min_sig_signals = min_nt_signals
            elif signal_group == 'aa_based':
                max_gap_size = max_aa_gap
                min_sig_signals = min_aa_signals
            elif signal_group == 'merged':
                continue
            min_signal_frac = min_sig_frac
        
            log_progress(f"{signal_group} evaluation of {self.name}", msglevel = 1, loglevel='DEBUG')
            for target, signal_array in self.signals[signal_group].items():
                log_progress(f"{np.count_nonzero(signal_array)} {'positions' if signal_group == 'nt_based' else 'proteins'} based on {target}", msglevel = 2, loglevel='DEBUG')
                # don't perform evaluation if less than min_sig_signals
                if np.count_nonzero(signal_array) < min_sig_signals: 
                    continue
                i_pos = -1
                i_gap = -1
                pos_len = 0
                gap_size = 0
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
                                    region_start = self.cdss[i_pos].location.nofuzzy_start
                                    region_end = self.cdss[i_gap - max_gap_size + 1].location.nofuzzy_end
                                candidate_region = Region(
                                    file_path = self.file_path,
                                    record = self.get_record()[region_start : region_end],
                                    start = region_start,
                                    end = region_end,
                                    signal_source = target,
                                    signal = signal_array[i_pos : i_gap - max_gap_size + 1],
                                    signal_group = signal_group,
                                    rno = len(self.regions[signal_group]) + 1
                                    )
                                
                                self.regions[signal_group].append(candidate_region)
                                log_progress(f"{len(self.regions[signal_group])}: {candidate_region}", msglevel = 3, loglevel='DEBUG')
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
                        z_pos = i - gap_size # last positive signal
                        if signal_group == 'nt_based':
                            region_start = i_pos
                            region_end = z_pos
                        elif signal_group == 'aa_based':
                            region_start = self.cdss[i_pos].location.nofuzzy_start
                            region_end = self.cdss[z_pos].location.nofuzzy_end
                        candidate_region = Region(
                            file_path = self.file_path,
                            record = self.get_record()[region_start : region_end],
                            start = region_start,
                            end = region_end,
                            signal_source = target,
                            signal = signal_array[i_pos : z_pos],
                            signal_group = signal_group,
                            rno = len(self.regions[signal_group]) + 1,
                            )
                        self.regions[signal_group].append(candidate_region)
                        log_progress(f"{len(self.regions[signal_group])}: {candidate_region}", msglevel = 3, loglevel='DEBUG')

    def merge_regions(self) -> None:
        """
        Merge regions to dereplicate them
        1. Get all regions
        2. Sort them by start position
        3. Check if the next region overlaps with the previous one
        4. If the next region overlaps completely with the previous one, remove the shorter one - as it is easier to trim rather than expand automatically or manually
        """

        log_progress(f"merging overlapping regions of {self.name}", msglevel = 1, loglevel='DEBUG')
        
        # get all regions
        regions = sorted(self.regions['nt_based'] + self.regions['aa_based'], key = lambda x: [x.start, -x.len()])
        # get coverage of all regions across sequence
        regions_cov = np.zeros(self.len())
        merged_signal = np.zeros(self.len())
        if 'combined' in self.signals['nt_based'].keys():
            merged_signal += self.signals['nt_based']['combined']
        else:
            # for signal_group in self.signals['nt_based'].keys():
            for signal_array in self.signals['nt_based'].values():
                merged_signal += signal_array

        # overlay positions of regions
        for region in regions:
            regions_cov[region.start : region.end] += 1
        
        # find non-overlapping regions
        # based on this handy reply https://stackoverflow.com/a/27642744
        edges, = np.nonzero(np.diff((regions_cov==0)*1))
        edge_vec = [edges + 1]
        if regions_cov[0] != 0:
            edge_vec.insert(0, [0])
        if regions_cov[-1] != 0:
            edge_vec.append([len(regions_cov)])
        edges = np.concatenate(edge_vec)
        ranges = zip(edges[::2], edges[1::2])

        # create merged regions
        signal_group = 'merged'
        for region_start, region_end in ranges:
            candidate_region = Region(
                file_path = self.file_path,
                record = self.get_record()[region_start : region_end],
                start = region_start,
                end = region_end,
                signal_source = signal_group,
                signal = merged_signal[region_start : region_end],
                signal_group = signal_group,
                rno = len(self.regions[signal_group]) + 1
                )
            self.regions[signal_group].append(candidate_region)
            log_progress(f"{len(self.regions[signal_group])}: {candidate_region}", msglevel = 2, loglevel='DEBUG')

        log_progress(f"found {len(self.regions['merged'])} non-overlapping regions: {len(self.regions['nt_based'])} based on nt signal group and {len(self.regions['aa_based'])} based on aa signal group", msglevel = 1, loglevel='DEBUG')

class Region(Record):
    """
    Region class
    """

    def __init__(self, file_path : str, record : SeqRecord, start : int, end : int, signal : np.ndarray, signal_source : str, signal_group : str, rno : int, category : str = 'pp', status : str = 'candidate'):
        self.file_path = file_path
        self.record = record
        self.start = start
        self.end = end
        self.signal = signal
        self.signal_source = signal_source
        self.signal_group = signal_group
        self.rno = rno
        self.category = category
        self.status = status
        self.id = f"{self.record.id}_{self.category}{self.rno}"
        self.desc = f"{self.category}{self.rno} of {self.record.description}"
        self.name = ''
        self.header = ''
        self.checkv = {}

        self.update_name()
        self.update_header()

    ### default methods ###
    def __repr__(self):
        return f"{self.name} ({self.len()}; {self.get_sig_frac():0.2f})"

    def __str__(self):
        return self.__repr__()

    ### get methods ###
    def get_sig_frac(self):
        return np.count_nonzero(self.signal) / len(self.signal)

    ### custom methods ###
    def update_name(self):
        """
        Updates the region name
        """

        self.name = f"{self.record.id}|{self.signal_source}|{self.signal_group}|{self.rno}|{self.len()}|{self.start + 1}..{self.end}"
    
    def update_header(self):
        """
        Updates the header string
        """

        self.header = f"{self.record.id}|{self.start + 1}..{self.end}..{self.end - self.start}|{self.signal_source}|{self.signal_group}|{self.rno}|{self.category}|{self.status}"

    def update_checkv(self, d : dict):
        """
        Updates Region based on CheckV output data.
        """
        self.checkv = d
        # change status
        self.status = f"CheckV {d['checkv_quality']}"
        
        # update record, coords and signal if considered for further use
        if self.status.startswith('CheckV'):
            region_types = d['region_types']
            if 'host' in region_types: # extract viral region
                region_types = region_types.split(',')
                region_coords_bp = d['region_coords_bp'].split(',')
                start, end = map(int, region_coords_bp[region_types.index('viral')].split('-'))
                start -= 1
                # coords
                self.start += start
                self.end -= self.len() - end
                # record
                self.record = self.record[start : end]
                # signal
                self.signal = self.signal[start : end]
            # update name and header
            self.update_name()
            self.update_header()

    def to_fasta_nt(self) -> str:
        """
        Write a fasta file of the query records
        :return: FASTA string
        """

        return f">{self.header}\n{format_seq(self.record.seq)}\n"

    def to_genbank(self) -> SeqRecord:
        """
        Returns updated SeqRecord object
        """
        
        rec = self.get_record()
        rec.id = self.id
        rec.description = self.desc
        rec.annotations['molecule_type'] = 'DNA'

        return rec

    def to_list(self) -> str:
        """
        Returns a list of the region data
        :retrun: separator-separated string
        """

        return [
            self.file_path,
            self.id,
            self.desc,
            self.record.id,
            self.start + 1,
            self.end,
            self.len(),
            self.signal_source,
            self.signal_group,
            self.rno,
            self.category,
            self.status,
            self.name,
            self.header] + list(self.checkv.values())