#!/usr/bin/env python3
from collections import defaultdict
import os
import ont_fast5_api.fast5_interface
import pandas as pd
import numpy as np
import sys
import joblib
import datetime
import click
import multiprocessing
import pickle
import gzip
import subprocess
import time
#import PolyAcaller

"""
version 0.6.2 2021.07.27
add debug plot

version 0.6.1 2021.07.10
compress_basecall fix bug when fast5 is error

version 0.6 2021.06.23
We used med_mad funciton from bonito (https://github.com/nanoporetech/bonito/blob/master/bonito/fast5.py) to scale the raw
singal. Note: guppy basecalled results also contain scaling median and deviation.
We first find PolyA region, and if max_extend_base is not 0, we try to extend the polyA region. 


version 0.5 2021.06.11
  used for test. contain some find polya method
version 0.4 2021.06.11

version 0.3 2021.05.20

version 0.2 2020.12.31
version 0.1 wirte by Jia Jinbu 2020.09.12

Identify polyA region from raw nanopore signal.

see more by help.

                    genome mapping region (|)      unmapping region (-)  3' adapter (*)
read: ---***----|||||||||||||||||||||||||||||||-------------------------***************----
search region:                       ||||||||||-------------------------*****
                                     |--10nt--|                         |5nt|


If you want to read fast5pos:
Please use function: PolyAcaller.read_read_fast5pos(filein). It return a dict to store the origin absoult fast5 path of 
each read.

"""
    
def parallel_extract_polya(fast5_to_reads, threads=1, basecall_group="", debug_dir="", return_debug_reads=False, **kwargs):
    #print(datetime.datetime.now())
    """
    fast5_to_reads is a dataframe, require read_id, polyA_type (T | A), search_start_base, search_end_base, 
                file_fast5 column. Option column: read_core_id.
    
    return: A dataframe, include columns: 
            ["read_id" or "read_core_id", "polya_start_raw", "polya_end_raw", 
            "polya_start_base", "polya_end_base", "polya_length", 
            "polya_score", "polya_type"]
            If provide read_core_id in fast5_to_reads, return read_core_id, or return read_id column.
    """
    max_threads = multiprocessing.cpu_count()
    if threads > max_threads:
        threads = max_threads
    if threads > 1:
        with joblib.Parallel(threads) as pool:                
            if "raw_file_fast5" in fast5_to_reads.columns:
                res = pool(
                        joblib.delayed(extract_polya_from_reads)(file_fast5, adapter_data, basecall_group, raw_file_fast5, return_debug_reads, debug_dir, **kwargs)
                        for (file_fast5, raw_file_fast5), adapter_data in fast5_to_reads.groupby(["file_fast5", "raw_file_fast5"])
                    )
            else:
                res = pool(
                        joblib.delayed(extract_polya_from_reads)(file_fast5, adapter_data, basecall_group, return_debug_reads, debug_dir, **kwargs)
                        for file_fast5, adapter_data in fast5_to_reads.groupby("file_fast5")
                    )
    else:
        if "raw_file_fast5" in fast5_to_reads.columns:
            res = []
            for (file_fast5, raw_file_fast5), adapter_data in fast5_to_reads.groupby(["file_fast5", "raw_file_fast5"]):
                res.append(extract_polya_from_reads(file_fast5, adapter_data, basecall_group, raw_file_fast5,  return_debug_reads, debug_dir, **kwargs))
        else:
            res = []
            for file_fast5, adapter_data in fast5_to_reads.groupby("file_fast5"):
                res.append(extract_polya_from_reads(file_fast5, adapter_data, basecall_group, return_debug_reads, debug_dir, **kwargs))
    
    real_res = []
    all_reads = []
    for real_r, reads in res:
        real_res.append(real_r)
        all_reads.extend(reads)
    df = pd.concat(real_res)
    print("parallel_extract_polya result df (line 99)")
    print(len(df))
    return [df, all_reads]
    #print(datetime.datetime.now())

def get_header(filein):
    """
    Return the header line (the first line) of one file.
    String format with "\n"
    """
    header = ""
    with open(filein) as f:
        header = next(f)
    return header

def get_header_index(filein, select_column_names, sep="\t"):
    """
    Return the indexs (0-based) of select_column_names.
    If has duplicated column name, return the minmum index.
    
    Input:
    filein: has column name, tab-seperated
    select_column_names: column name list
    Return:
    a list, the indexs of select_column_names
    """
    select_indexs = []
    column_names = get_header(filein).rstrip("\n").split(sep)
    select_indexs = [column_names.index(column_name) for column_name in select_column_names]
    return select_indexs

def read_by_column_names(filein, select_column_names, sep="\t"):
    select_indexs = get_header_index(filein, select_column_names, sep)
    with open(filein) as f:
        next(f)
        for l in f:
            d = l.rstrip("\n").split(sep)
            yield([d[i] for i in select_indexs]) 
        
def read_read_fast5pos(file_fast5pos=None, file_sequencing_summary=None, fast5_dir=None, return_df=False):
    
    def _read_pickle(file_fast5pos):
        if file_fast5pos.endswith(".gz"):
            read2file = pickle.load(gzip.open(file_fast5pos, 'rb'))
        else:
            read2file = pickle.load(open(file_fast5pos, 'r'))
        file_fast5_pos_dir = os.path.dirname(file_fast5pos)
        for read_id, file_name in  read2file.items():
            file_path = os.path.join(file_fast5_pos_dir, file_name)
            read2file[read_id] = file_path
        return read2file
    
    def _read_fast5pos_and_index(file_fast5pos):
        read2file = {}
        index_file = file_fast5pos + ".index.txt"
        file_fast5_pos_dir = os.path.dirname(file_fast5pos)
        index2file = {}
        for l in open(index_file):
            file_name, index = l.rstrip("\n").split("\t")[:2]
            if not os.path.isabs(file_name):
                file_name = os.path.join(file_fast5_pos_dir, file_name)
            index2file[index] = file_name
        for l in open(file_fast5pos):
            read_id, index = l.rstrip("\n").split("\t")[:2]
            read2file[read_id] = index2file[index]
        return read2file
        
    def _read_sequencing_summary(file_sequencing_summary, fast5_dir):
        read2file = {}
        for filename, read_id in read_by_column_names(file_sequencing_summary, ["filename", "read_id"]):
            filename = os.path.join(fast5_dir, filename)
            read2file[read_id] = filename
        return read2file
    
    if file_fast5pos:
        if file_fast5pos.endswith(".pkl.gz") or file_fast5pos.endswith(".pkl"):
            read2file = _read_pickle(file_fast5pos)
        else:
            read2file = _read_fast5pos_and_index(file_fast5pos)
    else:
        read2file = _read_sequencing_summary(file_sequencing_summary, fast5_dir)
    if return_df:
        read2file = pd.DataFrame({"read_id": read2file.keys(), "file_fast5": read2file.values()})
    return read2file

def read_adapter_info(fileadapter, file_fast5pos=None, file_sequencing_summary=None, fast5_dir=None, raw_fast5dir=None, file_select_reads=""):
    
    """
    Input:
    fileadapter
    file_sequencing_summary (optional) generated by Guppy
    fast5_dir  (optional)
    
    Default read fileadapter as pandas dataframe.
    Require this column, but if not exist will be generated by other column except polyA_type:
    read_core_id, file_fast5, polyA_type, search_start_base, search_end_base
    
    read_id will be generated based on read_core_id.
    If not file_fast5, need file_sequencing_summary and fast5_dir
    If not search_start_base, search_end_base:
        require: r_align_start, f_align_end, genome_align_start, genome_align_end
    
    Output:
    A dataframe include read_core_id, read_id, file_fast5, polyA_type, search_start_base, search_end_base
    """
    
    PAD_MAP_LENGTH = 10
    PAD_PRIMER_LENGTH = 5
    
    d = pd.read_table(fileadapter)
    print(f"nb rows fileadapter: {fileadapter}")
    print(len(d))
    
    if file_select_reads:
        select_reads = pd.read_csv(file_select_reads, sep="\t").read_core_id
        d = d.loc[d.read_core_id.isin(select_reads), :]
    
    if "read_id" not in d.columns:
        d["read_id"] = d["read_core_id"].map(lambda x: x.split(",")[0])
    if "file_fast5" not in d.columns:
        read2fast5file = read_read_fast5pos(file_fast5pos, file_sequencing_summary, fast5_dir)
        d["file_fast5"] = d["read_id"].map(read2fast5file)
    
    #append raw file fast5
    if raw_fast5dir:
        tmp_fast5_2_raw_file = {}
        for tmp_file_fast5 in d.file_fast5.unique():
            tmp_raw_file_fast5 = os.path.join(raw_fast5dir, os.path.split(tmp_file_fast5)[1].split(".")[0] + ".fast5")
            tmp_fast5_2_raw_file[tmp_file_fast5] = tmp_raw_file_fast5
        d["raw_file_fast5"] = d["file_fast5"].map(tmp_fast5_2_raw_file)    
    
    if ("search_start_base" not in d.columns) or ("search_end_base" not in d.columns):
        d["polyA_end"] = d["r_align_start"] - 1
        d["polyT_start"] = d["f_align_end"] + 1
        start_base = d["genome_align_end"] - PAD_MAP_LENGTH + 1
        end_base = d["r_align_start"] + PAD_PRIMER_LENGTH - 1
        tmp_flag = d["polyA_type"] == "T"
        start_base[tmp_flag] = d["f_align_end"][tmp_flag] - PAD_PRIMER_LENGTH  + 1
        end_base[tmp_flag] = d["genome_align_start"][tmp_flag] + PAD_MAP_LENGTH - 1
        start_base[start_base < 1] = 1
        #not used: will be checked in findpolyA method
        #tmp_flag = end_base > d["read_length"]
        #end_base[tmp_flag] = d["read_length"][tmp_flag]
        d["search_start_base"] = start_base
        d["search_end_base"] = end_base
    return(d)

def max_subarray(A):
    """
    !!!The function in included in both polyCaller.py and 
    pacbio_find_polyA.py. They are same function, but haven't
    be put in a module to keep each script can be run independently.
    If you want to modify one of them, please modify them at the 
    same time. 

    Maximum subarray problem: select subarray with maxmium sum
    modified Kadane's algorithm (not greedy)
    
    return (index is 0-based), you can get the subarray by A[start_index:(end_index+1)]:
    [start_index, end_index, sum]
    
    if the maxmium sum is <= 0, then return [-1, -1, 0]
    """
        
    max_ending_here = max_so_far = 0
    max_start_index = startIndex = 0
    max_end_index = -1
    for i, x in enumerate(A):
        if 0 >= max_ending_here + x:
        #For greedy at left side : if 0 > max_ending_here + x:
            startIndex = i+1
            max_ending_here = 0
        else:
            max_ending_here += x
        if max_ending_here > max_so_far:
        #For greedy at right side : if max_ending_here >= max_so_far:
            max_so_far = max_ending_here
            max_start_index = startIndex
            max_end_index = i
    
    if max_so_far <= 0 or (max_start_index > max_end_index):
        return ((-1, -1, 0))
    else:
        return (max_start_index, max_end_index, max_so_far)

def find_continous_true(d):
    #learn from https://stackoverflow.com/questions/35610834/find-consecutive-ones-in-numpy-array
    #d can be bool or int (0, 1)
    #d = d.astype(int)
    #return 0-based
    try:
        d_ext = np.concatenate([[0], d, [0]])
        diffs = np.diff((d_ext==1).astype(int))
        starts = np.where(diffs == 1)[0]
        stops = np.where(diffs == -1)[0]
        lengths = stops - starts
        max_length_idx = np.argmax(lengths)
        return [lengths[max_length_idx], starts[max_length_idx], stops[max_length_idx]]
    except:
        return [0, 0, 0]

def polyA_finder(seq, base="A", match = 1, mis = -1.5):
    """
    !!!The function in included in both polyCaller.py and 
    pacbio_find_polyA.py. They are same function, but haven't
    be put in a module to keep each script can be run independently.
    If you want to modify one of them, please modify them at the 
    same time.
    """
    scores = [match if base == s else mis for s in seq]
    start_index, end_index, max_score = max_subarray(scores)
    return (start_index+1, end_index+1, max_score, seq[start_index:(end_index+1)])

def extract_polya_from_reads(file_fast5, adapter_data, basecall_group="", raw_file_fast5=None, return_debug_reads=False, debug_dir="", **kwargs):
    """
    require import ont_fast5_api.fast5_interface   import pandas as pd
    
    file_fast5: basecalled fast5 file, .pkl, .pkl.gz. 
    adapter_data: Must contain read_id, polyA_type (T | A), search_start_base, search_end_base
                  Also can contain read_core_id, but not is necessary.
    raw_file_fast5: If you don't provide basecalled fast5 file, you can support raw_file_fast5 
                    to get raw signal, but it is not nesseary. Not used for now.
                    Wait for improving.
    return_debug_reads: True or False.
    
    Return: If return_debug_reads: return [results, debug_reads] else return results.
            results is A dataframe, include columns: 
            ["read_id" or "read_core_id", "polya_start_raw", "polya_end_raw", 
            "polya_start_base", "polya_end_base", "polya_length", 
            "polya_score", "polya_type"]
            If the adapter_data didn't contain read_core_id, it return read_id, or return read_core_id.
            debug_reads is a list of read (Fast5Read) object.
    """          
    results = []
    read_core_ids = []
    debug_reads = []
    
    if debug_dir:
        file_pre = os.path.split(file_fast5)[1].split(".")[0]
        this_debug_dir = os.path.join(debug_dir, file_pre)
        if not os.path.exists(this_debug_dir):
            os.makedirs(this_debug_dir)
        this_debug_plot_dir = os.path.join(this_debug_dir, "plot")
        this_debug_adapter_dir = os.path.join(this_debug_dir, "adapter")
        this_debug_read_dir = os.path.join(this_debug_dir, "read")
        this_debug_result_dir = os.path.join(this_debug_dir, "result")
        for tmp_dir in [this_debug_plot_dir, this_debug_adapter_dir, this_debug_read_dir, this_debug_result_dir]:
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)
                
    if debug_dir:
        adapter_data.reset_index(drop=True).to_feather(os.path.join(this_debug_adapter_dir, "adapter.ftr"))
        other_info = {}
        other_info["file_fast5"] = file_fast5
        other_info["basecall_group"] = basecall_group
        other_info["raw_file_fast5"] = raw_file_fast5
        other_info["kwargs"] = kwargs
        with open(os.path.join(this_debug_adapter_dir, "other_info.pkl"), 'wb') as f:
            pickle.dump(other_info, f)
    
    is_fast5 = True
    if file_fast5.endswith(".pkl.gz"):
        is_fast5 = False
        read_infos = pickle.load(gzip.open(file_fast5, 'rb'))["reads"]
    elif file_fast5.endswith(".pkl"):
        is_fast5 = False
        read_infos = pickle.load(open(file_fast5))["reads"]
    else:
        is_fast5 = True
        IN = ont_fast5_api.fast5_interface.get_fast5_file(file_fast5, mode="r")
    
    if raw_file_fast5:
        RAW_IN = ont_fast5_api.fast5_interface.get_fast5_file(raw_file_fast5, mode="r")
    else:
        RAW_IN = None
    
    try:
        for na, row in adapter_data.iterrows():
            read_id = row["read_id"]
            try:
                read_core_id = row["read_core_id"]
            except:
                read_core_id = None
            if is_fast5:
                read = Fast5Read(IN, read_id, basecall_group, config_info = row)
            else:
                read = BasecallPickleRead(read_id, read_infos[read_id], 
                                          config_info = row, 
                                          RAW_IN = RAW_IN)
            result = read.find_polyA(row["polyA_type"], 
                               row["search_start_base"], 
                               row["search_end_base"], **kwargs)
            if read_core_id is not None:
                read_core_ids.append(read_core_id)
            results.append(result)
            debug_reads.append(read)
    except:
        raise
    finally:
        if is_fast5:
            IN.close()
        if raw_file_fast5:
            RAW_IN.close()
    
    column_names = read.find_polyA(return_result_column_name=True)

    results = pd.DataFrame(results, columns=column_names)
    if read_core_ids:
        results["read_core_id"] = read_core_ids
        #move read_core_id to the first column
        results = results[["read_core_id"] + column_names]
    else:
        results["read_id"] = adapter_data.read_id.values.copy()
        results = results[["read_id"] + column_names]
    if debug_dir:
        results.reset_index(drop=True).to_feather(os.path.join(this_debug_result_dir, "result.ftr"))
        for read in debug_reads:
            with open(os.path.join(this_debug_read_dir, read.read_core_id + ".read.pkl"), 'wb') as f:
                pickle.dump(read, f)
            read.plot_polyA(signal_scaled_method=1)
            read.savefig(os.path.join(this_debug_plot_dir, read.read_core_id + ".pdf"))
    
    if return_debug_reads:
        return [results, debug_reads]
    else:
        return [results, []]

def med_mad(x, factor=1.4826):
    """
    Calculate signal median and median absolute deviation
    From https://github.com/nanoporetech/bonito/blob/master/bonito/fast5.py
    """
    med = np.median(x)
    mad = np.median(np.absolute(x - med)) * factor
    return med, mad

def scale_raw_data(raw_data, start_idx=None, end_idx=None):
        
    if not start_idx or start_idx < 0:
        start_idx = 0
    if not end_idx or end_idx > len(raw_data):
        end_idx = len(raw_data)
    med, mad = med_mad(raw_data[start_idx:end_idx])
    return ((raw_data - med) / mad, med, mad)

class Fast5Read():
    
    """
    require package: import numpy as np \n import pandas as pd
    require function: max_subarray
    
    Note: After fast5 IO is closed, the read object generated from fast5 IO didn't work (return None).
        
    """
    
    def __init__(self, IN, read_id=None, basecall_group="", config_info=None, need_scale_data=True):
        self.load_config(config_info)
        self.read(IN, read_id, basecall_group)
        if need_scale_data:
            self.scale_raw_data()
    
    def load_config(self, config_info=None):
        config_options = ["read_core_id", "read_align_strand", "rna_strand", 
                         "primer_type", "genome_align_start", "genome_align_end",
                         "primer_score", "f_primer_type", "f_primer_start", 
                         "f_align_end", "r_primer_type", "r_primer_start",
                         "r_align_start", "polyA_type", "search_start_base", "search_end_base"]
        for s in config_options:
            try:
                setattr(self, s, config_info[s])
            except:
                pass
        try:
            if (not hasattr(self, "search_start_base")) or (not hasattr(self, "search_end_base")):
                self.generate_polya_search_start_end(force_change=True)
        except:
            pass
    
    def get_base_start_raw(self, pos):
        return self.event_table.start.iloc[pos-1]
    
    def get_base_length_raw(self, pos):
        return self.event_table.raw_length[pos-1]
    
    def get_base_end_raw(self, pos):
        return self.get_base_start_raw(pos) + self.get_base_length_raw(pos) - 1
        
    def scale_raw_data(self, use_genome_align_region=False):
        raw_data = self.raw_data
        
        if use_genome_align_region:
            start_base = self.genome_align_start
            end_base = self.genome_align_end
        else:
            start_base = 1
            end_base = len(self.event_table)
            
        scaled_data, med, mad = scale_raw_data(self.raw_data, 
                       self.get_base_start_raw(start_base), 
                       self.get_base_end_raw(end_base))
        self.scaled_data = scaled_data
        return scaled_data

    def generate_polya_search_start_end(self, force_change=False):
        PAD_MAP_LENGTH = 10
        PAD_PRIMER_LENGTH = 5
        
        if self.polyA_type == "A":
            start_base = self.genome_align_end - PAD_MAP_LENGTH + 1
            end_base = self.r_align_start + PAD_PRIMER_LENGTH - 1
            
        else:
            start_base = self.f_align_end - PAD_PRIMER_LENGTH  + 1
            end_base = self.genome_align_start + PAD_MAP_LENGTH - 1
        
        if start_base < 1:
            start_base = 1
        if force_change or not hasattr(self, "search_start_base"):
            self.search_start_base = start_base
        if force_change or not hasattr(self, "search_end_base"):
            self.search_end_base = end_base
        
                        
    def read(self, IN, read_id=None, basecall_group=""):
        
        """
        read.get_analysis_attributes(basecall_group) is a dict:
        {'component': 'basecall_1d',
         'model_type': 'flipflop',
         'name': 'ONT Guppy basecalling software.',
         'segmentation': 'Segmentation_000',
         'time_stamp': '2020-03-10T09:56:33Z',
         'version': '3.3.0+ef22818'}
         
        read.get_summary_data(segmentation_name) 
        {'segmentation': {'adapter_max': 0.0,
         'duration_template': 6465,
         'first_sample_template': 766,
         'has_complement': 0,
         'has_template': 1,
         'med_abs_dev_template': 8.23161506652832,
         'median_template': 86.69466400146484,
         'num_events_template': 3232,
         'pt_detect_success': 0,
         'pt_median': 0.0,
        'pt_sd': 0.0}}

        event_table
        ------------------------------------
        |start|base_index|base|raw_length  |
        |766  |     1    |  G | 12         |
        |778  |     2    |  C | 2          |
        ------------------------------------
        """
        
        if not read_id:
            read_id = IN.get_read_ids()[0]
        read = IN.get_read(read_id)
        if not basecall_group:
            basecall_group = read.get_latest_analysis("Basecall_1D") #'Basecall_1D_000' 'Basecall_1D_001'
        
        self.read_id = read_id
        self.read = read
                
        self.io = IN
        self.basecall_group = basecall_group
        
        template_name = basecall_group + "/BaseCalled_template"  #'Basecall_1D_000/BaseCalled_template'
        self.basecall_group_attributes = read.get_analysis_attributes(basecall_group) 
        segmentation_name = self.basecall_group_attributes['segmentation'] #'Segmentation_000'
        #if self.basecall_group_attributes['model_type'] != 'flipflop':
        #    raise ValueError('model type is not flipflop')
        self.raw_data = read.get_raw_data()
        #raw_data: array([805, 496, 514, ..., 521, 531, 643], dtype=int16)
        self.basecall_summary = read.get_summary_data(basecall_group)['basecall_1d_template']
        self.stride = stride = self.basecall_summary['block_stride']
        self.scaling_median = self.basecall_summary['scaling_median']
        #guppy2用的scaling_med_abs_dev，4用的scaling_mean_abs_dev
        #segmentation里的也存有类似的值，但不太一样。
        try:
            self.scaling_mad = self.basecall_summary['scaling_mean_abs_dev']
        except:
            self.scaling_mad = self.basecall_summary['scaling_med_abs_dev']
        self.skip_prob = self.basecall_summary['skip_prob']
        self.stay_prob = self.basecall_summary['stay_prob']
        self.step_prob = self.basecall_summary['step_prob']
        self.strand_scroe = self.basecall_summary['strand_score']
        self.mean_qscore = self.basecall_summary['mean_qscore']
        self.fastq = read.get_analysis_dataset(group_name=template_name, dataset_name='Fastq')
        self.seq = self.fastq.split("\n")[1]
        self.segmentation_summary = read.get_summary_data(segmentation_name)['segmentation']
        self.start = start = self.segmentation_summary['first_sample_template']
        self.median = self.segmentation_summary['median_template']
        self.mad = self.segmentation_summary['med_abs_dev_template']
        
        self.move = move = read.get_analysis_dataset(group_name=template_name, dataset_name='Move')
        #move: array([1, 0, 0, ..., 0, 0, 0], dtype=uint8)
        self.event_length = event_length = len(self.move)
        #2020.12.14 change: self.end = end = start + (event_length - 1) * stride            
        self.end = end = start + event_length * stride - 1
        
        #generate event table
        ls_move_raw_start = (np.arange(event_length)[move==1])*stride + start
        #2020.12.14 change: ls_move_raw_start_with_end = np.append(ls_move_raw_start, end)
        ls_move_raw_start_with_end = np.append(ls_move_raw_start, end + 1)
        #2020.12.14 change: ls_event_length = ls_move_raw_start_with_end[1:] - ls_move_raw_start
        ls_raw_length = ls_move_raw_start_with_end[1:] - ls_move_raw_start

        # WARNING: Trimming the barcode changes the size of the sequence and crashes here
        #print(len(ls_move_raw_start), "\n", len(self.seq), "\n", len(ls_raw_length))

        self.event_table = pd.DataFrame({"start": ls_move_raw_start, 
                                    "base": list(self.seq), 
                                    "raw_length": ls_raw_length})
        print(self.event_table)
        ##2020.12.14 change:                             "event_length": ls_event_length})
        self.samples_per_nt = np.mean(ls_raw_length[ls_raw_length <=np.quantile(ls_raw_length, 0.95)])
        ##2020.12.14 add:
        self.samples_per_nt_median = np.median(ls_raw_length)


    
    def self_find_polyA(self, **kwargs):
        return self.find_polyA(self.polyA_type, self.search_start_base, self.search_end_base, **kwargs)
        
    def find_polyA(self, find_base="A", start_base=None, end_base=None, 
                   min_polya_length=15, match_score=1, mismatch_score=-1,
                   gap_score_ratio=-5,
                   long_event_base_ratio=15,
                   extend_deta=0.2,
                   extend_point_deta=0.6,
                   max_extend_base=1000,
                   extend_to_genome_max_num=20,
                   extend_to_adapter_max_num=0,
                   extend_mis_score = -5,
                   extend_gap_score = -5,
                   extend_up_score = -5,
                   return_result_column_name=False):
        #only need self.event_table and self.samples_per_nt
        #but modify self value
        #if you want to use old method before 2020.12.15, please set mismatch_score = -1.5 and gap_score_ratio = 0
        #not greedy
        
        if return_result_column_name:
            return ["polya_start_raw", "polya_end_raw", 
                    "polya_start_base", "polya_end_base", "polya_length", 
                    "polya_score", "polya_type",
                    "polya_scale_median", "length_continous_a",
                    "longest_continous_a_start_base", "longest_continous_a_end_base", 
                    "longest_continous_a_start_raw", "longest_continous_a_end_raw", 
                    "longest_continous_a_length", "init_polya_start",
                     "init_polya_end",
                     "init_polya_start_base",
                     "init_polya_end_base",
                     "init_polya_length"]
                     
        #if not provide start_base and end_base, search total read
        #2021.06以后的版本是启发式的。因为把以前的seed向外拓展，这就导致，可能从次优的seed拓展更好
        if not start_base or start_base < 1:
            start_base = 1
        seq_length = len(self.event_table.index)
        if not end_base or end_base > seq_length:
            end_base = seq_length
        #if not provide find_base, try find A and find T, 
        #then return the best result
        
        self.search_start_base = start_base
        self.search_end_base = end_base
        self.search_start_raw = self.get_base_start_raw(start_base)
        self.search_end_raw = self.get_base_end_raw(end_base)
        
        if not find_base:
            polya_result = self.find_polyA("A", start_base, end_base, 
                                min_polya_length, match_score, mismatch_score,
                                           gap_score_ratio,
                                long_event_base_ratio)
            polyt_result = self.find_polyA("T", start_base, end_base, 
                                min_polya_length, match_score, mismatch_score, 
                                           gap_score_ratio,
                                long_event_base_ratio)
            result = polya_result if polya_result[-3] >= polyt_result[-3] else polyt_result
            (self.polya_start, self.polya_end, self.polya_start_base, 
                 self.polya_end_base, self.polya_length, self.polya_score,
                 self.polya_type) = result[:7]
            return result
        #start_base bigger than end_base, wrong
        if start_base > end_base:
            (self.polya_start, self.polya_end, self.polya_start_base, 
                 self.polya_end_base, self.polya_length, self.polya_score,
                 self.polya_type) = result = [0, 0, 0, 0, 0, 0, find_base]
            result.extend([0] * 12)
            return result

        e = self.event_table.iloc[(start_base-1):end_base, :]
        samples_per_nt = self.samples_per_nt

        base_is_right = e.base.values == find_base

        #convert long event base near A to A
        #near_A = np.logical_or(np.insert(base_is_right[:-1], 0, False),
        #              np.append(base_is_right[1:], False))
        #base_is_right[np.logical_and(e.raw_length.values >= samples_per_nt * long_event_base_ratio,
        #               near_A)] = True
    
        #Note: This might raise flase-positive (pause of Polymerase?)
        need_convert_base = np.logical_and(e.raw_length.values >= samples_per_nt * long_event_base_ratio,
                       np.logical_not(base_is_right))
        base_is_right[need_convert_base] = True
    
        #call polya
        if gap_score_ratio:
            #method1
            #scores = np.where(base_is_right, match_score * e.raw_length.values, gap_score_ratio * samples_per_nt)
            scores = np.where(base_is_right, match_score, mismatch_score) *  e.raw_length.values
            open_gap = np.logical_and(np.logical_not(base_is_right), 
                              np.insert(base_is_right[:-1], 0, base_is_right[0]))
            scores[open_gap] = gap_score_ratio * samples_per_nt
        else:
            scores = np.where(base_is_right, match_score, mismatch_score) \
                    * e.raw_length.values
        start_index, end_index, max_score = max_subarray(scores) 
        
        #calculate polya position
        if start_index != -1: # -1 indicate not find polyA
            self.polya_start = polya_start = init_polya_start = e.start.values[start_index]
            #2020.12.14 change self.polya_end = e.start.values[end_index] + e.event_length.values[end_index]
            self.polya_end = polya_end = init_polya_end = e.start.values[end_index] + e.raw_length.values[end_index] - 1
            self.polya_start_base = polya_start_base = init_polya_start_base = start_base + start_index
            self.polya_end_base = polya_end_base = init_polya_end_base = start_base + end_index
            self.polya_length = polya_length = init_polya_length = round((self.polya_end - self.polya_start + 1)/samples_per_nt, 2)
            
            scaled_data = self.scaled_data
            polya_scale_median  = np.median(scaled_data[(self.polya_start-1):self.polya_end])
            
            #最多向基因组比对区域向上再搜索20bp，接头序列不再向上搜索
            #如果该碱基和polya_scale_median相差大于0.2，分为-5*length，否则为length
            #如果该碱基和polya_scale_median相差小于0.2，
            #但该碱基中存在比polya_scale_median相差超过0.6的点，该碱基得分为length，
            #但会再末端加一个gap，空罚分为-5*samples_per_nt
            #从左到右扫描，求最大值
            #PolyA中最后一个碱基，也添加在罚分中。
            #最后获得最大不为0的累积分所在的位置。
            if find_base == "A":
                extend_to_up_max_num, extend_to_down_max_num = extend_to_genome_max_num, extend_to_adapter_max_num
            else:
                extend_to_up_max_num, extend_to_down_max_num = extend_to_adapter_max_num, extend_to_genome_max_num            
            
            def get_extend_base_num(tmp_event_table_region, to_down=True):
                tmp_scores = []
                tmp_gaps = []
                for i, (na, row) in enumerate(tmp_event_table_region.iterrows()):
                    tmp_s = row["start"]
                    tmp_l = row["raw_length"]
                    tmp_e = tmp_s + tmp_l - 1
                    tmp_signals = scaled_data[(tmp_s - 1): tmp_e]
                    tmp_signal = np.median(tmp_signals)
                    tmp_signal_diff = tmp_signal - polya_scale_median
                    tmp_add_gap = False
                    tmp_gap_score = 0
                    if np.abs(tmp_signal_diff) <= extend_deta:
                        tmp_score = tmp_l
                        if i == 0:
                            if to_down:
                                tmp_point_high_num = (np.abs(tmp_signals[-50:] - polya_scale_median) > extend_point_deta).sum()
                            else:
                                tmp_point_high_num = (np.abs(tmp_signals[:50] - polya_scale_median) > extend_point_deta).sum()
                        else:
                            tmp_point_high_num = (np.abs(tmp_signals - polya_scale_median) > extend_point_deta).sum()
                        if tmp_point_high_num:
                            tmp_add_gap = True
                            tmp_gap_score = extend_gap_score * samples_per_nt
                        if i == 0 and not to_down:
                            tmp_add_gap = True
                            tmp_gap_score = extend_up_score * samples_per_nt                      
                    else:
                        tmp_score = extend_mis_score * tmp_l
                        if i == 0 and not to_down:
                            tmp_score = extend_up_score * samples_per_nt
                    tmp_scores.append(tmp_score)
                    if tmp_add_gap:
                        tmp_scores.append(tmp_gap_score)
                        tmp_gaps.append(i)
                tmp_scores = np.cumsum(tmp_scores)
                tmp_max_idx = np.argmax(tmp_scores)
                tmp_max_score = tmp_scores[tmp_max_idx]
                if tmp_max_score > 0:
                    for tmp_gap_idx in tmp_gaps:
                        if tmp_gap_idx < tmp_max_idx:
                            tmp_max_idx -= 1
                        else:
                            break
                else:
                    tmp_max_idx = 0
                return tmp_max_idx
            
            if max_extend_base > 0 and (polya_length >= 5 or ((end_index - start_index + 1) >= 5)):
                #对于较短的polyA，4.0guppy很容易产生大量持续短的A。
                #end_base是开始定义的搜索终止位置
                tmp_end_search_base = min(end_base + extend_to_down_max_num, polya_end_base + max_extend_base) #1-based
                down_extend_base_num = get_extend_base_num(self.event_table.iloc[(polya_end_base-1):tmp_end_search_base,:], to_down=True)
                
                end_index += down_extend_base_num
                self.polya_end_base = polya_end_base = polya_end_base + down_extend_base_num
                self.polya_end = self.event_table.start.values[polya_end_base - 1] + self.event_table.raw_length.values[polya_end_base - 1] -1
            
                tmp_start_search_base = max(0, start_base - extend_to_up_max_num, polya_start_base - max_extend_base)
                up_extend_base_num = get_extend_base_num(self.event_table.iloc[(polya_start_base-1):(tmp_start_search_base-2):-1, :], to_down=False)

                start_index -= up_extend_base_num
                self.polya_start_base = polya_start_base = polya_start_base - up_extend_base_num
                self.polya_start = polya_start = self.event_table.start.values[polya_start_base - 1]
                
                self.polya_length =  polya_length = round((self.polya_end - self.polya_start + 1)/samples_per_nt, 2)
            
            polya_base_seq = self.event_table.base.values[(polya_start_base-1):(polya_end_base)]
            length_continous_a, start_a_base_idx, end_a_base_idex = find_continous_true(polya_base_seq == find_base)
            if length_continous_a:
                longest_continous_a_start_base = start_a_base_idx + polya_start_base
                longest_continous_a_end_base = end_a_base_idex + polya_start_base - 1
                longest_continous_a_start_raw = self.event_table.start.values[longest_continous_a_start_base - 1]
                longest_continous_a_end_raw = self.event_table.start.values[longest_continous_a_end_base - 1] + self.event_table.raw_length.values[longest_continous_a_end_base - 1] -1
                longest_continous_a_length = round((longest_continous_a_end_raw - longest_continous_a_start_raw + 1)/samples_per_nt, 2)
            else:
                longest_continous_a_start_base = 0
                longest_continous_a_end_base = 0
                longest_continous_a_start_raw = 0
                longest_continous_a_end_raw = 0
                longest_continous_a_length = 0
        else:
            self.polya_start = 0
            self.polya_end = 0
            self.polya_start_base = 0
            self.polya_end_base = 0
            self.polya_length  = 0
            
            init_polya_start = 0
            init_polya_end = 0
            init_polya_start_base = 0
            init_polya_end_base = 0
            init_polya_length = 0

            polya_scale_median = 0
            length_continous_a = 0
            longest_continous_a_start_base = 0
            longest_continous_a_end_base = 0
            longest_continous_a_start_raw = 0
            longest_continous_a_end_raw = 0
            longest_continous_a_length = 0
        self.polya_score = max_score
        self.polya_type = find_base
        
        return([self.polya_start, self.polya_end, self.polya_start_base, 
             self.polya_end_base, self.polya_length, self.polya_score,
             self.polya_type, 
             polya_scale_median, length_continous_a, longest_continous_a_start_base,
             longest_continous_a_end_base, longest_continous_a_start_raw,
             longest_continous_a_end_raw, longest_continous_a_length, 
             init_polya_start,
             init_polya_end,
             init_polya_start_base,
             init_polya_end_base,
             init_polya_length
             ])
    
    def plot_polyA(self,  extend_xlim_left=0, extend_xlim_right=0, **kwargs):
        if "xlim" not in kwargs:
            kwargs["xlim"] = [self.search_start_raw - extend_xlim_left, self.search_end_raw + 1 + extend_xlim_right]
        if "plot_base" not in kwargs:
            kwargs["plot_base"] = True
        if "plot_base_line" not in kwargs:
            kwargs["plot_base_line"] = True
        self.plot(**kwargs)
        self.ax.set_title(self.read_id + " {:.1f}".format(self.polya_length) + " {:.2f}".format(self.samples_per_nt))
    
    def plot(self, figsize = None, plot_base=False, plot_base_line=False,
             plot_base_median_line=False,
             plot_genome_line = False,
             plot_adapter_line = False,
             xlim=None, ylim=None, signal_scaled_method=0, fig=None, ax=None):
        import matplotlib.pyplot as plt
        if signal_scaled_method == 0:
            raw_data = self.raw_data
        else:
            raw_data = self.scaled_data            
        event_table = self.event_table #don't change
        end = self.end
        start = self.start
        
        if plot_base_median_line:
            try:
                polya_median_singal_value = np.median(raw_data[(self.get_base_start_raw(self.polya_start_base)-1):self.get_base_end_raw(self.polya_end_base)])
            except:
                polya_median_singal_value = np.median(raw_data)
            polya_median_diff_thre = 0.2
            if signal_scaled_method == 0:
                try:
                    polya_median_diff_thre = 0.2 * np.std(raw_data[(polya_next_start-1):(polya_next_start+polya_next_long-1)])
                except:
                    polya_median_diff_thre = 0.2 * np.std(raw_data)
            
        plot_polya = True
        try:
            polya_start = self.polya_start
            polya_end = self.polya_end
            polya_type = self.polya_type
            polya_length = self.polya_length
        except:
            plot_polya = False
        
        if ax is None:
            if figsize:
                fig, ax = plt.subplots(figsize=figszie)
            else:
                fig, ax = plt.subplots()
            
        ax.plot(np.arange(len(raw_data)), raw_data, zorder=2)
        ax.axvspan(start, end+1, facecolor='g', alpha=0.25)
        ax.axvspan(polya_start, polya_end+1, facecolor='r', alpha=0.25)                    

        if plot_base:
            starts = event_table["start"]
            ends = np.append(event_table["start"][1:], end)
            bases = event_table["base"]
            for start, end, base in zip(starts, ends, bases):
                y_pos = raw_data[start:end].max() * 1.02
                x_pos = (start+end-1)/2
                if xlim:
                    if x_pos < xlim[0] or x_pos > xlim[1]:
                        continue
                ax.text(x_pos, y_pos, base, clip_on=True, horizontalalignment="center",
                       verticalalignment='center') #base + str(end-start)
                if plot_base_line:
                    ax.axvline(start, color="grey", zorder=1)
                if plot_base_median_line:
                    base_median_signal = np.median(raw_data[start:end])
                    color = "black"
                    if np.abs(base_median_signal -  polya_median_singal_value) <= polya_median_diff_thre:
                        color = "red"
                    ax.plot([start, end], [base_median_signal, base_median_signal], color=color, zorder=3)
            if plot_base_line:
                ax.axvline(end+1, color="grey", zorder=1)
        
        try:
            ax.axvline(self.get_base_start_raw(self.genome_align_start), color="blue", zorder=2)
            ax.axvline(self.get_base_start_raw(self.genome_align_end) + 1, color="blue", zorder=2)
            #ax.axvspan(self.get_base_start_raw(self.genome_align_start), self.get_base_end_raw(self.genome_align_end)+1,
            #                    facecolor='blue', alpha=0.25)
            ax.axvline(self.get_base_start_raw(self.r_align_start), color="red", zorder=2)
            ax.axvline(self.get_base_end_raw(self.f_align_end+1), color="red", zorder=2)
            ax.set_title(self.read_id)
        except:
            pass
        
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
            
        self.ax = ax
        self.fig = fig
        
    def savefig(self, filename):
        self.fig.savefig(filename)

class BasecallPickleRead(Fast5Read):
    
    def __init__(self, read_id, read_info_list, config_info=None, RAW_IN=None, need_scale_data=True):
        """
        read_info_list is:
        [
            read.stride, 
            read.basecall_summary["mean_qscore"], 
            read.start, 
            read.samples_per_nt,
            read.samples_per_nt_median,
            read.seq, 
            read.event_table.raw_length.values,
            read.scaling_median,
            read.scaling_mad,
            read.skip_prob,
            read.stay_prob,
            read.step_prob,
            read.strand_scroe,
            read.median,
            read.mad]
        """
        self.read_info_from_list(read_id, read_info_list)
        self.load_config(config_info)
        if RAW_IN:
            self.raw_data = RAW_IN.get_read(read_id).get_raw_data()
            if need_scale_data:
                self.scale_raw_data()
                
    def read_info_from_list(self, read_id, read_info_list):
        self.read_id = read_id
        ( self.stride, self.mean_qscore, 
        self.start, self.samples_per_nt, 
        self.samples_per_nt_median, self.seq, raw_length) = read_info_list[:7]
        try:
            if len(read_info_list) > 7:
                (self.scaling_median,
                self.scaling_mad,
                self.skip_prob,
                self.stay_prob,
                self.step_prob,
                self.strand_scroe,
                self.median,
                self.mad) = read_info_list[7:]
        except:
            pass
        raw_start = np.insert(np.cumsum(raw_length[:-1]) + self.start, 0, self.start)
        self.event_table = pd.DataFrame({"start": raw_start, 
                                    "base": list(self.seq), 
                                    "raw_length": raw_length})
        self.end = self.event_table.start.values[-1] + self.event_table.raw_length.values[-1]

def fast5_compress_raw_data(file_fast5, fileout):
    #60M to 35M, 意义有限
    #而且gzip直接压缩fast5文件，也是60M to 38M
    tmp_name = os.path.split(file_fast5)[1]
    with open(os.path.join(os.path.split(fileout)[0], tmp_name + ".init.txt"), 'w') as o:
        o.write("good")
    print(file_fast5)
    try:        
        with ont_fast5_api.fast5_interface.get_fast5_file(file_fast5, mode="r") as IN:
            reads = {}
            all_read_ids = IN.get_read_ids()
            wrong_read_ids = []
            all_right_ids = []
            for i, read_id in enumerate(all_read_ids):
                try:
                    read = IN.get_read(read_id)
                    raw_data = read.get_raw_data()
                    reads[read_id] = raw_data
                    all_right_ids.append(read_id)
                except:
                    wrong_read_ids.append(read_id)        
        with gzip.open(fileout, 'wb') as f:
            pickle.dump(reads, f)
        return [all_right_ids, wrong_read_ids]
    except:
        return [[], []]

def fast5_compress_basecall(file_fast5, fileout):
    tmp_name = os.path.split(file_fast5)[1]
    with open(os.path.split(fileout)[0] + "/" + tmp_name + ".init.txt", 'w') as o:
        o.write("good")
    print(file_fast5)
    try:        
        with ont_fast5_api.fast5_interface.get_fast5_file(file_fast5, mode="r") as IN:
            reads = {}
            all_read_ids = IN.get_read_ids()
            high_qulality_reads = []
            low_qulality_reads = []
            base_info = {}
            wrong_reads = []
            for i, read_id in enumerate(all_read_ids):                
                try:
                    read = Fast5Read(IN, read_id, need_scale_data=False)
                    if i == 0:
                        base_info["basecall_group_attributes"] = read.basecall_group_attributes
                    reads[read_id] = [
                        read.stride, 
                        read.basecall_summary["mean_qscore"], 
                        read.start, 
                        read.samples_per_nt,
                        read.samples_per_nt_median,
                        read.seq, 
                        read.event_table.raw_length.values,
                        read.scaling_median,
                        read.scaling_mad,
                        read.skip_prob,
                        read.stay_prob,
                        read.step_prob,
                        read.strand_scroe,
                        read.median,
                        read.mad
                        #,
                        #read.mean_qscore                
                        ]
                    if read.basecall_summary["mean_qscore"] >= 7:
                        high_qulality_reads.append([read_id, read.seq])
                    else:
                        low_qulality_reads.append([read_id, read.seq])
                except:
                    wrong_reads.append(read_id)        
        with gzip.open(fileout, 'wb') as f:
            pickle.dump({"base_info": base_info, "reads": reads, "reads_key": [
                "stride",
                "mean_qscore",
                "start",
                "samples_per_nt",
                "samples_per_nt_median",
                "seq",
                "event_table_raw_length",
                "scaling_median",
                "scaling_mad",
                "skip_prob",
                "stay_prob",
                "step_prob",
                "strand_scroe",
                "median",
                "mad"
            ]}, f)
        return [high_qulality_reads, low_qulality_reads, wrong_reads]
    except:
        return [[], [], []]

def fast5_compress_basecall_dir(fast5_dirs, fileout_dir, threads=1, read2file=None):
    
    def iter_fast5_dir_inout(fast5_dirs, fileout_dir):        
        if not os.path.exists(fileout_dir):
            os.makedirs(fileout_dir)
        
        for fast5_dir in fast5_dirs:
            for filein_name in os.listdir(fast5_dir):
                filein_path = os.path.join(fast5_dir, filein_name)
                file_pre, file_last = os.path.splitext(filein_name)
                if file_last != ".fast5":
                    continue
                fileout_name = file_pre + ".basecalled.pkl.gz"
                fileout_path = os.path.join(fileout_dir, fileout_name)
                yield((filein_path, fileout_path))
    
    print(time.asctime( time.localtime(time.time()) ))
    #print(datetime.datetime.now())
    with joblib.Parallel(threads) as pool:    
        
        filein_fileout_pair = list(iter_fast5_dir_inout(fast5_dirs, fileout_dir))
        
        res = pool(
                joblib.delayed(fast5_compress_basecall)(filein, fileout)
                for filein, fileout in filein_fileout_pair
            )
        
        fileout_high_seq = os.path.join(fileout_dir, "all.qscore_more7.fa")
        fileout_low_seq = os.path.join(fileout_dir, "all.qscore_less7.fa")
        fileout_high_seq_gz = os.path.join(fileout_dir, "all.qscore_more7.fa.gz")
        fileout_low_seq_gz = os.path.join(fileout_dir, "all.qscore_less7.fa.gz")
        fileout_read2file = os.path.join(fileout_dir, "read2file.pkl.gz")
        fileout_error_log = os.path.join(fileout_dir, "error.log.txt")
        read_id2fileout = {}
        
        #gzip.open is slower than subprocess(gzip) 60s vs 15s
        with open(fileout_high_seq, 'wt') as o1, open(fileout_low_seq, 'wt') as o2, open(fileout_error_log, 'w') as o3:
            for i, ((filein, fileout), (high_qulality_reads, low_qulality_reads, wrong_reads)) in enumerate(zip(filein_fileout_pair, res)):
                if not high_qulality_reads and not low_qulality_reads:
                    o3.write(filein + "\n")
                if high_qulality_reads:
                    for read_id, seq in high_qulality_reads:
                        fileout_name = os.path.split(fileout)[1]
                        real_fileout_dir = os.path.relpath(fileout_dir, os.path.split(fileout_read2file)[0])
                        if real_fileout_dir == ".": real_fileout_dir = ""
                        read_id2fileout[read_id] = os.path.join(real_fileout_dir, fileout_name)
                        o1.write(f">{read_id}\n{seq}\n")
                if low_qulality_reads:
                    for read_id, seq in low_qulality_reads:
                        fileout_name = os.path.split(fileout)[1]
                        real_fileout_dir = os.path.relpath(fileout_dir, os.path.split(fileout_read2file)[0])
                        if real_fileout_dir == ".": real_fileout_dir = ""
                        read_id2fileout[read_id] = os.path.join(real_fileout_dir, fileout_name)
                        o2.write(f">{read_id}\n{seq}\n")
                for read_id in wrong_reads:
                    o3.write(filein + "\t" + read_id + "\n")
                
        subprocess.run(["gzip", "-c", fileout_high_seq], stdout=open(fileout_high_seq_gz, "wb"))
        subprocess.run(["gzip", "-c", fileout_low_seq], stdout=open(fileout_low_seq_gz, "wb"))
        os.remove(fileout_high_seq)
        os.remove(fileout_low_seq)
        
        """"
        with gzip.open(fileout_high_seq, 'wt') as o1, gzip.open(fileout_low_seq, 'wt') as o2:
            for i, ((filein, fileout), (high_qulality_reads, low_qulality_reads)) in enumerate(zip(filein_fileout_pair, res)):
                for read_id, seq in high_qulality_reads:
                    read_id2fileout[read_id] = fileout
                    o1.write(f">{read_id}\n{seq}\n")
                for read_id, seq in low_qulality_reads:
                    read_id2fileout[read_id] = os.path.split(fileout)[1]
                    o2.write(f">{read_id}\n{seq}\n")
        """
        print(time.asctime( time.localtime(time.time()) ))
        with gzip.open(fileout_read2file, 'wb') as f:
            pickle.dump(read_id2fileout, f)
        print(time.asctime( time.localtime(time.time()) ))

@click.command()
@click.option('-i', '--filein', help=('Input the file of adapter information of each read.'
                                         'Option format1: The file can be generated by bin/adapterFinder.py.'
                                         'In this case, the file contain read_core_id, polyA_type,'
                                         'r_align_start, f_align_end, genome_align_start, genome_align_end.'
                                         'Option format2: You also can provide one file directly contain column '
                                         'read_core_id, polyA_type, search_start_base, search_end_base'
                                         'polyA_type is A or T. If format is not format2, the script will '
                                         'convert them to format1. You also can directly provide file_fast5 column,'
                                         'In this case, you do not need set summary and fast5dir option'), 
                    required=False, default="")
@click.option('-o', '--fileout', required=False, default="", help=("""

                    Output a gzip pickle file of a dict:

                    \b
                    key: read_id
                    value:
                    [
                        stride
                        mean_qscore 
                        read.start
                        read.samples_per_nt
                        read.samples_per_nt_median
                        read.seq 
                        event_table.raw_length.values
                    ]                    
                    """))
@click.option('--fast5dirs', help=('Input the directorys of fast5 files, separated by ,'), required=False, default="")
@click.option("--outdir", required=False, default="", help=("Output directory"))
@click.option('-t', '--threads', required=False, default=10, help='Number of threads to use. (default: 10)')
def compress_basecall_main(filein, fileout, fast5dirs, outdir, threads):
    """
    The compress result didn't contain raw data, you can obtain them from origin guppy file.
    
    import ont_fast5_api.fast5_interface
    IN = ont_fast5_api.fast5_interface.get_fast5_file(file_fast5, mode="r")
    IN.get_read(read_id).get_raw_data()
    """
    max_threads = multiprocessing.cpu_count()
    fast5dirs = fast5dirs.split(",")
    if threads > max_threads:
        threads = max_threads
        
    if filein:
        fast5_compress_basecall(filein, fileout)
    else:
        fast5_compress_basecall_dir(fast5dirs, outdir, threads)


@click.command()
@click.option('-i', '--filein', help=('Input the file of adapter information of each read.'
                                         'Option format1: The file can be generated by bin/adapterFinder.py.'
                                         'In this case, the file contain read_core_id, polyA_type,'
                                         'r_align_start, f_align_end, genome_align_start, genome_align_end.'
                                         'Option format2: You also can provide one file directly contain column '
                                         'read_core_id, polyA_type, search_start_base, search_end_base'
                                         'polyA_type is A or T. If format is not format2, the script will '
                                         'convert them to format1. You also can directly provide file_fast5 column,'
                                         'In this case, you do not need set summary and fast5dir option'), 
                    required=False, default="")
@click.option('-o', '--fileout', required=False, default="", help=("""

                    Output a gzip pickle file of a dict:

                    \b
                    key: read_id
                    value:
                    raw_data                 
                    """))
@click.option('--fast5dirs', help=('Input the directorys of fast5 files, separated by ,'), required=False, default="")
@click.option("--outdir", required=False, default="", help=("Output directory"))
@click.option('-t', '--threads', required=False, default=10, help='Number of threads to use. (default: 10)')
def compress_raw_data_main(filein, fileout, fast5dirs, outdir, threads):
    """
    The compress result didn't contain raw data, you can obtain them from origin guppy file.
    
    import ont_fast5_api.fast5_interface
    IN = ont_fast5_api.fast5_interface.get_fast5_file(file_fast5, mode="r")
    IN.get_read(read_id).get_raw_data()
    """
    max_threads = multiprocessing.cpu_count()
    fast5dirs = fast5dirs.split(",")
    if threads > max_threads:
        threads = max_threads
        
    if filein:
        fast5_compress_raw_data(filein, fileout)
    else:
        fast5_compress_raw_data_dir(fast5dirs, outdir, threads)



@click.command()
@click.option('-i', '--inadapter', help=('Input the file of adapter information of each read.'
                                         'Option format1: The file can be generated by bin/adapterFinder.py.'
                                         'In this case, the file contain read_core_id, polyA_type,'
                                         'r_align_start, f_align_end, genome_align_start, genome_align_end.'
                                         'Option format2: You also can provide one file directly contain column '
                                         'read_core_id, polyA_type, search_start_base, search_end_base'
                                         'polyA_type is A or T. If format is not format2, the script will '
                                         'convert them to format1. You also can directly provide file_fast5 column,'
                                         'In this case, you do not need set summary and fast5dir option'), 
                    required=True, type=click.Path(exists=True))
@click.option('-p', '--fast5pos', help=('Input the fast5 pos of each read.'
                                        'Two differnt format. One is pos txt, and require a index file ".index.txt"'
                                        'And one is read2file.pkl.gz which store a dict, key is read, value is filename (not contain directory path)'
                                        'If the file is .endswith(".pkl.gz") or .endswith(".pkl"), will be processed as the second format.'
                                        'If you do not provide, you must provide -s/-f or provide file_fast5 column in adapter file'),
                    required=False, default="")  
@click.option('-s', '--summary', help=('Input sequencing_summary.txt generated by basecalling software. '
                                       'If you provide -p, the option will be not used. If you provide -s, you must provide -f. If you do not provide -p and -s, you must provide file_fast5 column in adapter file'), 
                    required=False, default="")          
@click.option('-f', '--fast5dir', help=('Input the directory of fast5 files. '
                                        'If you provide -p, the option will be not used. If you provide -f, you must provide -s. If you do not provide -p and -s, you must provide file_fast5 column in adapter file'),
                    required=False, default="")
@click.option('-o', '--out', help=("""Output File of polyA results.

                    \b
                    Output: Tab-seperated fromat
                    read_core_id        e909fc71-7798-4be9-81ea-11b8d5602005,chr1,1185559,1186502
                    polya_start_raw     1185
                    polya_end_raw       1924
                    polya_start_base    73
                    polya_end_base      77
                    polya_length        83.32693538067818
                    polya_score         740
                    polya_type          T
                    polya_start_raw, polya_end_raw: The raw signal event index of potential polyA region
                    polya_start_base, polya_end_base: The basecalled base position
                    polya_type: A or T
                    All is 1-based.
                    """), required=True)
@click.option('-t', '--threads', required=False, default=10, help='Number of threads to use. (default: 10)')
@click.option('-b', '--basecall_group', help=('Basecall group, default is the latest.'
                                        'If only contain one basecall result, it should be Basecall_1D_000'
                                        'If contain two, it should be Basecall_1D_000, Basecall_1D_001, and Basecall_1D_001'
                                        'is the latest.'),
                    required=False, default="")
@click.option('-r', '--raw_fast5dir', help=('Input the directory of raw fast5 files to extract raw signal data. '
                                         'If you provide basecalled result in pkl fromat, you need to provide the directory of raw fast5 file.'
                                        ),
                    required=False, default="")
@click.option('-d', '--debug_dir', help=('Output debug dir'
                                        ),
                    required=False, default="")
@click.option('--file_select_reads', help=('Input only select this reads'
                                        ),
                    required=False, default="")
def main(inadapter, fast5pos, summary, fast5dir, out, threads, basecall_group, raw_fast5dir, debug_dir, file_select_reads):
    """
    Identify polyA region from raw nanopore signal.
    
    \b
    require package:
    ont_fast5_api
    pandas
    numpy
    matplotlib #only if plot
    joblib
    
    """
    fast5_to_reads = read_adapter_info(inadapter, fast5pos, summary, fast5dir, raw_fast5dir, file_select_reads)
    #df, all_reads = PolyAcaller.parallel_extract_polya(fast5_to_reads, threads, basecall_group, debug_dir)
    df, all_reads = parallel_extract_polya(fast5_to_reads, threads, basecall_group, debug_dir)
    #if not use PolyAcaller.parallel_extract_polya, but parallel_extract_polya.
    #it will error when use pickle.dump #local variable 'BasecallPickleRead' referenced before assignment
    df.to_csv(out, sep="\t", index=False)

if __name__ == '__main__':
    main()
    