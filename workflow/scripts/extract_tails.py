#!/usr/bin/env python3
import pandas as pd
import numpy as np
import swifter
import click
import pyfastx
import sys
from Bio import SeqIO
from itertools import zip_longest
import edlib
import regex
import gzip
from tqdm import tqdm
from colorama import Fore, Style, init
from collections import Counter
init(autoreset=True)


equalities=[("M", "A"), ("M", "C"),("R", "A"), ("R", "A"), ("W", "A"), ("W", "A"), ("S", "C"), ("S", "C"), ("Y", "C"), ("Y", "C"), 
("K", "G"), ("K", "G"), ("V", "A"), ("V", "C"), ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"),
("B", "C"), ("B", "G"), ("B", "T"), ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C")]
    

FWD=["F-R", "F-UR", "F-UUR", "F-UF", "F-F", "F-N", "N-R", "N-N", "N-UR", "N-UUR", "UF-R", "UF-N", "UF-UF", "UF-UR", "UF-UUR", "UR-R", "UUR-R", "UUR-UR"]
REV=["R-F", "R-R", "R-N", "R-UF", "R-UUR", "R-UR", "N-F", "N-UF", "UF-F", "UR-F", "UR-N", "UR-UF", "UR-UR", "UR-UUR", "UUR-F", "UUR-N", "UUR-UF", "UUR-UUR"]
 
@click.command()
@click.option('-i', '--inadapter', help='Input adapter file', required=True, type=click.Path(exists=True))                         
@click.option('-s', '--inseq', help='Input read fastq file', required=True, type=click.Path(exists=True))                             
@click.option('-o', '--out', help='Output adapter information of each read', required=True)
@click.option('-v', '--verbose', is_flag=True, help="Print sequences with colors for polytail, adapter and delimiter")
@click.option('-d', '--debug', is_flag=True, help="for developping purposes, prints additional information")

def main(inadapter, inseq, out, constant_seq="CTGAC", umi_seq="NNNNNNNNNN", adapt_seq="CTGTAGGCACCATCAAT", verbose=False, debug=False):

    fields = ['read_core_id','mRNA','init_polya_start_base', 'init_polya_end_base','primer_type']
    dtypes= {'read_core_id': str,
    'mRNA': str,
    'init_polya_start_base': int,
    'init_polya_end_base': int,
    'primer_type': str}
    
    df = pd.read_csv(inadapter, delimiter = "\t", usecols=fields,dtype=dtypes)
    print("ok")
    
    
    df['readname'] = df['read_core_id'].str.split(",",n=1).str[0]
    #fq = SeqIO.to_dict(SeqIO.parse(gzip.open(inseq, "rt"),'fastq'))
    readnames=df['readname'].to_list
    fq={}

    with gzip.open(inseq, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq")):
            fq[record.id] = record.seq

    tqdm.pandas()

    df['read_seqs']=df.swifter.apply(lambda x: get_seq(fq, x.readname, x.primer_type in REV ), axis=1)


    df[['polytail', 'additional_tail', 'delimiter', 'dist_delim', 'umi', 'adapter', 'dist_adapter']] = df.apply(get_three_primes_parts_row, constant_seq=constant_seq, umi_seq=umi_seq, adapt_seq=adapt_seq, debug=debug, axis=1 )
    df = pd.concat([df, df.apply(lambda col: get_composition(col["polytail"], "A_tail"), axis=1, result_type="expand")], axis = 1)
    df = pd.concat([df, df.apply(lambda col: get_composition(col["additional_tail"], "add_tail"), axis=1, result_type="expand")], axis = 1)
    df.drop('read_seqs', axis=1, inplace=True)
    df = df.replace(r'^\s*$', np.nan, regex=True)
    df.to_csv(out, encoding='utf-8', index=False, sep='\t', na_rep="NA")


def get_three_primes_parts_row(row, constant_seq, umi_seq, adapt_seq, debug):

    three_p_motive=constant_seq+umi_seq+adapt_seq
    read_seq=row['read_seqs']
    primer_type=row['primer_type']

    if row['init_polya_start_base'] == row['init_polya_end_base']:
        
        return pd.Series([np.nan, np.nan, np.nan,np.nan , np.nan, np.nan, np.nan])

    polytail=np.nan
    additional_tail=np.nan
    if primer_type in FWD:
        gene = read_seq[:row['init_polya_start_base']-1]
        polytail = read_seq[row['init_polya_start_base']-1:row['init_polya_end_base']]

        three_p_seq=read_seq[row['init_polya_end_base']:]
        if len(three_p_seq)>200:
            three_p_seq=three_p_seq[:200]


        align = get_umi_alignment(three_p_motive, three_p_seq)
        start_adapt=align['locations'][0][0]
        additional_tail=read_seq[row['init_polya_end_base']:row['init_polya_end_base']+start_adapt]

    elif primer_type in REV:

        gene = read_seq[:len(read_seq) -row['init_polya_end_base']]
        polytail = read_seq[len(read_seq) -row['init_polya_end_base']: len(read_seq) -row['init_polya_start_base']+1]
        three_p_seq=read_seq[len(read_seq) - row['init_polya_start_base']+1:]
        if len(three_p_seq)>200:
            three_p_seq=three_p_seq[:200]

        align = get_umi_alignment(three_p_motive, three_p_seq)
        start_adapt=align['locations'][0][0]
        additional_tail=read_seq[len(read_seq) - row['init_polya_start_base']+1:len(read_seq) - row['init_polya_start_base']+1+start_adapt]  
            
    else :
        raise Exception(f"Unknown primer type: {primer_type}\n{FWD}\n{REV}")

    umi = align["umi"]  
    three_p_seq_without_addtail=three_p_seq[len(additional_tail):]

    delimiter, adapter= three_p_seq.partition(umi)[0][len(additional_tail):], three_p_seq.partition(umi)[2][:len(adapt_seq)] 
    #delimiter, adapter= three_p_seq_without_addtail.partition(umi)[0][len(additional_tail):], three_p_seq_without_addtail.partition(umi)[2][:len(adapt_seq)] 
    if delimiter=="":
        return pd.Series([np.nan, np.nan, np.nan,np.nan , np.nan, np.nan, np.nan])

    
    dist_adapter=edlib.align(adapter, adapt_seq,task="path", mode='HW')["editDistance"]
    dist_delim=edlib.align(delimiter, constant_seq,task="path", mode='HW')["editDistance"] 
    reconstructed_seq=gene+polytail+additional_tail+delimiter + umi + adapter

    

    if debug:
        print("")
        print("###########")
        print(row['read_core_id'])
        print(three_p_seq.partition(umi)[0])
        print(three_p_seq.partition(umi)[2])
        print(primer_type)
        print(three_p_seq)
        print(three_p_seq_without_addtail)
        print("----")
        print(three_p_motive)
        print(align)
        print(read_seq)
        print("gene")
        print(gene)
        print("polytail")
        print(row['init_polya_start_base'], row['init_polya_end_base'])
        print(polytail)
        print("additional_tail")
        print(additional_tail)
        print("delimiter")
        print(delimiter)
        print(delimiter=="")
        print("umi")
        print(umi)
        print("adapter")
        print(adapter)
        print(dist_adapter)
        print(dist_delim)
        input("press enter")

    try:
        assert(reconstructed_seq in read_seq)
    except AssertionError:
        print("")
        print("###########")
        print(row['read_core_id'])
        print(three_p_seq.partition(umi)[0])
        print(three_p_seq.partition(umi)[2])
        print(primer_type)
        print(three_p_seq)
        print(three_p_seq_without_addtail)
        print("----")
        print(three_p_motive)
        print(align)
        print(read_seq)
        print("gene")
        print(gene)
        print("polytail")
        print(row['init_polya_start_base'], row['init_polya_end_base'])
        print(polytail)
        print("additional_tail")
        print(additional_tail)
        print("delimiter")
        print(delimiter)
        print(delimiter=="")
        print("umi")
        print(umi)
        print("adapter")
        print(adapter)
        print(dist_adapter)
        print(dist_delim)

        sys.exit("Assertion error for reconstructed sequence")

    return pd.Series([polytail, additional_tail, delimiter, dist_delim, umi, adapter, dist_adapter])

"""
def get_three_primes_parts(inadapter, inseq, out, constant_seq="CTGAC", umi_seq="NNNNNNNNNN", adapt_seq="CTGTAGGCACCATCAAT", verbose=False, debug=False):
    three_p_motive=constant_seq+umi_seq+adapt_seq

    df = pd.read_csv(inadapter, delimiter = "\t")

    read_core_id_list, mRNA_list, polytail_list, additional_tail_list,delimiter_list, dist_delim_list, umi_list, adapter_list, dist_adapter_list = [], [], [], [], [], [], [], [], []
    Count=0
    CountFwd=0
    CountRev=0
    nb_proc=0

    for readid, mRNA, polyA_start, polyA_end, primer_type in tqdm(zip(df['read_core_id'], df['mRNA'], df['init_polya_start_base'], df['init_polya_end_base'], df['primer_type']), total = df.shape[0]):
        #input("wait input\n\n")
        if polyA_start==polyA_end:
            continue
        Count+=1

        readname = readid.split(",")[0]
        read_seq=get_seq(inseq, readname, primer_type in REV )

        if primer_type in FWD:
            CountFwd+=1
            gene = read_seq[:polyA_start-1]
            polytail = read_seq[polyA_start-1:polyA_end]

            three_p_seq=read_seq[polyA_end:]
            align = get_umi_alignment(three_p_motive, three_p_seq)
            start_adapt=align['locations'][0][0]
            additional_tail=read_seq[polyA_end:polyA_end+start_adapt]

        elif primer_type in REV:
            CountRev+=1
            gene = read_seq[:len(read_seq) -polyA_end]
            polytail = read_seq[len(read_seq) -polyA_end: len(read_seq) -polyA_start+1]
            three_p_seq=read_seq[len(read_seq) - polyA_start+1:]            
            align = get_umi_alignment(three_p_motive, three_p_seq)
            start_adapt=align['locations'][0][0]
            additional_tail=read_seq[len(read_seq) - polyA_start+1:len(read_seq) - polyA_start+1+start_adapt]            

       
        else :
            raise Exception(f"Unknown primer type: {primer_type}")

        umi = align["umi"]
        delimiter, adapter= three_p_seq.partition(umi)[0][len(additional_tail):], three_p_seq.partition(umi)[2][:len(adapt_seq)]            
        dist_adapter=edlib.align(adapter, adapt_seq,task="path", mode='HW')["editDistance"]
        dist_delim=edlib.align(delimiter, constant_seq,task="path", mode='HW')["editDistance"] 
        reconstructed_seq=gene+polytail+additional_tail+delimiter + umi + adapter
        if debug:
            print("")
            print("###########")
            print(read_seq)
            print("gene")
            print(gene)
            print("polytail")
            print(polytail)
            print("additional_tail")
            print(additional_tail)
            print("delimiter")
            print(delimiter)
            print("umi")
            print(umi)
            print("adapter")
            print(adapter)
            print(dist_adapter)
            print(dist_delim)

        assert(reconstructed_seq in read_seq)

        read_core_id_list.append(readid)
        mRNA_list.append(mRNA)
        polytail_list.append(polytail)
        additional_tail_list.append(additional_tail)
        delimiter_list.append(delimiter)
        dist_delim_list.append(dist_delim)
        umi_list.append(umi)
        adapter_list.append(adapter)
        dist_adapter_list.append(dist_adapter)        
        
        if verbose:
            print(f">{readname},-")
            print(f"{Fore.GREEN}{gene}{Fore.LIGHTYELLOW_EX}{polytail}{Fore.GREEN}{delimiter}{Fore.BLUE}{umi}{Style.RESET_ALL}{Fore.GREEN}{adapter}")
            
    print(Counter(dist_adapter_list))
    print(Counter(dist_delim_list))
    res = pd.DataFrame(
        {'read_core_id':read_core_id_list,
          'gene': mRNA_list, 
          'polytail': polytail_list,
          'additional_tail' : additional_tail_list,
          'delimiter': delimiter_list,
          'umi': umi_list, 
          'adapter': adapter_list,
          'dist_adapter_list': dist_adapter_list
          } 
    )

    print(Count, "mRNAs")
    print(CountFwd, "fwd")
    print(CountRev, "rev")
    print(f"{len(umi_list)} identified 3', with {len(set(umi_list))} unique UMIs")


    print("Adding base composition")
    res = pd.concat([res, res.apply(lambda col: get_composition(col["polytail"], "A-tail"), axis=1, result_type="expand")], axis = 1)
    res = pd.concat([res, res.apply(lambda col: get_composition(col["additional_tail"], "add-tail"), axis=1, result_type="expand")], axis = 1)


    res.to_csv(out, encoding='utf-8', index=False)
"""

def get_umi_alignment(seq_to_find, seq_to_look, task="path", mode='HW'):
    wildcard=[]
    seq = seq_to_find
    for c in 'actgACTG':
        seq = seq.replace(c, "")
    wildcard = set(''.join(seq))

    result=edlib.align(seq_to_find, seq_to_look, task=task, mode=mode, additionalEqualities=equalities)
    align = edlib.getNiceAlignment(result, seq_to_find, seq_to_look)


    umi=""

    for q, t in zip(align["query_aligned"], align["target_aligned"]):
        if q not in wildcard:
            continue
        if t == "-":
            umi += "N"
        else:
            umi += t
    result["umi"] = umi

    return(result)


def get_seq(fq, readname, antisense) -> str:

    seq=fq[readname]
    if antisense:
       return str(seq.reverse_complement())
    else:
        return str(seq)

def get_composition(seq: str, col_prefix) -> pd.Series:
    """
    JP
    Get the composition and the percentage of ATGC in a DNA sequence. 
    Ignore all nucleotides other than ATGC (like N).
    """
    nucl_count = {col_prefix+'_A': np.nan, col_prefix+'_T': np.nan, col_prefix+'_G': np.nan, col_prefix+'_C': np.nan}
    nucl_perc = {col_prefix+'_pct_A': np.nan, col_prefix+'_pct_T': np.nan, col_prefix+'_pct_G': np.nan, col_prefix+'_pct_C': np.nan}

    if seq and not pd.isna(seq):

        seq = regex.sub('[^ATCG]', '', seq)
        seq_len=len(seq)
    
        
        if seq_len>0:
            res = Counter(seq)
            #print(res['A'])
            nucl_count = {col_prefix+'_A': res['A'], col_prefix+'_T': res['T'], col_prefix+'_G': res['G'], col_prefix+'_C': res['C']}
            nucl_perc = {col_prefix+'_pct_A': res['A']/seq_len, col_prefix+'_pct_T': res['T']/seq_len, col_prefix+'_pct_G': res['G']/seq_len, col_prefix+'_pct_C': res['C']/seq_len}
    
    return_val = pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])
    
    return return_val



if __name__ == '__main__':
    main()






