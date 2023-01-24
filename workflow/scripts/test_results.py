import sys
import os
import tempfile
import pandas as pd
import numpy as np
import multiprocessing
import gzip
import pyfastx

FWD=["F-R", "F-UR", "F-UUR", "F-UF", "F-F", "F-N", "N-R", "N-N", "N-UR", "N-UUR", "UF-R", "UF-N", "UF-UF", "UF-UR", "UF-UUR", "UR-R", "UUR-R", "UUR-UR"]
REV=["R-F", "R-R", "R-N", "R-UF", "R-UUR", "R-UR", "N-F", "N-UF", "UF-F", "UR-F", "UR-N", "UR-UF", "UR-UR", "UR-UUR", "UUR-F", "UUR-N", "UUR-UF", "UUR-UUR"]
 

def get_fasta_seq_np(df, fasta):
    fa=pyfastx.Fasta(fasta)
    seq_length=df['polya_end_base'] - df[['polya_start_base']]
    if ( seq_length == 0 ):
        return 'NO TAIL'
    elif ( seq_length > 0) :
        return fa.fetch(df['read'], (df['polya_start_base'], df['polya_end_base']))
    else:
        return 'REVERSE TAIL!!'

def get_fasta_seq(df, fasta):
    fa=pyfastx.Fasta(fasta)
    seq_length=df['polya_end_base'] - df[['polya_start_base']]
    if ( seq_length == 0 ):
        return 'NO TAIL'
    elif ( seq_length > 0) :
        return fa.fetch(df['read'], (df['polya_start_base'], df['polya_end_base']))
    else:
        return 'REVERSE TAIL!!'

def get_adapter_interval(row):
    if row['primer_type'] in FWD:
        adapter_interval= (row['polya_end_base'], row['r_align_start']-1)
    elif row['primer_type'] in REV:
        adapter_interval= (row['f_align_end'], row['polya_start_base']-1)
    else:
        sys.exit(f"ERROR: Primer type unknown: {row['primer_type']}")
    return adapter_interval


def read_merged_file_np(file_in, fasta, sep='\t', header=True):
    
    
    data=pd.read_csv(file_in, sep=sep)
    data["read"] = data["read_core_id"].str.split(',', expand=True)[0]
    data['polyA_length'] = data['polya_end_base'] - data['polya_start_base']
    data['polyA_interval'] = list(zip(data['polya_start_base'], data['polya_end_base']-1))
    data['mRNA_interval'] = list(zip(data['genome_align_start'], data['genome_align_end']-1))
    data['adapter_interval'] = data.apply(get_adapter_interval, axis=1) 


    l_errors=[]
    adapters_length=[]
    fa=pyfastx.Fasta(fasta)
    
    for index, row in data.iterrows():
        print("")
        print("row neuve!", index)
        row["comments"]=[]
        if row['type'] != "polya":
            row["comments"].append("-1: not polya")
            l_errors.append(-1)
            continue

        for i, v in list(zip(row.index.values, row.values)):
            print(f"{i} : {v}")

        print("initial_values. DONE")
        print(fa[row['read']])
        print("/read")

        


        if row['adapter_interval'][0] >= row['adapter_interval'][1]:
            row['comments'].append("0: adapter interval impossible")
            l_errors.append(0)
            #input("Press Enter to continue...")
            continue

        if row['polyA_length']==0:
            row['comments'].append("1: polyA length =0")
            l_errors.append(1)
            if row['low_accuracy_3end_mapped']== True:
                row['comments'].append("13: low_accuracy_3end_mapped")
                l_errors.append(13)

                
            #input("Press Enter to continue...")
            #input("Press Enter to continue...")
            continue
        elif row['polyA_length'] ==1:
            row['comments'].append("2: polyA length =1")
            l_errors.append(2)
            if row['low_accuracy_3end_mapped']== True:
                row['comments'].append("3: low_accuracy_3end_mapped")
                l_errors.append(23)

            #input("Press Enter to continue...")
            
            continue
        row, read_seq, polyA_seq, adapter_seq= get_read_parts(row, fasta)

        
        #print(read_seq, polyA_seq, adapter_seq)
        
        print("READ_SEQ:")
        print(read_seq)
        print(row['polyA_length'])
        print(row['polyA_interval'])
        

        print("---")
        print("ADAPTER", row['adapter_interval'])
        if pd.isna(adapter_seq):
            print("0")
            adapters_length.append(0)
        else:
            print(len(adapter_seq))
            adapters_length.append(len(adapter_seq))
            if adapter_seq.startswith("CTGAC"):
                print("ADAPTER STARTS WITH GCTAC")
        print(adapter_seq)
        print("---")
        if pd.isnull(polyA_seq):
            print("NO POLY A")
        else:

            print("POLYA:", row['polyA_interval'], row['polyA_length'])
            print(fa.fetch(row['read'], (row['polyA_interval'][0]-5, row['polyA_interval'][0]-1) ), polyA_seq, fa.fetch(row['read'], (row['polyA_interval'][1]+1, row['polyA_interval'][1]+5)))
        #input("Press Enter to continue...")

    print("over")
    from collections import Counter
    count = Counter(l_errors)
    print(count)
    count = Counter(adapters_length)
    print(count)


def get_read_parts(row, fasta):
    fa=pyfastx.Fasta(fasta)  

    if row['primer_type'] in FWD:
        read_seq=fa[row['read']].seq
        if row['polyA_length'] >1:
            polyA_seq=fa.fetch(row['read'], row['polyA_interval'])
            adapter_seq=fa.fetch(row['read'], row['adapter_interval'] )
        else:
            polyA_seq=pd.NA
            adapter_seq=pd.NA
 
    elif row['primer_type'] in REV:

        read_seq=fa[row['read']].antisense

        if row['polyA_length'] >1:
            print(fa.fetch(row['read'], row['polyA_interval']))
            polyA_seq = fa.fetch(row['read'], row['polyA_interval'], strand='-')            
            adapter_seq=fa.fetch(row['read'], row['adapter_interval'], strand='-')
            
        else:
            polyA_seq=pd.NA
            adapter_seq=pd.NA
       
    else:
        sys.exit(f"ERROR: Primer type unknown: {row['primer_type']}")
    
    return (row, read_seq, polyA_seq, adapter_seq)


def revcom(seq):
    """
    Return the reverse complement sequence of origin sequence.
    
    The origin sequence can be in uppercase and lowercase letters.
    But the output is in uppercase letters.
    All letters not in `ATCGatcg` will be converted to `N`. 
    """

    def complement(seq):
        seq = seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
        def _com(base):
            try:
                return basecomplement[base]
            except:
                return "N"
        letters = list(seq)
        letters = [_com(base) for base in letters]
        return ''.join(letters)
    
            
    return complement(seq[::-1])
    



res="/home/jpeter/DATA/FLEPseq/JP_TestRACE/3_PolyA/barcode01.nanopore.read_info.result.merged.txt"
fasta="/home/jpeter/DATA/FLEPseq/JP_TestRACE/00_intermediate_files/barcode01.nanopore.fasta"
read_merged_file_np(res, fasta)


