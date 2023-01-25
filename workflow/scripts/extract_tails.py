#!/usr/bin/env python3
import pandas as pd
import click
import pysam
import pyfastx
import collections
import regex
from tqdm import tqdm
from colorama import Fore, Style, init
from collections import Counter
init(autoreset=True)

FWD=["F-R", "F-UR", "F-UUR", "F-UF", "F-F", "F-N", "N-R", "N-N", "N-UR", "N-UUR", "UF-R", "UF-N", "UF-UF", "UF-UR", "UF-UUR", "UR-R", "UUR-R", "UUR-UR"]
REV=["R-F", "R-R", "R-N", "R-UF", "R-UUR", "R-UR", "N-F", "N-UF", "UF-F", "UR-F", "UR-N", "UR-UF", "UR-UR", "UR-UUR", "UUR-F", "UUR-N", "UUR-UF", "UUR-UUR"]
 

@click.command()
@click.option('-i', '--inadapter', help='Input adapter file', required=True, type=click.Path(exists=True))
#@click.option('-s', '--inseq', help='Input read fasta file', required=True, type=click.Path(exists=True))                             
@click.option('-s', '--inseq', help='Input read fastq file', required=True, type=click.Path(exists=True))                             
@click.option('-o', '--out', help='Output adapter information of each read', required=True)
@click.option('-v', '--verbose', is_flag=True, help="Print sequences with colors for polytail, adapter and delimiter")


def main(inadapter, inseq, out, verbose):  
    out2=out+"2"
    df = pd.read_csv(inadapter, delimiter = "\t")
    df.head()
    
    #pysam.faidx(inseq)
    
    #fasta = pyfastx.Fastq(inseq)
    String_AAAA_mm = '(%s){e<=0}(%s){e<=2}'% ("CTG", "TAGGCACCATCAAT") ## Search for adaptor in the sequence
    read_core_id_list, mRNA_list, polytail_list, tail_list, adapter_list = [], [], [], [], []
    Count=0
    CountFor=0
    CountRev=0
    CountForA=0
    CountRevA=0
    for readid, mRNA, start, end, primer_type in tqdm(zip(df['read_core_id'], df['mRNA'], df['init_polya_start_base'], df['init_polya_end_base'], df['primer_type']), total = df.shape[0]):
        readname = readid.split(",")[0]
        seq = get_seq(inseq, readname)
        Count += 1
        if primer_type in FWD: # gene_AAAAA + Adapter
            CountFor += 1
            gene = seq[:int(start-1)]
            polytail = seq[int(start-1):int(end)]
            adapter_to_end = seq[int(end):]
            if regex.compile(String_AAAA_mm).search(adapter_to_end):
                research = regex.search(String_AAAA_mm, adapter_to_end, regex.BESTMATCH)
                tail = adapter_to_end[:research.start()]
                adapter =adapter_to_end[research.start():research.start()+17]
                CountForA += 1
                read_core_id_list.append(readid)
                polytail_list.append(polytail)
                adapter_list.append(adapter)
                mRNA_list.append(mRNA)
                tail_list.append(tail)
                if verbose:
                    print(f">{readname},-")
                    print(f"{gene}{Fore.LIGHTYELLOW_EX}{polytail}{Fore.GREEN}{tail}{Fore.BLUE}{adapter}{Style.RESET_ALL}")

        else: # Adapter + 10N + GTCAG + TTTTTT_gene
            CountRev += 1
            adapter_to_end = revcom(seq[:int(start-1)])
            polytail = revcom(seq[int(start-1):int(end)])
            gene = revcom(seq[int(end):])
            if regex.compile(String_AAAA_mm).search(adapter_to_end):
                research = regex.search(String_AAAA_mm, adapter_to_end, regex.BESTMATCH)
                tail = adapter_to_end[:research.start()]
                adapter =adapter_to_end[research.start():research.start()+17]
                CountRevA += 1
                read_core_id_list.append(readid)
                polytail_list.append(polytail)
                adapter_list.append(adapter)
                mRNA_list.append(mRNA)
                tail_list.append(tail)
                if verbose:
                    print(f">{readname},+")
                    #print(f"{Fore.BLUE}{adapter}{Fore.MAGENTA}{random}{Fore.GREEN}{delimiter}{Fore.LIGHTYELLOW_EX}{polytail}{Style.RESET_ALL}{gene}")
                    print(f"{gene}{Fore.LIGHTYELLOW_EX}{polytail}{Fore.BLUE}{adapter}{Style.RESET_ALL}")
        #input("press enter")


    print(f"{Fore.BLUE}{Count} total reads")
    print(f"{Fore.RED}{CountFor} forward reads")
    print(f"{Fore.RED}{CountForA} indentified forward reads")
    print(f"{Fore.GREEN}{CountRev} reverse reads")
    print(f"{Fore.GREEN}{CountRevA} indentified reverse reads")

    res = pd.DataFrame(
        {'read_core_id': read_core_id_list,
          'mRNA': mRNA_list, 
          'polytail': polytail_list,
          'tail': tail_list,
          'adapter': adapter_list
          } 
    )

    

    print("Adding base composition")
    import time
    start_time = time.time()
    res1 = pd.concat([res, res.apply(lambda row: get_composition(row["polytail"]), axis=1, result_type="expand")], axis = 1)
    res1 = pd.concat([res, res.apply(lambda row: get_composition(row["tail"]), axis=1, result_type="expand")], axis = 1)
    res1.to_csv(out, encoding='utf-8', index=False)
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    res2 = pd.concat([res, res.apply(lambda row: get_composition2(row["polytail"]), axis=1, result_type="expand")], axis = 1)
    res2 = pd.concat([res, res.apply(lambda row: get_composition2(row["tail"]), axis=1, result_type="expand")], axis = 1)
    res2.to_csv(out2, encoding='utf-8', index=False)
    print("--- %s seconds ---" % (time.time() - start_time))


def get_seq(fastq, readname) -> str:
    
    fq=pyfastx.Fastq(fastq)
    seq = fq[readname].seq
    return seq

def get_composition2(seq: str) -> pd.Series:
    """
    HZ
    Get the composition and the percentage of ATGC in a DNA sequence. 
    Ignore all nucleotides other than ATGC (like N).
    """
    nucl_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    nucl_perc = {'%A': 0, '%T': 0, '%G': 0, '%C': 0}
    if seq and not pd.isna(seq):
        for nucl in seq:
            if nucl in nucl_count.keys():
                nucl_count[nucl] += 1
    total_nucl = sum(nucl_count.values())
    if total_nucl == 0:
        return pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])
    else:
        for nucl, count in nucl_count.items():
            pct = count * 100.0 / total_nucl
            nucl_perc["%" + nucl] = round(pct, 2)
    return pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])

def get_composition(seq: str) -> pd.Series:
    """
    JP
    Get the composition and the percentage of ATGC in a DNA sequence. 
    Ignore all nucleotides other than ATGC (like N).
    """
    nucl_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    nucl_perc = {'%A': 0, '%T': 0, '%G': 0, '%C': 0}
    seq = regex.sub('[^ATCG]', '', seq)
    seq_len=len(seq)
    if seq and not pd.isna(seq):
        if seq_len>0:
            res = Counter(seq)
            #print(res)
            nucl_count = {'A': res['A'], 'T': res['T'], 'G': res['G'], 'C': res['C']}
            nucl_perc = {'%A': res['A']/seq_len, '%T': res['T']/seq_len, '%G': res['G']/seq_len, '%C': res['C']/seq_len}
    
            return pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])
        else:
            return pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])




def revcom(seq):
    """
    !!!The function in included in both adapterFinder.py and 
    pacbio_find_polyA.py and extract_read_info.py. 
    They are same function, but haven't be put in a module to 
    keep each script can be run independently. If you want to 
    modify one of them, please modify them at the same time.
    
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
    


if __name__ == '__main__':
    main()






