# Nano3RACE-Seq

## Description

Snakemake pipeline for polyA detection via:

    - FLEP-seq (Full-Length Elongating and Polyadenylated RNA sequencing) 
    - Nano3RACE

This snakemake regroups the major steps described in the FLEPSeq2 github repository (https://github.com/ZhaiLab-SUSTech/FLEPSeq). 
We added an extra step (extract_tails.py) to extract the polyA tails and to help filter the PCR duplicates.  

Some minor changes have been done to the original FLEPSeq2 code:
- Handling directly FASTQ files without needing to convert them to FASTA

Steps of the workflow
```mermaid
graph TD
    A(.nanopore.fast5) --> |Basecalling with Guppy| B(.fastq)
    B --> |Mapping with minimap2| C(.bam)
    C --> |adapterFinder.py| D(.adapter_results.txt)
    D --> |PolyA Caller| E(results)
    D --> |Downstream Analysis|E(results)
    D --> |3' terminal non-adenosine analysis|E(results)

```
