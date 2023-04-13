# Nano3RACE-Seq

## Description

Snakemake pipeline for polyA detection via Nano3'RACE

This snakemake regroups the major steps described in the FLEPSeq2 github repository (https://github.com/ZhaiLab-SUSTech/FLEPSeq). 
We added an extra step (extract_tails.py) to extract the additionnal tail after the polyA, and analyse its composition. We also added an unique molecule identifier (UMI) sequence in the librairies to allow the deduplication of the PCR duplicates. This UMI consists in a random nucleotide sequence of 10 bp (NNNNNNNNNN)

Some minor changes have been done to the original FLEPSeq2 code:
- Handling directly FASTQ files without needing to convert them to FASTA
- Detecting the UMI to allow deduplication step. The NanoRACE construction contains one unique molecule identifier (UMI) allowing us to identify PCR duplicates. A column will be added in the final result file (dedup_state), and will indicate wether this transcript has the best read sequencing quality ("best") or not ("duplicate"), allowing you to discard duplicates, should you want to.

```mermaid
graph TD
    A(.fast5) --> |Basecalling with Guppy| B(.fastq)
    B --> |Mapping with minimap2| C(.bam)
    C --> |adapterFinder.py| D(.adapter_results.txt)
    D --> |PolyA Caller| E(results)
    D --> |Downstream Analysis|E(results)
    D --> |3' terminal non-adenosine analysis|E(results)
    E --> |plot_tail.R| F(result.merged.parts.csv)

```
## Data Preparation

the fastq files of each sample (genotype) have to be concatenated in the output directory:
```bash
mkdir -p /home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR2/1_Runs

cat barcode12/*.fastq.gz > /home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR2/1_Runs/barcode12.fastq.gz
cat barcode13/*.fastq.gz > /home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR2/1_Runs/barcode13.fastq.gz
cat barcode14/*.fastq.gz > /home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR2/1_Runs/barcode14.fastq.gz
```
A barcode correspondance file (tab separated) is required and looks as follows:

|  |  |
| ----------- | ------------|
| barcode12   | WT          |
| barcode13   | mut1        |
| barcode14   | mut2        |



The first line will be the reference for statistical tests.

## Getting started

Clone the repository

```bash
git clone https://github.com/jackson-peter/NanoRACE.git
```

A configuration file (.yaml) is required to run the workflow. An example is included in the FLEPseq2/config folder. You can modify it to suit your data.

```yaml
### EXPERIMENT SPECIFIC 
basecalled_dir: "/ssd_workspace/jpeter/ssData/Guppy_basecalling/20230131_Nano/workspace"
outdir: "/home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR3"
barcode_corr: "/home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR3/barcode_corr.tsv"
sequencing_summary: "/home/jpeter/DATA/NanoRACE/Test_RACE/OUTDIR3/1_Runs/sequencing_summary.txt"

### OUTPUT FILES ORGANIZATION
runs_dir: "1_Runs"
mapping_dir: "2_Mapping"
polyA_dir: "3_PolyA"
tail_dir: "4_Tail"

### ANNOTATIONS
introns_exons: "/home/jpeter/DATA/NanoRACE/Test_RACE/NanoRACE_RUN03_31012023_exon_intron_posV2.bed"
select_introns: "/home/jpeter/DATA/ReferenceGenomes/Athaliana/TAIR10/select_introns.txt"
reference_genome: "/home/jpeter/DATA/NanoRACE/Test_RACE/NanoRACE_RUN03_31012023_mapping.fa" # Reference genome in fasta

### MINIMAP MAPPING ADDITIONAL PARAMETERS
minimap2_add_opts: "--secondary=no -G 5000"

```

## Usage

To run the whole pipeline, just execute the runFLEPseq.sh file containing the snakemakes commands
```bash
conda activate snakemake # or mamba activate snakemake if it is installed
# Go in the 'workflow' folder where the Snakefile is 
cd FLEPseq2/workflow
bash runFLEPseq.sh
```

## Requirements

The only requirement to run the snakefiles is:
- snakemake (>=7)

The other requirements will be handled directly by snakemake with the --use-conda option. This will automatically set up an environement containing all the requirements for the pipeline's execution.


## Usefull links

- Research article describing FLEPSeq2: https://www.nature.com/articles/s41596-021-00581-7


