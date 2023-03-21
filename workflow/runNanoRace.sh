
snakemake --profile slurm --use-conda --configfile ../config/config_Nano.yaml -s Snakefile_NanoRace
snakemake --profile slurm --use-conda --configfile ../config/config_Nano.yaml -s SnakefileDownstream

