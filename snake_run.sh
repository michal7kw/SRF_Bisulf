#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/log.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/log.out"

# Load the appropriate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Create and activate the environment
conda env create -f environment.yaml
conda activate bisulfite_seq

# Run Snakemake for all samples
# snakemake --unlock
# snakemake \
#     --snakefile Snakefile \
#     --latency-wait 60 \
#     --configfile config.yaml \
#     -p \
#     --printshellcmds \
#     --keep-going \
#     --verbose \
#     --rerun-incomplete \
#     --cores all

# Run Snakemake for a single sample (e.g., NEU1)
snakemake --unlock
snakemake \
    --snakefile Snakefile \
    --latency-wait 60 \
    --configfile config.yaml \
    --config debug_sample=NEU1 \
    -p \
    --printshellcmds \
    --verbose \
    --rerun-incomplete \
    --cores all

#    --keep-going \

# Dry run Snakemake for a single sample (e.g., NEU1)
# snakemake \
#     --snakefile Snakefile \
#     --configfile config.yaml \
#     --config debug_sample=NEU1 \
#     -np \
#     --rerun-incomplete \
#     --cores all

# snakemake --unlock
# snakemake \
#     --snakefile Snakefile \
#     --latency-wait 60 \
#     --configfile config.yaml \
#     --config debug_sample=NEU1 \
#     -p \
#     --printshellcmds \
#     --keep-going \
#     --verbose \
#     --rerun-incomplete \
#     --use-conda \
#     --cores all

# snakemake --unlock
# snakemake \
#     --snakefile Snakefile \
#     --latency-wait 60 \
#     --configfile config.yaml \
#     --config debug_sample=NEU1 \
#     -p \
#     --printshellcmds \
#     --keep-going \
#     --verbose \
#     --rerun-incomplete \
#     --use-conda \
#     --cores all \
#     --until multiqc