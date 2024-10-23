#!/bin/bash
#SBATCH --job-name=Bisulf
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/bisulf.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/bisulf.out"

# # Load the appropriate conda environment
# source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# # Create and activate the environment
# conda env create -f environment.yaml
# conda activate bisulfite_seq

# Load the conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate bisulfite_seq

# Check if the genome file exists
GENOME_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/DATA/genome.fa"
if [ ! -f "$GENOME_FILE" ]; then
    echo "Error: Genome file not found at $GENOME_FILE"
    exit 1
fi

# Unlock the working directory
snakemake --unlock

# Run Snakemake
snakemake \
    --snakefile Snakefile_mod \
    --latency-wait 60 \
    --configfile config.yaml \
    --config debug_sample=NEU1 \
    -p \
    --printshellcmds \
    --verbose \
    --rerun-incomplete \
    --use-conda \
    --cores all \
    --jobs 100 \
    --max-jobs-per-second 1 \
    --keep-going \
    --restart-times 3

# snakemake \
#     --snakefile Snakefile \
#     --latency-wait 60 \
#     --configfile config.yaml \
#     --config debug_sample=NEU1 \
#     -p \
#     --printshellcmds \
#     --verbose \
#     --rerun-incomplete \
#     --use-conda \
#     --cores 32
    # --cluster "sbatch --parsable --job-name=Bisulf_{rule} --mem={resources.mem_mb} --time={resources.time} --cpus-per-task={threads}"

