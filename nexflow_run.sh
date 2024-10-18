#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/log.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/log.out"

# Load the appropriate conda environment (if needed)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Change to the working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf

# Run Nextflow
nextflow run nf-core/methylseq \
    --input /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/samplesheet.csv \
    --outdir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/output \
    --genome GRCm38 \
    -profile conda
