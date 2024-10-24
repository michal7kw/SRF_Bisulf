#!/bin/bash
#SBATCH --job-name=Bisulf_preparation
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/preparation.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/logs/preparation.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate bisulfite_seq

bismark_genome_preparation --verbose /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/DATA