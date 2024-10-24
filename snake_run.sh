#!/bin/bash
#SBATCH --job-name=Bisulf
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
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

# Unlock the working directory
snakemake --unlock

# Run Snakemake with default resource handling
snakemake \
    --snakefile Snakefile \
    --latency-wait 60 \
    --configfile config.yaml \
    --config debug_sample=NEU1 \
    -p \
    --printshellcmds \
    --verbose \
    --use-conda \
    --cores all \
    --jobs 100 \
    --max-jobs-per-second 1 \
    --keep-going \
#   --rerun-incomplete \

############### Test commands ###############
# # Create necessary directories
# mkdir -p logs

# # Create a log file with timestamp
# LOG_FILE="logs/snakemake_test_$(date +%Y%m%d_%H%M%S).log"

# # Function to log commands with timestamps
# log_command() {
#     echo -e "\n=== $(date '+%Y-%m-%d %H:%M:%S') === Running: $1 ===" | tee -a "$LOG_FILE"
#     eval "$1" 2>&1 | tee -a "$LOG_FILE"
#     echo -e "\n=== $(date '+%Y-%m-%d %H:%M:%S') === Finished: $1 ===\n" | tee -a "$LOG_FILE"
# }

# # Function to check prerequisites
# check_prerequisites() {
#     # Check if genome file exists
#     GENOME_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Bisulf/DATA/genome.fa"
#     if [ ! -f "$GENOME_FILE" ]; then
#         echo "Error: Genome file not found at $GENOME_FILE" | tee -a "$LOG_FILE"
#         exit 1
#     fi

#     # Check if conda environment is active
#     if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "bisulfite_seq" ]; then
#         echo "Error: bisulfite_seq conda environment is not active" | tee -a "$LOG_FILE"
#         exit 1
#     fi
# }

# # Run checks
# echo -e "\n=== $(date '+%Y-%m-%d %H:%M:%S') === Running prerequisite checks ===" | tee -a "$LOG_FILE"
# check_prerequisites


# # Run each command with logging
# echo -e "\n=== STARTING SNAKEMAKE TESTS ===" | tee -a "$LOG_FILE"
# echo "Test started at: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG_FILE"

# # 1. Check syntax and best practices
# log_command "snakemake --lint"

# # 2. Create conda environment if needed
# log_command "snakemake --use-conda --conda-create-envs-only"

# # 3. Do a dry run with detailed output
# log_command "snakemake -n -p --config debug_sample=NEU1"

# # 4. Generate workflow visualization
# log_command "snakemake --dag | dot -Tpng > workflow_$(date +%Y%m%d_%H%M%S).png"

# # 5. Check input files
# log_command "snakemake --list-input-changes"

# # 6. Test with touch (creates empty output files without running commands)
# log_command "snakemake -n --touch"

# # Optional: Test specific rules
# # log_command "snakemake -n -p --until fastqc_raw trim_galore --config debug_sample=NEU1"

# # Create a summary of the test results
# echo -e "\n=== TEST SUMMARY ===" | tee -a "$LOG_FILE"
# echo "Test completed at: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG_FILE"
# echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"

############### End of test commands ###############
# snakemake -n
# snakemake -n -p
# snakemake --dag | dot -Tpng > dag.png
# snakemake --lint
# snakemake --list-input-changes
# snakemake --use-conda --conda-create-envs-only