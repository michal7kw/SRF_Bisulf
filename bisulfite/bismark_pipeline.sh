#!/bin/bash
#SBATCH --job-name=metWGBSNPC2
#SBATCH --account=moscato.giosue
#SBATCH --mem=100GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=30
#SBATCH --nodes=1 # not really useful for not mpi jobs
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=moscato.giosue@hsr.it
#SBATCH --error="/beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/metWGBSNPC2.err"
#SBATCH --output="/beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/metWGBSNPC2.out"

echo "my job strart now" > myjob.log;
date >> myjob.log;

mkdir -p /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark
cd /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bismark

#bismark /beegfs/scratch/ric.cosr/moscato.giosue/Genomes/mm10/bismark -1 /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/NPC2_R1_paired.fastq.gz -2 /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/NPC2_R2_paired.fastq.gz -o /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark -p 8  --multicore 8

### Deduplication
#deduplicate_bismark --bam [options] <filenames>
#deduplicate_bismark --bam /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark/NPC1_R1_paired_bismark_bt2_pe.bam

#### Methylation extraction
#bismark_methylation_extractor [options] <filenames>
#A typical command to extract context-dependent (CpG/CHG/CHH) methylation could look like this:
#bismark_methylation_extractor --cytosine_report --gzip --bedGraph  --genome_folder /beegfs/scratch/ric.cosr/moscato.giosue/Genomes/mm10/bismark --paired-end --multicore 8 /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark/NPC2_R1_paired_bismark_bt2_pe.deduplicated.bam 


#### Sample report
#mkdir -p /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark/report
bismark2report #--dir /beegfs/scratch/ric.cosr/moscato.giosue/MIRKO_WGBS2/fastq/trimmed/trimmomatic/mapped/bismark/report 


#### Summary report
bismark2summary 


date >> myjob.log;
echo "all done!!" >> myjob.log