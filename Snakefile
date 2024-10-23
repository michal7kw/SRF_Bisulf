# Import necessary libraries
import os
from snakemake.utils import min_version
import pandas as pd

# Set minimum Snakemake version to ensure compatibility
min_version("6.0")

# Load configuration from config.yaml file
configfile: "config.yaml"

# Read sample information from CSV file
samples = pd.read_csv("samplesheet.csv")

if config.get("debug_sample"):
    SAMPLES = [config["debug_sample"]]
else:
    SAMPLES = samples['sample'].tolist()

# Define output directory from config file
OUTPUT_DIR = config["output_dir"]

# Specify the conda environment to use for all rules
conda: "environment.yaml"

# Default resources for all rules
default_resources = {
    "mem_mb": 8000,
    "time": "2:00:00"
}

# Default resources function
def get_default_resources(wildcards, attempt):
    """Default resources for all rules"""
    return {
        "mem_mb": attempt * 8000,
        "time": f"{attempt*2}:00:00",
        "threads": 1
    }

# Default rule resources
rule_resources = {
    "trim_galore": {
        "mem_mb": lambda wildcards, attempt: attempt * 32000,
        "time": lambda wildcards, attempt: f"{attempt*4}:00:00"
    },
    "bismark_align": {
        "mem_mb": lambda wildcards, attempt: attempt * 64000,
        "time": lambda wildcards, attempt: f"{attempt*8}:00:00"
    },
    "deduplicate": {
        "mem_mb": lambda wildcards, attempt: attempt * 32000,
        "time": lambda wildcards, attempt: f"{attempt*4}:00:00"
    },
    "genome_index": {
        "mem_mb": 16000,
        "time": "2:00:00"
    }
}

# Final target rule that defines all expected output files
rule all:
    input:
        # FastQC reports
        expand(os.path.join(OUTPUT_DIR, "fastqc", "{sample}_{read}_fastqc.html"), sample=SAMPLES, read=["R1", "R2"]),
        # Trimmed fastq files
        expand(os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R{read}_val_{read}.fq.gz"), sample=SAMPLES, read=[1,2]),
        # Aligned BAM files
        expand(os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam"), sample=SAMPLES),
        # Deduplicated BAM files
        expand(os.path.join(OUTPUT_DIR, "dedup", "{sample}_bismark_bt2_pe.deduplicated.bam"), sample=SAMPLES),
        # Methylation extraction results
        expand(os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.CpG_report.txt.gz"), sample=SAMPLES),
        # Bismark reports for each sample
        expand(os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html"), sample=SAMPLES),
        # Bismark summary report
        os.path.join(OUTPUT_DIR, "bismark_summary", "bismark_summary_report.html"),
        # Qualimap reports
        expand(os.path.join(OUTPUT_DIR, "qualimap", "{sample}", "qualimapReport.html"), sample=SAMPLES),
        # Preseq complexity estimates
        expand(os.path.join(OUTPUT_DIR, "preseq", "{sample}_complexity.txt"), sample=SAMPLES),
        # MultiQC report
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html"),
        # FastQC on trimmed reads
        expand(os.path.join(OUTPUT_DIR, "fastqc_trimmed", "{sample}_R{read}_val_{read}_fastqc.html"), 
               sample=SAMPLES, read=[1,2]),
        # M-bias plots
        expand(os.path.join(OUTPUT_DIR, "mbias", "{sample}_bismark_bt2_pe.M-bias.txt"), 
               sample=SAMPLES),
        # Coverage statistics
        expand(os.path.join(OUTPUT_DIR, "coverage_stats", "{sample}_coverage_stats.txt"), 
               sample=SAMPLES),
        # Methylation statistics
        expand(os.path.join(OUTPUT_DIR, "methylation_stats", "{sample}_methylation_stats.txt"), 
               sample=SAMPLES)

# Rule to run FastQC on raw reads
rule fastqc_raw:
    input:
        r1 = lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_1'].iloc[0],
        r2 = lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_2'].iloc[0]
    output:
        html = [os.path.join(OUTPUT_DIR, "fastqc", "{sample}_R1_fastqc.html"),
                os.path.join(OUTPUT_DIR, "fastqc", "{sample}_R2_fastqc.html")]
    log:
        os.path.join(OUTPUT_DIR, "logs", "fastqc", "{sample}.log")
    threads: 2
    resources:
        **default_resources
    message:
        "Running FastQC on raw reads for sample {wildcards.sample}"
    shell:
        "fastqc --outdir {OUTPUT_DIR}/fastqc --threads {threads} {input.r1} {input.r2} > {log} 2>&1"

# Rule to trim reads using Trim Galore
# Issue 1: trim_galore rule output naming doesn't match input patterns
rule trim_galore:
    input:
        r1 = lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_1'].iloc[0],
        r2 = lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_2'].iloc[0]
    output:
        r1 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R1_val_1.fq.gz"),
        r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_val_2.fq.gz")
    log:
        os.path.join(OUTPUT_DIR, "logs", "trim_galore", "{sample}.log")
    threads: lambda wildcards, attempt: attempt * 8
    resources:
        **rule_resources["trim_galore"]
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} \
            --quality 20 \
            --phred33 \
            --stringency 3 \
            --length 20 \
            --output_dir $(dirname {output.r1}) \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

rule bismark_genome_preparation:
    input:
        genome = config["genome_dir"]
    output:
        directory(os.path.join(config["genome_dir"], "Bisulfite_Genome"))
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark_genome_preparation.log")
    resources:
        mem_mb = 32000,
        time = "4:00:00"
    message:
        "Preparing genome for Bismark alignment"
    shell:
        """
        bismark_genome_preparation --verbose {input.genome} 2>&1 | tee {log}
        """

# Update bismark_align rule to match trim_galore output files
rule bismark_align:
    input:
        r1 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R1_val_1.fq.gz"),
        r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_val_2.fq.gz"),
        genome_prep = rules.bismark_genome_preparation.output,
        index_check = os.path.join(OUTPUT_DIR, "genome_index_check.done")  # Add this line
    output:
        bam = os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam")
    params:
        genome = config["genome_dir"],
        output_dir = os.path.join(OUTPUT_DIR, "aligned")
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark", "{sample}.log")
    threads: lambda wildcards, attempt: attempt * 8
    resources:
        **rule_resources["bismark_align"]
    message:
        "Aligning reads for sample {wildcards.sample}"
    shell:
        """
        bismark --genome {params.genome} \
            -1 {input.r1} -2 {input.r2} \
            --output_dir {params.output_dir} \
            --basename {wildcards.sample} \
            --multicore {threads} \
            --bowtie2 \
            --temp_dir /tmp/ \
            2>&1 | tee {log}
        """
# --non_directional


# Rule to deduplicate aligned reads
rule deduplicate:
    input:
        bam = rules.bismark_align.output.bam
    output:
        dedup = os.path.join(OUTPUT_DIR, "dedup", "{sample}_bismark_bt2_pe.deduplicated.bam")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "dedup")
    log:
        os.path.join(OUTPUT_DIR, "logs", "deduplicate", "{sample}.log")
    resources:
        **rule_resources["deduplicate"]
    message:
        "Deduplicating reads for sample {wildcards.sample}"
    shell:
        """
        deduplicate_bismark --bam {input.bam} --output_dir {params.output_dir} 2>&1 | tee {log}
        """

# Rule to extract methylation information
rule methylation_extraction:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        bedgraph = os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz"),
        cg_report = os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.CpG_report.txt.gz")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "methylation_extraction")
    log:
        os.path.join(OUTPUT_DIR, "logs", "methylation_extraction", "{sample}.log")
    threads: 4
    resources:
        mem_mb = 32000
    message:
        "Extracting methylation for sample {wildcards.sample}"
    shell:
        """
        bismark_methylation_extractor --bedGraph --gzip --cytosine_report \
            --multicore {threads} \
            --buffer_size 10G \
            --comprehensive \
            --merge_non_CpG \
            --output {params.output_dir} {input.bam} 2>&1 | tee {log}
        """

# Rule to generate Bismark report for each sample
rule bismark_report:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        report = os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "bismark_reports")
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark_report", "{sample}.log")
    resources:
        **default_resources
    message:
        "Generating Bismark report for sample {wildcards.sample}"
    shell:
        """
        bismark2report --alignment_report {input.bam} \
        --output {params.output_dir} 2>&1 | tee {log}
        """

# Rule to generate Bismark summary report for all samples
rule bismark_summary:
    input:
        expand(os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html"), sample=SAMPLES)
    output:
        summary = os.path.join(OUTPUT_DIR, "bismark_summary", "bismark_summary_report.html")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "bismark_summary")
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark_summary", "summary.log")
    resources:
        **default_resources
    message:
        "Generating Bismark summary report"
    shell:
        """
        bismark2summary --output_dir {params.output_dir} 2>&1 | tee {log}
        """

# Rule to perform alignment QC using Qualimap
rule qualimap:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        report = os.path.join(OUTPUT_DIR, "qualimap", "{sample}", "qualimapReport.html")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "qualimap", "{sample}")
    log:
        os.path.join(OUTPUT_DIR, "logs", "qualimap", "{sample}.log")
    threads: 4
    resources:
        mem_mb = 16000
    message:
        "Running Qualimap for sample {wildcards.sample}"
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt {threads} -c -outdir {params.output_dir} \
        --java-mem-size={resources.mem_mb}M 2>&1 | tee {log}
        """

# Rule to estimate library complexity using Preseq
rule preseq:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        complexity = os.path.join(OUTPUT_DIR, "preseq", "{sample}_complexity.txt")
    log:
        os.path.join(OUTPUT_DIR, "logs", "preseq", "{sample}.log")
    resources:
        mem_mb = 8000
    message:
        "Estimating library complexity for sample {wildcards.sample}"
    shell:
        """
        preseq lc_extrap -bam {input.bam} -o {output.complexity} 2>&1 | tee {log}
        """

# Rule to generate a comprehensive MultiQC report
rule multiqc:
    input:
        expand(os.path.join(OUTPUT_DIR, "fastqc", "{sample}_{read}_fastqc.html"), sample=SAMPLES, read=["R1", "R2"]),
        expand(os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html"), sample=SAMPLES),
        os.path.join(OUTPUT_DIR, "bismark_summary", "bismark_summary_report.html"),
        expand(os.path.join(OUTPUT_DIR, "qualimap", "{sample}", "qualimapReport.html"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "preseq", "{sample}_complexity.txt"), sample=SAMPLES)
    output:
        report = os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "multiqc")
    log:
        os.path.join(OUTPUT_DIR, "logs", "multiqc", "multiqc.log")
    resources:
        **default_resources
    message:
        "Generating MultiQC report"
    shell:
        """
        multiqc {OUTPUT_DIR} -o {params.output_dir} 2>&1 | tee {log}
        """

# Add a rule to check input files exist
rule check_inputs:
    input:
        lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_1'].iloc[0],
        lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_2'].iloc[0]
    output:
        touch(os.path.join(OUTPUT_DIR, "checks", "{sample}.done"))
    run:
        for f in input:
            if not os.path.exists(f):
                raise ValueError(f"Input file not found: {f}")

rule fastqc_trimmed:
    input:
        r1 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R1_val_1.fq.gz"),
        r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_val_2.fq.gz")
    output:
        html = [os.path.join(OUTPUT_DIR, "fastqc_trimmed", "{sample}_R1_val_1_fastqc.html"),
                os.path.join(OUTPUT_DIR, "fastqc_trimmed", "{sample}_R2_val_2_fastqc.html")]
    log:
        os.path.join(OUTPUT_DIR, "logs", "fastqc_trimmed", "{sample}.log")
    threads: 2
    resources:
        **default_resources
    message:
        "Running FastQC on trimmed reads for sample {wildcards.sample}"
    shell:
        "fastqc --outdir {OUTPUT_DIR}/fastqc_trimmed --threads {threads} {input.r1} {input.r2} > {log} 2>&1"

rule mbias_plots:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        mbias = os.path.join(OUTPUT_DIR, "mbias", "{sample}_bismark_bt2_pe.M-bias.txt")
    params:
        output_dir = os.path.join(OUTPUT_DIR, "mbias")
    log:
        os.path.join(OUTPUT_DIR, "logs", "mbias", "{sample}.log")
    threads: 4
    resources:
        mem_mb = 16000,
        time = "2:00:00"
    message:
        "Generating M-bias plots for sample {wildcards.sample}"
    shell:
        """
        bismark_methylation_extractor --mbias_only \
            --multicore {threads} \
            --output {params.output_dir} \
            {input.bam} 2>&1 | tee {log}
        """

rule coverage_stats:
    input:
        bedgraph = os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz")
    output:
        coverage = os.path.join(OUTPUT_DIR, "coverage_stats", "{sample}_coverage_stats.txt")
    log:
        os.path.join(OUTPUT_DIR, "logs", "coverage_stats", "{sample}.log")
    resources:
        **default_resources
    message:
        "Calculating coverage statistics for sample {wildcards.sample}"
    shell:
        """
        zcat {input.bedgraph} | awk '{{sum+=$4; sumsq+=$4*$4}} END {{
            mean=sum/NR;
            variance=sumsq/NR - (mean*mean);
            stddev=sqrt(variance);
            print "Total positions:", NR;
            print "Mean coverage:", mean;
            print "Standard deviation:", stddev;
            print "Variance:", variance;
        }}' > {output.coverage} 2> {log}
        """

rule methylation_stats:
    input:
        cg_report = os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.CpG_report.txt.gz")
    output:
        stats = os.path.join(OUTPUT_DIR, "methylation_stats", "{sample}_methylation_stats.txt")
    log:
        os.path.join(OUTPUT_DIR, "logs", "methylation_stats", "{sample}.log")
    resources:
        **default_resources
    message:
        "Calculating methylation statistics for sample {wildcards.sample}"
    shell:
        """
        zcat {input.cg_report} | awk '
            BEGIN {{meth=0; unmeth=0; total=0}}
            {{
                meth+=$4;
                unmeth+=$5;
                total++;
            }}
            END {{
                print "Total CpG sites:", total;
                print "Methylated calls:", meth;
                print "Unmethylated calls:", unmeth;
                print "Global methylation rate:", (meth/(meth+unmeth))*100 "%";
            }}' > {output.stats} 2> {log}
        """

rule check_genome_index:
    input:
        genome_dir = config["genome_dir"]
    output:
        touch(os.path.join(OUTPUT_DIR, "genome_index_check.done"))
    resources:
        **default_resources
    message:
        "Checking Bismark genome index status"
    run:
        index_path = os.path.join(input.genome_dir, "Bisulfite_Genome")
        if not os.path.exists(index_path):
            raise ValueError(f"Bismark genome index not found in {index_path}. Please run bismark_genome_preparation first.")
        
        required_files = ["Bisulfite_Genome/CT_conversion/genome.1.bt2",
                         "Bisulfite_Genome/GA_conversion/genome.1.bt2"]
        
        for f in required_files:
            full_path = os.path.join(input.genome_dir, f)
            if not os.path.exists(full_path):
                raise ValueError(f"Required index file not found: {full_path}")

