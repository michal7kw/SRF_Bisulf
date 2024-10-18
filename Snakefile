# Import necessary libraries
import os
from snakemake.utils import min_version
import pandas as pd

# Set minimum Snakemake version to ensure compatibility
min_version("6.0")

# Load configuration from config.yaml file
configfile: "config.yaml"

# Read sample information from CSV file
# This allows for easy management of multiple samples
samples = pd.read_csv("samplesheet.csv")

if config.get("debug_sample"):
    SAMPLES = [config["debug_sample"]]
else:
    SAMPLES = samples['sample'].tolist()

# Define output directory from config file
OUTPUT_DIR = config["output_dir"]

# Specify the conda environment to use for all rules
conda: "environment.yaml"

# Final target rule that defines all expected output files
rule all:
    input:
        # Trimmed fastq files
        expand(os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R{read}_paired_trimmed.fq.gz"), sample=SAMPLES, read=[1,2]),
        # Aligned BAM files
        expand(os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam"), sample=SAMPLES),
        # Deduplicated BAM files
        expand(os.path.join(OUTPUT_DIR, "dedup", "{sample}_bismark_bt2_pe.deduplicated.bam"), sample=SAMPLES),
        # Methylation extraction results
        expand(os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz"), sample=SAMPLES),
        # Bismark reports for each sample
        expand(os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html"), sample=SAMPLES),
        # Bismark summary report
        os.path.join(OUTPUT_DIR, "bismark_summary", "bismark_summary_report.html"),
        # Qualimap reports
        expand(os.path.join(OUTPUT_DIR, "qualimap", "{sample}", "qualimapReport.html"), sample=SAMPLES),
        # Preseq complexity estimates
        expand(os.path.join(OUTPUT_DIR, "preseq", "{sample}_complexity.txt"), sample=SAMPLES),
        # MultiQC report
        os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html")

# Rule to trim reads using Trim Galore
rule trim_galore:
    input:
        r1 = lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_1'].iloc[0],
        r2 = lambda wildcards: samples.loc[samples['sample'] == wildcards.sample, 'fastq_2'].iloc[0]
    output:
        r1 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R1_paired_trimmed.fq.gz"),
        r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_paired_trimmed.fq.gz")
    log:
        os.path.join(OUTPUT_DIR, "logs", "trim_galore", "{sample}.log")
    threads: 32
    resources:
        mem_mb = 32000
    message:
        "Trimming reads for sample {wildcards.sample}"
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} \
            --output_dir $(dirname {output.r1}) \
            {input.r1} {input.r2} \
            > {log} 2>&1

        echo "Trimming completed for {wildcards.sample}"
        """
        # Rename output files to match expected names
        # mv $(dirname {output.r1})/*_val_1.fq.gz {output.r1}
        # mv $(dirname {output.r1})/*_val_2.fq.gz {output.r2}

rule bismark_align:
    input:
        r1 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R1_paired_trimmed.fq.gz"),
        # r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_paired_val_2.fq.gz")
        r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_paired_trimmed.fq.gz")
    output:
        bam = os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam")
    params:
        genome = config["genome_dir"],
        output_dir = os.path.join(OUTPUT_DIR, "aligned")
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark", "{sample}.log")
    message:
        "Aligning reads for sample {wildcards.sample}"
    shell:
        """
        bismark --genome {params.genome} -1 {input.r1} -2 {input.r2} \
        --output_dir {params.output_dir} \
        --basename {wildcards.sample} \
        2>&1 | tee {log}

        # Ensure the output file has the correct name
        if [ -f {params.output_dir}/{wildcards.sample}_pe.bam ]; then
            mv {params.output_dir}/{wildcards.sample}_pe.bam {output.bam}
        elif [ -f {params.output_dir}/{wildcards.sample}_bismark_bt2_pe.bam ]; then
            mv {params.output_dir}/{wildcards.sample}_bismark_bt2_pe.bam {output.bam}
        else
            echo "Expected output BAM file not found" >&2
            exit 1
        fi

        echo "Alignment completed for {wildcards.sample}"
        """

# rule bismark_align:
#     input:
#         r1 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R1_paired_trimmed.fq.gz"),
#         r2 = os.path.join(OUTPUT_DIR, "trimmed", "{sample}_R2_paired_trimmed.fq.gz")
#     output:
#         bam = os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam")
#     params:
#         genome = config["genome_dir"]
#     log:
#         os.path.join(OUTPUT_DIR, "logs", "bismark", "{sample}.log")
#     message:
#         "Aligning reads for sample {wildcards.sample}"
#     shell:
#         """
#         bismark --genome {params.genome} -1 {input.r1} -2 {input.r2} \
#         --output_dir {OUTPUT_DIR}/aligned 2>&1 | tee {log} && \
#         mv {OUTPUT_DIR}/aligned/*_pe.bam {output.bam}
#         echo "Alignment completed for {wildcards.sample}"
#         """

# # Rule to align trimmed reads using Bismark
# rule bismark_align:
#     input:
#         r1 = rules.trim_galore.output.r1,
#         r2 = rules.trim_galore.output.r2
#     output:
#         bam = os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam")
#     params:
#         genome = config["genome_dir"]
#     log:
#         os.path.join(OUTPUT_DIR, "logs", "bismark", "{sample}.log")
#     message:
#         "Aligning reads for sample {wildcards.sample}"
#     shell:
#         """
#         bismark --genome {params.genome} -1 {input.r1} -2 {input.r2} \
#         --output_dir {OUTPUT_DIR}/aligned 2>&1 | tee {log} && \
#         mv {OUTPUT_DIR}/aligned/*_pe.bam {output.bam}
#         echo "Alignment completed for {wildcards.sample}"
#         """

# Rule to deduplicate aligned reads
rule deduplicate:
    input:
        bam = rules.bismark_align.output.bam
    output:
        dedup = os.path.join(OUTPUT_DIR, "dedup", "{sample}_bismark_bt2_pe.deduplicated.bam")
    log:
        os.path.join(OUTPUT_DIR, "logs", "deduplicate", "{sample}.log")
    message:
        "Deduplicating reads for sample {wildcards.sample}"
    shell:
        """
        deduplicate_bismark --bam {input.bam} --output_dir {OUTPUT_DIR}/dedup 2>&1 | tee {log}
        echo "Deduplication completed for {wildcards.sample}"
        """

# Rule to extract methylation information
rule methylation_extraction:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        bedgraph = os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz")
    log:
        os.path.join(OUTPUT_DIR, "logs", "methylation_extraction", "{sample}.log")
    message:
        "Extracting methylation for sample {wildcards.sample}"
    shell:
        """
        bismark_methylation_extractor --bedGraph --gzip \
        --output {OUTPUT_DIR}/methylation_extraction {input.bam} 2>&1 | tee {log}
        echo "Methylation extraction completed for {wildcards.sample}"
        """

# Rule to generate Bismark report for each sample
rule bismark_report:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        report = os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html")
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark_report", "{sample}.log")
    message:
        "Generating Bismark report for sample {wildcards.sample}"
    shell:
        """
        bismark2report --alignment_report {input.bam} \
        --output {OUTPUT_DIR}/bismark_reports 2>&1 | tee {log}
        echo "Bismark report generated for {wildcards.sample}"
        """

# Rule to generate Bismark summary report for all samples
rule bismark_summary:
    input:
        expand(os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html"), sample=SAMPLES)
    output:
        summary = os.path.join(OUTPUT_DIR, "bismark_summary", "bismark_summary_report.html")
    log:
        os.path.join(OUTPUT_DIR, "logs", "bismark_summary", "summary.log")
    message:
        "Generating Bismark summary report"
    shell:
        """
        bismark2summary --output_dir {OUTPUT_DIR}/bismark_summary 2>&1 | tee {log}
        echo "Bismark summary report generated"
        """

# Rule to perform alignment QC using Qualimap
rule qualimap:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        report = os.path.join(OUTPUT_DIR, "qualimap", "{sample}", "qualimapReport.html")
    log:
        os.path.join(OUTPUT_DIR, "logs", "qualimap", "{sample}.log")
    message:
        "Running Qualimap for sample {wildcards.sample}"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {OUTPUT_DIR}/qualimap/{wildcards.sample} 2>&1 | tee {log}
        echo "Qualimap analysis completed for {wildcards.sample}"
        """

# Rule to estimate library complexity using Preseq
rule preseq:
    input:
        bam = rules.deduplicate.output.dedup
    output:
        complexity = os.path.join(OUTPUT_DIR, "preseq", "{sample}_complexity.txt")
    log:
        os.path.join(OUTPUT_DIR, "logs", "preseq", "{sample}.log")
    message:
        "Estimating library complexity for sample {wildcards.sample}"
    shell:
        """
        preseq lc_extrap -bam {input.bam} -o {output.complexity} 2>&1 | tee {log}
        echo "Library complexity estimation completed for {wildcards.sample}"
        """

# Rule to generate a comprehensive MultiQC report
rule multiqc:
    input:
        expand(os.path.join(OUTPUT_DIR, "aligned", "{sample}_bismark_bt2_pe.bam"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "methylation_extraction", "{sample}_bismark_bt2_pe.bedGraph.gz"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "bismark_reports", "{sample}_bismark_bt2_PE_report.html"), sample=SAMPLES),
        os.path.join(OUTPUT_DIR, "bismark_summary", "bismark_summary_report.html"),
        expand(os.path.join(OUTPUT_DIR, "qualimap", "{sample}", "qualimapReport.html"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "preseq", "{sample}_complexity.txt"), sample=SAMPLES)
    output:
        report = os.path.join(OUTPUT_DIR, "multiqc", "multiqc_report.html")
    log:
        os.path.join(OUTPUT_DIR, "logs", "multiqc", "multiqc.log")
    message:
        "Generating MultiQC report"
    shell:
        """
        multiqc {OUTPUT_DIR} -o {OUTPUT_DIR}/multiqc 2>&1 | tee {log}
        echo "MultiQC report generated"
        """
