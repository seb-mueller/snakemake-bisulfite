import pandas as pd
import os

# Usage:
# source activate srna_mapping
# snakemake --use-conda --conda-prefix ~/.myconda -n

# to make shell rule to work we need to determine the base path of Snakefile since we expect
# the scripts directory there as well
SRCDIR = srcdir("")

configfile: "config.yaml"
missmatches =  config['MAPPING']['missmatches']
REFERENCE   =  config['MAPPING']['reference']
refbase     =  os.path.basename(REFERENCE)
mode        =  config['MAPPING']['mode']

CONTEXT = ['CpG','CHG','CHH']
PAIRS = ['R1', 'R2']

# Get sample names from samples.csv
SAMPLES = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
# SAMPLES = glob_wildcards("data/{S}_R1.fastq.gz").S
# print("My samples:")
print(SAMPLES)
print(CONTEXT)
print(SAMPLES)
# threads     =  config['THREADS']
#software requirements
#fastqc
#trim_galore (cutadapt)
#samtools
#bismark
#bowtie1


#bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# ruleorder: fastqc_raw > trim > alignment > methylation_extractor > bismark2bed

rule all:
    input:
        expand("trimmed/{sample}_{pair}_trimmed.fq.gz", sample=SAMPLES.index, pair = PAIRS),
        expand("logs/fastqc/raw/{sample}_{pair}_fastqc.html", sample=SAMPLES.index, pair = PAIRS),
        # expand("logs/fastqc/raw/{sample}_R2_fastqc.zip", sample=SAMPLES.index),

        # expand("mapped/{sample}/{sample}_trimmed_bismark.bam", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.bam", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam.bai",sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam",sample=SAMPLES),

        # expand("mapped/{sample}/{sample}_trimmed_bismark.nucleotide_stats.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplication_report.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark_SE_report.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.M-bias.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted_splitting_report.txt", sample=SAMPLES),

        # expand("mapped/{sample}/{context}_context_{sample}_trimmed_bismark.deduplicated.sorted.txt.gz", sample=SAMPLES, context=CONTEXT),

        # expand("mapped/{sample}/{sample}_{context}.gz.bismark.cov.gz", sample=SAMPLES,context=CONTEXT),
        # expand("mapped/{sample}/{sample}_{context}.gz", sample=SAMPLES,context=CONTEXT),

rule fastqc_raw:
    """Create fastqc report"""
    input:
        read="data/{sample}_{pair}.fastq.gz",
        # read1="data/{sample}_{pair}.fastq.gz",
        # read2="data/{sample}_R2.fastq.gz",
    output:
        qual="logs/fastqc/raw/{sample}_{pair}_fastqc.html",
        zip ="logs/fastqc/raw/{sample}_{pair}_fastqc.zip",
        # qual2="logs/fastqc/raw/{sample}_R2_fastqc.html",
        # zip1= "logs/fastqc/raw/{sample}_R1_fastqc.zip",
        # zip2= "logs/fastqc/raw/{sample}_R2_fastqc.zip",
    # log: "logs"
    shell:
        """
        fastqc {input.read} -f fastq --outdir logs/fastqc/raw
        """

# https://cutadapt.readthedocs.io/en/stable/guide.html#bisulfite-sequencing-rrbs
rule trim:
    input:
        read1="data/{sample}_" + PAIRS[0] + ".fastq.gz",
        read2="data/{sample}_" + PAIRS[1] + ".fastq.gz",
    # params:
    #     " -a " +     config['FILTER']['cutadapt']['adapter'] +
    #     " -q " + str(config['FILTER']['cutadapt']['quality-filter']) +
    #     " -m " + str(config['FILTER']['cutadapt']['minimum-length']) +
    #     " -M " + str(config['FILTER']['cutadapt']['maximum-length'])
    log:
        "logs/cutadapt/{sample}_{pair}.log"
    output:
        read="trimmed/{sample}_{pair}_trimmed.fq.gz",
        # l1="trimmed/{sample}.fastq.gz_trimming_report.txt",
    shell:
        """
        trim_galore --paired --trim1 {input.read1} {input.read2} --output_dir trimmed
        """

# trim_galore [options] <filename(s)>
# -o/--output_dir <DIR>   If specified all output will be written to this directory instead of the current directory.
    # """--length <INT>"""
# --paired                This option performs length trimming of quality/adapter/RRBS trimmed reads for
# -t/--trim1              Trims 1 bp off every read from its 3' end. This may be needed for FastQ files that
# trim_galore --paired bsseq_sample1_R1.fastq.gz bsseq_sample1_R2.fastq.gz
# results in:
# bsseq_sample1_R1_val_1.fq.gz
# bsseq_sample1_R2_val_2.fq.gz

# rule trim_galor:
#     """--gzip"""
#     input:
#         read1="data/{sample}_{PAIRS[0]}.fastq.gz",
#         read2="data/{sample}_{PAIRS[1]}.fastq.gz",
#     output:
#         read="trimmed/{sample}_R1_trimmed.fq.gz",
#         l1="trimmed/{sample}.fastq.gz_trimming_report.txt",
#     shell:
#         """
#         trim_galore {input.read} --output_dir cleaned
#         """

rule bismark:
    input:
        read="trimmed/{sample}_trimmed.fq.gz",
    output:
        "mapped/{sample}/{sample}_trimmed_bismark.bam",
        se="mapped/{sample}/{sample}_trimmed_bismark_SE_report.txt",
        nuc="mapped/{sample}/{sample}_trimmed_bismark.nucleotide_stats.txt",
    params:
        ref=REFERENCE,
    shell:
        """
        bismark --bowtie1 --multicore 2 -n 1 -l 28 --gzip --nucleotide_coverage {params.ref} {input.read} --output_dir mapped/{wildcards.sample}
        """

rule deduplicate:
    input:
        bam="mapped/{sample}/{sample}_trimmed_bismark.bam",
    output:
        bam="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.bam",
        dup="mapped/{sample}/{sample}_trimmed_bismark.deduplication_report.txt",
    shell:
        """
        deduplicate_bismark --bam {input.bam}
        """

rule sort:
    input:
        bam="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.bam",
    output:
        sort="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam",
        bai="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam.bai",
    shell:
        """
        samtools sort {input.bam} > {output.sort}
        samtools index {output.sort}
        """

rule methylation_extractor:
    input:
        bam="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam",
    output:
        split="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted_splitting_report.txt",
        mbias="mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.M-bias.txt",
        cg="mapped/{sample}/CpG_context_{sample}_trimmed_bismark.deduplicated.sorted.txt.gz",
        chg="mapped/{sample}/CHG_context_{sample}_trimmed_bismark.deduplicated.sorted.txt.gz",
        chh="mapped/{sample}/CHH_context_{sample}_trimmed_bismark.deduplicated.sorted.txt.gz",
    params:
        ref=REFERENCE,
    shell:
        """
        bismark_methylation_extractor --gzip --multicore 4 --comprehensive --report --genome_folder {params.ref} --buffer_size 8G -s {input.bam} --output mapped/{wildcards.sample}
        """

rule bismark2bed:
    input:
        extract="mapped/{sample}/{context}_context_{sample}_trimmed_bismark.deduplicated.sorted.txt.gz",
    output:
        bed="mapped/{sample}/{sample}_{context}.gz",
        cgcov="mapped/{sample}/{sample}_{context}.gz.bismark.cov.gz",
    shell:
        """
        bismark2bedGraph --CX {input.extract} -o {wildcards.sample}_{wildcards.context} --dir mapped/{wildcards.sample}
        """
