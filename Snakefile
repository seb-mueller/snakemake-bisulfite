import pandas as pd
import os

# Usage:
# source activate srna_mapping
# snakemake --use-conda --conda-prefix ~/.myconda -n

# to make shell rule to work we need to determine the base path of Snakefile since we expect
# the scripts directory there as well
# TODO
# -bismark crushes if -p 1. need to capture if no --cores or --cores 1
SRCDIR = srcdir("")


configfile: "config.yaml"
# missmatches =  config['MAPPING']['missmatches']
REFERENCE   =  config['MAPPING']['reference']
refbase     =  config['MAPPING']['reference_short']
# mode        =  config['MAPPING']['mode']

CONTEXT = ['CpG','CHG','CHH']
PAIRS = ['R1', 'R2']

# Get sample names from samples.csv
SAMPLES = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
# SAMPLES = glob_wildcards("data/{S}_R1.fastq.gz").S
# print("My samples:")
print(REFERENCE)
print("refbase")
print(refbase)
print(CONTEXT)
print(SAMPLES)
print(PAIRS)
#software requirements
#fastqc
#trim_galore (cutadapt)
#samtools
#bismark
#bowtie1


#bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# ruleorder: fastqc_raw > trim > bismark >deduplicate > sort > methylation_extractor > bismark2bedGraph
# to make shell rule to work we need to determine the base path of Snakefile since we expect

rule all:
    input:
        # trimming
        expand("trimmed/{sample}_{pair}_trim.fq.gz", sample=SAMPLES.index, pair = PAIRS),
        # QC
        expand("logs/fastqc/raw/{sample}_{pair}_fastqc.html", sample=SAMPLES.index, pair = PAIRS),
        # bismark mapping
        expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.bam", sample=SAMPLES.index),
        # deduplication
        expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.bam", sample=SAMPLES.index),
        # sorting
        expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.bam", sample=SAMPLES.index),
        expand("mapped/bws/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.bw", sample=SAMPLES.index),
        # "mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bw",
        # extracting methyl for all 3 contexts
        expand("methylation_extracted/CHG_context_{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.txt.gz", sample=SAMPLES.index),
        expand("methylation_extracted/CpG_context_{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.txt.gz", sample=SAMPLES.index),
        expand("methylation_extracted/CHH_context_{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.txt.gz", sample=SAMPLES.index),
        # CHH="mapped/CHH_context_{sample}_MappedOn_{refbase}_trimmed_bismark_pe.deduplicated.sorted.bam",
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam.bai",sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.bam",sample=SAMPLES),

        # expand("mapped/{sample}/{sample}_trimmed_bismark.nucleotide_stats.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplication_report.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark_SE_report.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted.M-bias.txt", sample=SAMPLES),
        # expand("mapped/{sample}/{sample}_trimmed_bismark.deduplicated.sorted_splitting_report.txt", sample=SAMPLES),

        # expand("mapped/{sample}/{context}_context_{sample}_trimmed_bismark.deduplicated.sorted.txt.gz", sample=SAMPLES, context=CONTEXT),

        # expand("methylation_extracted/{sample}_{context}.gz.bismark.cov.gz", sample=SAMPLES,context=CONTEXT),
        # bedGraph="methylation_extracted/{sample}_{context}.gz",
        # cov=     "mapped/{sample}/{sample}_{context}.gz.bismark.cov.gz",

        expand("coverage/{sample}_MappedOn_" + refbase + "_{context}", sample=SAMPLES.index, context=CONTEXT),
        expand("coverage/{sample}_MappedOn_" + refbase + "_{context}.gz.bismark.cov", sample=SAMPLES.index, context=CONTEXT),
        expand("coverage/{sample}_MappedOn_" + refbase + "_{context}.CX_report.txt", sample=SAMPLES.index, context=CONTEXT),
        expand("coverage/{sample}_MappedOn_" + refbase + "_{context}.bw", sample=SAMPLES.index, context=CONTEXT),

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
        # trim_galore alway seems to output fq indpendent of fq|fastq input
        read="trimmed/{sample}_{pair}_trim.fq.gz",
        # l1="trimmed/{sample}.fastq.gz_trimming_report.txt",
    shell:
        """
        trim_galore --paired --trim1 {input.read1} {input.read2} --output_dir trimmed
        rename 's/val_[12]/trim/g' trimmed/*
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

rule bismark:
    input:
        # read="trimmed/{sample}_trimmed.fq.gz",
        read1="trimmed/{sample}_" + PAIRS[0] + "_trim.fq.gz",
        read2="trimmed/{sample}_" + PAIRS[1] + "_trim.fq.gz",
    output:
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.bam",
        # se="mapped/{sample}/{sample}_trimmed_bismark_SE_report.txt",
        # nuc="mapped/{sample}/{sample}_trimmed_bismark.nucleotide_stats.txt",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}.log",
    params:
        ref=REFERENCE,
    threads: 4
    shell:
        # USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}
        """
        bismark --bowtie2 -p {threads} --nucleotide_coverage {params.ref} -1 {input.read1} -2 {input.read2} --basename {wildcards.sample}_MappedOn_{refbase}_trim_bismark --output_dir mapped 2> {log}
        """
        # mapped/{wildcards.sample}


rule deduplicate:
    input:
        # bam="mapped/{sample}/{sample}_trim_bismark.bam",
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.bam",
    output:
        # bam="mapped/{sample}/{sample}_trim_bismark.deduplicated.bam",
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.bam",
        # dup="mapped/{sample}/{sample}_trim_bismark.deduplication_report.txt",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}.deduplication.log",
    shell:
        """
        deduplicate_bismark --paired --bam {input.bam} --output_dir mapped 2> {log}
        """

rule sort:
    input:
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.bam",
    output:
        sort="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam",
        index="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam.bai",
    shell:
        """
        samtools sort {input.bam} > {output.sort}
        samtools index {output.sort}
        """

rule calccoverage:
    """compute coverage using deeptools into bigWig(bw) file"""
    input:
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam",
        bai="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam.bai",
    output:
        "mapped/bws/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bw",
    params:
        binsize="50"
    threads: 8
    shell:
        """
        bamCoverage -b {input.bam} -o {output}  --normalizeUsing CPM --binSize {params.binsize} -p {threads}
        """

#%_fusedICv2_all.cov2cyt: %_fusedICv2_pe.bam
#	#not doing --bedgraph option since the filenames are messed up, doing this manualy
#	#p stands for paired end
#	${bismark_dir}/bismark_methylation_extractor -p --no_overlap --ample_memory --ignore_r2 2 --report --multicore $(threads) $?
#	${bismark_dir}/bismark_methylation_extractor -s --ample_memory --report --multicore $(threads) $*_fusedICv2_unmapped_reads.bam

#	$(eval files_CHG = CHG_OB_$*_fusedICv2_pe.txt CHG_OT_$*_fusedICv2_pe.txt CHG_OB_$*_fusedICv2_unmapped_reads.txt CHG_OT_$*_fusedICv2_unmapped_reads.txt)
#	$(eval files_CHH = CHH_OB_$*_fusedICv2_pe.txt CHH_OT_$*_fusedICv2_pe.txt CHH_OB_$*_fusedICv2_unmapped_reads.txt CHH_OT_$*_fusedICv2_unmapped_reads.txt)
#	$(eval files_CpG = CpG_OB_$*_fusedICv2_pe.txt CpG_OT_$*_fusedICv2_pe.txt CpG_OB_$*_fusedICv2_unmapped_reads.txt CpG_OT_$*_fusedICv2_unmapped_reads.txt)
#	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_CHG.bedgraph $(files_CHG)
#	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_CHH.bedgraph $(files_CHH)
#	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_CpG.bedgraph $(files_CpG)
#	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_all.bedgraph $(files_CpG) $(files_CHH) $(files_CHG)

rule methylation_extractor:
    # USAGE: bismark_methylation_extractor [options] <filenames>
    input:
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam",
    output:
        CHH="methylation_extracted/CHH_context_{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.txt.gz",
        CHG="methylation_extracted/CHG_context_{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.txt.gz",
        CpG="methylation_extracted/CpG_context_{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.txt.gz",
        # "methylation_extracted/{CONTEXT}_context_{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.txt.gz",
        # split="mapped/{sample}/{sample}_trim_bismark.deduplicated.sorted_splitting_report.txt",
        # mbias="mapped/{sample}/{sample}_trim_bismark.deduplicated.sorted.M-bias.txt",
        # cg="mapped/{sample}/CpG_context_{sample}_trim_bismark.deduplicated.sorted.txt.gz",
        # chg="mapped/{sample}/CHG_context_{sample}_trim_bismark.deduplicated.sorted.txt.gz",
        # chh="mapped/{sample}/CHH_context_{sample}_trim_bismark.deduplicated.sorted.txt.gz",
    params:
        ref=REFERENCE,
    threads: 4
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}.methylation_extractor.log",
    shell:
        """
        bismark_methylation_extractor --gzip --paired-end --multicore {threads} --comprehensive --report --genome_folder {params.ref} --buffer_size 8G -s {input.bam} --output methylation_extracted 2> {log}
        """

rule bismark2bedGraph:
    # bedGraph file as well as a coverage file which are both sorted by chromosomal position.
    input:
        # extract="mapped/{sample}/{context}_context_{sample}_trim_bismark.deduplicated.sorted.txt.gz",
        methylex="methylation_extracted/{context}_context_{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.txt.gz",
    output:
        bedGraph="coverage/{sample}_MappedOn_{refbase}_{context}",
        cov=     "coverage/{sample}_MappedOn_{refbase}_{context}.gz.bismark.cov",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.bismark2bedGraph.log",
    shell:
        # gzip outfiles even without --gzip option. Need to unzip manual for later bigwiggle rule
        """
        bismark2bedGraph --CX {input.methylex} -o {wildcards.sample}_MappedOn_{refbase}_{wildcards.context} --dir coverage 2> {log}
        gunzip {output.bedGraph}.gz
        gunzip {output.cov}.gz
        """

rule coverage2cytosine:
    # generates a cytosine methylation report for a genome of interest and a sorted methylation input file produced
    # by the script bismark2bedGraph or the bismark_methylation_extractor when '--bedGraph' was specified. The input files
    # (.cov or .cov.gz) are expected to use 1-based genomic coordinates. By default, the output uses 1-based chromosomal coordinates
    # and reports CpG positions only (for both strands individually and not merged in any way).
    # The input file needs to have been generated with the script bismark2bedGraph (the file is called *.cov, or .cov.gz) or
    # otherwise be sorted by position and exactly in the following format:
    # USAGE: coverage2cytosine [options] --genome_folder <path> -o <output> [input]
    input:
        bedGraph="coverage/{sample}_MappedOn_{refbase}_{context}",
    output:
        cov2cyt="coverage/{sample}_MappedOn_{refbase}_{context}.CX_report.txt",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.coverage2cytosine.log",
    params:
        ref=REFERENCE,
    shell:
        # genome_folder needs to contain fa rather than fasta!
        """
        coverage2cytosine --CX --genome_folder {params.ref} -o {input.bedGraph} --dir . {input.bedGraph} 2> {log}
        """

rule bedGraphToBigWig:
    # usage:
    #    bedGraphToBigWig in.bedGraph chrom.sizes out.bw
    # where in.bedGraph is a four column file in the format:
    #       <chrom> <start> <end> <value>
    # and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>
    # and out.bw is the output indexed big wig file.
    input:
        bedGraph="coverage/{sample}_MappedOn_{refbase}_{context}",
    output:
        bw="coverage/{sample}_MappedOn_{refbase}_{context}.bw",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.bedGraphToBigWig.log",
    params:
        ref=REFERENCE,
    shell:
        # fai files needs to be present to serve as chrom.sizes and can be created by:
        # samtools faidx input.fa
        # cut -f1,2 input.fa.fai > sizes.genome
        """
        bedGraphToBigWig {input.bedGraph} {params.ref}/{refbase}.fa.fai {input.bedGraph}.bw 2> {log}
        """


	# ${bismark_dir}coverage2cytosine --genome_folder ${index_bismarck} --CX -o $*_fusedICv2_all.cov2cyt $*_fusedICv2_all.bedgraph.gz.bismark.cov.gz

	# gunzip $*_fusedICv2_CHG.bedgraph.gz $*_fusedICv2_CpG.bedgraph.gz $*_fusedICv2_CHH.bedgraph.gz
	# ${bedGraphToBigWig} $*_fusedICv2_CpG.bedgraph  ${chromsizes} $*_fusedICv2_CpG.bw
	# ${bedGraphToBigWig} $*_fusedICv2_CHH.bedgraph  ${chromsizes} $*_fusedICv2_CHH.bw
	# ${bedGraphToBigWig} $*_fusedICv2_CHG.bedgraph  ${chromsizes} $*_fusedICv2_CHG.bw
	# rm $*_fusedICv2_CHG.bedgraph $*_fusedICv2_CpG.bedgraph $*_fusedICv2_CHH.bedgraph
