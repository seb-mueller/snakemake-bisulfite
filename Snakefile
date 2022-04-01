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
REFERENCE                      = config['MAPPING']['reference']
refbase                        = config['MAPPING']['reference_short']
extra_params_bismark           = config['MAPPING']['extra_params_bismark']
extra_params_trim              = config['FILTER']['extra_params_trim']
BINSIZE                        = config['REPORT']['binsize']
extra_params_bamCoverage       = config['REPORT']['extra_params_bamCoverage']
extra_params_meth_extractor    = config['REPORT']['extra_params_meth_extractor']
extra_params_bismark2bedGraph  = config['REPORT']['extra_params_bismark2bedGraph']
extra_params_coverage2cytosine = config['REPORT']['extra_params_coverage2cytosine']

# input file parameters
r1_suffix = config["Fastq"]["suffix_R1"]
r2_suffix = config["Fastq"]["suffix_R2"]
fastq_ext = config["Fastq"]["file_extension"]


CONTEXT = ['CpG','CHG','CHH']
PAIRS = [r1_suffix, r2_suffix]

# Get sample names from samples.csv
SAMPLES = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
# an alternative approach would be globing for filenames (but I like specifying files more):
# SAMPLES = glob_wildcards("data/{S}_R1.fastq.gz").S
print(REFERENCE)
print("refbase")
print(refbase)
print(CONTEXT)
print(SAMPLES)
print(PAIRS)

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
        # expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.bam", sample=SAMPLES.index),
        # deduplication
        # expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.bam", sample=SAMPLES.index),
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
        # this rule is buggy, CX_report.txt seems only to contain 0s, deactivated for now.
        # expand("coverage/{sample}_MappedOn_" + refbase + "_{context}.CX_report.txt", sample=SAMPLES.index, context=CONTEXT),
        expand("coverage/bws/{sample}_MappedOn_" + refbase + "_{context}.bw", sample=SAMPLES.index, context=CONTEXT),


# this rule is suitable if only mapping is required
# for example for mapping against lambda phage or mitochondrion to test for conversion
# the converstion rates are reported in the log/bismark folder
# call: snakemake map
rule map:
    input:
        expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.bam", sample=SAMPLES.index),

rule trim:
    input:
        expand("trimmed/{sample}_{pair}_trim.fq.gz", sample=SAMPLES.index, pair = PAIRS),

rule sort:
    input:
        expand("mapped/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.bam", sample=SAMPLES.index),

rule bamCoverage:
    input:
        expand("mapped/bws/{sample}_MappedOn_" + refbase +"_trim_bismark_pe.deduplicated.sorted.bw", sample=SAMPLES.index),

rule bw:
    input:
        expand("coverage/bws/{sample}_MappedOn_" + refbase + "_{context}.bw", sample=SAMPLES.index, context=CONTEXT),

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
    conda: "environment.yaml"
    shell:
        """
        fastqc {input.read} -f fastq --outdir logs/fastqc/raw
        """

# https://cutadapt.readthedocs.io/en/stable/guide.html#bisulfite-sequencing-rrbs
rule trim_individual:
    input:
        read1="data/{sample}_" + PAIRS[0] + ".fastq.gz",
        read2="data/{sample}_" + PAIRS[1] + ".fastq.gz",
    log:
        log1="{sample}_" + PAIRS[0] + ".fastq.gz_trimming_report.txt",
        log2="{sample}_" + PAIRS[1] + ".fastq.gz_trimming_report.txt",
    benchmark:
        "{sample}_" + PAIRS[0] + ".fastq.gz_trimming_report.benchmark",
    output:
        # trim_galore alway seems to output fq indpendent of fq|fastq input
        read1="trimmed/{sample}_" + PAIRS[0] + "_trim.fq.gz",
        read2="trimmed/{sample}_" + PAIRS[1] + "_trim.fq.gz",
    conda: "environment.yaml"
    shell:
        """
        trim_galore {extra_params_trim} --paired --trim1 {input.read1} {input.read2} --output_dir trimmed
        rename 's/val_[12]/trim/g' trimmed/{wildcards.sample}*
        mkdir -p logs/trim
        mv  trimmed/{log.log1} logs/trim/
        mv  trimmed/{log.log2} logs/trim/
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
        bam=temp("mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.bam"),
        # se="mapped/{sample}/{sample}_trimmed_bismark_SE_report.txt",
        # nuc="mapped/{sample}/{sample}_trimmed_bismark.nucleotide_stats.txt",
    log:
        log1="logs/bismark/{sample}_MappedOn_{refbase}.log",
        log2="logs/bismark/{sample}_MappedOn_{refbase}_trim_bismark_PE_report.txt",
        log3="logs/bismark/{sample}_MappedOn_{refbase}_trim_bismark_pe.nucleotide_stats.txt",
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}_bismark.benchmark",
    threads: 4
    conda: "environment.yaml"
    shell:
        # USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}
        """
        bismark {extra_params_bismark} --bowtie2 -p {threads} --nucleotide_coverage {REFERENCE} -1 {input.read1} -2 {input.read2} --basename {wildcards.sample}_MappedOn_{refbase}_trim_bismark --output_dir mapped 2> {log.log1}
        mv mapped/{wildcards.sample}_MappedOn_{refbase}_trim_bismark_PE_report.txt {log.log2}
        mv mapped/{wildcards.sample}_MappedOn_{refbase}_trim_bismark_pe.nucleotide_stats.txt {log.log3}
        """


rule deduplicate:
    input:
        # bam="mapped/{sample}/{sample}_trim_bismark.bam",
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.bam",
    output:
        # bam="mapped/{sample}/{sample}_trim_bismark.deduplicated.bam",
        bam=temp("mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.bam"),
        # dup="mapped/{sample}/{sample}_trim_bismark.deduplication_report.txt",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}.deduplication.log",
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}.deduplication.benchmark",
    conda: "environment.yaml"
    shell:
        """
        deduplicate_bismark --paired --bam {input.bam} --output_dir mapped 2> {log}
        """

rule sort_individual:
    input:
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.bam",
    output:
        sort="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam",
        index="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam.bai",
    conda: "environment.yaml"
    shell:
        """
        samtools sort {input.bam} > {output.sort}
        samtools index {output.sort}
        """

rule bamCoverage_individual:
    """compute coverage using deeptools into bigWig(bw) file"""
    input:
        bam="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam",
        bai="mapped/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bam.bai",
    output:
        "mapped/bws/{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.bw",
    conda: "environment.yaml"
    params:
        binsize=BINSIZE,
        extra=extra_params_bamCoverage,
    threads: 8
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}.bamCoverage_individual.log",
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}.bamCoverage_individual.benchmark",
    shell:
        """
        bamCoverage {params.extra} -b {input.bam} -o {output}  --binSize {params.binsize} -p {threads}
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
    conda: "environment.yaml"
    params:
        ref=REFERENCE,
        extra=extra_params_meth_extractor,
    threads: 4
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}.methylation_extractor.log",
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}.methylation_extractor.benchmark",
    shell:
        """
        bismark_methylation_extractor {params.extra} --gzip --paired-end --multicore {threads}  --genome_folder {params.ref} -s {input.bam} --output methylation_extracted 2> {log}
        """

rule bismark2bedGraph:
    # bedGraph file as well as a coverage file which are both sorted by chromosomal position.
    input:
        # extract="mapped/{sample}/{context}_context_{sample}_trim_bismark.deduplicated.sorted.txt.gz",
        methylex="methylation_extracted/{context}_context_{sample}_MappedOn_{refbase}_trim_bismark_pe.deduplicated.sorted.txt.gz",
    output:
        bedGraph="coverage/{sample}_MappedOn_{refbase}_{context}",
        cov=     "coverage/{sample}_MappedOn_{refbase}_{context}.gz.bismark.cov",
    params:
        extra=extra_params_bismark2bedGraph,
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.bismark2bedGraph.log",
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.bismark2bedGraph.benchmark",
    conda: "environment.yaml"
    threads: 4
    shell:
        # gzip outfiles even without --gzip option. Need to unzip manual for later bigwiggle rule
        """
        ulimit -n 10240
        bismark2bedGraph {params.extra} --CX {input.methylex} -o {wildcards.sample}_MappedOn_{refbase}_{wildcards.context} --dir coverage 2> {log}
        pigz -p {threads} -d {output.bedGraph}.gz
        pigz -p {threads} -d {output.cov}.gz
        """

# this rule is buggy, CX_report.txt seems only to contain 0s, deactivated for now.
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
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.coverage2cytosine.benchmark",
    params:
        ref=REFERENCE,
        extra=extra_params_coverage2cytosine,
    conda: "environment.yaml"
    shell:
        # genome_folder needs to contain fa rather than fasta!
        """
        coverage2cytosine {params.extra} --genome_folder {params.ref} -o {input.bedGraph} --dir . {input.bedGraph} 2> {log}
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
        bw="coverage/bws/{sample}_MappedOn_{refbase}_{context}.bw",
    log:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.bedGraphToBigWig.log",
    benchmark:
        "logs/bismark/{sample}_MappedOn_{refbase}_{context}.bedGraphToBigWig.benchmark",
    params:
        ref=REFERENCE,
    conda: "environment.yaml"
    shell:
        # fai files needs to be present to serve as chrom.sizes and can be created by:
        # samtools faidx input.fa
        # cut -f1,2 input.fa.fai > sizes.genome
        """
        ulimit -n 10240
        sed 1d {input.bedGraph} | LC_COLLATE=C sort -k1,1 -k2,2n > {input.bedGraph}_sorted
        bedGraphToBigWig {input.bedGraph}_sorted {params.ref}/{refbase}.fa.fai {output} 2> {log}
        rm {input.bedGraph}_sorted
        """

onstart:
    shell("ulimit -n 10240")

onsuccess:
    shell("ulimit -n ")

onerror:
    shell("ulimit -n ")
    print("\n !!! ERROR in index creation workflow! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
	# ${bismark_dir}coverage2cytosine --genome_folder ${index_bismarck} --CX -o $*_fusedICv2_all.cov2cyt $*_fusedICv2_all.bedgraph.gz.bismark.cov.gz

	# gunzip $*_fusedICv2_CHG.bedgraph.gz $*_fusedICv2_CpG.bedgraph.gz $*_fusedICv2_CHH.bedgraph.gz
	# ${bedGraphToBigWig} $*_fusedICv2_CpG.bedgraph  ${chromsizes} $*_fusedICv2_CpG.bw
	# ${bedGraphToBigWig} $*_fusedICv2_CHH.bedgraph  ${chromsizes} $*_fusedICv2_CHH.bw
	# ${bedGraphToBigWig} $*_fusedICv2_CHG.bedgraph  ${chromsizes} $*_fusedICv2_CHG.bw
	# rm $*_fusedICv2_CHG.bedgraph $*_fusedICv2_CpG.bedgraph $*_fusedICv2_CHH.bedgraph
