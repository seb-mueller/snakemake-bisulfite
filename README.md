# snakemake-bisulfite
Automatic pipeline for genome wide DNA methylome processing data based on snakemake

# Automated workflow for WGBS and EM-Seq (NEB) data
Snakemake workflow for processing BS-seq libaries produced by Illumina bisulfite sequencing kits.

Supported protocols:
- [**WGBS**](https://en.wikipedia.org/wiki/Bisulfite_sequencing)
- [NEBNext Enzymatic Methyl-seq **(EM-seq)**](https://international.neb.com/about-neb/news-and-press-releases/new-england-biolabs-to-present-latest-innovations-for-ngs-sample-preparation-at-agbt-2017)

Note, both WGBS and EM-seq produce the same kind of data, so no adjustment for postprocessing needed.

# Requirments
- demultiplex fastq files in located in `data` directory. They need to be in the form `{sample}_R1.fastq.gz`
- `Snakefile` shipped with this repository.
- `config.yaml` shipped with this repository. It contains all parameters and settings to customize the processing of the current dataset.
-  `samples.csv` listing all samples in the `data` directory withoug the `_R1.fastq.gz` suffix. The first line is the header i.e. the work `library`. An example is shipped with this repository which can be used as a template.
- Optionall: `environment.yaml` to create the software environment if conda is used.
- Installation of [snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- If conda is not used, `bowtie`, `fastqc`, `samtools` and `deeptools` need to be in the PATH.

    The above files can be downloaded as a whole by cloning the repository (which requires git):

```
git clone https://github.com/seb-mueller/snakemake-bisulfite
```
Or individually for example the `Snakemake` file using `wget`:

```
wget https://raw.githubusercontent.com/seb-mueller/snakemake_sRNAseq/master/Snakefile
```

# creating conda environment
```
conda env create --file environment.yaml --name bsseq_pipeline
```

# activate

```
conda activate bsseq_pipeline
```
To `deactivate` the environment, run:

```
conda deactivate
```
# Genome preparation:

The pipeline requires a bismark index file which might have to be created at first.
This best done having the reference genome as fasta file located in a directory by itself (there shouldn't be any other fasta files since bismark works in mysterious ways and you can't specify a conrete file, just the directory).

Say our genome is located here: `ref/genome/mygenome.fa` with no other fa files in it.
Run bismark indexer which will create a folder within that:

```
bismark_genome_preparation ref/genome 
```
This will create `ref/genome/Bisulfite_Genome`.
To let the pipeline know where the index is located, change the `config.yaml`:

```
  reference: "ref/genome/"
  reference_short: "mygenome"
```

Also, the fasta-index is required at some point and should be generated in the same directory as follows:

```
samtools faidx mygenome.fa
```


# Usage:

Navigate in a Unix shell to the base directory contains the files listed above plus the `data` directory including the data like int this example:

```

.
├── bsseq.makefile
├── config.yaml
├── data
│   ├── bsseq_sample1_R1.fastq.gz
│   ├── bsseq_sample1_R2.fastq.gz
│   ├── bsseq_sample2_R1.fastq.gz
│   └── bsseq_sample2_R2.fastq.gz
├── environment.yaml
├── README.md
├── samples.csv
└── Snakefile

```

The `sample.csv` could look something like this:
```
library,optional_info
bsseq_sample1,WT
bsseq_sample2,Mutant
```

Then just run snakmake in base directory:


```
snakemake
```
## useful parameters:
- `--cores` max number of threads
- `-n` dryrun
- `-p` print commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun postmapping` forces rerun of a given rule (e.g. `postmapping`)


# Output:

`trimmed`, `log` and `mapped` directory with trimming and mapping results.

Below is the structure of all generated files once the pipeline is finished:
```
.
├── config.yaml
├── coverage
│   ├── bseq_sample2_MappedOn_chloroplast_CHH.gz
│   ├── bseq_sample2_MappedOn_chloroplast_CHH.gz.bismark.cov.gz
│   ├── bsseq_sample1_MappedOn_chloroplast_CHG
│   ├── bsseq_sample1_MappedOn_chloroplast_CHG.bw
│   ├── bsseq_sample1_MappedOn_chloroplast_CHG.CX_report.txt
│   ├── bsseq_sample1_MappedOn_chloroplast_CHG.gz.bismark.cov
│   ├── bsseq_sample1_MappedOn_chloroplast_CHH
│   ├── bsseq_sample1_MappedOn_chloroplast_CHH.bw
│   ├── bsseq_sample1_MappedOn_chloroplast_CHH.CX_report.txt
│   ├── bsseq_sample1_MappedOn_chloroplast_CHH.gz.bismark.cov
│   ├── bsseq_sample1_MappedOn_chloroplast_CpG
│   ├── bsseq_sample1_MappedOn_chloroplast_CpG.bw
│   ├── bsseq_sample1_MappedOn_chloroplast_CpG.CX_report.txt
│   ├── bsseq_sample1_MappedOn_chloroplast_CpG.gz.bismark.cov
│   ├── bsseq_sample2_MappedOn_chloroplast_CHG
│   ├── bsseq_sample2_MappedOn_chloroplast_CHG.bw
│   ├── bsseq_sample2_MappedOn_chloroplast_CHG.CX_report.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_CHG.gz.bismark.cov
│   ├── bsseq_sample2_MappedOn_chloroplast_CHH
│   ├── bsseq_sample2_MappedOn_chloroplast_CHH.bw
│   ├── bsseq_sample2_MappedOn_chloroplast_CHH.CX_report.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_CHH.gz.bismark.cov
│   ├── bsseq_sample2_MappedOn_chloroplast_CpG
│   ├── bsseq_sample2_MappedOn_chloroplast_CpG.bw
│   ├── bsseq_sample2_MappedOn_chloroplast_CpG.CX_report.txt
│   └── bsseq_sample2_MappedOn_chloroplast_CpG.gz.bismark.cov
├── data
│   ├── bsseq_sample1_R1.fastq.gz
│   ├── bsseq_sample1_R2.fastq.gz
│   ├── bsseq_sample2_R1.fastq.gz
│   ├── bsseq_sample2_R2.fastq.gz
│   └── index
│       ├── Bisulfite_Genome
│       │   ├── CT_conversion
│       │   │   ├── BS_CT.1.bt2
│       │   │   ├── BS_CT.2.bt2
│       │   │   ├── BS_CT.3.bt2
│       │   │   ├── BS_CT.4.bt2
│       │   │   ├── BS_CT.rev.1.bt2
│       │   │   ├── BS_CT.rev.2.bt2
│       │   │   └── genome_mfa.CT_conversion.fa
│       │   └── GA_conversion
│       │       ├── BS_GA.1.bt2
│       │       ├── BS_GA.2.bt2
│       │       ├── BS_GA.3.bt2
│       │       ├── BS_GA.4.bt2
│       │       ├── BS_GA.rev.1.bt2
│       │       ├── BS_GA.rev.2.bt2
│       │       └── genome_mfa.GA_conversion.fa
│       ├── chloroplast.fa
│       ├── chloroplast.fa.fai
│       └── genomic_nucleotide_frequencies.txt
├── environment.yaml
├── logs
│   ├── bismark
│   │   ├── bsseq_sample1.
│   │   ├── bsseq_sample1_CHG_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample1_CHH_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample1_CpG_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample1.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CHG.bedGraphToBigWig.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CHG.bismark2bedGraph.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CHG.coverage2cytosine.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CHH.bedGraphToBigWig.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CHH.bismark2bedGraph.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CHH.coverage2cytosine.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CpG.bedGraphToBigWig.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CpG.bismark2bedGraph.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast_CpG.coverage2cytosine.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast.deduplication.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast.log
│   │   ├── bsseq_sample1_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample2_CHG_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample2_CHH_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample2_CpG_MappedOn_chloroplast.methylation_extractor.log
│   │   ├── bsseq_sample2.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CHG.bedGraphToBigWig.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CHG.bismark2bedGraph.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CHG.coverage2cytosine.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CHH.bedGraphToBigWig.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CHH.bismark2bedGraph.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CHH.coverage2cytosine.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CpG.bedGraphToBigWig.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CpG.bismark2bedGraph.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast_CpG.coverage2cytosine.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast.deduplication.log
│   │   ├── bsseq_sample2_MappedOn_chloroplast.log
│   │   └── bsseq_sample2_MappedOn_chloroplast.methylation_extractor.log
│   ├── cutadapt
│   │   ├── bsseq_sample1_R1.log
│   │   ├── bsseq_sample1_R2.log
│   │   ├── bsseq_sample2_R1.log
│   │   └── bsseq_sample2_R2.log
│   └── fastqc
│       └── raw
│           ├── bsseq_sample1_R1_fastqc.html
│           ├── bsseq_sample1_R1_fastqc.zip
│           ├── bsseq_sample1_R2_fastqc.html
│           ├── bsseq_sample1_R2_fastqc.zip
│           ├── bsseq_sample2_R1_fastqc.html
│           ├── bsseq_sample2_R1_fastqc.zip
│           ├── bsseq_sample2_R2_fastqc.html
│           └── bsseq_sample2_R2_fastqc.zip
├── mapped
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.bam
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.bam
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bam
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bam.bai
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bw
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplication_report.txt
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.nucleotide_stats.txt
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_PE_report.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.bam
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.bam
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bam
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bam.bai
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bw
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplication_report.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.nucleotide_stats.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_PE_report.txt
│   ├── bws
│   │   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bw
│   │   └── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.bw
│   └── test
├── methylation_extracted
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.M-bias.txt
│   ├── bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted_splitting_report.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.M-bias.txt
│   ├── bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted_splitting_report.txt
│   ├── CHG_context_bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt.gz
│   ├── CHG_context_bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt
│   ├── CHG_context_bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt.gz
│   ├── CHH_context_bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt.gz
│   ├── CHH_context_bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt
│   ├── CHH_context_bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt.gz
│   ├── CpG_context_bsseq_sample1_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt.gz
│   ├── CpG_context_bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt
│   └── CpG_context_bsseq_sample2_MappedOn_chloroplast_trim_bismark_pe.deduplicated.sorted.txt.gz
├── README.md
├── samples.csv
├── Snakefile
└── trimmed
    ├── bsseq_sample1_R1.fastq.gz_trimming_report.txt
    ├── bsseq_sample1_R1_trim.fq.gz
    ├── bsseq_sample1_R2.fastq.gz_trimming_report.txt
    ├── bsseq_sample1_R2_trim.fq.gz
    ├── bsseq_sample2_R1.fastq.gz_trimming_report.txt
    ├── bsseq_sample2_R1_trim.fq.gz
    ├── bsseq_sample2_R1_val_1.fq.gz
    ├── bsseq_sample2_R2.fastq.gz_trimming_report.txt
    ├── bsseq_sample2_R2_trim.fq.gz
    └── bsseq_sample2_R2_val_2.fq.gz
```

# Update pipeline/environment:
```
git pull
conda env update --file environment.yaml --name bsseq_pipeline
```
# Post processing

Results can be imported in various software packages such as [https://github.com/al2na/methylKit](https://github.com/al2na/methylKit) using the R-snippet below. 
Note, here all 3 context are considered seperately which can by useful for example for anyalysing plant DNA-metyhylation.
The imported files are `*.gz.bismark.cov` from the `coverage` folder:

```r
dir_base <- "/mypath..." # set this as base path manually!

meta.data <- read_csv(file.path( dir_base, "samples.csv"))
config    <- read_yaml(file.path(dir_base, "config.yaml"))
genome_suffix <- config$MAPPING$reference_short
dir_files <- file.path(dir_base, "coverage")

contexts <- c("CpG", "CHG", "CHH")
for (context in contexts) {
  fls_cov_short <- paste(meta.data$library, "_MappedOn_", genome_suffix, "_", context, ".gz.bismark.cov", sep = "")
  fls_cov <- file.path(dir_files, fls_cov_short)
  myobj <- methRead(as.list(fls_cov), as.list(as.character(meta.data$library)), "assembly", dbtype = NA,
                    pipeline = "bismarkCoverage", header = FALSE, skip = 0, sep = "\t",
                    context = context, resolution = "base",
                    treatment = as.integer(as.factor(meta.data$condition)), dbdir = getwd(),
                    mincov = 1)
  # ...
}
```
