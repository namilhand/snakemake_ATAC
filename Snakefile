# Snakemake workflow for mapping ChIP-seq libraries to a reference genome

# Usage (snakemake --cores should reflect available cores):
# conda env create --file environment.yaml --name ChIPseq_mapping
# conda activate ChIPseq_mapping
# snakemake -p --cores 48
# conda deactivate

import pandas as pd
import os

# To make the per_base_coverage rule work with a shell script invoked using the "shell" directive,
# we need to determine the base path of Snakefile since we expect the scripts directory to be there as well
SRCDIR = srcdir("")

# Specify config file parameters
configfile: "config.yaml"
sample        = config["SAMPLES"]
reference     = config["MAPPING"]["reference"]
refbase       = os.path.basename(reference)
genomeBinName = config["COVERAGE"]["genomeBinName"]

# Specify the desired end target file(s)
rule all:
    input:
        expand("logs/fastqc/raw/{sample}",
               sample = sample),
        expand("results/01_trimmed/{sample}_{read}.tr.fastq.gz",
               sample = sample,
               read = [1,2]),
        expand("logs/fastqc/trimmed/{sample}",
               sample = sample),
        expand("results/02_bowtie/{sample}_MappedOn_{refbase}.bam",
                sample = sample,
                refbase = refbase),
        expand("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.bam",
               sample = sample,
               refbase = refbase),
        expand("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
               sample = sample,
               refbase = refbase),
        expand("results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai",
               sample = sample,
               refbase = refbase),
        expand("results/03_bamCoverage/bw/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bw",
               sample = sample,
               refbase = refbase),
        expand("results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bedgraph",
               sample = sample,
               refbase = refbase),
        expand("results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_binSize{genomeBinName}.bedgraph",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName),
        expand("results/03_bamCoverage/bw/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_binSize{genomeBinName}.bw",
               sample = sample,
               refbase = refbase,
               genomeBinName = genomeBinName)
        # expand("mapped/both/tsv/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.tsv",
        #        sample = sample,
        #        refbase = refbase,
        #        genomeBinName = genomeBinName),
        # expand("mapped/both/pb/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.perbase",
        #        sample = sample,
        #        refbase = refbase)

# Run fastqc on single-end raw data
rule fastqc_raw:
    """Create fastqc report"""
    input:
        read1 = "raw/{sample}_1.fastq.gz",
        read2 = "raw/{sample}_2.fastq.gz"
    output:
        # html = "logs/fastqc/raw/{sample}_fastqc.html",
        # zip  = "logs/fastqc/raw/{sample}_fastqc.zip"
        directory("logs/fastqc/raw/{sample}")
    params:
        " --extract" +
        " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
    log:
        "logs/fastqc/raw/{sample}.log"
    shell:
        "mkdir -p {output}; fastqc -o {output} {params} {input.read1} {input.read2}"

# Trim off adapters
rule cutadapt:
    """Remove adapters"""
    input:
        read1 = "raw/{sample}_1.fastq.gz",
        read2 = "raw/{sample}_2.fastq.gz"
    output:
        tr1 = "results/01_trimmed/{sample}_1.tr.fastq.gz",
        tr2 = "results/01_trimmed/{sample}_2.tr.fastq.gz",
        qc    = "qc/cutadapt/{sample}_cutadapt.qc.txt"
    params:
        adapter=config["FILTER"]["cutadapt"]["adapter"], 
        min_overlap=config["FILTER"]["cutadapt"]["minimum-overlap"], 
    log:
        "logs/cutadapt/{sample}_trimmed.log"
    threads: config["THREADS"]
    shell:
        r"""
        cutadapt -a {params.adapter} -A {params.adapter} \
        -O {params.min_overlap} \
        {input.read1} {input.read2} \
        --cores={threads} \
        -o {output.tr1} -p {output.tr2} > {output.qc} 2> {log}
        """


# Run fastqc on trimmed data
rule fastqc_trimmed:
    """Create fastqc report"""
    input:
        tr1="results/01_trimmed/{sample}_1.tr.fastq.gz",
        tr2="results/01_trimmed/{sample}_2.tr.fastq.gz"
    output:
        # html = "logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.html",
        # zip  = "logs/fastqc/trimmed/{sample}_dedup_trimmed_fastqc.zip"
        directory("logs/fastqc/trimmed/{sample}")
    params:
        " --extract" +
        " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
    log:
        "logs/fastqc/trimmed/{sample}_trimmed.log"
    shell:
        "mkdir -p {output}; fastqc -o {output} {params} {input.tr1} {input.tr2}"

# Align to reference genome
# Only primary reads with MAPQ > MAPQmaxi are retained
rule bowtie2:
    """Map reads using bowtie2 and filter alignments using samtools"""
    input:
        # fastq = "data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz",
        tr1="results/01_trimmed/{sample}_1.tr.fastq.gz",
        tr2="results/01_trimmed/{sample}_2.tr.fastq.gz"
    output:
        "results/02_bowtie/{sample}_MappedOn_{refbase}.bam"
    params:
        MAPQmaxi = config["MAPPING"]["MAPQmaxi"]
    threads: config["THREADS"]
    log:
        "logs/bowtie2/{sample}_MappedOn_{refbase}.log"
    shell:
        # -F 2308 excludes unmapped reads,
        # as well as secondary and supplementary alignments
        # Exclude alignments with MAPQ < config["MAPPING"]["MAPQmaxi"]
        "(bowtie2 --very-sensitive --no-mixed"
        " -X 2000"
        " --threads {threads}"
        " -x {reference} -1 {input.tr1} -2 {input.tr2}"
        " | samtools view -bh -@ {threads} -F 2308 -q {params.MAPQmaxi} -o {output} - ) 2> {log}"

# filtering alignments
## -f 2 to retain only properly mapped read pairs
## alignments to organellar DNA are excluded.
rule filter_bam:
    input:"results/02_bowtie/{sample}_MappedOn_{refbase}.bam"
    output:"results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.bam"
    threads: config["THREADS"]
    params: 
        sortMemory = config["MAPPING"]["sortMemory"],
    log:
        "logs/samtools/{sample}_MappedOn_{refbase}_nuclear_sort.log",
    shell:
        r"""
        (samtools view -h -f 2 -@ {threads} {input} | 
        grep -v -e "ChrM" -e "ChrC" | 
        samtools view -u - |
        samtools sort -@ {threads} -m {params.sortMemory} -O bam -o {output} -) 2> {log}
        """

rule markdup:
    output:
        bam="results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
        metric="logs/markdup/{sample}_MappedOn_{refbase}_nuclear_sort.md.txt"
        # index="results/02_bowtie2/nuclear/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    input: "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.bam"
    shell:
        r"""
        picard MarkDuplicates -I {input} \
        -O {output.bam} \
        -M {output.metric} \
        --REMOVE_DUPLICATES true;
        """
rule postmapping:
    """bam.bai samtools flagstat idxstats"""
    input:
        "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam"
    output:
        "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    log:
        flagstat   = "logs/samtools/stats/{sample}_MappedOn_{refbase}_nuclear_sort_md_flagstat.log",
        idxstats   = "logs/samtools/stats/{sample}_MappedOn_{refbase}_nuclear_sort_md_idxstats.log"
    shell:
        """
        samtools index    {input}
        samtools flagstat {input} > {log.flagstat}
        samtools idxstats {input} > {log.idxstats}
        """
rule calc_coverage:
    """Calculate library-size-normalized coverage"""
    input:
        BAM   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
        BAMidx   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    output:
        BW   = "results/03_bamCoverage/bw/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bw",
        BG   = "results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.bedgraph"
    params:
        normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
        ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
        binSize                = config["COVERAGE"]["binSize"]
    log:
        "logs/bamCoverage/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_1bp-resolution.log"
    threads: config["THREADS"]  
    shell:
        "(bamCoverage -b {input.BAM} -o {output.BW}"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads"
        " --binSize {params.binSize} -p {threads}; "
        "bamCoverage -b {input.BAM} -o {output.BG} -of bedgraph"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads"
        " --binSize {params.binSize} -p {threads}) 2> {log}"
rule calc_coverage_genome:
    """Calculate library-size-normalized coverage in adjacent windows"""
    input:
        BAM   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam",
        BAMidx   = "results/02_bowtie/filtered/{sample}_MappedOn_{refbase}_nuclear_sort.md.bam.bai"
    output:
        bg="results/03_bamCoverage/bg/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_binSize{genomeBinName}.bedgraph",
        bw="results/03_bamCoverage/bw/{sample}_MappedOn_{refbase}_nuclear_sort_md_norm_binSize{genomeBinName}.bw"
    params:
        normalizeUsing         = config["COVERAGE"]["normalizeUsing"],
        ignoreForNormalization = config["COVERAGE"]["ignoreForNormalization"],
        genomeBinSize          = config["COVERAGE"]["genomeBinSize"]
    log:
        "logs/bamCoverage/{sample}_MappedOn_{refbase}_nuclear_both_sort_norm_binSize{genomeBinName}.log"
    threads: config["THREADS"]  
    shell:
        "(bamCoverage -b {input.BAM} -o {output.bw}"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads"
        " --binSize {params.genomeBinSize} -p {threads};"
        "bamCoverage -b {input.BAM} -o {output.bg} -of bedgraph"
        " --normalizeUsing {params.normalizeUsing}"
        " --ignoreForNormalization {params.ignoreForNormalization}"
        " --exactScaling"
        " --extendReads"
        " --binSize {params.genomeBinSize} -p {threads}) 2> {log}"

# Use R script genomeBin_bedgraphToTSV.R to convert *{genomeBinName}.bedgraph files into TSV files
# These TSV files can be imported into R for calculating and plotting log2(ChIP/control) chromosome-scale profiles
rule bedgraphToTSV:
    """Convert *{genomeBinName}.bedgraph files into TSV files"""
    input:
        uniqueBGgenome = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.bedgraph",
        bothBGgenome   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.bedgraph"
    output:
        uniqueTSVgenome = "mapped/unique/tsv/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.tsv",
        bothTSVgenome   = "mapped/both/tsv/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.tsv"
    params:
        genomeBinSize = config["COVERAGE"]["genomeBinSize"]
    log:
        unique = "logs/genomeBin_bedgraphToTSV/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_binSize{genomeBinName}.log",
        both   = "logs/genomeBin_bedgraphToTSV/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_binSize{genomeBinName}.log"
    threads: config["THREADS"]
    shell:
        "(scripts/genomeBin_bedgraphToTSV.R"
        " {wildcards.sample}"
        " {refbase}"
        " unique"
        " {params.genomeBinSize}) 2> {log.unique}; "
        "(scripts/genomeBin_bedgraphToTSV.R"
        " {wildcards.sample}"
        " {refbase}"
        " both"
        " {params.genomeBinSize}) 2> {log.both}"

# Convert bedgraph to per-base 1-based coverage file
rule per_base_coverage:
    """Convert bedgraph to per-base 1-based coverage file"""
    input:
        uniqueBG = "mapped/unique/bg/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.bedgraph",
        bothBG   = "mapped/both/bg/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.bedgraph"
    output:
        uniquePB = "mapped/unique/pb/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm.perbase",
        bothPB   = "mapped/both/pb/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm.perbase"
    log:
        unique = "logs/perBaseCoverage/{sample}_MappedOn_{refbase}_lowXM_unique_sort_norm_pb.log",
        both = "logs/perBaseCoverage/{sample}_MappedOn_{refbase}_lowXM_both_sort_norm_pb.log"
    shell:
        "(bash scripts/perbase_1based_coverage.sh {input.bothBG} {output.bothPB} ) 2> {log.both}; "
        "(bash scripts/perbase_1based_coverage.sh {input.uniqueBG} {output.uniquePB} ) 2> {log.unique}"
