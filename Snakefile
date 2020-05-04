#!/usr/bin/env python3

import os
import re
import sys
import tempfile
import pandas as pd
from snake_setup import set_config, load_samples, get_grouping, load_regions, load_vcf_paths, load_genomes, get_allele_groupings

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

if not config:
    configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'workdir':           workflow.basedir,
    'threads':           1          ,
    'data':
        {'samples':      ''         ,
         'phased_vcf':   None       ,},
    'genome':
        {'table':        None        ,
         'genes':        None        ,
         'ctcf':         None        ,
         'ctcf_orient':  None        ,},
    'protocol':
        {'regions':      ''          ,
         're1':          'enzyme'    ,
         're1_seq':      ''          ,
         're2':          None        ,
         're2_seq':      None        ,
         'arima':        False       ,},
    'hicup':
        {'shortest' :    150          ,
         'longest' :     850          ,
         'nofill' :      False        ,},
    'HiCcompare':
        {'fdr' :         0.05         ,
         'logFC' :       1            ,},
    'compareMatrices':
        {'operation' :   'log2ratio'  ,
         'colourmap' :   'bwr'        ,
         'vMin' :        -3           ,
         'vMax' :        3            ,},
    'binsize':           [5000, 10000],
    'fastq_screen':      None,
    'tmpdir':            tempfile.gettempdir(),
    'known_sites':       None
}

config = set_config(config, default_config)

workdir: config['workdir']
THREADS = config['threads']
RE1 = config['protocol']['re1']
RE1_SEQ = config['protocol']['re1_seq']
RE2 = config['protocol']['re2']
RE2_SEQ = config['protocol']['re2_seq']
READS = ['R1', 'R2']
BINS = config['binsize']
BASE_BIN = BINS[0]

# Read path to samples in pandas
samples = load_samples(config['data']['samples'])

# Extract groups and replicates.
ORIGINAL_SAMPLES, ORIGINAL_GROUPS, CELL_TYPES = get_grouping(samples)

REGIONS = load_regions(config['protocol']['regions'])

if config['data']['phased_vcf']:
    PHASED_VCFS = load_vcf_paths(config['data']['phased_vcf'], samples)
    workdir: config['workdir'] + '/allele'
    GROUPS, SAMPLES = get_allele_groupings(ORIGINAL_SAMPLES)
    ALLELE_SPECIFIC = True
else:
    SAMPLES = ORIGINAL_SAMPLES
    GROUPS = ORIGINAL_GROUPS
    ALLELE_SPECIFIC = False
    if not config['known_sites']:
        sys.stderr.write(
            f'\033[31mNo configuration provided for "known_sites" and '
             'no default available.\n')
        sys.exit(1)

GENOMES = load_genomes(config['genome']['table'])


wildcard_constraints:
    cell_type = r'[^-\.\/]+',
    pre_group = r'[^-\.\/g]+',
    pre_sample = r'[^-\.\/g]+-\d+',
    region = r'[^-\.\/]+',
    allele = r'[12]',
    rep = r'\d+',
    read = r'R[12]',
    bin = r'\d+',
    mode = r'SNP|INDEL'

if ALLELE_SPECIFIC:
    wildcard_constraints:
        group = r'[^-\.\/]+_g\d+',
        sample = r'[^-\.\/]+_g\d+-\d+'
else:
    wildcard_constraints:
        group = r'[^-\.\/g]+',
        sample = r'[^-\.\/g]+-\d+'

preQC_mode = ['qc/multiqc', 'qc/multiBamQC', 'qc/filterQC/ditag_length.png']
HiC_mode = [expand('matrices/{region}/{bin}/plots/matrices/{all}-{region}-{bin}.png',
                region=REGIONS.index, bin=BINS, all=SAMPLES+list(GROUPS)),
            expand('matrices/{region}/{bin}/{method}/{all}-{region}-{bin}.{ext}',
                region=REGIONS.index, bin=BINS, all=SAMPLES+list(GROUPS),
                method=['norm', 'ice'], ext=['h5', 'gz']),
            expand('matrices/{region}/{bin}/raw/{sample}-{region}-{bin}.{ext}',
                region=REGIONS.index, bin=BINS,
                sample=SAMPLES, ext=['h5', 'gz']),
            expand('qc/hicrep/{region}-{bin}-hicrep.png',
                region=REGIONS.index, bin=BINS),
            expand('matrices/{region}/{bin}/plots/{region}-{bin}-{group1}-vs-{group2}.png',
                region=REGIONS.index, bin=BINS,
                group1 = list(GROUPS), group2 = list(GROUPS)),
            expand('matrices/{region}/{bin}/plots/{group}-{region}-{bin}.png',
                region=REGIONS.index, bin=BINS, group=list(GROUPS)),
            expand('diffhic/bams/{sample}.bam', sample=SAMPLES),
            expand('diffhic/genome/{cell_type}-custom.fa',
                cell_type=list(CELL_TYPES))]
phase_gatk = [expand('allele/hapcut2/{cell_type}-phased.vcf',
                cell_type=list(CELL_TYPES))] if not ALLELE_SPECIFIC else []
phase_bcf = [expand('qc/variant_quality/{cell_type}-{region}-bcftoolsStats.txt',
                region=REGIONS.index,
                cell_type=list(CELL_TYPES))] if not ALLELE_SPECIFIC else []

rule all:
    input:
        preQC_mode,
        HiC_mode,
        phase_gatk,

if ALLELE_SPECIFIC:
    rule maskPhased:
        input:
            genome = lambda wc: GENOMES[wc.cell_type],
            vcf = lambda wc: PHASED_VCFS[wc.cell_type]
        output:
            'genome/masked/{cell_type}.fa'
        log:
            'logs/maskPhased/{cell_type}.log'
        conda:
            f'{ENVS}/bedtools.yaml'
        shell:
            'bedtools maskfasta -fullHeader '
            '-fi <(zcat -f {input.genome}) '
            '-bed {input.vcf} -fo {output} 2> {log}'


rule vcf2SNPsplit:
    input:
        vcf = lambda wc: PHASED_VCFS[wc.cell_type]
    output:
        'snpsplit/{cell_type}-snpsplit.txt'
    log:
        'logs/vcf2SNPsplit/{cell_type}.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/reformat_snpsplit.sh {input} > {output} 2> {log}'


rule bgzipGenome:
    input:
        lambda wc: GENOMES[wc.cell_type]
    output:
        'genome/{cell_type}.fa.gz'
    log:
        'logs/bgzipGenome/{cell_type}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip > {output}) 2> {log}'


rule indexGenome:
    input:
        rules.bgzipGenome.output
    output:
        f'{rules.bgzipGenome.output}.fai'
    log:
        'logs/indexGenome/{cell_type}-indexGenome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule getChromSizes:
    input:
        rules.indexGenome.output
    output:
        'genome/chrom_sizes/{cell_type}.chrom.sizes'
    log:
        'logs/getChromSizes/{cell_type}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


def bowtie2BuildInput(wildcards):
    if ALLELE_SPECIFIC:
        return rules.maskPhased.output
    else:
        return rules.bgzipGenome.output


rule bowtie2Build:
    input:
        bowtie2BuildInput
    output:
        expand('genome/index/{{cell_type}}.{n}.bt2',
               n=['1', '2', '3', '4', 'rev.1', 'rev.2'])
    params:
        basename = lambda wc: f'genome/index/{wc.cell_type}'
    log:
        'logs/bowtie2Build/{cell_type}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params.basename} &> {log}'

rule reformatGenes:
    input:
        config['genome']['genes']
    output:
        f'genome/genes.bed'
    log:
        f'logs/reformatGenes.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'zcat -f {input} | {SCRIPTS}/reformatGenes.py > {output} 2> {log}'


rule fastQC:
    input:
        lambda wc: samples.xs(wc.single, level=3)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule reformatFastQC:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/reformatFastQC/{single}.raw.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/modify_fastqc.sh {input} {output} '
        '{wildcards.single} &> {log}'

rule hicupTruncate:
    input:
        lambda wc: samples.xs(wc.pre_sample, level=2)['path']
    output:
        truncated = ['fastq/truncated/{pre_sample}-R1.trunc.fastq.gz',
                     'fastq/truncated/{pre_sample}-R2.trunc.fastq.gz'],
        summary = 'qc/hicup/{pre_sample}-truncate-summary.txt'
    params:
        re1_seq = RE1_SEQ,
        fill = '--nofill' if config['hicup']['nofill'] else ''
    threads:
        2 if THREADS > 2 else THREADS
    log:
        'logs/hicupTruncate/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicupTruncate.py {params.fill} '
        '--output {output.truncated} '
        '--summary {output.summary} '
        '--re1 {params.re1_seq} '
        '--threads {threads} {input} &> {log}'


rule cutadapt:
    input:
        rules.hicupTruncate.output.truncated
    output:
        trimmed = ['fastq/trimmed/{pre_sample}-R1.trim.fastq.gz',
                   'fastq/trimmed/{pre_sample}-R2.trim.fastq.gz'],
        qc = 'qc/cutadapt/unmod/{pre_sample}.cutadapt.txt'
    group:
        'cutadapt'
    params:
        others = '--minimum-length 20 --quality-cutoff 20 --discard-trimmed '
        '--gc-content 46 --overlap 6 --error-rate 0.1'
    log:
        'logs/cutadapt/{pre_sample}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        THREADS
    shell:
        'cutadapt '
        '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA '
        '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT '
        '{params.others} --cores {threads} '
        '-o {output.trimmed[0]} -p {output.trimmed[1]} {input} '
        '> {output.qc} 2> {log}'


rule reformatCutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{pre_sample}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/reformatCutadapt/{pre_sample}.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        'awk -v sample={wildcards.pre_sample} -f {SCRIPTS}/modify_cutadapt.awk '
        '{input} > {output} 2> {log}'


if config['fastq_screen'] is not None:
    rule fastQScreen:
        input:
            'fastq/trimmed/{single}.trim.fastq.gz'
        output:
            txt = 'qc/fastq_screen/{single}.fastq_screen.txt',
            png = 'qc/fastq_screen/{single}.fastq_screen.png'
        params:
            fastq_screen_config = config['fastq_screen'],
            subset = 100000,
            aligner = 'bowtie2'
        log:
            'logs/fastq_screen/{single}.log'
        threads:
            8
        wrapper:
            "0.49.0/bio/fastq_screen"


rule fastQCTrimmed:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{single}.trim_fastqc.html',
        zip = 'qc/fastqc/{single}.trim_fastqc.zip'
    log:
        'logs/fastqc_trimmed/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule hicupDigest:
    input:
        rules.bgzipGenome.output
    output:
        f'genome/digest/{{cell_type}}-{RE1}-hicup-digest.txt.gz'
    params:
        re1 = RE1,
        re1_seq = RE1_SEQ,
        re2 = f'--re2_name {RE2}' if RE2 else '',
        re2_seq = '--re2 {RE2_SEQ}' if RE2_SEQ else '',
        arima = '--arima' if config['protocol']['arima'] else ''
    log:
        'logs/hicupDigest/{cell_type}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_digester.py '
        '--output {output} --genome {wildcards.cell_type} '
        '--re1 {params.re1_seq} --re1_name {params.re1} '
        '{params.arima} {params.re2_seq} {params.re2} {input} &> {log}'


def hicup_map_index(wildcards):
    """ Retrieve bowtie2 index associated with sample. """

    for cell_type, samples in CELL_TYPES.items():
        if wildcards.pre_sample in samples:
            type = cell_type

    return expand('genome/index/{cell_type}.{n}.bt2',
        cell_type=type, n=['1', '2', '3', '4', 'rev.1', 'rev.2'])


def hicup_map_basename(wildcards):
    """ Retrieve bowtie2 index basename associated with sample. """

    for cell_type, samples in CELL_TYPES.items():
        if wildcards.pre_sample in samples:
            type = cell_type

    return expand('genome/index/{cell_type}', cell_type=type)


rule hicupMap:
    input:
        reads = rules.cutadapt.output.trimmed,
        bt2_index = hicup_map_index
    output:
        mapped = 'mapped/{pre_sample}.pair.bam',
        summary = 'qc/hicup/{pre_sample}-map-summary.txt'
    params:
        basename = hicup_map_basename
    threads:
        12
    log:
        'logs/hicupMap/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_mapper.py '
        '--output {output.mapped} '
        '--summary {output.summary} '
        '--index {params.basename} '
        '--threads {threads} {input.reads} &> {log}'


rule digest:
    input:
        rules.bgzipGenome.output
    output:
        f'genome/digest/{{cell_type}}-{RE1}-pyHiCtools-digest.txt'
    params:
        re1_seq = RE1_SEQ
    log:
        f'logs/digest/{{cell_type}}-{RE1}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools digest --restriction {params.re1_seq} <(zcat -f {input}) '
        '> {output} 2> {log}'


rule sampleReads:
    input:
        rules.hicupMap.output.mapped
    output:
        pipe('mapped/subsampled/{pre_sample}-subsample.bam')
    group:
        'filterQC'
    params:
        seed = '42',
        frac = '20'
    log:
        'logs/sampleReads/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools view -@ {threads} -s {params.seed}.{params.frac} {input} '
        '> {output} 2> {log}'


def processHiC_digest(wildcards):
    """ Retrieve pyHiCTools digest file associated with sample. """

    for cell_type, samples in CELL_TYPES.items():
        if wildcards.pre_sample in samples:
            type = cell_type

    return expand('genome/digest/{cell_type}-{re1}-pyHiCtools-digest.txt',
        cell_type=type, re1=RE1)


rule processHiC:
    input:
        dedup_nmsort = rules.sampleReads.output,
        digest = processHiC_digest
    output:
        pipe('mapped/subsampled/{pre_sample}.processed.sam')
    group:
        'filterQC'
    log:
        'logs/process/{pre_sample}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools process --digest {input.digest} {input.dedup_nmsort} '
        '> {output} 2> {log}'


rule extractHiC:
    input:
        rules.processHiC.output
    output:
        'mapped/subsampled/{pre_sample}-processed.txt'
    group:
        'filterQC'
    log:
        'logs/extractHiC/{pre_sample}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools extract '
        '--sample {wildcards.pre_sample} {input} '
        '> {output} 2> {log}'


rule plotQC:
    input:
        expand('mapped/subsampled/{pre_sample}-processed.txt',
            pre_sample=ORIGINAL_SAMPLES)
    output:
        expand('qc/filterQC/{fig}',
               fig=['trans_stats.csv', 'insert_size_frequency.png',
                    'ditag_length.png'])
    params:
        outdir = 'qc/filterQC'
    log:
        'logs/plot_filter/plot_subsample.log'
    conda:
        f'{ENVS}/ggplot2.yaml'
    shell:
        '{SCRIPTS}/figures/plot_filter.R {params.outdir} {input} 2> {log}'


def getHicupDigest(wildcards):
    """ Retrieve hicup filter digest file associated with pre_sample. """

    for cell_type, samples in CELL_TYPES.items():
        if wildcards.pre_sample in samples:
            type = cell_type

    return expand('genome/digest/{cell_type}-{re1}-hicup-digest.txt.gz',
        cell_type=type, re1=RE1)


rule hicupFilter:
    input:
        bam = rules.hicupMap.output.mapped,
        digest = getHicupDigest
    output:
        filtered = 'mapped/{pre_sample}.filt.bam',
        summary = 'qc/hicup/{pre_sample}-filter-summary.txt',
        rejects = directory('qc/hicup/{pre_sample}-ditag_rejects')
    params:
        shortest = config['hicup']['shortest'],
        longest = config['hicup']['longest']
    log:
        'logs/hicupFilter/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_filter.py '
        '--output {output.filtered} '
        '--digest {input.digest} '
        '--outdir {output.rejects} '
        '--summary {output.summary} '
        '--shortest {params.shortest} '
        '--longest {params.longest} {input.bam} &> {log}'


rule hicupDeduplicate:
    input:
        rules.hicupFilter.output.filtered
    output:
        deduped = 'mapped/{pre_sample}.dedup.bam',
        summary = 'qc/hicup/{pre_sample}-deduplicate-summary.txt',
    log:
        'logs/hicupDeduplicate/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_deduplicator.py '
        '--output {output.deduped} '
        '--summary {output.summary} {input} &> {log}'


rule mergeHicupQC:
    input:
        truncater = rules.hicupTruncate.output.summary,
        mapper = rules.hicupMap.output.summary,
        filter = rules.hicupFilter.output.summary,
        deduplicator = rules.hicupDeduplicate.output.summary
    output:
        'qc/hicup/HiCUP_summary_report-{pre_sample}.txt'
    log:
        'logs/mergeHicupQC/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/combine_summary.py '
        '--truncater {input.truncater} '
        '--mapper {input.mapper}  '
        '--filter {input.filter}  '
        '--deduplicator {input.deduplicator} '
        '> {output} 2> {log}'


# Currently need to remove PG header due to bug in hicup_deduplicator (0.7.2)
rule remove_PG_header:
    input:
        rules.hicupDeduplicate.output.deduped
    output:
        pipe('mapped/{sample}.dedup-no-PG.bam')
    log:
        'logs/remove_PG_header/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        "samtools reheader <(samtools view -H {input} | grep -v '^@PG') "
        "{input} > {output} 2> {log}"


def SNPsplit_input(wildcards):
    """ Retrieve cell type associated with sample. """

    for cell_type, samples in CELL_TYPES.items():
        if wildcards.pre_sample in samples:
            type = cell_type

    return f'snpsplit/{type}-snpsplit.txt'


rule SNPsplit:
    input:
        bam = rules.hicupDeduplicate.output.deduped,
        snps = SNPsplit_input
    output:
        expand('snpsplit/{{pre_sample}}.dedup.{ext}',
            ext = ['G1_G1.bam', 'G1_G2.bam', 'G1_UA.bam', 'G2_G2.bam',
                   'G2_UA.bam', 'SNPsplit_report.txt', 'SNPsplit_sort.txt',
                   'UA_UA.bam', 'allele_flagged.bam'])
    params:
        outdir = 'snpsplit/'
    log:
        'logs/SNPsplit/SNPsplit-{pre_sample}.log'
    conda:
        f'{ENVS}/snpsplit.yaml'
    shell:
        'SNPsplit {input.bam} --snp_file {input.snps} '
        '--hic --output_dir {params.outdir} &> {log}'


rule mergeSNPsplit:
    input:
        'snpsplit/{pre_group}-{rep}.dedup.G{allele}_G{allele}.bam',
        'snpsplit/{pre_group}-{rep}.dedup.G{allele}_UA.bam'
    output:
        'snpsplit/merged/{pre_group}_g{allele}-{rep}.dedup.bam'
    log:
        'logs/mergeSNPsplit/{pre_group}_g{allele}-{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools merge -n {output} {input} &> {log}'


rule sortBam:
    input:
        rules.hicupDeduplicate.output.deduped
    output:
        'mapped/{pre_sample}.sort.bam'
    params:
        mem = '1G'
    threads:
        THREADS
    log:
        'logs/sortBam/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -m {params.mem} {input} '
        '> {output} 2> {log}'


rule indexBam:
    input:
        rules.sortBam.output
    output:
        f'{rules.sortBam.output}.bai'
    threads:
        THREADS
    log:
        'logs/indexBam/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools index -@ 6 {input} &> {log}'


rule samtoolsStats:
    input:
        rules.sortBam.output
    output:
        'qc/samtools/stats/{pre_sample}.stats.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtoolsStats/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'


rule samtoolsIdxstats:
    input:
        bam = rules.sortBam.output,
        index = rules.indexBam.output
    output:
        'qc/samtools/idxstats/{pre_sample}.idxstats.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtoolsIdxstats/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input.bam} > {output} 2> {log}'


rule samtoolsFlagstat:
    input:
        rules.sortBam.output
    output:
        'qc/samtools/flagstat/{pre_sample}.flagstat.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtoolsFlagstat/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'


rule bamQC:
    input:
        bam = rules.sortBam.output,
        regions = config['protocol']['regions']
    output:
        directory('qc/bamQC/{pre_sample}')
    resources:
        mem_mb = 3000
    log:
        'logs/bamQC/{pre_sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap bamqc '
        '-bam {input.bam} -gff {input.regions} -outdir {output} '
        '-nt 6 --outside-stats --java-mem-size={resources.mem_mb}M '
        '--paint-chromosome-limits &> {log}'


rule multibamQCConfig:
    input:
        expand('qc/bamQC/{pre_sample}', pre_sample=ORIGINAL_SAMPLES)
    output:
        'qc/bamQC/multibamQC.config'
    group:
        'multiBamQC'
    log:
        'logs/multibamQC_config/multibamQC_config.log'
    shell:
        '{SCRIPTS}/multibamqc_config.py {input} > {output} 2> {log}'


rule multibamQC:
    input:
        rules.multibamQCConfig.output
    output:
        directory('qc/multiBamQC')
    group:
        'multiBamQC'
    log:
        'logs/multibamQC/multibamQC.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    shell:
        'qualimap multi-bamqc --data {input} -outdir {output} &> {log}'


def split_input(wildcards):
    if ALLELE_SPECIFIC:
        return 'snpsplit/merged/{sample}.dedup.bam'
    else:
        return f'mapped/{wildcards.sample}.dedup.bam'


rule splitPairedReads:
    input:
        split_input
    output:
        'mapped/split/{sample}-{read}.hic.bam'
    params:
        read = READS,
        flag = lambda wc: '0x40' if wc.read == READS[0] else '0x80'
    log:
        'logs/splitPairedReads/{sample}-{read}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -@ {threads} -f {params.flag} -b {input} '
        '> {output} 2> {log}'


rule buildBaseMatrix:
    input:
        expand('mapped/split/{{sample}}-{read}.hic.bam', read=READS)
    output:
        hic = f'matrices/{{region}}/base/raw/{{sample}}-{{region}}.{BASE_BIN}.h5',
        bam = 'matrices/{region}/{sample}-{region}.bam',
        qc = directory(f'qc/hicexplorer/{{sample}}-{{region}}.{BASE_BIN}_QC')
    params:
        bin = BASE_BIN,
        region = REGIONS.index,
        chr = lambda wildcards: REGIONS['chr'][wildcards.region],
        start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
        end = lambda wildcards: REGIONS['end'][wildcards.region]
    log:
        'logs/buildBaseMatrix/{sample}-{region}.log'
    threads:
        THREADS
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicBuildMatrix --samFiles {input} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--binSize {params.bin} '
        '--outFileName {output.hic} '
        '--outBam {output.bam} '
        '--QCfolder {output.qc} '
        '--skipDuplicationCheck '
        '--threads {threads} '
        '&> {log} || mkdir -p {output.qc}; touch {output.hic} {output.bam}'


rule mergeBins:
    input:
        f'matrices/{{region}}/base/raw/{{sample}}-{{region}}.{BASE_BIN}.h5'
    output:
        'matrices/{region}/{bin}/raw/{sample}-{region}-{bin}.h5'
    params:
        bin = config['binsize'],
        nbins = lambda wildcards: int(int(wildcards.bin) / BASE_BIN)
    log:
        'logs/mergeBins/{sample}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log} || touch {output}'


rule readCountNormalise:
    input:
        expand('matrices/{{region}}/{{bin}}/raw/{sample}-{{region}}-{{bin}}.h5',
               sample=SAMPLES)
    output:
        expand('matrices/{{region}}/{{bin}}/norm/{sample}-{{region}}-{{bin}}.h5',
               sample=SAMPLES)
    params:
        method = 'smallest'
    log:
        'logs/readCountNormalise/{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicNormalize --matrices {input} '
        '--outFileName {output} --normalize {params.method} '
        '&> {log} || touch {output} '


rule sumReplicates:
    input:
        lambda wildcards: expand(
            'matrices/{{region}}/{{bin}}/norm/{group}-{rep}-{{region}}-{{bin}}.h5',
            group=wildcards.group, rep=GROUPS[wildcards.group])
    output:
        f'matrices/{{region}}/{{bin}}/norm/{{group}}-{{region}}-{{bin}}.h5'
    log:
        'logs/sumReplicates/{group}-{bin}-{region}.h5'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicSumMatrices --matrices {input} --outFileName {output} '
        '&> {log} || touch {output}'


rule IceMatrix:
    input:
        'matrices/{region}/{bin}/norm/{all}-{region}-{bin}.h5'
    output:
        plot = 'qc/IceMatrix/{all}-{region}-{bin}-diagnosic_plot.png',
        matrix = 'matrices/{region}/{bin}/ice/{all}-{region}-{bin}.h5'
    params:
        iternum = 1000,
        upper_threshold = 5
    log:
        'logs/IceMatrix/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '{SCRIPTS}/hic_correct.sh -p {output.plot} '
        '-o {output.matrix} -u {params.upper_threshold} '
        '-i {params.iternum} -@ 1 {input} &> {log} || touch {output}'


rule distanceNormalise:
    input:
        rules.IceMatrix.output.matrix
    output:
        'matrices/{region}/{bin}/ice/obs_exp/{all}-{region}-{bin}.h5'
    params:
        method = 'obs_exp'
    log:
        'logs/distanceNormalise/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicTransform -m {input} --method {params.method} -o {output} '
        '&> {log} || touch {output}'


rule plotMatrix:
    input:
        rules.distanceNormalise.output
    output:
        'matrices/{region}/{bin}/plots/matrices/{all}-{region}-{bin}.png'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        title = '{all}-{region}-{bin}-norm-ice-obs_exp',
        dpi = 600,
        colour = 'YlGn'
    log:
        'logs/plotMatrix/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicPlotMatrix --matrix {input} '
        '--outFileName {output} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--colorMap {params.colour} '
        '--title {params.title} '
        '--vMin 0 --vMax 2 --dpi {params.dpi} '
        '&> {log} || touch {output}'


rule TadInsulation:
    input:
        rules.IceMatrix.output.matrix
    output:
        expand(
            'matrices/{{region}}/{{bin}}/tads/{{all}}-{{region}}-{{bin}}{ext}',
            ext = ['_boundaries.bed', '_boundaries.gff', '_domains.bed',
                   '_score.bedgraph', '_zscore_matrix.h5']),
        score = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}_tad_score.bm'
    params:
        method = 'fdr',
        bin = lambda wc: wc.bin,
        region = lambda wc: wc.region,
        all = lambda wc: wc.all,
        min_depth = lambda wc: int(wc.bin) * 3,
        max_depth = lambda wc: int(wc.bin) * 10,
        prefix = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}'
    threads:
        THREADS
    log:
        'logs/TadInsulation/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicFindTADs --matrix {input} '
        '--minDepth {params.min_depth} --maxDepth {params.max_depth} '
        '--step {wildcards.bin} --outPrefix {params.prefix} '
        '--correctForMultipleTesting {params.method} '
        '--numberOfProcessors {threads} &> {log} || touch {output}'


rule detectLoops:
    input:
        rules.IceMatrix.output.matrix
    output:
        'matrices/{region}/{bin}/loops/{all}-{region}-{bin}.bedgraph'
    params:
        minLoop = 5000,
        maxLoop = 1000000,
        windowSize = 10,
        peakWidth = 6,
        pValuePre = 0.05,
        pValue = 0.05,
        peakInter = 5
    log:
        'logs/detectLoops/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    threads:
        THREADS
    shell:
        'hicDetectLoops --matrix {input} --outFileName {output} '
        '--minLoopDistance {params.minLoop} '
        '--maxLoopDistance {params.maxLoop} '
        '--windowSize {params.windowSize} '
        '--peakWidth {params.peakWidth} '
        '--pValuePreselection {params.pValuePre} '
        '--pValue {params.pValue} '
        '--peakInteractionsThreshold {params.peakInter} '
        '--threads {threads} '
        '&> {log} || touch {output} && touch {output} '


# Deprecetad rule - now we build comparison matrices from HiCcompare output see homerToH5
rule compareMatrices:
    input:
        expand('matrices/{{region}}/{{bin}}/ice/{group}-{{region}}-{{bin}}.h5',
               group = list(GROUPS))
    output:
        expand('matrices/{{region}}/{{bin}}/ice/{group1}-vs-{group2}-{{region}}-{{bin}}.h5',
               group1 = list(GROUPS), group2 = list(GROUPS))
    log:
        'logs/compareMatrices/{region}-{bin}.log'
    params:
        operation = config['compareMatrices']['operation']
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '{SCRIPTS}/compareMatrices.py --operation {params.operation} '
        '{input} &> {log}'


rule reformatHomer:
    input:
        'matrices/{region}/{bin}/{method}/{all}-{region}-{bin}.h5'
    output:
        'matrices/{region}/{bin}/{method}/{all}-{region}-{bin}.gz'
    log:
        'logs/reformatHomer/{all}-{region}-{bin}-{method}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat h5 --outputFormat homer &> {log} || touch {output}'


rule reformatNxN3p:
    input:
        'matrices/{region}/{bin}/raw/{sample}-{region}-{bin}.gz'
    output:
        'matrices/{region}/{bin}/raw/{sample}-{region}-{bin}.nxnp3.tsv'
    params:
        region = REGIONS.index,
    log:
        'logs/reformat_hicrep/{sample}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/reformat_hicrep.sh '
        '-r {wildcards.region} -b {wildcards.bin} {input}'
        '> {output} 2> {log}'


rule HiCRep:
    input:
        expand(
            'matrices/{{region}}/{{bin}}/raw/{sample}-{{region}}-{{bin}}.nxnp3.tsv',
            sample = SAMPLES)
    output:
        'qc/hicrep/{region}-{bin}-hicrep.png'
    params:
        bin = BINS,
        region = REGIONS.index,
        start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
        end = lambda wildcards: REGIONS['end'][wildcards.region]
    log:
        'logs/HiCRep/{region}-{bin}.log'
    conda:
        f'{ENVS}/hicrep.yaml'
    shell:
        '{SCRIPTS}/figures/plot_hicrep.R {output} {wildcards.bin} '
        '{params.start} {params.end} {input} &> {log}'


rule reformatNxN:
    input:
        'matrices/{region}/{bin}/ice/{all}-{region}-{bin}.gz'
    output:
        'matrices/{region}/{bin}/ice/{all}-{region}-{bin}.nxn.tsv'
    log:
        'logs/reformatNxN/{all}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/reformatNxN.sh {input} > {output} 2> {log} || touch {output}'


rule OnTAD:
    input:
        rules.reformatNxN.output
    output:
        bed = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.bed',
        tad = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.tad'
    params:
        bin = BINS,
        region = REGIONS.index,
        chr = lambda wc: re.sub('chr', '', str(REGIONS['chr'][wc.region])),
        length = lambda wc: REGIONS['length'][wc.region],
        outprefix = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad'
    log:
        'logs/OnTAD/{all}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/OnTAD {input} -o {params.outprefix} -bedout chr{params.chr} '
        '{params.length} {wildcards.bin} &> {log} || touch {output}'


rule reformatLinks:
    input:
        rules.OnTAD.output.bed
    output:
        'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.links'
    params:
        region = REGIONS.index,
        start = lambda wildcards: REGIONS['start'][wildcards.region]
    log:
        'logs/reformatLinks/{all}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/reformat_ontad.sh {params.start} {wildcards.bin} {input} '
        '> {output} 2> {log} || touch {output}'


rule createConfig:
    input:
        #obsexp_matrix = 'matrices/{region}/{bin}/ice/obs_exp/{group}-{region}-{bin}.h5',
        matrix = 'matrices/{region}/{bin}/ice/{group}-{region}-{bin}.h5',
        loops = 'matrices/{region}/{bin}/loops/{group}-{region}-{bin}.bedgraph',
        insulations = 'matrices/{region}/{bin}/tads/{group}-{region}-{bin}_tad_score.bm',
        tads = 'matrices/{region}/{bin}/tads/{group}-{region}-{bin}-ontad.links',
        genes = rules.reformatGenes.output
    output:
        'matrices/{region}/{bin}/plots/configs/{group}-{region}-{bin}.ini'
    conda:
        f'{ENVS}/python3.yaml'
    params:
        ctcf_orientation = config['genome']['ctcf_orient'],
        ctcf = config['genome']['ctcf'],
        depth = lambda wc: int(REGIONS['length'][wc.region] / 2),
    log:
        'logs/createConfig/{group}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} --flip '
        '--insulations {input.insulations} '
        '--loops {input.loops} '
        '--tads {input.tads} '
        '--ctcfs {params.ctcf} '
        '--ctcf_orientation {params.ctcf_orientation} '
        '--genes {input.genes} '
        '--depth {params.depth} > {output} 2> {log}'


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        'matrices/{region}/{bin}/plots/{group}-{region}-{bin}.png'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        title = lambda wc: f'"{wc.group} : {wc.region} at {wc.bin} bin size"',
        dpi = 600
    log:
        'logs/plotHiC/{group}-{region}-{bin}.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    threads:
        12 # Need to ensure it is run 1 at a time!
    shell:
        'pyGenomeTracks --tracks {input} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


rule mergeBamByReplicate:
    input:
        lambda wildcards: expand(
            'matrices/{{region}}/{group}-{rep}-{{region}}.bam',
            group = wildcards.group, rep = GROUPS[wildcards.group])
    output:
        'matrices/{region}/{group}-{region}.bam'
    log:
        'logs/mergeBamByReplicate/{group}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {input} 2> {log}'


rule reformatPre:
    input:
        'matrices/{region}/{all}-{region}.bam'
    output:
        'matrices/{region}/base/raw/{all}-{region}.pre.tsv'
    log:
        'logs/reformatPre/{all}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        '(samtools view -@ {threads} {input} '
        '| awk -f {SCRIPTS}/bam_to_pre.awk > {output}) 2> {log} '


def getChromSizes(wildcards):
    """ Retrieve chromSizes file associated with group or sample. """

    for cell_type, samples in CELL_TYPES.items():
        if ALLELE_SPECIFIC:
            groups, samples = get_allele_groupings(samples)
            groups = list(groups) # Get keys from dictionary as list
        else:
            groups = [sample.split('-')[0] for sample in samples]
        all = samples + groups
        if wildcards.all in all:
            type = cell_type

    return expand('genome/chrom_sizes/{cell_type}.chrom.sizes', cell_type=type)


rule juicerPre:
    input:
        tsv = 'matrices/{region}/base/raw/{all}-{region}.pre.tsv',
        chrom_sizes = getChromSizes
    output:
        'matrices/{region}/{all}-{region}.hic'
    log:
        'logs/juicerPre/{all}-{region}.log'
    params:
        chr = lambda wildcards: REGIONS['chr'][wildcards.region],
        resolutions = ','.join([str(bin) for bin in BINS])
    conda:
        f'{ENVS}/openjdk.yaml'
    shell:
        'java -jar {SCRIPTS}/juicer_tools_1.14.08.jar pre '
        '-c {params.chr} -r {params.resolutions} '
        '{input.tsv} {output} {input.chrom_sizes} 2> {log}'


# Reform for HiCcompare input
rule straw:
    input:
        rules.juicerPre.output
    output:
        'matrices/{region}/{bin}/{all}-{region}-{bin}-sutm.txt'
    log:
        'logs/straw/{all}-{region}-{bin}.log'
    params:
        # Strip 'chr' as juicer removes by default
        chr = lambda wc: re.sub('chr', '', str(REGIONS['chr'][wc.region])),
        start = lambda wildcards: REGIONS['start'][wildcards.region],
        end = lambda wildcards: REGIONS['end'][wildcards.region],
    conda:
        f'{ENVS}/hic-straw.yaml'
    shell:
        '{SCRIPTS}/run-straw.py NONE {input} '
        '{params.chr}:{params.start}:{params.end} '
        '{params.chr}:{params.start}:{params.end} '
        'BP {wildcards.bin} {output} &> {log}'


rule HiCcompare:
    input:
        expand('matrices/{{region}}/{{bin}}/{group}-{{region}}-{{bin}}-sutm.txt',
               group = list(GROUPS))
    output:
        matrices = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}.homer',
            group1 = list(GROUPS), group2 = list(GROUPS)),
        links = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}.links',
            group1 = list(GROUPS), group2 = list(GROUPS)),
        loess = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-loess.png',
            group1 = list(GROUPS), group2 = list(GROUPS)),
        filter = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-filter_params.png',
            group1 = list(GROUPS), group2 = list(GROUPS)),
        compare = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-hicCompare.png',
            group1 = list(GROUPS), group2 = list(GROUPS))
    group:
        'HiCcompare'
    log:
        'logs/HiCcompare/{region}-{bin}.log'
    params:
        dir = lambda wc: f'HiCcompare/{wc.region}/{wc.bin}',
        chr = lambda wc: REGIONS['chr'][wc.region]
    threads:
        12
    conda:
        f'{ENVS}/HiCcompare.yaml'
    shell:
        '({SCRIPTS}/HiCcompare.R {params.dir} {params.chr} {wildcards.bin} '
        '{input} || touch {output}) &> {log}'


rule homerToH5:
    input:
        'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.homer'
    output:
        'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.h5'
    log:
        'logs/homerToH5/{group1}-vs-{group2}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '(hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 || touch {output})  &> {log}'


rule reformatInteract:
    input:
        'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.links'
    output:
        up = 'HiCcompare/{region}/{bin}/{group1}-vs-{group2}-up.interact',
        down = 'HiCcompare/{region}/{bin}/{group1}-vs-{group2}-down.interact'
    group:
        'HiCcompare'
    log:
        'logs/reformatInteract/{region}-{bin}-{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/reformatInteract.py '
        '--up {output.up} --down {output.down} {input} &> {log}'


rule filterHiCcompare:
    input:
        'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.links'
    output:
        up = f'HiCcompare/{{region}}/{{bin}}/{{group1}}-vs-{{group2}}-up.links',
        down = f'HiCcompare/{{region}}/{{bin}}/{{group1}}-vs-{{group2}}-down.links'
    params:
        p_value = config['HiCcompare']['fdr'],
        log_fc = config['HiCcompare']['logFC'],
    group:
        'HiCcompare'
    log:
        'logs/filterHiCcompare/{region}-{bin}-{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filterHiCcompare.py '
        '--up {output.up} --down {output.down} '
        '--p_value {params.p_value} --log_fc {params.log_fc} {input} &> {log}'


rule createCompareConfig:
    input:
        matrix = 'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.h5',
        links = ['HiCcompare/{region}/{bin}/{group1}-vs-{group2}-up.links',
                 'HiCcompare/{region}/{bin}/{group1}-vs-{group2}-down.links'],
        genes = rules.reformatGenes.output
    output:
        'matrices/{region}/{bin}/plots/configs/{region}-{bin}-{group1}-vs-{group2}-compare.ini',
    params:
        ctcf_orientation = config['genome']['ctcf_orient'],
        ctcf = config['genome']['ctcf'],
        depth = lambda wc: int(REGIONS['length'][wc.region] / 2),
        colourmap = config['compareMatrices']['colourmap'],
        vMin = config['compareMatrices']['vMin'],
        vMax = config['compareMatrices']['vMax']
    log:
        'logs/createCompareConfig/{region}-{bin}-{group1}-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} --compare '
        '--genes {input.genes} --loops {input.links} '
        '--ctcfs {params.ctcf} --ctcf_orientation {params.ctcf_orientation} '
        '--depth {params.depth} --colourmap {params.colourmap} '
        '--vMin {params.vMin} --vMax {params.vMax} > {output} 2> {log}'


# Deprecetated rule
rule mergeConfigs:
    input:
        links_up = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-up.links',
            group1 = list(GROUPS), group2 = list(GROUPS)),
        links_down = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-down.links',
            group1 = list(GROUPS), group2 = list(GROUPS)),
        configs = expand(
            'matrices/{{region}}/{{bin}}/plots/configs/{group}-{{region}}-{{bin}}.ini',
            group = list(GROUPS))
    output:
        expand(
            'matrices/{{region}}/{{bin}}/plots/configs/{{region}}-{{bin}}-{group1}-vs-{group2}.ini',
            group1 = list(GROUPS), group2 = list(GROUPS))
    group:
        'HiCcompare'
    log:
        'logs/mergeConfigs/{region}-{bin}.log'
    params:
        prefix = lambda wc: f'matrices/{wc.region}/{wc.bin}/plots/configs/{wc.region}-{wc.bin}-'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeConfigsHiCcompare.py '
            '--up {input.links_up} '
            '--down {input.links_down} '
            '--configs {input.configs} '
            '--prefix {params.prefix} &> {log}'


rule plotCompare:
    input:
        rules.createCompareConfig.output
    output:
        'matrices/{region}/{bin}/plots/{region}-{bin}-{group1}-vs-{group2}.png'
    params:
        # Ensure region is 1 bin smaller because it is trimmed by HiCcompare
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        dpi = 600
    log:
        'logs/plotAnalysis/{region}-{bin}-{group1}-vs-{group2}.log'
    threads:
        12 # Need to ensure it is run 1 at a time!
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'export NUMEXPR_MAX_THREADS=1; pyGenomeTracks --tracks {input} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--outFileName {output} '
        '--title "{wildcards.region} at {wildcards.bin} bin size" '
        '--dpi {params.dpi} &> {log}'


rule modifyBam:
    input:
        rules.buildBaseMatrix.output.bam
    output:
        'diffhic/bams/{region}/{sample}-{region}.bam'
    params:
        region = REGIONS.index,
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region],
        end = lambda wc: REGIONS['end'][wc.region]
    log:
        'logs/modifyBam/{sample}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        '{SCRIPTS}/modify_bams.sh '
        '-r {wildcards.region} -c {params.chr} '
        '-s {params.start} -e {params.end} {input} '
        '> {output} 2> {log} || touch {output}'


rule mergeBamByRegion:
    input:
        expand('diffhic/bams/{region}/{{sample}}-{region}.bam',
            region = REGIONS.index)
    output:
        'diffhic/bams/{sample}.bam'
    log:
        'logs/mergeBamByRegion/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools merge -n {output} {input} 2> {log}'


rule createCustomGenome:
    input:
        genome = rules.bgzipGenome.output,
        index = rules.indexGenome.output,
        regions = config['protocol']['regions']
    output:
        'diffhic/genome/{cell_type}-custom.fa'
    log:
        'logs/createCustomGenome/{cell_type}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        '{SCRIPTS}/create_custom_genome.sh {input.regions} {input.genome} '
        '> {output} 2> {log}'


if not ALLELE_SPECIFIC:

    # Need to add IndexFeature for each input feature provided in config
    rule mergeBamByCellType:
        input:
            lambda wildcards: expand('mapped/{pre_sample}.pair.bam',
                pre_sample = CELL_TYPES[wildcards.cell_type])
        output:
            pipe('mapped/merged_by_cell/{cell_type}.merged.bam')
        group:
            'mergeCellType'
        log:
            'logs/mergeBamByCellType/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            max(1, (THREADS - 1) / 2)
        shell:
            'samtools merge -@ {threads} -O bam,level=0 '
            '-n - {input} > {output} 2> {log}'


    rule fixmate:
        input:
            rules.mergeBamByCellType.output
        output:
            pipe('mapped/merged_by_cell/{cell_type}.fixed.bam')
        group:
            'mergeCellType'
        log:
            'logs/fixmate/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools fixmate -O bam,level=0 '
            '-pmr {input} - > {output} 2> {log}'


    rule addReadGroup:
        input:
            rules.fixmate.output
        output:
            temp('mapped/merged_by_cell/{cell_type}.fixed-RG.bam')
        group:
            'mergeCellType'
        log:
            'logs/addReadGroup/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            max(1, (THREADS - 1) / 2)
        shell:
            'samtools addreplacerg -@ {threads} -O bam,level=-1 '
            '-r "ID:1\tPL:.\tPU:.\tLB:.\tSM:{wildcards.cell_type}" '
            '{input} > {output} 2> {log}'


    rule sortMergedBam:
        input:
            rules.addReadGroup.output
        output:
            pipe('mapped/merged_by_cell/{cell_type}.sort.bam')
        group:
            'dedup'
        params:
            mem = '1G'
        threads:
            THREADS - 2 if THREADS > 2 else 1
        log:
            'logs/sortMergedBam/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools sort -@ {threads} -O bam,level=0 '
            '-m {params.mem} {input} > {output} 2> {log}'


    rule deduplicate:
        input:
            rules.sortMergedBam.output
        output:
            bam = 'mapped/merged_by_cell/{cell_type}.dedup.bam',
            qc = 'qc/deduplicate/{cell_type}.txt'
        group:
            'dedup'
        threads:
            2 if THREADS - 2 > 1 else 1
        log:
            'logs/deduplicate/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools markdup -@ {threads} -O bam,level=-1 '
            '-rsf {output.qc} {input} {output.bam} &> {log}'


    rule indexMergedBam:
        input:
            rules.deduplicate.output.bam
        output:
            f'{rules.deduplicate.output.bam}.bai'
        threads:
            THREADS
        log:
            'logs/samtools/index/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools index -@ {threads} {input} &> {log}'


    rule createSequenceDictionary:
        input:
            rules.bgzipGenome.output
        output:
            'genome/{cell_type}.dict'
        params:
            tmp = config['tmpdir']
        log:
            'logs/gatk/createSequenceDictionary/{cell_type}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard CreateSequenceDictionary R={input} O={output} '
            'TMP_DIR={params.tmp} &> {log}'


    def known_sites(input_known):
        input_known_string = ""
        for known in input_known:
            input_known_string += f' --known-sites {known}'
        return input_known_string


    rule baseRecalibrator:
        input:
            bam = rules.deduplicate.output.bam,
            bam_index = rules.indexMergedBam.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            recal_table = 'gatk/baseRecalibrator/{cell_type}.recal.table'
        params:
            tmp = config['tmpdir'],
            known = known_sites(config['known_sites']),
            extra = ''
        log:
            'logs/gatk/baseRecalibrator/{cell_type}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk BaseRecalibrator {params.extra} {params.known} '
             '--input {input.bam} --reference {input.ref} '
             '--output {output.recal_table} '
             '--sequence-dictionary {input.ref_dict} '
             '--tmp-dir {params.tmp} &> {log}'


    rule applyBQSR:
        input:
            bam = rules.deduplicate.output.bam,
            bam_index = rules.indexMergedBam.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output,
            recal_table = rules.baseRecalibrator.output
        output:
            'mapped/merged_by_cell/{cell_type}.recalibrated.bam'
        params:
            tmp = config['tmpdir'],
            extra = ''
        log:
            'logs/gatk/applyBQSR/{cell_type}.log'
        threads:
            THREADS
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk ApplyBQSR {params.extra} '
            '--input {input.bam} --reference {input.ref} '
            '--bqsr-recal-file {input.recal_table} --output {output} '
            '--tmp-dir {params.tmp} &> {log}'


    rule haplotypeCaller:
        input:
            bam = rules.applyBQSR.output,
            bam_index = rules.indexMergedBam.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'gatk/{cell_type}-g.vcf.gz'
        params:
            tmp = config['tmpdir'],
            intervals = config['protocol']['regions'],
            java_opts = '-Xmx6G',
            extra = ''
        log:
            'logs/gatk/haplotypeCaller/{cell_type}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk --java-options {params.java_opts} HaplotypeCaller '
            '{params.extra} --input {input.bam} --output {output} '
            '--reference {input.ref} --intervals {params.intervals} '
            '--tmp-dir {params.tmp} -ERC GVCF &> {log}'


    # Possibly should combine gvcfs
    rule genotypeGVCFs:
        input:
            gvcf = rules.haplotypeCaller.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'gatk/{cell_type}.vcf.gz'
        params:
            tmp = config['tmpdir'],
            java_opts = '-Xmx4G',
            extra = '',  # optional
        log:
            'logs/gatk/genotypeGVCFs/{cell_type}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk --java-options {params.java_opts} GenotypeGVCFs '
             '--reference {input.ref} --variant {input.gvcf} '
             '--tmp-dir {params.tmp} --output {output} &> {log}'


    rule selectVariants:
        input:
            vcf = rules.genotypeGVCFs.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'gatk/{cell_type}-{mode}.vcf.gz'
        params:
            tmp = config['tmpdir']
        log:
            'logs/gatk/selectVariants/{cell_type}-{mode}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk SelectVariants --reference {input.ref} '
            '--variant {input.vcf} --select-type-to-include {wildcards.mode} '
            '--tmp-dir {params.tmp} --output {output} &> {log}'


    rule variantRecalibratorSNPs:
        input:
            vcf = 'gatk/{cell_type}-SNP.vcf.gz',
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            recal = 'gatk/{cell_type}-SNP.vcf.recal',
            tranches = 'gatk/{cell_type}-SNP.vcf.tranches'
        params:
            tmp = config['tmpdir'],
            truth = '--resource:hapmap,known=false,training=true,truth=true,'
            'prior=15.0 /media/stephen/Data/HiC-subsample/data/hapmap_3.3.hg38.vcf.gz',
            training1 = '--resource:omni,known=false,training=true,truth=false,'
            'prior=12.0 /media/stephen/Data/HiC-subsample/data/1000G_omni2.5.hg38.vcf.gz',
            training2 = '--resource:1000G,known=false,training=true,truth=false,'
            'prior=10.0 /media/stephen/Data/HiC-subsample/data/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
            known = '--resource:dbsnp,known=true,training=false,truth=false,'
            'prior=2.0 /media/stephen/Data/HiC-subsample/data/dbsnp_146.hg38.vcf.gz',
            max_gaussians = 3,
            java_opts = '-Xmx4G',
            extra = '',  # optional
        log:
            'logs/gatk/VariantRecalibrator/{cell_type}-SNP.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk --java-options {params.java_opts} VariantRecalibrator '
            '--reference {input.ref} --variant {input.vcf} --mode SNP '
            '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR '
            '--output {output.recal} --tranches-file {output.tranches} '
            '--max-gaussians {params.max_gaussians} '
            '{params.known} {params.truth} {params.training1} '
            '{params.training2} {params.extra} '
            '--tmp-dir {params.tmp} &> {log}'


    rule variantRecalibratorINDELS:
        input:
            vcf = 'gatk/{cell_type}-INDEL.vcf.gz',
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            recal = 'gatk/{cell_type}-INDEL.vcf.recal',
            tranches = 'gatk/{cell_type}-INDEL.vcf.tranches'
        params:
            tmp = config['tmpdir'],
            mills = '--resource:mills,known=false,training=true,truth=true,'
            'prior=12.0 /media/stephen/Data/HiC-subsample/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
            known = '--resource:dbsnp,known=true,training=false,truth=false,'
            'prior=2.0 /media/stephen/Data/HiC-subsample/data/dbsnp_146.hg38.vcf.gz',
            max_gaussians = 3,
            java_opts = '-Xmx4G',
            extra = '',  # optional
        log:
            'logs/gatk/VariantRecalibrator/{cell_type}-INDEL.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk --java-options {params.java_opts} VariantRecalibrator '
            '--reference {input.ref} --variant {input.vcf} --mode INDEL '
            '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR '
            '--output {output.recal} --tranches-file {output.tranches} '
            '--max-gaussians {params.max_gaussians} '
            '{params.known} {params.mills} {params.extra} '
            '--tmp-dir {params.tmp} &> {log}'

    rule applyVQSR:
        input:
            vcf = rules.selectVariants.output,
            tranches = 'gatk/{cell_type}-{mode}.vcf.tranches',
            recal = 'gatk/{cell_type}-{mode}.vcf.recal',
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'gatk/{cell_type}-{mode}.filt.vcf.gz'
        params:
            tmp = config['tmpdir']
        log:
            'logs/gatk/applyVQSR/{cell_type}-{mode}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk ApplyVQSR --reference {input.ref} --variant {input.vcf} '
            '--tranches-file {input.tranches} --mode {wildcards.mode} '
            '--exclude-filtered --recal-file {input.recal} '
            '--truth-sensitivity-filter-level 99.0 '
            '--tmp-dir {params.tmp} --output {output} &> {log}'


    rule mergeVCFs:
        input:
            SNP = 'gatk/{cell_type}-SNP.filt.vcf.gz',
            INDEL = 'gatk/{cell_type}-INDEL.filt.vcf.gz',
        output:
            'gatk/{cell_type}-all.filt.vcf.gz'
        params:
            tmp = config['tmpdir']
        log:
            'logs/picard/merge_vcfs/{cell_type}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard MergeVcfs INPUT={input.SNP} INPUT={input.INDEL} '
            'TMP_DIR={params.tmp} OUTPUT={output} &> {log}'


    rule splitVCFS:
        input:
            rules.mergeVCFs.output
        output:
            'gatk/{cell_type}-all-{region}.filt.vcf'
        params:
            region = REGIONS.index,
            chr = lambda wildcards: REGIONS['chr'][wildcards.region],
            start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
            end = lambda wildcards: REGIONS['end'][wildcards.region]
        log:
            'logs/splitVCFS/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view --regions {params.chr}:{params.start}-{params.end} '
            '{input} > {output} 2> {log}'


    rule extractHAIRS:
        input:
            vcf = rules.splitVCFS.output,
            bam = rules.deduplicate.output.bam
        output:
            'allele/hapcut2/{region}/{cell_type}-{region}.fragments'
        group:
            'hapcut2'
        log:
            'logs/extractHAIRS/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/hapcut2.yaml'
        shell:
            'extractHAIRS --hic 1 --bam {input.bam} '
            '--VCF {input.vcf} --out {output} &> {log}'


    rule hapCut2:
        input:
            fragments = rules.extractHAIRS.output,
            vcf = rules.splitVCFS.output
        output:
            block = 'allele/hapcut2/{region}/{cell_type}-{region}',
            vcf = 'allele/hapcut2/{region}/{cell_type}-{region}.phased.VCF'
        group:
            'hapcut2'
        log:
            'logs/hapCut2/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/hapcut2.yaml'
        shell:
            'HAPCUT2 --hic 1 --fragments {input.fragments} '
            '--VCF {input.vcf} --outvcf 1 --out {output.block} &> {log}'


    rule bgzipPhased:
        input:
            rules.hapCut2.output.vcf
        output:
            f'{rules.hapCut2.output.vcf}.gz'
        log:
            'logs/bgzip_phased/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/tabix.yaml'
        shell:
            'bgzip -c {input} > {output} 2> {log}'


    rule indexPhased:
        input:
            rules.bgzipPhased.output
        output:
            f'{rules.bgzipPhased.output}.tbi'
        log:
            'logs/index_phased/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/tabix.yaml'
        shell:
            'tabix {input} &> {log}'


    rule extractBestBlock:
        input:
            rules.hapCut2.output.block
        output:
            'allele/hapcut2/{region}/{cell_type}-{region}.tsv'
        log:
            'logs/extractBestPhase/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/coreutils.yaml'
        shell:
            '{SCRIPTS}/extract_best_hapcut2.sh {input} > {output} 2> {log}'


    rule extractVCF:
        input:
            block = rules.extractBestBlock.output,
            vcf = rules.bgzipPhased.output,
            vcf_index = rules.indexPhased.output
        output:
            'allele/hapcut2/{region}/{cell_type}-{region}-best.vcf'
        log:
            'logs/extractVCF/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -R {input.block} {input.vcf} '
            '> {output} 2> {log} || touch {output}'


    rule mergeVCFsbyRegion:
        input:
            vcfs = expand(
                'allele/hapcut2/{region}/{{cell_type}}-{region}-best.vcf',
                region = REGIONS.index),
        output:
            'allele/hapcut2/{cell_type}-phased.vcf'
        log:
            'logs/mergeVCFsbyRegion/{cell_type}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools concat {input.vcfs} > {output} 2> {log}'


    rule mpileup:
        input:
            bam = rules.deduplicate.output.bam,
            genome = rules.bgzipGenome.output,
        output:
            pipe('allele/vcfs/{region}/{cell_type}-{region}-mpileup.bcf')
        group:
            'bcftoolsVariants'
        log:
            'logs/mpileup/{cell_type}-{region}.log'
        threads:
            (THREADS - 4) * 0.5
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools mpileup -q 15 --ignore-RG --count-orphans '
            '--max-depth 100000 --output-type u -f {input.genome} '
            '--threads {threads} {input.bam} > {output} 2> {log} '


    rule callVariants:
        input:
            rules.mpileup.output
        output:
            pipe('allele/vcfs/{region}/{cell_type}-{region}-calls.bcf')
        group:
            'bcftoolsVariants'
        log:
            'logs/callVariants/{cell_type}-{region}.log'
        threads:
            (THREADS - 4) * 0.5
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools call --multiallelic-caller --variants-only '
            '--output-type u --threads {threads} '
            '{input} > {output} 2> {log}'


    rule filterVariants:
        input:
            rules.callVariants.output
        output:
            'allele/vcfs/{region}/{cell_type}-{region}-filt.bcf'
        group:
            'bcftoolsVariants'
        log:
            'logs/filterVariants/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -i "%QUAL>=20" --output-type u '
            '{input} > {output} 2> {log}'


    rule bcftoolsStats:
        input:
            rules.filterVariants.output
        output:
            'qc/variant_quality/{cell_type}-{region}-bcftoolsStats.txt'
        group:
            'bcftoolsVariants'
        log:
            'logs/bcftoolsStats/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools stats {input} > {output} 2> {log}'



    rule run_hapcompass:
        input:
            vcf = rules.splitVCFS.output,
            bam = rules.deduplicate.output.bam
        output:
            expand(
                'allele/hapcompass/{{cell_type}}-{{region}}_{ext}',
                ext = ['MWER_solution.txt', 'reduced_representation.vcf',
                       'reduced_representation.sam', 'frags.txt',
                       'phasedSolution.txt', 'reads.sam',])
        params:
            dir = 'allele/hapcompass'
        log:
            'logs/run_hapcompass/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/openjdk.yaml'
        resources:
            mem_mb = 120000
        shell:
            'java -Xmx120g -jar {SCRIPTS}/hapcompass.jar '
            '--bam {input.bam} --vcf {input.vcf} --debug '
            '--output {params.dir}/{wildcards.cell_type}-{wildcards.region} '
            '&> {log} || touch {output} '


    rule extract_best_hapcompass_phasing:
        input:
            'allele/hapcompass/{cell_type}-{region}_MWER_solution.txt'
        output:
            'allele/hapcompass/{cell_type}-{region}-best.txt'
        log:
            'logs/extract_best_hapcompass_phasing/{cell_type}-{region}.log'
        shell:
            '{SCRIPTS}/extract_best_hapcompass.sh {input} > {output} 2> {log}'


    rule reformat_hapcompass:
        input:
            mwer = rules.extract_best_hapcompass_phasing.output,
            vcf = rules.splitVCFS.output
        output:
            'allele/hapcompass/{cell_type}-{region}-best.txt.vcf'
        log:
            'logs/reformat_vcf/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/openjdk.yaml'
        shell:
            'java -jar {SCRIPTS}/hc2vcf.jar {input.mwer} {input.vcf} 2 true '
            '2> {log} || touch {output}'


    rule compress_hapcompass:
        input:
             rules.reformat_hapcompass.output
        output:
            f'{rules.reformat_hapcompass.output}.gz'
        log:
            'logs/compress_hapcompass/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -O z {input} > {output} 2> {log} || touch {output}'


    rule index_hapcompass:
        input:
            rules.compress_hapcompass.output
        output:
            f'{rules.compress_hapcompass.output}.csi'
        log:
            'logs/index_hapcompass/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools index -f {input} 2> {log} || touch {output}'


    rule merge_phased_regions:
        input:
            vcfs = expand(
                'allele/hapcompass/{{cell_type}}-{region}-best.txt.vcf.gz',
                region = REGIONS.index),
            indexes = expand(
                'allele/hapcompass/{{cell_type}}-{region}-best.txt.vcf.gz.csi',
                region = REGIONS.index)
        output:
            'allele/hapcompass/{cell_type}-phased.vcf.gz'
        log:
            'logs/merge_phased_regions/{cell_type}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools concat --naive {input.vcfs} > {output} 2> {log}'


rule multiqc:
    input:
        [expand('qc/fastqc/{sample}-{read}.raw_fastqc.zip',
                sample=ORIGINAL_SAMPLES, read=READS),
         expand('qc/cutadapt/{sample}.cutadapt.txt', sample=ORIGINAL_SAMPLES),
         expand('qc/fastqc/{sample}-{read}.trim_fastqc.zip',
                sample=ORIGINAL_SAMPLES, read=READS),
         expand('qc/hicup/HiCUP_summary_report-{sample}.txt', sample=ORIGINAL_SAMPLES),
         expand('qc/samtools/stats/{sample}.stats.txt', sample=ORIGINAL_SAMPLES),
         expand('qc/samtools/idxstats/{sample}.idxstats.txt', sample=ORIGINAL_SAMPLES),
         expand('qc/samtools/flagstat/{sample}.flagstat.txt', sample=ORIGINAL_SAMPLES),
         expand('qc/bamQC/{sample}', sample=ORIGINAL_SAMPLES),
         expand('qc/hicexplorer/{sample}-{region}.{bin}_QC',
                sample=SAMPLES, region=REGIONS.index, bin=BASE_BIN),
         expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
                sample=ORIGINAL_SAMPLES, read=READS) if config['fastq_screen'] else []]
    output:
        directory('qc/multiqc')
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc --outdir {output} --force '
        '--config {BASE}/config/multiqc_config.yaml {input} '
        '&> {log}'
