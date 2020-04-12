#!/usr/bin/env python3

import sys
import pandas as pd
from snake_setup import set_config, load_samples, get_grouping, load_regions, load_vcf_paths, get_allele_groupings

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
    'threads':           1,
    'data':
        {'samples':      ''      ,
         'phased_vcf':   None    ,},
    'genome':
        {'build':        'genome',
         'sequence':     ''      ,
         'genes':        None    ,
         'ctcf':         None    ,
         'ctcf_orient':  None    ,},
    'protocol':
        {'regions':      ''      ,
         're1':          'enzyme',
         're1_seq':      ''      ,
         're2':          None    ,
         're2_seq':      None    ,
         'arima':        False   ,},
    'hicup':
        {'shortest' :    150     ,
         'longest' :     850     ,},
    'HiCcompare':
        {'fdr' :         0.01    ,
         'logFC' :       1       ,},
    'binsize':      [5000, 10000],
    'fastq_screen': None,
}

config = set_config(config, default_config)

workdir: config['workdir']
THREADS = config['threads']
BUILD = config['genome']['build']
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
    workdir: config['workdir'] + 'allele_specific'
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

# Need to prevent '.' also!
wildcard_constraints:
    cell_type = r'[^-\/]+',
    pre_group = r'[^-\/g]+',
    pre_sample = r'[^-\/g]+-\d+',
    region = r'[^-\/]+',
    allele = r'[12]',
    rep = r'\d+',
    read = r'R[12]',
    bin = r'\d+',
    mode = r'SNP|INDEL'

if ALLELE_SPECIFIC:
    wildcard_constraints:
        group = r'[^-\/]+_g\d+',
        sample = r'[^-\/]+_g\d+-\d+'
else:
    wildcard_constraints:
        group = r'[^-\/g]+',
        sample = r'[^-\/g]+-\d+'


rule all:
    input:
        ['qc/multiqc', 'qc/multibamqc',
         'qc/filter_qc/insert_size_frequency.png',
         expand('matrices/{region}/{bin}/plots/{all}-{region}-{bin}.png',
                region=REGIONS.index, bin=BINS, all=SAMPLES+list(GROUPS)),
         expand('matrices/{region}/{bin}/{method}/{all}-{region}-{bin}.{ext}',
                region=REGIONS.index, bin=BINS, all=SAMPLES+list(GROUPS),
                method=['norm', 'ice'], ext=['h5', 'gz',]),
         expand('qc/hicrep/{region}-{bin}-hicrep.png',
                region=REGIONS.index, bin=BINS),
         expand('matrices/{region}/{bin}/plots/{region}-{bin}-{group1}-vs-{group2}.png',
                region=REGIONS.index, bin=BINS,
                group1 = list(GROUPS), group2 = list(GROUPS))],
         [expand('qc/variant_quality/{cell_type}-{region}-bcftools_stats.txt',
                region=REGIONS.index, cell_type=list(CELL_TYPES)),
          expand('allele/hapcut2/{region}/{cell_type}-{region}.phased.VCF',
                region=REGIONS.index, cell_type=list(CELL_TYPES)),
          expand('allele/hapcompass/{cell_type}-phased.vcf.gz',
                cell_type=list(CELL_TYPES))] if not ALLELE_SPECIFIC else []


rule mask_genome:
    input:
        genome = config['genome']['sequence'],
        vcf = lambda wc: PHASED_VCFS[wc.cell_type]
    output:
        f'allele/genome/masked/{BUILD}-{{cell_type}}.fa'
    log:
        'logs/mask_genome/{cell_type}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -fullHeader '
        '-fi <(zcat -f {input.genome}) '
        '-bed {input.vcf} -fo {output} 2> {log}'


rule reformat_SNPsplit:
    input:
        vcf = lambda wc: PHASED_VCFS[wc.cell_type]
    output:
        'allele/snpsplit/{cell_type}-snpsplit.txt'
    log:
        'logs/reformat_SNPsplit/{cell_type}.log'
    conda:
        f'{ENVS}/gawk.yaml'
    shell:
        '{SCRIPTS}/reformat_snpsplit.sh {input} > {output} 2> {log}'


rule bgzip_genome:
    input:
        config['genome']['sequence']
    output:
        f'genome/{BUILD}.fa.gz'
    log:
        f'logs/bgzip_genome/{BUILD}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip > {output}) 2> {log}'


rule index_genome:
    input:
        rules.bgzip_genome.output
    output:
        f'{rules.bgzip_genome.output}.fai'
    log:
        'logs/index_genome/index_genome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule get_chrom_sizes:
    input:
        rules.index_genome.output
    output:
        f'genome/chrom_sizes/{BUILD}.chrom.sizes'
    log:
        f'logs/get_chrom_sizes/{BUILD}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


if ALLELE_SPECIFIC:
    rule bowtie2Build:
        input:
            rules.mask_genome.output
        output:
            expand('genome/index/{build}-{{cell_type}}.{n}.bt2',
                   build=BUILD, n=['1', '2', '3', '4', 'rev.1', 'rev.2'])
        params:
            basename = f'genome/index/{BUILD}'
        log:
            f'logs/bowtie2Build/{BUILD}-{{cell_type}}.log'
        conda:
            f'{ENVS}/bowtie2.yaml'
        threads:
            THREADS
        shell:
            'bowtie2-build --threads {threads} {input} {params.basename} &> {log}'
else:
    rule bowtie2Build:
        input:
            rules.bgzip_genome.output
        output:
            expand('genome/index/{build}.{n}.bt2',
                   build=BUILD, n=['1', '2', '3', '4', 'rev.1', 'rev.2'])
        params:
            basename = f'genome/index/{BUILD}'
        log:
            f'logs/bowtie2Build/{BUILD}.log'
        conda:
            f'{ENVS}/bowtie2.yaml'
        threads:
            THREADS
        shell:
            'bowtie2-build --threads {threads} {input} {params.basename} &> {log}'

rule process_gff3:
    input:
        config['genome']['genes']
    output:
        f'genome/{BUILD}-genes.bed'
    log:
        f'logs/process_gff3/{BUILD}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'zcat -f {input} | {SCRIPTS}/process_gff3.py > {output} 2> {log}'


rule fastqc:
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


rule modify_fastqc:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/modify_fastqc/{single}.raw.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/modify_fastqc.sh {input} {output} '
        '{wildcards.single} &> {log}'


rule hicup_truncate:
    input:
        lambda wc: samples.xs(wc.pre_sample, level=2)['path']
    output:
        truncated = ['fastq/truncated/{pre_sample}-R1.trunc.fastq.gz',
                     'fastq/truncated/{pre_sample}-R2.trunc.fastq.gz'],
        summary = 'qc/hicup/{pre_sample}-truncate-summary.txt'
    params:
        re1_seq = RE1_SEQ
    threads:
        2 if THREADS > 2 else THREADS
    log:
        'logs/hicup_truncate/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_truncater.py '
        '--output {output.truncated} '
        '--summary {output.summary} '
        '--re1 {params.re1_seq} '
        '--threads {threads} {input} &> {log}'


rule cutadapt:
    input:
        rules.hicup_truncate.output.truncated
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


rule modify_cutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{pre_sample}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/modify_cutadapt/{pre_sample}.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        'awk -v sample={wildcards.pre_sample} -f {SCRIPTS}/modify_cutadapt.awk '
        '{input} > {output} 2> {log}'


if config['fastq_screen'] is not None:
    rule fastq_screen:
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


rule fastqc_trimmed:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{single}.trim_fastqc.html',
        zip = 'qc/fastqc/{single}.trim_fastqc.zip'
    log:
        'logs/fastqc_trimmed/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule hicup_digest:
    input:
        rules.bgzip_genome.output
    output:
        f'genome/digest/{BUILD}-{RE1}-hicup-digest.txt.gz'
    params:
        genome = BUILD,
        re1 = RE1,
        re1_seq = RE1_SEQ,
        re2 = f'--re2_name {RE2}' if RE2 else '',
        re2_seq = '--re2 {RE2_SEQ}' if RE2_SEQ else '',
        arima = '--arima' if config['protocol']['arima'] else ''
    log:
        f'logs/hicup_digest/{BUILD}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_digester.py '
        '--output {output} --genome {params.genome} '
        '--re1 {params.re1_seq} --re1_name {params.re1} '
        '{params.arima} {params.re2_seq} {params.re2} {input} &> {log}'


def hicup_map_index(wildcards):
    """ Retrieve cell type associated with sample. """

    if ALLELE_SPECIFIC:
        for cell_type, samples in CELL_TYPES.items():
            if wildcards.pre_sample in samples:
                type = cell_type

        return expand('genome/index/{build}-{cell_type}.{n}.bt2',
            build=BUILD, cell_type=type,
            n=['1', '2', '3', '4', 'rev.1', 'rev.2'])
    else:
        return expand('genome/index/{build}.{n}.bt2',
            build=BUILD, n=['1', '2', '3', '4', 'rev.1', 'rev.2'])


rule hicup_map:
    input:
        reads = rules.cutadapt.output.trimmed,
        bt2_index = hicup_map_index
    output:
        mapped = 'mapped/{pre_sample}.pair.bam',
        summary = 'qc/hicup/{pre_sample}-map-summary.txt'
    params:
        basename = f'genome/index/{BUILD}'
    threads:
        12
    log:
        'logs/hicup_mapper/{pre_sample}.log'
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
        rules.bgzip_genome.output
    output:
        f'genome/digest/{BUILD}-{RE1}-pyHiCtools-digest.txt'
    params:
        re1_seq = RE1_SEQ
    log:
        f'logs/digest/{BUILD}-{RE1}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools digest --restriction {params.re1_seq} <(zcat -f {input}) '
        '> {output} 2> {log}'


rule subsample_reads:
    input:
        rules.hicup_map.output.mapped
    output:
        pipe('mapped/subsampled/{pre_sample}-subsample.bam')
    group:
        'filter_qc'
    params:
        seed = '42',
        frac = '20'
    log:
        'logs/subsample_reads/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools view -@ {threads} -s {params.seed}.{params.frac} {input} '
        '> {output} 2> {log}'


rule pyHiCTools_process:
    input:
        dedup_nmsort = rules.subsample_reads.output,
        digest = rules.digest.output
    output:
        pipe('mapped/subsampled/{pre_sample}.processed.sam')
    group:
        'filter_qc'
    log:
        'logs/process/{pre_sample}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools process --digest {input.digest} {input.dedup_nmsort} '
        '> {output} 2> {log}'


rule extract_hic_stats:
    input:
        rules.pyHiCTools_process.output
    output:
        'mapped/subsampled/{pre_sample}-processed.txt'
    group:
        'filter_qc'
    log:
        'logs/extract_hic_stats/{pre_sample}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools extract '
        '--sample {wildcards.pre_sample} {input} '
        '> {output} 2> {log}'


rule plot_filter_QC:
    input:
        expand('mapped/subsampled/{pre_sample}-processed.txt', pre_sample=ORIGINAL_SAMPLES)
    output:
        expand('qc/filter_qc/{fig}',
               fig=['trans_stats.csv', 'insert_size_frequency.png',
                    'ditag_length.png'])
    params:
        outdir = 'qc/filter_qc'
    log:
        'logs/plot_filter/plot_subsample.log'
    conda:
        f'{ENVS}/ggplot2.yaml'
    shell:
        '{SCRIPTS}/figures/plot_filter.R {params.outdir} {input} 2> {log}'


rule hicup_filter:
    input:
        bam = rules.hicup_map.output.mapped,
        digest = rules.hicup_digest.output
    output:
        filtered = 'mapped/{pre_sample}.filt.bam',
        summary = 'qc/hicup/{pre_sample}-filter-summary.txt',
        rejects = directory('qc/hicup/{pre_sample}-ditag_rejects')
    params:
        shortest = config['hicup']['shortest'],
        longest = config['hicup']['longest']
    log:
        'logs/hicup_filter/{pre_sample}.log'
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


rule hicup_deduplicate:
    input:
        rules.hicup_filter.output.filtered
    output:
        deduped = 'mapped/{pre_sample}.dedup.bam',
        summary = 'qc/hicup/{pre_sample}-deduplicate-summary.txt',
    log:
        'logs/hicup_deduplicate/{pre_sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_deduplicator.py '
        '--output {output.deduped} '
        '--summary {output.summary} {input} &> {log}'


rule merge_hicup_summary:
    input:
        truncater = rules.hicup_truncate.output.summary,
        mapper = rules.hicup_map.output.summary,
        filter = rules.hicup_filter.output.summary,
        deduplicator = rules.hicup_deduplicate.output.summary
    output:
        'qc/hicup/HiCUP_summary_report-{pre_sample}.txt'
    log:
        'logs/merge_hicup_summary/{pre_sample}.log'
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
        rules.hicup_deduplicate.output.deduped
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

    return f'allele/snpsplit/{type}-snpsplit.txt'


rule SNPsplit:
    input:
        bam = rules.hicup_deduplicate.output.deduped,
        snps = SNPsplit_input
    output:
        expand('allele/snpsplit/{{pre_sample}}.matepairs.{ext}',
            ext = ['G1_G1.bam', 'G1_G2.bam', 'G1_UA.bam', 'G2_G2.bam',
                   'G2_UA.bam', 'SNPsplit_report.txt', 'SNPsplit_sort.txt',
                   'UA_UA.bam', 'allele_flagged.bam'])
    params:
        outdir = 'allele/snpsplit/'
    log:
        'logs/SNPsplit/SNPsplit-{pre_sample}.log'
    conda:
        f'{ENVS}/snpsplit.yaml'
    shell:
        'SNPsplit {input.bam} --snp_file {input.snps} '
        '--hic --output_dir {params.outdir} &> {log}'


rule merge_SNPsplit:
    input:
        'allele/snpsplit/{pre_group}-{rep}.matepairs.G{allele}_G{allele}.bam',
        'allele/snpsplit/{pre_group}-{rep}.matepairs.G{allele}_UA.bam'
    output:
        'allele/snpsplit/merged/{pre_group}_g{allele}-{rep}.matepairs.bam'
    log:
        'logs/allele/merge_SNPsplit/{pre_group}_g{allele}-{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools merge -n {output} {input} &> {log}'


rule coordinate_sort:
    input:
        rules.hicup_deduplicate.output.deduped
    output:
        'mapped/{pre_sample}.sort.bam'
    params:
        mem = '1G'
    threads:
        THREADS
    log:
        'logs/coordinate_sort/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -m {params.mem} {input} '
        '> {output} 2> {log}'


rule index_bam:
    input:
        rules.coordinate_sort.output
    output:
        f'{rules.coordinate_sort.output}.bai'
    threads:
        THREADS
    log:
        'logs/index_bam/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools index -@ 6 {input} &> {log}'


rule samtools_stats:
    input:
        rules.coordinate_sort.output
    output:
        'qc/samtools/stats/{pre_sample}.stats.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtools_stats/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'


rule samtools_idxstats:
    input:
        bam = rules.coordinate_sort.output,
        index = rules.index_bam.output
    output:
        'qc/samtools/idxstats/{pre_sample}.idxstats.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtools_idxstats/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input.bam} > {output} 2> {log}'


rule samtools_flagstat:
    input:
        rules.coordinate_sort.output
    output:
        'qc/samtools/flagstat/{pre_sample}.flagstat.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtools_flagstat/{pre_sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'


rule bamqc:
    input:
        bam = rules.coordinate_sort.output,
        regions = config['protocol']['regions']
    output:
        directory('qc/bamqc/{pre_sample}')
    resources:
        mem_mb = 3000
    log:
        'logs/bamqc/{pre_sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap bamqc '
        '-bam {input.bam} -gff {input.regions} -outdir {output} '
        '-nt 6 --outside-stats --java-mem-size={resources.mem_mb}M '
        '--paint-chromosome-limits &> {log}'


rule multibamqc_config:
    input:
        expand('qc/bamqc/{pre_sample}', pre_sample=ORIGINAL_SAMPLES)
    output:
        'qc/bamqc/multibamqc.config'
    group:
        'multiBAM_QC'
    log:
        'logs/multibamqc_config/multibamqc_config.log'
    shell:
        '{SCRIPTS}/multibamqc_config.py {input} > {output} 2> {log}'


rule multibamqc:
    input:
        rules.multibamqc_config.output
    output:
        directory('qc/multibamqc')
    group:
        'multiBAM_QC'
    log:
        'logs/multibamqc/multibamqc.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    shell:
        'qualimap multi-bamqc --data {input} -outdir {output} &> {log}'


def split_input(wildcards):
    if ALLELE_SPECIFIC:
        return 'allele/snpsplit/merged/{sample}.matepairs.bam'
    else:
        return f'mapped/{wildcards.sample}.dedup.bam'


rule split_paired_bam:
    input:
        split_input
    output:
        'mapped/split/{sample}-{read}.hic.bam'
    params:
        read = READS,
        flag = lambda wc: '0x40' if wc.read == READS[0] else '0x80'
    log:
        'logs/split/{sample}-{read}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -@ {threads} -f {params.flag} -b {input} '
        '> {output} 2> {log}'


rule build_base_matrix:
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
        'logs/build_base_matrix/{sample}-{region}.log'
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


rule normalise_matrices:
    input:
        expand('matrices/{{region}}/base/raw/{sample}-{{region}}.{bin}.h5',
               bin=BASE_BIN, sample=SAMPLES)
    output:
        expand('matrices/{{region}}/base/norm/{sample}-{{region}}-{bin}.h5',
               bin=BASE_BIN, sample=SAMPLES)
    params:
        method = 'smallest'
    log:
        'logs/normalise_matrices/{region}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicNormalize --matrices {input} '
        '--outFileName {output} --normalize {params.method} '
        '&> {log} || touch {output} '


rule merge_replicates:
    input:
        lambda wildcards: expand(
            'matrices/{{region}}/base/norm/{group}-{rep}-{{region}}-{bin}.h5',
            bin=BASE_BIN, group=wildcards.group,
            rep=GROUPS[wildcards.group])
    output:
        f'matrices/{{region}}/base/norm/{{group}}-{{region}}-{BASE_BIN}.h5'
    log:
        'logs/sum_matrix/{group}-{region}.h5'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicSumMatrices --matrices {input} --outFileName {output} '
        '&> {log} || touch {output}'


rule merge_bins:
    input:
        f'matrices/{{region}}/base/norm/{{all}}-{{region}}-{BASE_BIN}.h5'
    output:
        'matrices/{region}/{bin}/norm/{all}-{region}-{bin}.h5'
    params:
        bin = config['binsize'],
        nbins = lambda wildcards: int(int(wildcards.bin) / BASE_BIN)
    log:
        'logs/merge_matrix/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log} || touch {output}'


rule ICE_correction:
    input:
        rules.merge_bins.output
    output:
        plot = 'qc/ICE_correction/{all}-{region}-{bin}-diagnosic_plot.png',
        matrix = 'matrices/{region}/{bin}/ice/{all}-{region}-{bin}.h5'
    params:
        iternum = 1000,
        upper_threshold = 5
    log:
        'logs/ICE_correction/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '{SCRIPTS}/hic_correct.sh -p {output.plot} '
        '-o {output.matrix} -u {params.upper_threshold} '
        '-i {params.iternum} -@ 1 {input} &> {log} || touch {output}'


rule distance_normalise:
    input:
        rules.ICE_correction.output.matrix
    output:
        'matrices/{region}/{bin}/ice/obs_exp/{all}-{region}-{bin}.h5'
    params:
        method = 'obs_exp'
    log:
        'logs/distance_normalise/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicTransform -m {input} --method {params.method} -o {output} '
        '&> {log} || touch {output}'


rule plot_hic:
    input:
        rules.distance_normalise.output
    output:
        'matrices/{region}/{bin}/plots/{all}-{region}-{bin}.png'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        title = '{all}-{region}-{bin}-norm-ice-obs_exp',
        dpi = 600,
        colour = 'YlGn'
    log:
        'logs/plot_hic/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicPlotMatrix --matrix {input} '
        '--outFileName {output} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--colorMap {params.colour} '
        '--title {params.title} '
        '--vMin 0 --vMax 2 --dpi {params.dpi} &> {log}'


rule compute_tad_insulation:
    input:
        rules.ICE_correction.output
    output:
        expand(
            'matrices/{{region}}/{{bin}}/tads/{{all}}-{{region}}-{{bin}}_{ext}',
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
        'logs/compute_tad_insulation/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicFindTADs --matrix {input} '
        '--minDepth {params.min_depth} --maxDepth {params.max_depth} '
        '--step {wildcards.bin} --outPrefix {params.prefix} '
        '--correctForMultipleTesting {params.method} '
        '--numberOfProcessors {threads} &> {log} || touch {output}'


rule detect_loops:
    input:
        rules.ICE_correction.output.matrix
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
        'logs/detect_loops/{all}-{region}-{bin}.log'
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


rule reformat_homer:
    input:
        'matrices/{region}/{bin}/{method}/{all}-{region}-{bin}.h5'
    output:
        'matrices/{region}/{bin}/{method}/{all}-{region}-{bin}.gz'
    log:
        'logs/reformat_homer/{all}-{region}-{bin}-{method}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat h5 --outputFormat homer &> {log} || touch {output}'


rule reformat_nxn3p:
    input:
        'matrices/{region}/{bin}/norm/{sample}-{region}-{bin}.gz'
    output:
        'matrices/{region}/{bin}/norm/{sample}-{region}-{bin}.nxnp3.tsv'
    params:
        region = REGIONS.index,
    log:
        'logs/reformat_hicrep/{sample}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/reformat_hicrep.sh '
        '-r {wildcards.region} -b {wildcards.bin} {input}'
        '> {output} 2> {log}'


rule run_hicrep:
    input:
        expand(
            'matrices/{{region}}/{{bin}}/norm/{sample}-{{region}}-{{bin}}.nxnp3.tsv',
            sample = SAMPLES)
    output:
        'qc/hicrep/{region}-{bin}-hicrep.png'
    params:
        bin = BINS,
        region = REGIONS.index,
        start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
        end = lambda wildcards: REGIONS['end'][wildcards.region]
    log:
        'logs/run_hicrep/{region}-{bin}.log'
    conda:
        f'{ENVS}/hicrep.yaml'
    shell:
        '{SCRIPTS}/figures/plot_hicrep.R {output} {wildcards.bin} '
        '{params.start} {params.end} {input} &> {log}'


rule reformat_nxn:
    input:
        'matrices/{region}/{bin}/ice/{all}-{region}-{bin}.gz'
    output:
        'matrices/{region}/{bin}/ice/{all}-{region}-{bin}.nxn.tsv'
    log:
        'logs/reformat_nxn/{all}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/reformat_nxn.sh {input} > {output} 2> {log} || touch {output}'


rule run_OnTAD:
    input:
        rules.reformat_nxn.output
    output:
        bed = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.bed',
        tad = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.tad'
    params:
        bin = BINS,
        region = REGIONS.index,
        chr = lambda wildcards: REGIONS['chr'][wildcards.region],
        length = lambda wildcards: REGIONS['length'][wildcards.region],
        outprefix = 'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad'
    log:
        'logs/run_OnTAD/{all}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/OnTAD {input} -o {params.outprefix} -bedout chr{params.chr} '
        '{params.length} {wildcards.bin} &> {log} || touch {output}'


rule reformat_OnTAD:
    input:
        rules.run_OnTAD.output.bed
    output:
        'matrices/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.links'
    params:
        region = REGIONS.index,
        start = lambda wildcards: REGIONS['start'][wildcards.region]
    log:
        'logs/reformat_OnTAD/{all}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/reformat_ontad.sh {params.start} {input} '
        '> {output} 2> {log} || touch {output}'


rule generate_config:
    input:
        matrix = 'matrices/{region}/{bin}/ice/obs_exp/{group}-{region}-{bin}.h5',
        loops = 'matrices/{region}/{bin}/loops/{group}-{region}-{bin}.bedgraph',
        insulations = 'matrices/{region}/{bin}/tads/{group}-{region}-{bin}_tad_score.bm',
        tads = 'matrices/{region}/{bin}/tads/{group}-{region}-{bin}-ontad.links',
        genes = rules.process_gff3.output
    output:
        'matrices/{region}/{bin}/plots/configs/{group}-{region}-{bin}.ini'
    conda:
        f'{ENVS}/python3.yaml'
    params:
        ctcf_orientation = config['genome']['ctcf_orient'],
        ctcf = config['genome']['ctcf'],
        depth = lambda wc: int(REGIONS['length'][wc.region] / 1.5),
    log:
        'logs/generate_config/{group}-{region}-{bin}.log'
    shell:
        '{SCRIPTS}/generate_config.py '
        '-i {input.insulations} '
        '-l {input.loops} '
        '-t {input.tads} '
        '-m {input.matrix} '
        '-s {wildcards.group} '
        '-c {params.ctcf} '
        '-r {params.ctcf_orientation} '
        '-g {input.genes} '
        '-d {params.depth} > {output} 2> {log}'


rule merge_bam_replicate:
    input:
        lambda wildcards: expand(
            'matrices/{{region}}/{group}-{rep}-{{region}}.bam',
            group = wildcards.group, rep = GROUPS[wildcards.group])
    output:
        'matrices/{region}/{group}-{region}.bam'
    log:
        'logs/merge_bam_replicate/{group}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {input} 2> {log}'


rule bam_to_pre:
    input:
        'matrices/{region}/{all}-{region}.bam'
    output:
        'matrices/{region}/base/raw/{all}-{region}.pre.tsv'
    log:
        'logs/bam_to_pre/{all}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        '(samtools view -@ {threads} {input} '
        '| awk -f {SCRIPTS}/bam_to_pre.awk > {output}) 2> {log} '


rule juicer_pre:
    input:
        tsv = 'matrices/{region}/base/raw/{all}-{region}.pre.tsv',
        chrom_sizes = rules.get_chrom_sizes.output
    output:
        'matrices/{region}/{all}-{region}.hic'
    log:
        'logs/juicer_pre/{all}-{region}.log'
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
        rules.juicer_pre.output
    output:
        'matrices/{region}/{bin}/{all}-{region}-{bin}-sutm.txt'
    log:
        'logs/straw/{all}-{region}-{bin}.log'
    params:
        # Strip 'chr' as juicer removes by default
        chr = lambda wildcards: REGIONS['chr'][wildcards.region].lstrip('chr'),
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
        dir = lambda wildcards: f'HiCcompare/{wildcards.region}/{wildcards.bin}',
        chr = lambda wildcards: REGIONS['chr'][wildcards.region]
    conda:
        f'{ENVS}/HiCcompare.yaml'
    shell:
        '({SCRIPTS}/hicCompare.R '
        'HiCcompare/{wildcards.region}/{wildcards.bin} '
        '{params.chr} {wildcards.bin} {input} || touch {output}) &> {log}'


rule links2interact:
    input:
        'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.links'
    output:
        up = 'HiCcompare/{region}/{bin}/{group1}-vs-{group2}-up.interact',
        down = 'HiCcompare/{region}/{bin}/{group1}-vs-{group2}-down.interact'
    group:
        'HiCcompare'
    log:
        'logs/links2interact/{region}-{bin}-{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/links2interact.py '
        '--up {output.up} --down {output.down} {input} &> {log}'


rule filter_links:
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
        'logs/filter_links/{region}-{bin}-{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filter_links.py '
        '--up {output.up} --down {output.down} '
        '--p_value {params.p_value} --log_fc {params.log_fc} {input} &> {log}'


rule merge_config:
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
        'logs/merge_config/{region}-{bin}.log'
    params:
        prefix = lambda wc: f'matrices/{wc.region}/{wc.bin}/plots/configs/{wc.region}-{wc.bin}-'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/merge_config_hic_compare.py '
            '--up {input.links_up} '
            '--down {input.links_down} '
            '--configs {input.configs} '
            '--prefix {params.prefix} &> {log}'


rule plot_tads:
    input:
        'matrices/{region}/{bin}/plots/configs/{region}-{bin}-{group1}-vs-{group2}.ini'
    output:
        'matrices/{region}/{bin}/plots/{region}-{bin}-{group1}-vs-{group2}.png'
    params:
        chr = lambda wildcards: REGIONS['chr'][wildcards.region],
        start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
        end = lambda wildcards: REGIONS['end'][wildcards.region],
        dpi = 600
    log:
        'logs/plot_tads/{region}-{bin}-{group1}-vs-{group2}-plot_tads.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'pyGenomeTracks --tracks {input} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--outFileName {output} '
        '--title "{wildcards.region} at {wildcards.bin} bin size" '
        '--dpi {params.dpi} &> {log}'


if not ALLELE_SPECIFIC:

    # Need to add IndexFeature for each input feature provided in config
    rule merge_bam_cell_type_gatk:
        input:
            lambda wildcards: expand('mapped/{pre_sample}.pair.bam',
                pre_sample = CELL_TYPES[wildcards.cell_type])
        output:
            'mapped/merged_by_cell/{cell_type}.bam'
        log:
            'logs/merge_bam_cell_type_gatk/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            THREADS
        shell:
            'samtools merge -@ {threads} {output} {input} 2> {log}'


    rule coordinate_sort_gatk:
        input:
            rules.merge_bam_cell_type_gatk.output
        output:
            'mapped/merged_by_cell/{cell_type}.sort.bam'
        params:
            mem = '1G'
        threads:
            THREADS
        log:
            'logs/coordinate_sort_gatk/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools sort -@ {threads} -m {params.mem} {input} '
            '> {output} 2> {log}'


    rule MarkDuplicates:
        input:
            rules.coordinate_sort_gatk.output
        output:
            bam = 'mapped/merged_by_cell/{cell_type}.dedup.bam',
            metrics = 'qc/picard_dedup/{cell_type}.metrics.txt'
        log:
            'logs/picard/MarkDuplicates/{cell_type}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard MarkDuplicates REMOVE_DUPLICATES=true '
            'INPUT={input} OUTPUT={output.bam} '
            'METRICS_FILE={output.metrics} &> {log}'


    rule AddReadGroup:
        input:
            rules.MarkDuplicates.output.bam
        output:
            'mapped/merged_by_cell/{cell_type}.dedup-RG.bam'
        params:
            lib = 'lib',
            platform = 'platform',
            unit = 'unit'
        log:
            'logs/picard/AddReadGroup/{cell_type}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard AddOrReplaceReadGroups INPUT={input} '
            'OUTPUT={output} RGLB={params.lib} '
            'RGPL={params.platform} RGPU={params.unit} '
            'RGSM={wildcards.cell_type} &> {log}'


    rule index_gatk:
        input:
            rules.AddReadGroup.output
        output:
            f'{rules.AddReadGroup.output}.bai'
        threads:
            THREADS
        log:
            'logs/samtools/index/{cell_type}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools index -@ {threads} {input} &> {log}'


    rule CreateSequenceDictionary:
        input:
            rules.bgzip_genome.output
        output:
            f'{rules.bgzip_genome.output}.dict'
        log:
            f'logs/gatk/CreateSequenceDictionary/{BUILD}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard CreateSequenceDictionary R={input} O={output} 2> {log}'


    def known_sites(input_known):
        input_known_string = ""
        for known in input_known:
            input_known_string += f' --known-sites {known}'
        return input_known_string


    rule BaseRecalibrator:
        input:
            bam = rules.AddReadGroup.output,
            bam_index = rules.index_gatk.output,
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            recal_table = 'gatk/recalibation/{cell_type}.recal.table'
        params:
            intervals = config['protocol']['regions'],
            known = known_sites(config['known_sites']),
            extra = ''
        log:
            'logs/gatk/BaseRecalibrator/{cell_type}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk BaseRecalibrator {params.extra} {params.known} '
             '--input {input.bam} --reference {input.ref} '
             '--output {output.recal_table} --intervals {params.intervals} '
             '--sequence-dictionary {input.ref_dict} &> {log}'


    rule ApplyBQSR:
        input:
            bam = rules.AddReadGroup.output,
            bam_index = rules.index_gatk.output,
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output,
            recal_table = rules.BaseRecalibrator.output
        output:
            'mapped/merged_by_cell/{cell_type}.recalibrated.bam'
        params:
            intervals = config['protocol']['regions'],
            extra = ''
        log:
            'logs/gatk/ApplyBQSR/{cell_type}.log'
        threads:
            THREADS
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk ApplyBQSR {params.extra} '
            '--input {input.bam} --reference {input.ref} '
            '--bqsr-recal-file {input.recal_table} --output {output} '
            '--intervals {params.intervals} &> {log}'


    rule HaplotypeCaller:
        input:
            bam = rules.ApplyBQSR.output,
            bam_index = rules.index_gatk.output,
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            'gatk/{cell_type}-g.vcf.gz'
        params:
            intervals = config['protocol']['regions'],
            java_opts = '-Xmx6G',
            extra = ''
        log:
            'logs/gatk/HaplotypeCaller/{cell_type}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk --java-options {params.java_opts} HaplotypeCaller '
            '{params.extra} --input {input.bam} --output {output} '
            '--reference {input.ref}  --intervals {params.intervals} '
            '-ERC GVCF &> {log}'


    # Possibly should combine gvcfs
    rule GenotypeGVCFs:
        input:
            gvcf = rules.HaplotypeCaller.output,
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            'gatk/{cell_type}.vcf.gz'
        params:
            java_opts = '-Xmx4G',
            extra = '',  # optional
        log:
            'logs/gatk/GenotypeGVCFs/{cell_type}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk --java-options {params.java_opts} GenotypeGVCFs '
             '--reference {input.ref} --variant {input.gvcf} '
             '--output {output} &> {log}'


    rule SelectVariants:
        input:
            vcf = rules.GenotypeGVCFs.output,
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            'gatk/{cell_type}-{mode}.vcf.gz'
        log:
            'logs/gatk/SelectVariants/{cell_type}-{mode}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk SelectVariants --reference {input.ref} '
            '--variant {input.vcf} --select-type-to-include {wildcards.mode} '
            '--output {output} &> {log}'


    rule VariantRecalibrator_SNPs:
        input:
            vcf = 'gatk/{cell_type}-SNP.vcf.gz',
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            recal = 'gatk/{cell_type}-SNP.vcf.recal',
            tranches = 'gatk/{cell_type}-SNP.vcf.tranches'
        params:
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
            '{params.training2} {params.extra} &> {log}'


    rule VariantRecalibrator_INDELs:
        input:
            vcf = 'gatk/{cell_type}-INDEL.vcf.gz',
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            recal = 'gatk/{cell_type}-INDEL.vcf.recal',
            tranches = 'gatk/{cell_type}-INDEL.vcf.tranches'
        params:
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
            '{params.known} {params.mills} {params.extra} &> {log}'

    rule ApplyVQSR:
        input:
            vcf = rules.SelectVariants.output,
            tranches = 'gatk/{cell_type}-{mode}.vcf.tranches',
            recal = 'gatk/{cell_type}-{mode}.vcf.recal',
            ref = rules.bgzip_genome.output,
            ref_index = rules.index_genome.output,
            ref_dict = rules.CreateSequenceDictionary.output
        output:
            'gatk/{cell_type}-{mode}.filt.vcf.gz'
        log:
            'logs/gatk/ApplyVQSR/{cell_type}-{mode}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk ApplyVQSR --reference {input.ref} --variant {input.vcf} '
            '--tranches-file {input.tranches} --mode {wildcards.mode} '
            '--exclude-filtered --recal-file {input.recal} '
            '--truth-sensitivity-filter-level 99.0 '
            '--output {output} &> {log}'


    rule MergeVCFs:
        input:
            SNP = 'gatk/{cell_type}-SNP.filt.vcf.gz',
            INDEL = 'gatk/{cell_type}-INDEL.filt.vcf.gz',
        output:
            'gatk/{cell_type}-all.filt.vcf.gz'
        log:
            'logs/picard/merge_vcfs/{cell_type}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard MergeVcfs INPUT={input.SNP} INPUT={input.INDEL} '
            'OUTPUT={output} &> {log}'


    rule SplitVCFS:
        input:
            rules.MergeVCFs.output
        output:
            'gatk/{cell_type}-all-{region}.filt.vcf'
        params:
            region = REGIONS.index,
            chr = lambda wildcards: REGIONS['chr'][wildcards.region],
            start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
            end = lambda wildcards: REGIONS['end'][wildcards.region]
        log:
            'logs/SplitVCFS/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view --regions {params.chr}:{params.start}-{params.end} '
            '{input} > {output} 2> {log}'


    rule merge_bam_cell_type:
        input:
            lambda wildcards: expand(
                'matrices/{{region}}/{sample}-{{region}}.bam',
                sample = CELL_TYPES[wildcards.cell_type]),
        output:
            'matrices/{region}/merged_by_cell/{cell_type}-{region}.bam'
        log:
            'logs/merge_bam_cell_type/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            THREADS
        shell:
            'samtools merge -@ {threads} {output} {input} 2> {log}'


    rule sort:
        input:
            rules.merge_bam_cell_type.output
        output:
            'matrices/{region}/merged_by_cell/{cell_type}-{region}.sort.bam'
        params:
            mem = '1G'
        threads:
            THREADS
        log:
            'logs/sort/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools sort -@ {threads} -m {params.mem} {input} > {output} 2> {log}'


    rule mpileup:
        input:
            bam = rules.sort.output,
            genome = rules.bgzip_genome.output,
        output:
            pipe('allele/vcfs/{region}/{cell_type}-{region}-mpileup.bcf')
        group:
            'variant_calling'
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


    rule call_variants:
        input:
            rules.mpileup.output
        output:
            pipe('allele/vcfs/{region}/{cell_type}-{region}-calls.bcf')
        group:
            'variant_calling'
        log:
            'logs/call_variants/{cell_type}-{region}.log'
        threads:
            (THREADS - 4) * 0.5
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools call --multiallelic-caller --variants-only '
            '--output-type u --threads {threads} '
            '{input} > {output} 2> {log}'


    rule filter_variants:
        input:
            rules.call_variants.output
        output:
            'allele/vcfs/{region}/{cell_type}-{region}-filt.bcf'
        group:
            'variant_calling'
        log:
            'logs/filter_variants/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -i "%QUAL>=20" --output-type u '
            '{input} > {output} 2> {log}'


    rule bcftools_stats:
        input:
            rules.filter_variants.output
        output:
            'qc/variant_quality/{cell_type}-{region}-bcftools_stats.txt'
        group:
            'variant_calling'
        log:
            'logs/bcftools_stats/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools stats {input} > {output} 2> {log}'


    rule sort_vcf:
        input:
            rules.filter_variants.output
        output:
            'allele/vcfs/{region}/{cell_type}-{region}.sorted.vcf'
        group:
            'variant_calling'
        log:
            'logs/sort_vcf/{cell_type}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools sort --output-type v {input} > {output} 2> {log}'

    rule extractHAIRS:
        input:
            vcf = rules.SplitVCFS.output,
            bam = rules.AddReadGroup.output
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
            #vcf = rules.sort_vcf.output
            vcf = rules.SplitVCFS.output
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


    rule run_hapcompass:
        input:
            vcf = rules.SplitVCFS.output,
            bam = rules.AddReadGroup.output
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
            vcf = rules.sort_vcf.output
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
         expand('qc/bamqc/{sample}', sample=ORIGINAL_SAMPLES),
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
