#!/usr/bin/env python3

import pandas as pd
from set_config import set_config, load_samples, get_grouping

wildcard_constraints:
    group = '[^-\/g]+',
    sample = '[^-\/g]+-\d+',
    group_AS = '[^-\/]+_g\d+',
    sample_AS = '[^-\/]+_g\d+-\d+',
    rep = '\d+',
    read = 'R[12]',
    bin = '\d+'

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'workdir' :     workflow.basedir,
    'samples' :     '',
    'build' :       'genome',
    'genome' :      '',
    'fastq_screen': None,
    're1' :         'enzyme',
    're1_seq' :     '',
    'regions' :     '',
    'ctcf' :        None,
    'ctcf_orient' : None,
    'genes' :       None,
    'reads' :       ['R1', 'R2'],
    'threads' :     1,
    'arima' :       False,
    're2' :         None,
    're2_seq' :     None
}
config = set_config(config, default_config)

workdir : config['workdir']

BUILD = config['build']
THREADS = config['threads']
RE1 = config['re1']
RE1_SEQ = config['re1_seq']
RE2 = config['re2']
RE2_SEQ = config['re2_seq']
READS = config['reads']
BINS = config['binsize']
BASE_BIN = BINS[0]

# Read path to samples in pandas
samples = load_samples(config['samples'])

# Extract groups and replicates.
SAMPLES, GROUPS = get_grouping(samples)

# Load capture regions
REGIONS = pd.read_table(config['regions'],
    names = ['chr', 'start', 'end', 'region'],
    index_col = 'region',
    dtype = {'start' : int, 'end' : int})
REGIONS['length'] = REGIONS['end'] - REGIONS['start']

rule all:
    input:
        ['qc/multiqc', 'qc/filter_qc/insert_size_frequency.png']


rule bgzip_genome:
    input:
        config['genome']
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


rule bowtie2Build:
    input:
        rules.bgzip_genome.output
    output:
        expand('genome/index/{build}.{n}.bt2',
            build = BUILD,
            n = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
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


rule fastqc:
    input:
        lambda wc: samples.xs(wc.single, level = 2)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.50.4/bio/fastqc'

# Modify FASTQ filename to match {sample}-{read} for multiQC
rule modify_fastqc:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    params:
        name = lambda wc: f'{wc.single}'
    log:
        'logs/modify_fastqc/{single}.raw.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/modify_fastqc.sh {input} {output} {params.name} &> {log}'


rule hicup_truncate:
    input:
        lambda wc: samples.xs(wc.sample, level = 1)['path']
    output:
        truncated = ['fastq/truncated/{sample}-R1.trunc.fastq.gz',
                     'fastq/truncated/{sample}-R2.trunc.fastq.gz'],
        summary = 'qc/hicup/{sample}-truncate-summary.txt'
    params:
        re1_seq = RE1_SEQ
    threads:
        2 if THREADS > 2 else THREADS
    log:
        'logs/hicup_truncate/{sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_truncater.py '
            '--output {output.truncated} '
            '--summary {output.summary} '
            '--re1 {params.re1_seq} '
            '--threads {threads} {input} '
        '&> {log}'


rule cutadapt:
    input:
        #lambda wc: DATA.xs(wc.sample, level = 1)['path']
        rules.hicup_truncate.output.truncated
    output:
        trimmed = ['fastq/trimmed/{sample}-R1.trim.fastq.gz',
                   'fastq/trimmed/{sample}-R2.trim.fastq.gz'],
        qc = 'qc/cutadapt/unmod/{sample}.cutadapt.txt'
    group:
        'cutadapt'
    params:
        others = '--minimum-length 20 --quality-cutoff 20 --discard-trimmed '
                 '--gc-content 46 --overlap 6 --error-rate 0.1'
    log:
        'logs/cutadapt/{sample}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        THREADS
    shell:
        'cutadapt '
            '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA '
            '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT '
            '{params.others} --cores {THREADS} '
            '-o {output.trimmed[0]} -p {output.trimmed[1]} '
            '{input} > {output.qc} '
        '2> {log}'


rule modify_cutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{sample}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/modify_cutadapt/{sample}.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        'awk -v sample={wildcards.sample} '
            '-f {SCRIPTS}/modify_cutadapt.awk {input} > {output} 2> {log}'


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
        "0.50.4/bio/fastq_screen"


rule fastqc_trimmed:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{single}.trim_fastqc.html',
        zip = 'qc/fastqc/{single}.trim_fastqc.zip'
    log:
        'logs/fastqc_trimmed/{single}.log'
    wrapper:
        '0.50.4/bio/fastqc'


rule hicup_digest:
    input:
        '/media/stephen/Data/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    output:
        f'genome/digest/{BUILD}-{RE1}-hicup-digest.txt.gz'
    params:
        genome = BUILD,
        re1 = RE1,
        re1_seq = RE1_SEQ,
        re2 = f'--re2_name {RE2}' if RE2 else '',
        re2_seq = '--re2 {RE2_SEQ}' if RE2_SEQ else '',
        arima = '--arima' if config['arima'] else ''
    log:
        f'logs/hicup_digest/{BUILD}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_digester.py '
            '--output {output} '
            '--genome {params.genome} '
            '--re1 {params.re1_seq} '
            '--re1_name {params.re1} '
            '{params.arima} {params.re2_seq} {params.re2} '
            '{input} '
        '&> {log}'


rule hicup_map:
    input:
        reads = rules.cutadapt.output.trimmed,
        bt2_index = rules.bowtie2Build.output
    output:
        mapped = 'mapped/{sample}.pair.bam',
        summary = 'qc/hicup/{sample}-map-summary.txt'
    threads:
        12
    log:
        'logs/hicup_mapper/{sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_mapper.py '
            '--output {output.mapped} '
            '--summary {output.summary} '
            '--index {input.bt2_index} '
            '--threads 1 {input.reads} '
        '&> {log}'


rule fixmate:
    input:
        rules.hicup_map.output.mapped
    output:
        'mapped/{sample}.fixed.bam'
    threads:
        THREADS
    log:
        'logs/fixmate/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools fixmate -@ {threads} {input} {output} &> {log}'


rule digest:
    input:
        '/media/stephen/Data/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    output:
        f'genome/digest/{BUILD}-{RE1}-pyHiCtools-digest.txt'
    params:
        re1_seq = config['re1_seq']
    log:
        f'logs/digest/{BUILD}-{RE1}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools digest '
            '--restriction {params.re1_seq} '
            '<(zcat -f {input}) > {output} '
        '2> {log}'


rule subsample_reads:
    input:
        rules.fixmate.output
    output:
        pipe('mapped/subsampled/{sample}-subsample.bam')
    group:
        'filter_qc'
    params:
        seed = '42',
        frac = '20'
    log:
        'logs/subsample_reads/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools view -@ {threads} -s {params.seed}.{params.frac} '
            '{input} > {output} '
        '2> {log}'


rule pyHiCTools_process:
    input:
        dedup_nmsort = rules.subsample_reads.output,
        digest = rules.digest.output
    output:
        pipe('mapped/subsampled/{sample}.processed.sam')
    group:
        'filter_qc'
    log:
        'logs/process/{sample}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools process '
            '--digest {input.digest} {input.dedup_nmsort} '
            '> {output} '
        '2> {log}'


rule extract_hic_stats:
    input:
        rules.pyHiCTools_process.output
    output:
        'mapped/subsampled/{sample}-processed.txt'
    group:
        'filter_qc'
    log:
        'logs/extract_hic_stats/{sample}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools extract '
            '--sample {wildcards.sample} {input} > {output} '
        '2> {log}'


rule plot_filter_QC:
    input:
        expand('mapped/subsampled/{sample}-processed.txt', sample = SAMPLES)
    output:
        expand('qc/filter_qc/{fig}',
            fig = ['trans_stats.csv', 'insert_size_frequency.png',
                   'ditag_length.png'])
    params:
        outdir = 'qc/filter_qc'
    log:
        'logs/plot_filter/plot_subsample.log'
    conda:
        f'{ENVS}/ggplot2.yaml'
    shell:
        '{SCRIPTS}/figures/plot_filter.R {params.outdir} {input} '
        '2> {log}'


rule hicup_filter:
    input:
        bam = rules.fixmate.output,
        digest = rules.hicup_digest.output
    output:
        filtered = 'mapped/{sample}.filt.bam',
        summary = 'qc/hicup/{sample}-filter-summary.txt',
        rejects = directory('qc/hicup/{sample}-ditag_rejects')
    params:
        shortest = 150,
        longest = 800
    log:
        'logs/hicup_filter/{sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_filter.py '
            '--output {output.filtered} '
            '--digest {input.digest} '
            '--outdir {output.rejects} '
            '--summary {output.summary} '
            '--shortest {params.shortest} '
            '--longest {params.longest} {input.bam} '
        '&> {log}'


rule hicup_deduplicate:
    input:
        rules.hicup_filter.output.filtered
    output:
        deduped = 'mapped/{sample}.dedup.bam',
        summary = 'qc/hicup/{sample}-deduplicate-summary.txt',
    log:
        'logs/hicup_deduplicate/{sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicup_deduplicator.py '
            '--output {output.deduped} '
            '--summary {output.summary} {input} '
        '&> {log}'


rule merge_hicup_summary:
    input:
        truncater = rules.hicup_truncate.output.summary,
        mapper = rules.hicup_map.output.summary,
        filter = rules.hicup_filter.output.summary,
        deduplicator = rules.hicup_deduplicate.output.summary
    output:
        'qc/hicup/HiCUP_summary_report-{sample}.txt'
    log:
        'logs/merge_hicup_summary/{sample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/combine_summary.py '
            '--truncater {input.truncater} '
            '--mapper {input.mapper}  '
            '--filter {input.filter}  '
            '--deduplicator {input.deduplicator} > {output} '
        '2> {log}'


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
            "{input} > {output} "
        "2> {log}"



rule coordinate_sort:
    input:
        rules.hicup_deduplicate.output.deduped
    output:
        'mapped/{sample}.sort.bam'
    params:
        mem = '1G'
    threads:
        THREADS
    log:
        'logs/coordinate_sort/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -m {params.mem} {input} > {output} 2> {log}'


rule index_bam:
    input:
        rules.coordinate_sort.output
    output:
        f'{rules.coordinate_sort.output}.bai'
    threads:
        THREADS
    log:
        'logs/index_bam/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools index -@ 6 {input} &> {log}'


rule samtools_stats:
    input:
        rules.coordinate_sort.output
    output:
        'qc/samtools/stats/{sample}.stats.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtools_stats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'


rule samtools_idxstats:
    input:
        bam = rules.coordinate_sort.output,
        index = rules.index_bam.output
    output:
        'qc/samtools/idxstats/{sample}.idxstats.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtools_idxstats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input.bam} > {output} 2> {log}'


rule samtools_flagstat:
    input:
        rules.coordinate_sort.output
    output:
        'qc/samtools/flagstat/{sample}.flagstat.txt'
    group:
        'samtools_qc'
    log:
        'logs/samtools_flagstat/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'


rule bamqc:
    input:
        bam = rules.coordinate_sort.output,
        regions = '/home/stephen/phd/scripts/capture_regions.bed'
    output:
        directory('qc/bamqc/{sample}')
    resources:
        mem_mb = 3000
    log:
        'logs/bamqc/{sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap bamqc '
            '-bam {input.bam} -gff {input.regions} -outdir {output} '
            '-nt 6 --outside-stats --java-mem-size={resources.mem_mb}M '
            '--paint-chromosome-limits '
        '&> {log}'



rule split_paired_bam:
    input:
        rules.hicup_deduplicate.output.deduped
    output:
        'mapped/split/{sample}-{read}.hic.bam'
    params:
        read = READS,
        flag = lambda wc : '0x40' if wc.read == READS[0] else '0x80'
    log:
        'logs/split/{sample}-{read}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -@ {threads} -f {params.flag} -b {input} '
            '> {output} '
        '2> {log}'


rule build_base_matrix:
    input:
        expand('mapped/split/{{sample}}-{read}.hic.bam', read = READS)
    output:
        matrix = f'matrices/{{region}}/base/{{sample}}-{{region}}.{BASE_BIN}.h5',
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
        'hicBuildMatrix '
            '--samFiles {input} '
            '--region {params.chr}:{params.start}-{params.end} '
            '--binSize {params.bin} '
            '--outFileName {output.matrix} '
            '--outBam {output.bam} '
            '--QCfolder {output.qc} '
            '--skipDuplicationCheck '
            '--threads {threads} '
        '&> {log} || mkdir -p {output.qc}; touch {output.matrix} {output.bam}'


rule multiqc:
    input:
        expand('qc/fastqc/{sample}-{read}.raw_fastqc.zip',
            sample = SAMPLES, read = READS),
        expand('qc/cutadapt/{sample}.cutadapt.txt', sample=SAMPLES),
        expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
            sample = SAMPLES, read = READS),
        expand('qc/fastqc/{sample}-{read}.trim_fastqc.zip',
            sample = SAMPLES, read = READS),
        expand('qc/hicup/HiCUP_summary_report-{sample}.txt', sample=SAMPLES),
        expand('qc/samtools/stats/{sample}.stats.txt', sample=SAMPLES),
        expand('qc/samtools/idxstats/{sample}.idxstats.txt', sample=SAMPLES),
        expand('qc/samtools/flagstat/{sample}.flagstat.txt', sample=SAMPLES),
        expand('qc/bamqc/{sample}', sample=SAMPLES),
        expand('qc/hicexplorer/{sample}-{region}.{bin}_QC',
            sample=SAMPLES, region=REGIONS.index, bin=BASE_BIN)
    output:
        directory('qc/multiqc')
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc --outdir {output} '
            '--force --config {BASE}/config/multiqc_config.yaml {input} '
        '&> {log}'
