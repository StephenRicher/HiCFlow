#!/usr/bin/env python3

import pandas as pd
from snake_setup import set_config, load_samples, get_grouping, load_regions

wildcard_constraints:
    group = r'[^-\/g]+',
    sample = r'[^-\/g]+-\d+',
    group_AS = r'[^-\/]+_g\d+',
    sample_AS = r'[^-\/]+_g\d+-\d+',
    rep = r'\d+',
    read = r'R[12]',
    bin = r'\d+'

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'workdir':      workflow.basedir,
    'samples':      '',
    'build':        'genome',
    'genome':       '',
    'fastq_screen': None,
    're1':          'enzyme',
    're1_seq':      '',
    'regions':      '',
    'ctcf':         None,
    'ctcf_orient':  None,
    'genes':        None,
    'reads':        ['R1', 'R2'],
    'threads':      1,
    'arima':        False,
    're2':          None,
    're2_seq':      None
}
config = set_config(config, default_config)

workdir: config['workdir']
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

REGIONS = load_regions(config['regions'])

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
               group1 = list(GROUPS), group2 = list(GROUPS)),
        expand('allele/vcfs/{region}/{group}-{region}-filt.bcf',
               region=REGIONS.index, group=list(GROUPS))
        ]


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


rule get_chrom_sizes:
    input:
        rules.index_genome.output
    output:
        f'genome/chrom_sizes/{BUILD}.chrom.sizes'
    log:
        f'logs/get_chrom_sizes/{BUILD}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


rule process_gff3:
    input:
        config['genes']
    output:
        f'genome/{BUILD}-genes.bed'
    log:
        f'logs/process_gff3/{BUILD}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'zcat -f {input} | {SCRIPTS}/process_gff3.py > {output} 2> {log}'


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


rule fastqc:
    input:
        lambda wc: samples.xs(wc.single, level=2)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.50.4/bio/fastqc'


rule modify_fastqc:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    log:
        'logs/modify_fastqc/{single}.raw.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/modify_fastqc.sh {input} {output} '
        '{wildcards.single} &> {log}'


rule hicup_truncate:
    input:
        lambda wc: samples.xs(wc.sample, level=1)['path']
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
        '--threads {threads} {input} &> {log}'


rule cutadapt:
    input:
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
        '-o {output.trimmed[0]} -p {output.trimmed[1]} {input} '
        '> {output.qc} 2> {log}'


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
        'awk -v sample={wildcards.sample} -f {SCRIPTS}/modify_cutadapt.awk '
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
        rules.bgzip_genome.output
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
        '--output {output} --genome {params.genome} '
        '--re1 {params.re1_seq} --re1_name {params.re1} '
        '{params.arima} {params.re2_seq} {params.re2} {input} &> {log}'


rule hicup_map:
    input:
        reads = rules.cutadapt.output.trimmed,
        bt2_index = rules.bowtie2Build.output
    output:
        mapped = 'mapped/{sample}.pair.bam',
        summary = 'qc/hicup/{sample}-map-summary.txt'
    params:
        basename = f'genome/index/{BUILD}'
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
        '--index {params.basename} '
        '--threads 1 {input.reads} &> {log}'


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
        rules.bgzip_genome.output
    output:
        f'genome/digest/{BUILD}-{RE1}-pyHiCtools-digest.txt'
    params:
        re1_seq = config['re1_seq']
    log:
        f'logs/digest/{BUILD}-{RE1}.log'
    conda:
        f'{ENVS}/pyHiCTools.yaml'
    shell:
        'pyHiCTools digest --restriction {params.re1_seq} <(zcat -f {input}) '
        '> {output} 2> {log}'


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
        'samtools view -@ {threads} -s {params.seed}.{params.frac} {input} '
        '> {output} 2> {log}'


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
        'pyHiCTools process --digest {input.digest} {input.dedup_nmsort} '
        '> {output} 2> {log}'


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
        '--sample {wildcards.sample} {input} '
        '> {output} 2> {log}'


rule plot_filter_QC:
    input:
        expand('mapped/subsampled/{sample}-processed.txt', sample=SAMPLES)
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
        '--longest {params.longest} {input.bam} &> {log}'


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
        '--summary {output.summary} {input} &> {log}'


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
        '--paint-chromosome-limits &> {log}'


rule multibamqc_config:
    input:
        expand('qc/bamqc/{sample}', sample=SAMPLES)
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


rule split_paired_bam:
    input:
        rules.hicup_deduplicate.output.deduped
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
        ctcf_orientation = config['ctcf_orientation'],
        ctcf = config['ctcf'],
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


rule merge_replicate_bams:
    input:
        lambda wildcards: expand(
            'matrices/{{region}}/{group}-{rep}-{{region}}.bam',
            group = wildcards.group, rep = GROUPS[wildcards.group])
    output:
        'matrices/{region}/{group}-{region}.bam'
    log:
        'logs/merge_replicate_bams/{group}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {input} 2> {log}'


rule index:
    input:
        rules.merge_replicate_bams.output
    output:
        f'{rules.merge_replicate_bams.output}.bai'
    log:
        'logs/index_split_bam/{group}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools index -@ {threads} {input} 2> {log}'


rule bam_to_pre:
    input:
        'matrices/{region}/{all}-{{region}}.bam'
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
        chr = lambda wildcards: REGIONS['chr'][wildcards.region],
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


p_thresh = 0.05
fc_thresh = 0
rule filter_links:
    input:
        'HiCcompare/{region}/{bin}/{group1}-vs-{group2}.links'
    output:
        up = f'HiCcompare/{{region}}/{{bin}}/{{group1}}-vs-{{group2}}-p{p_thresh}-fc{fc_thresh}-up.links',
        down = f'HiCcompare/{{region}}/{{bin}}/{{group1}}-vs-{{group2}}-p{p_thresh}-fc{fc_thresh}-down.links',
    group:
        'HiCcompare'
    log:
        'logs/filter_links/{region}-{bin}-{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filter_links.py '
        '--up {output.up} --down {output.down} '
        '--p_value {p_thresh} --log_fc {fc_thresh} {input} &> {log}'


rule merge_config:
    input:
        links_up = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-p{p_thresh}-fc{fc_thresh}-up.links',
            p_thresh = p_thresh, fc_thresh = fc_thresh,
            group1 = list(GROUPS), group2 = list(GROUPS)),
        links_down = expand(
            'HiCcompare/{{region}}/{{bin}}/{group1}-vs-{group2}-p{p_thresh}-fc{fc_thresh}-down.links',
            p_thresh = p_thresh, fc_thresh = fc_thresh,
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


rule mpileup:
    input:
        bam_index = rules.index.output,
        merged_bam = rules.merge_replicate_bams.output,
        genome = rules.bgzip_genome.output,
        genome_index = rules.index_genome.output
    output:
        pipe('allele/vcfs/{region}/{group}-{region}-mpileup.bcf')
    params:
        chr = lambda wildcards: REGIONS['chr'][wildcards.region],
        start = lambda wildcards: REGIONS['start'][wildcards.region] + 1,
        end = lambda wildcards: REGIONS['end'][wildcards.region]
    group:
        'variant_calling'
    log:
        'logs/mpileup/{group}-{region}.log'
    threads:
        (THREADS - 4) * 0.5
    conda:
        f'{ENVS}/bcftools.yaml'
    shell:
        'bcftools mpileup -q 15 --ignore-RG --count-orphans '
        '--max-depth 100000 --output-type u -f {input.genome} '
        '--regions {params.chr}:{params.start}-{params.end} '
        '--threads {threads} {input.merged_bam} > {output} 2> {log} '


rule call_variants:
    input:
        rules.mpileup.output
    output:
        'allele/vcfs/{region}/{group}-{region}-calls.bcf'
    group:
        'variant_calling'
    log:
        'logs/call_variants/{group}-{region}.log'
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
        'allele/vcfs/{region}/{group}-{region}-filt.bcf'
    group:
        'variant_calling'
    log:
        'logs/filter_variants/{group}-{region}.log'
    conda:
        f'{ENVS}/bcftools.yaml'
    shell:
        'bcftools view -i "%QUAL>=20" --output-type u '
        '{input} > {output} 2> {log}'


multiqc_input = (
     [expand('qc/fastqc/{sample}-{read}.raw_fastqc.zip',
             sample=SAMPLES, read=READS),
      expand('qc/cutadapt/{sample}.cutadapt.txt', sample=SAMPLES),
      expand('qc/fastqc/{sample}-{read}.trim_fastqc.zip',
             sample=SAMPLES, read=READS),
      expand('qc/hicup/HiCUP_summary_report-{sample}.txt', sample=SAMPLES),
      expand('qc/samtools/stats/{sample}.stats.txt', sample=SAMPLES),
      expand('qc/samtools/idxstats/{sample}.idxstats.txt', sample=SAMPLES),
      expand('qc/samtools/flagstat/{sample}.flagstat.txt', sample=SAMPLES),
      expand('qc/bamqc/{sample}', sample=SAMPLES),
      expand('qc/hicexplorer/{sample}-{region}.{bin}_QC',
             sample=SAMPLES, region=REGIONS.index, bin=BASE_BIN)]
)
if config['fastq_screen'] is not None:
    multiqc_input = multiqc_input.append(
        expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
               sample=SAMPLES, read=READS))

rule multiqc:
    input:
        [expand('qc/fastqc/{sample}-{read}.raw_fastqc.zip',
                sample=SAMPLES, read=READS),
         expand('qc/cutadapt/{sample}.cutadapt.txt', sample=SAMPLES),
         expand('qc/fastqc/{sample}-{read}.trim_fastqc.zip',
                sample=SAMPLES, read=READS),
         expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
                sample=SAMPLES, read=READS,
                proxy=[] if config['fastq_screen'] else [None]),
         expand('qc/hicup/HiCUP_summary_report-{sample}.txt', sample=SAMPLES),
         expand('qc/samtools/stats/{sample}.stats.txt', sample=SAMPLES),
         expand('qc/samtools/idxstats/{sample}.idxstats.txt', sample=SAMPLES),
         expand('qc/samtools/flagstat/{sample}.flagstat.txt', sample=SAMPLES),
         expand('qc/bamqc/{sample}', sample=SAMPLES),
         expand('qc/hicexplorer/{sample}-{region}.{bin}_QC',
                sample=SAMPLES, region=REGIONS.index, bin=BASE_BIN)]
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
