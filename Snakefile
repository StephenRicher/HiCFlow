#!/usr/bin/env python3

import os
import re
import sys
import math
import tempfile
from snake_setup import set_config, load_regions, load_coords, filterRegions, HiCSamples

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
    'tmpdir':            tempfile.gettempdir(),
    'threads':           workflow.cores,
    'data':              ''          ,
    'phased_vcf':        None        ,
    'genome':            ''          ,
    'bigWig':            {}          ,
    'bed':               {}          ,
    'regions':           ''          ,
    'cutadapt':
        {'forwardAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
         'reverseAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
         'overlap':         3                                 ,
         'errorRate':       0.1                               ,
         'minimumLength':   0                                 ,
         'qualityCutoff':  '0,0'                              ,
         'GCcontent':       50                                ,},
    'restrictionSeqs':      ''      ,
    'HiCParams':
        {'minBins':              50  ,
         'minDistance':          300  ,
         'maxLibraryInsertSize': 1000 ,
         'minMappingQuality':    15   ,
         'removeSelfLigation':   True ,
         'keepSelfCircles':      False,
         'skipDuplicationCheck': False,
         'nofill':               False,},
    'HiCcompare':
        {'fdr' :         0.05         ,
         'logFC' :       1            ,
         'multi' :       False        ,},
    'compareMatrices':
        {'vMin'       : -2.5          ,
         'vMax'       :  2.5          ,
         'size'       :  1            ,
         'maxDistance':  1000000      ,},
    'gatk':
        {'hapmap'      : None         ,
         'omni'        : None         ,
         'G1K'         : None         ,
         'dbsnp'       : None         ,
         'mills'       : None         ,
         'all_known'   : None         ,
         'trustPoly'   : False        ,
         'downSample'  : None         ,
         'scatterCount': 100},
    'resolution':
        {'base':          5000        ,
         'bins':         [5000, 10000],},
    'plot_coordinates':  None,
    'fastq_screen':      None,
    'phase':             True,
    'createValidBam':    False,
    'runQC':             True,
    'runHiCRep':         True,
    'plotRep':           True,
    'colourmap':         'Purples',
    'multiQCconfig':     None,
    'groupJobs':         False,
}

config = set_config(config, default_config)

workdir: config['workdir']
THREADS = config['threads']
BASE_BIN = config['resolution']['base']
ALLELE_SPECIFIC = True if config['phased_vcf'] else False

HiC = HiCSamples(config['data'], config['restrictionSeqs'], ALLELE_SPECIFIC)
REGIONS = load_regions(config['regions'])

# Remove region-binSize combinations with too few bins
regionBin = filterRegions(REGIONS, config['resolution']['bins'], nbins=config['HiCParams']['minBins'])

# Turn of phasing if allele specific mode is running
config['phase'] = False if ALLELE_SPECIFIC else config['phase']

if config['phase']:
    if config['gatk']['all_known']:
        PHASE_MODE = 'GATK'
    else:
        PHASE_MODE = 'BCFTOOLS'
else:
    PHASE_MODE = None


wildcard_constraints:
    cellType = rf'{"|".join(HiC.cellTypes())}',
    preGroup = rf'{"|".join(HiC.originalGroups())}',
    preSample = rf'{"|".join(HiC.originalSamples())}',
    region = rf'{"|".join(REGIONS.index)}',
    allele = r'[12]',
    rep = r'\d+',
    read = r'R[12]',
    bin = r'\d+',
    mode = r'SNP|INDEL',
    set = r'logFC|sig|fdr',
    compare = r'HiCcompare|multiHiCcompare',
    group = rf'{"|".join(HiC.groups())}',
    group1 = rf'{"|".join(HiC.groups())}',
    group2 = rf'{"|".join(HiC.groups())}',
    sample = rf'{"|".join(HiC.samples())}',
    all = rf'{"|".join(HiC.samples() + list(HiC.groups()))}'

# Generate dictionary of plot coordinates, may be multple per region
COORDS = load_coords([config['plot_coordinates'], config['regions']])

tools = (
    ['HiCcompare', 'multiHiCcompare'] if config['HiCcompare']['multi']
    else ['HiCcompare'])
HiC_mode = ([
    [expand('plots/{region}/{bin}/{tool}/{set}/{compare}-{region}-{coords}-{bin}-{set}.png',
        region=region, coords=COORDS[region], set=['logFC', 'sig'],
        compare=HiC.groupCompares(), bin=regionBin[region],
        tool=tools) for region in regionBin],
    [expand('plots/{region}/{bin}/pyGenomeTracks/{group}-{region}-{coords}-{bin}.png',
        region=region, coords=COORDS[region],
        group=HiC.groups(), bin=regionBin[region]) for region in regionBin],
    [expand('plots/{region}/{bin}/obs_exp/{all}-{region}-{bin}.png',
        all=(HiC.all() if config['plotRep'] else list(HiC.groups())),
        region=region, bin=regionBin[region]) for region in regionBin],
    'qc/hicup/.tmp.aggregatehicupTruncate'])

rule all:
    input:
        HiC_mode,
        (expand('phasedVCFs/{cellType}-phased.vcf', cellType=HiC.cellTypes())
         if config['phase'] else []),
        (expand('dat/gatk/.tmp.{cellType}-applyBQSR', cellType=HiC.cellTypes())
         if (config['phase'] and PHASE_MODE == 'GATK') else []),
        (['qc/multiqc', 'qc/filterQC/ditag_length.png',
         'qc/fastqc/.tmp.aggregateFastqc'] if config['runQC'] else []),
        ([expand('qc/hicrep/{region}-{bin}-hicrep.png', region=region,
            bin=regionBin[region]) for region in regionBin]
         if config['runHiCRep'] else []),
        (expand('dat/mapped/{sample}-validHiC.bam', sample=HiC.samples())
         if (config['createValidBam'] and regionBin) else [])


if ALLELE_SPECIFIC:
    rule maskPhased:
        input:
            genome = lambda wc: config['genome'][wc.cellType],
            vcf = lambda wc: config['phased_vcf'][wc.cellType]
        output:
            'dat/genome/masked/{cellType}.fa'
        group:
            'prepareGenome'
        log:
            'logs/maskPhased/{cellType}.log'
        conda:
            f'{ENVS}/bedtools.yaml'
        shell:
            'bedtools maskfasta -fullHeader '
            '-fi <(zcat -f {input.genome}) '
            '-bed {input.vcf} -fo {output} 2> {log}'


rule vcf2SNPsplit:
    input:
        lambda wc: config['phased_vcf'][wc.cellType]
    output:
        'snpsplit/{cellType}-snpsplit.txt'
    group:
        'prepareGenome'
    log:
        'logs/vcf2SNPsplit/{cellType}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/reformatSNPsplit.py {input} > {output} 2> {log}'


rule bgzipGenome:
    input:
        lambda wc: config['genome'][wc.cellType]
    output:
        'dat/genome/{cellType}.fa.gz'
    group:
        'prepareGenome'
    log:
        'logs/bgzipGenome/{cellType}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip > {output}) 2> {log}'


rule indexGenome:
    input:
        rules.bgzipGenome.output
    output:
        f'{rules.bgzipGenome.output}.fai'
    group:
        'prepareGenome'
    log:
        'logs/indexGenome/{cellType}-indexGenome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule getChromSizes:
    input:
        rules.indexGenome.output
    output:
        'dat/genome/chrom_sizes/{cellType}.chrom.sizes'
    group:
        'prepareGenome'
    log:
        'logs/getChromSizes/{cellType}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


rule findRestSites:
    input:
        rules.bgzipGenome.output
    output:
        'dat/genome/{cellType}-{re}-restSites.bed'
    params:
        reSeq = lambda wc: HiC.restrictionNames(removeCut=True)[wc.re]
    group:
        'prepareGenome'
    log:
        'logs/findRestSites/{cellType}-{re}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicFindRestSite --fasta {input} --searchPattern {params.reSeq} '
        '--outFile {output} &> {log}'


def bowtie2BuildInput(wildcards):
    if ALLELE_SPECIFIC:
        return rules.maskPhased.output
    else:
        return rules.bgzipGenome.output


rule bowtie2Build:
    input:
        bowtie2BuildInput
    output:
        expand('dat/genome/index/{{cellType}}.{n}.bt2',
               n=['1', '2', '3', '4', 'rev.1', 'rev.2'])
    params:
        basename = lambda wc: f'dat/genome/index/{wc.cellType}'
    group:
        'prepareGenome'
    log:
        'logs/bowtie2Build/{cellType}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params.basename} &> {log}'


rule fastQC:
    input:
        lambda wc: HiC.path(wc.preSample, [wc.read])
    output:
        html = 'qc/fastqc/{preSample}-{read}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{preSample}-{read}.raw.fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/fastqc/{preSample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule reformatFastQC:
    input:
        rules.fastQC.output.zip
    output:
        'qc/fastqc/{preSample}-{read}.raw_fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/reformatFastQC/{preSample}-{read}.raw.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyFastQC.py {input} {output} '
        '{wildcards.preSample}-{wildcards.read} &> {log}'


rule fastQCTrimmed:
    input:
        'dat/fastq/trimmed/{preSample}-{read}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{preSample}-{read}.trim_fastqc.html',
        zip = 'qc/fastqc/{preSample}-{read}.trim_fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/fastqc_trimmed/{preSample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule aggregateFastqc:
    input:
        expand('qc/fastqc/{sample}-{read}.{type}_fastqc.zip',
            sample=HiC.originalSamples(),
            read=['R1', 'R2'], type=['trim', 'raw'])
    output:
        touch('qc/fastqc/.tmp.aggregateFastqc')
    group:
        'fastqc' if config['groupJobs'] else 'aggregateTarget'


if config['fastq_screen'] is not None:

    rule fastQScreen:
        input:
            'dat/fastq/trimmed/{preSample}-{read}.trim.fastq.gz'
        output:
            txt = 'qc/fastq_screen/{preSample}-{read}.fastq_screen.txt',
            png = 'qc/fastq_screen/{preSample}-{read}.fastq_screen.png'
        params:
            fastq_screen_config = config['fastq_screen'],
            subset = 100000,
            aligner = 'bowtie2'
        log:
            'logs/fastq_screen/{preSample}-{read}.log'
        threads:
            THREADS
        wrapper:
            "0.60.0/bio/fastq_screen"


rule cutadapt:
    input:
        lambda wc: HiC.path(wc.preSample, ['R1', 'R2'])
    output:
        trimmed = [temp('dat/fastq/trimmed/{preSample}-R1.trim.fastq.gz'),
                   temp('dat/fastq/trimmed/{preSample}-R2.trim.fastq.gz')],
        qc = 'qc/cutadapt/unmod/{preSample}.cutadapt.txt'
    group:
        'cutadapt'
    params:
        forwardAdapter = config['cutadapt']['forwardAdapter'],
        reverseAdapter = config['cutadapt']['reverseAdapter'],
        overlap = config['cutadapt']['overlap'],
        errorRate = config['cutadapt']['errorRate'],
        minimumLength = config['cutadapt']['minimumLength'],
        qualityCutoff = config['cutadapt']['qualityCutoff'],
        GCcontent = config['cutadapt']['GCcontent']
    log:
        'logs/cutadapt/{preSample}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        THREADS
    shell:
        'cutadapt -a {params.forwardAdapter} -A {params.reverseAdapter} '
        '--overlap {params.overlap} --error-rate {params.errorRate} '
        '--minimum-length {params.minimumLength} '
        '--quality-cutoff {params.qualityCutoff} '
        '--gc-content {params.GCcontent} --cores {threads} '
        '-o {output.trimmed[0]} -p {output.trimmed[1]} {input} '
        '> {output.qc} 2> {log}'


rule reformatCutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{preSample}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/reformatCutadapt/{preSample}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyCutadapt.py {wildcards.preSample} {input} '
        '> {output} 2> {log}'


rule hicupTruncate:
    input:
        rules.cutadapt.output.trimmed
    output:
        truncated = [temp('dat/fastq/truncated/{preSample}-R1.trunc.fastq.gz'),
                     temp('dat/fastq/truncated/{preSample}-R2.trunc.fastq.gz')],
        summary = 'qc/hicup/{preSample}-truncate-summary.txt'
    params:
        re1 = lambda wc: list(HiC.restrictionSeqs()[wc.preSample].values())[0],
        fill = '--nofill' if config['HiCParams']['nofill'] else ''
    group:
        'hicupTruncate'
    threads:
        2 if THREADS > 2 else THREADS
    log:
        'logs/hicupTruncate/{preSample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/hicup/hicupTruncate.py {params.fill} '
        '--output {output.truncated} '
        '--summary {output.summary} '
        '--re1 {params.re1} '
        '--threads {threads} {input} &> {log}'


rule aggregatehicupTruncate:
    input:
        expand('qc/hicup/{sample}-truncate-summary.txt',
            sample=HiC.originalSamples())
    output:
        touch('qc/hicup/.tmp.aggregatehicupTruncate')
    group:
        'hicupTruncate' if config['groupJobs'] else 'aggregateTarget'


def bowtie2Index(wc):
    """ Retrieve bowtie2 index associated with sample. """
    return expand('dat/genome/index/{cellType}.{n}.bt2',
        cellType=HiC.sample2Cell()[wc.preSample],
        n=['1', '2', '3', '4', 'rev.1', 'rev.2'])


def bowtie2Basename(wc):
    """ Retrieve bowtie2 index basename associated with sample. """
    cellType = HiC.sample2Cell()[wc.preSample]
    return f'dat/genome/index/{cellType}'


rule bowtie2:
    input:
        fastq = 'dat/fastq/truncated/{preSample}-{read}.trunc.fastq.gz',
        bt2_index = bowtie2Index
    output:
        sam = pipe('dat/mapped/{preSample}-{read}.sam'),
        qc = 'qc/bowtie2/{preSample}-{read}.bowtie2.txt'
    params:
        index = bowtie2Basename,
        cellType = lambda wc: HiC.sample2Cell()[wc.preSample],
        sensitivity = 'sensitive'
    group:
        'bowtie2'
    log:
        'logs/bowtie2/{preSample}-{read}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        THREADS - 2 if THREADS > 1 else 1
    shell:
        'bowtie2 -x {params.index} -U {input.fastq} '
        '--reorder --threads {threads} --{params.sensitivity} '
        '> {output.sam} 2> {log} && cp {log} {output.qc}'


rule addReadFlag:
    input:
        rules.bowtie2.output.sam
    output:
        pipe('dat/mapped/{preSample}-{read}-addFlag.sam')
    params:
        flag = lambda wc: '0x41' if wc.read == 'R1' else '0x81'
    group:
        'bowtie2'
    log:
        'logs/addReadFlag/{preSample}-{read}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        '{SCRIPTS}/addReadFlag.awk -v flag={params.flag} {input} '
        '> {output} 2> {log}'


rule sam2bam:
    input:
        rules.addReadFlag.output
    output:
        temp('dat/mapped/{preSample}-{read}-addFlag.bam')
    group:
        'bowtie2'
    log:
        'logs/sam2bam/{preSample}-{read}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools view -u {input} > {output} 2> {log}'


rule catBam:
    input:
        'dat/mapped/{preSample}-R1-addFlag.bam',
        'dat/mapped/{preSample}-R2-addFlag.bam'
    output:
        pipe('dat/mapped/{preSample}-merged.bam')
    group:
        'prepareBAM'
    log:
        'logs/catBam/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools cat {input} > {output} 2> {log}'


rule collateBam:
    input:
        rules.catBam.output
    output:
        pipe('dat/mapped/{preSample}-collate.bam')
    params:
        # Add '/' to path if not present
        tmpPrefix = os.path.join(config['tmpdir'], '')
    group:
        'prepareBAM'
    log:
        'logs/collateBam/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools collate -Ofu {input} {params.tmpPrefix} > {output} 2> {log}'


rule fixmateBam:
    input:
        rules.collateBam.output
    output:
        temp('dat/mapped/{preSample}.fixed.bam')
    group:
        'prepareBAM'
    log:
        'logs/fixmateBam/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS - 2 if THREADS > 3 else 1
    shell:
        'samtools fixmate -@ {threads} -mp {input} {output} 2> {log}'


# Input to SNPsplit
rule removeUnmapped:
    input:
        rules.fixmateBam.output
    output:
        temp('dat/mapped/{preSample}.hic.bam')
    group:
        'prepareBAM'
    log:
        'logs/removeUnmapped/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -@ {threads} -b -F 12 {input} > {output} 2> {log}'


def SNPsplitInput(wc):
    """ Retrieve cell type associated with sample. """
    cellType = HiC.sample2Cell()[wc.preSample]
    return f'snpsplit/{cellType}-snpsplit.txt'


rule SNPsplit:
    input:
        bam = rules.removeUnmapped.output,
        snps = SNPsplitInput
    output:
        temp(expand('dat/snpsplit/{{preSample}}.hic.{ext}',
            ext = ['G1_G1.bam', 'G1_UA.bam',
                   'G2_G2.bam', 'G2_UA.bam',
                   'G1_G2.bam', 'UA_UA.bam'])),
        expand('dat/snpsplit/{{preSample}}.hic.{ext}',
            ext = ['SNPsplit_report.txt', 'SNPsplit_sort.txt'])
    params:
        outdir = 'dat/snpsplit/'
    group:
        'SNPsplit'
    log:
        'logs/SNPsplit/SNPsplit-{preSample}.log'
    conda:
        f'{ENVS}/snpsplit.yaml'
    shell:
        'SNPsplit {input.bam} --snp_file {input.snps} '
        '--hic --output_dir {params.outdir} &> {log}'


rule mergeSNPsplit:
    input:
        'dat/snpsplit/{preGroup}-{rep}.hic.G{allele}_G{allele}.bam',
        'dat/snpsplit/{preGroup}-{rep}.hic.G{allele}_UA.bam'
    output:
        temp('dat/snpsplit/merged/{preGroup}_a{allele}-{rep}.hic.bam')
    group:
        'SNPsplit'
    log:
        'logs/mergeSNPsplit/{preGroup}_a{allele}-{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools merge -n {output} {input} &> {log}'


def splitInput(wc):
    if ALLELE_SPECIFIC:
        return 'dat/snpsplit/merged/{sample}.hic.bam'
    else:
        return 'dat/mapped/{sample}.fixed.bam'


rule splitPairedReads:
    input:
        splitInput
    output:
        'dat/mapped/split/{sample}-{read}.bam'
    params:
        flag = lambda wc: '0x40' if wc.read == 'R1' else '0x80'
    group:
        'prepareBAM'
    log:
        'logs/splitPairedReads/{sample}-{read}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -@ {threads} -f {params.flag} -b {input} '
        '| samtools sort > {output} 2> {log}'


def getRestSites(wc):
    """ Retrieve restSite files associated with sample wildcard """
    try:
        sample = wc.sample
    except AttributeError:
        sample = wc.preSample
    return expand('dat/genome/{cellType}-{re}-restSites.bed',
        cellType=HiC.sample2Cell()[sample],
        re=HiC.restrictionSeqs()[sample])


def getRestrictionSeqs(wc):
    sequences = HiC.restrictionSeqs(removeCut=True)[wc.sample].values()
    return ' '.join(sequences)


def getDanglingSequences(wc):
    sequences = HiC.restrictionSeqs(dangling=True)[wc.sample].values()
    return ' '.join(sequences)


rule buildBaseMatrix:
    input:
        bams = expand('dat/mapped/split/{{sample}}-{read}.bam', read=['R1', 'R2']),
        restSites = getRestSites
    output:
        hic = f'dat/matrix/{{region}}/base/raw/{{sample}}-{{region}}.{BASE_BIN}.h5',
        bam = 'dat/matrix/{region}/{sample}-{region}.bam',
        qc = directory(f'qc/hicexplorer/{{sample}}-{{region}}.{BASE_BIN}_QC')
    params:
        bin = BASE_BIN,
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        reSeqs = getRestrictionSeqs,
        danglingSequences = getDanglingSequences,
        maxLibraryInsertSize = config['HiCParams']['maxLibraryInsertSize'],
        minMappingQuality = config['HiCParams']['minMappingQuality'],
        minDistance = config['HiCParams']['minDistance'],
        removeSelfLigation = (
            'True' if config['HiCParams']['removeSelfLigation'] else 'False'),
        keepSelfCircles = (
            '--keepSelfCircles' if config['HiCParams']['keepSelfCircles'] else ''),
        skipDuplicationCheck = (
            '--skipDuplicationCheck' if config['HiCParams']['skipDuplicationCheck'] else '')
    log:
        'logs/buildBaseMatrix/{sample}-{region}.log'
    threads:
        max(2, THREADS)
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicBuildMatrix --samFiles {input.bams} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--restrictionCutFile {input.restSites} '
        '--restrictionSequence {params.reSeqs} '
        '--maxLibraryInsertSize {params.maxLibraryInsertSize} '
        '--minDistance {params.minDistance} '
        '--minMappingQuality {params.minMappingQuality} '
        '--removeSelfLigation {params.removeSelfLigation} '
        '--danglingSequence {params.danglingSequences} '
        '{params.keepSelfCircles} '
        '{params.skipDuplicationCheck} --binSize {params.bin} '
        '--outFileName {output.hic} --outBam {output.bam} '
        '--QCfolder {output.qc} --threads {threads} '
        '&> {log}  || mkdir -p {output.qc}; touch {output.hic} {output.bam}'


rule mergeValidHiC:
    input:
        expand('dat/matrix/{region}/{{sample}}-{region}.bam',
            region=regionBin.keys())
    output:
        'dat/mapped/{sample}-validHiC.bam'
    log:
        'logs/mergeValidHiC/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {input} 2> {log}'


def nonEmpty(wc, output, input):
    """ Build find command to remove non-empty input files at runtime. """
    findCmd = '$(find -size +0 \( '
    for i, file in enumerate(input):
        # Skip index files
        if file.endswith('.csi'):
            continue
        elif i > 0:
            findCmd += ' -o '
        findCmd += f"-path './{file}'"
    findCmd += ' \))'
    return findCmd


rule sumReplicates:
    input:
        lambda wc: expand(
            'dat/matrix/{{region}}/base/raw/{{group}}-{rep}-{{region}}.{bin}.h5',
            rep=HiC.groups()[wc.group], bin=BASE_BIN),
    output:
        f'dat/matrix/{{region}}/base/raw/{{group}}-{{region}}.{BASE_BIN}.h5'
    params:
        nonEmpty = nonEmpty
    log:
        'logs/sumReplicates/{group}-{region}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicSumMatrices --matrices {params.nonEmpty} --outFileName {output} '
        '&> {log} || touch {output}'


rule mergeBins:
    input:
        f'dat/matrix/{{region}}/base/raw/{{all}}-{{region}}.{BASE_BIN}.h5'
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}.h5'
    params:
        nbins = lambda wc: int(int(wc.bin) / BASE_BIN)
    group:
        'processHiC'
    log:
        'logs/mergeBins/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log} || touch {output}'


rule IceMatrix:
    input:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}.h5'
    output:
        plot = 'qc/IceMatrix/{all}-{region}-{bin}-diagnosic_plot.png',
        matrix = 'dat/matrix/{region}/{bin}/ice/{all}-{region}-{bin}.h5'
    params:
        iternum = 1000,
        upper_threshold = 5
    group:
        'processHiC'
    log:
        'logs/IceMatrix/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '{SCRIPTS}/hicCorrect.sh -p {output.plot} '
        '-o {output.matrix} -u {params.upper_threshold} '
        '-i {params.iternum} {input} &> {log} || touch {output}'


rule TadInsulation:
    input:
        rules.IceMatrix.output.matrix
    output:
        expand(
            'dat/matrix/{{region}}/{{bin}}/tads/{{all}}-{{region}}-{{bin}}{ext}',
            ext = ['_boundaries.bed', '_boundaries.gff', '_domains.bed',
                   '_score.bedgraph', '_zscore_matrix.h5']),
        score = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}_tad_score.bm'
    params:
        method = 'fdr',
        bin = lambda wc: wc.bin,
        region = lambda wc: wc.region,
        all = lambda wc: wc.all,
        min_depth = lambda wc: int(wc.bin) * 3,
        max_depth = lambda wc: int(wc.bin) * 10,
        prefix = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}'
    group:
        'processHiC'
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
        'dat/matrix/{region}/{bin}/loops/{all}-{region}-{bin}.bedgraph'
    params:
        maxLoop = 2000000,
        windowSize = 10,
        peakWidth = 6,
        pValuePre = 0.05,
        pValue = 0.05,
        peakInter = 5
    group:
        'processHiC'
    log:
        'logs/detectLoops/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    threads:
        THREADS
    shell:
        'hicDetectLoops --matrix {input} --outFileName {output} '
        '--maxLoopDistance {params.maxLoop} '
        '--windowSize {params.windowSize} '
        '--peakWidth {params.peakWidth} '
        '--pValuePreselection {params.pValuePre} '
        '--pValue {params.pValue} '
        '--peakInteractionsThreshold {params.peakInter} '
        '--threads 1 --threadsPerChromosome {threads} '
        '&> {log} || touch {output} && touch {output} '


rule hicPCA:
    input:
        rules.IceMatrix.output.matrix
    output:
        'dat/matrix/{region}/{bin}/PCA/{all}-{region}-{bin}.bedgraph'
    params:
        method = "dist_norm",
        format = 'bedgraph'
    group:
        'processHiC'
    log:
        'logs/hicPCA/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicPCA --matrix {input} --outputFileName {output} '
        '--format {params.format} '
        '--numberOfEigenvectors 1 --method {params.method} '
        '--ignoreMaskedBins &> {log} || touch {output} && touch {output} '


rule fixBedgraph:
    input:
        rules.hicPCA.output
    output:
        'dat/matrix/{region}/{bin}/PCA/{all}-{region}-{bin}-fix.bedgraph'
    params:
        pos = lambda wc: REGIONS['end'][wc.region]
    group:
        'processHiC'
    log:
        'logs/fixBedgraph/{all}-{region}-{bin}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/fixBedgraph.py {input} --pos {params.pos} '
        '> {output} 2> {log}'


rule reformatHomer:
    input:
        'dat/matrix/{region}/{bin}/{method}/{all}-{region}-{bin}.h5'
    output:
        'dat/matrix/{region}/{bin}/{method}/{all}-{region}-{bin}.gz'
    group:
        'processHiC'
    log:
        'logs/reformatHomer/{all}-{region}-{bin}-{method}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat h5 --outputFormat homer &> {log} || touch {output}'


rule reformatNxN:
    input:
        'dat/matrix/{region}/{bin}/ice/{all}-{region}-{bin}.gz'
    output:
        'dat/matrix/{region}/{bin}/ice/{all}-{region}-{bin}.nxn.tsv'
    group:
        'processHiC'
    log:
        'logs/reformatNxN/{region}/{bin}/{all}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/reformatNxN.py <(zcat {input}) '
        '> {output} 2> {log} || touch {output}'


rule OnTAD:
    input:
        rules.reformatNxN.output
    output:
        bed = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.bed',
        tad = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.tad'
    params:
        chr = lambda wc: re.sub('chr', '', str(REGIONS['chr'][wc.region])),
        length = lambda wc: REGIONS['length'][wc.region],
        outprefix = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad'
    group:
        'processHiC'
    log:
        'logs/OnTAD/{region}/{bin}/{all}.log'
    shell:
        '{SCRIPTS}/OnTAD {input} -o {params.outprefix} -bedout chr{params.chr} '
        '{params.length} {wildcards.bin} &> {log} || touch {output}'


rule reformatLinks:
    input:
        rules.OnTAD.output.bed
    output:
        'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad.links'
    params:
        start = lambda wc: REGIONS['start'][wc.region]
    group:
        'processHiC'
    log:
        'logs/reformatLinks/{region}/{bin}/{all}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/reformatLinks.py {params.start} {wildcards.bin} {input} '
        '> {output} 2> {log}'


def getTracks(wc):
    """ Build track command for generate config """
    command = ''
    for title, track in config['bigWig'].items():
        command += f'--bigWig {title},{track} '
    for title, track in config['bed'].items():
        command += f'--bed {title},{track} '
    return command


rule createConfig:
    input:
        matrix = 'dat/matrix/{region}/{bin}/ice/{group}-{region}-{bin}.h5',
        loops = 'dat/matrix/{region}/{bin}/loops/{group}-{region}-{bin}.bedgraph',
        insulations = 'dat/matrix/{region}/{bin}/tads/{group}-{region}-{bin}_tad_score.bm',
        tads = 'dat/matrix/{region}/{bin}/tads/{group}-{region}-{bin}-ontad.links',
        pca = 'dat/matrix/{region}/{bin}/PCA/{group}-{region}-{bin}-fix.bedgraph'
    output:
        'plots/{region}/{bin}/pyGenomeTracks/configs/{group}-{region}-{bin}.ini'
    params:
        tracks = getTracks,
        depth = lambda wc: int(REGIONS['length'][wc.region]),
        colourmap = config['colourmap']
    group:
        'processHiC'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/createConfig/{region}/{bin}/{group}.log'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} '#--flip '
        '--insulations {input.insulations} --log '
        '--loops {input.loops} --colourmap {params.colourmap} '
        '--bigWig PCA1,{input.pca} '
        '--tads {input.tads} {params.tracks} '
        '--depth {params.depth} > {output} 2> {log}'


def setRegion(wc):
    """ Replace underscores with : and - for valid --region argument. """
    region = list(wc.coord)
    # Find indices of all underscores
    inds = [i for i,c in enumerate(region) if c == '_']
    # Replace penultimate and last underscore with : and -
    region[inds[-1]] = '-'
    region[inds[-2]] = ':'
    return ''.join(region)


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        'plots/{region}/{bin}/pyGenomeTracks/{group}-{region}-{coord}-{bin}.png'
    params:
        region = setRegion,
        title = '"{group} : {region} at {bin} bin size"',
        dpi = 600
    group:
        'processHiC'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    log:
        'logs/plotHiC/{region}/{bin}/{group}-{coord}.log'
    threads:
        THREADS
    shell:
        'export NUMEXPR_MAX_THREADS=1; pyGenomeTracks --tracks {input} '
        '--region {params.region} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


rule distanceNormalise:
    input:
        rules.IceMatrix.output.matrix
    output:
        'dat/matrix/{region}/{bin}/ice/obs_exp/{all}-{region}-{bin}.h5'
    params:
        method = 'obs_exp'
    group:
        'plotObsExp'
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
        'plots/{region}/{bin}/obs_exp/{all}-{region}-{bin}.png'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        title = '"{all} : {region} at {bin} bin size"',
        dpi = 600,
        colour = 'YlGn'
    group:
        'plotObsExp'
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


rule reformatNxN3p:
    input:
        'dat/matrix/{region}/{bin}/raw/{sample}-{region}-{bin}.gz'
    output:
        'dat/matrix/{region}/{bin}/raw/{sample}-{region}-{bin}.nxnp3.tsv'
    group:
        'HiCRep'
    log:
        'logs/reformatNxN3p/{sample}-{region}-{bin}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/reformatNxN3p.py {wildcards.bin} {wildcards.region} '
        '<(zcat {input}) > {output} 2> {log}'


rule HiCRep:
    input:
        'dat/matrix/{region}/{bin}/raw/{sample1}-{region}-{bin}.nxnp3.tsv',
        'dat/matrix/{region}/{bin}/raw/{sample2}-{region}-{bin}.nxnp3.tsv',
    output:
        'qc/hicrep/data/{sample1}-vs-{sample2}-{region}-{bin}-hicrep.csv'
    params:
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'HiCRep'
    log:
        'logs/HiCRep/{sample1}-vs-{sample2}-{region}-{bin}.log'
    conda:
        f'{ENVS}/hicrep.yaml'
    shell:
        '{SCRIPTS}/runHiCRep.R {output} {wildcards.bin} '
        '{params.start} {params.end} {input} &> {log}'


rule plotHiCRep:
    input:
        expand('qc/hicrep/data/{compare}-{{region}}-{{bin}}-hicrep.csv',
            compare=HiC.sampleCompares())
    output:
        'qc/hicrep/{region}-{bin}-hicrep.png'
    params:
        dpi = 600,
    group:
        'HiCRep'
    log:
        'logs/plotHiCRep/{region}-{bin}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotHiCRep.py --out {output} --dpi {params.dpi} '
        '{input} &> {log}'


rule mergeBamByReplicate:
    input:
        lambda wc: expand(
            'dat/matrix/{{region}}/{{group}}-{rep}-{{region}}.bam',
            rep = HiC.groups()[wc.group]),
    output:
        'dat/matrix/{region}/{group}-{region}.bam'
    log:
        'logs/mergeBamByReplicate/{region}/{group}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {input} 2> {log}'


rule reformatPre:
    input:
        'dat/matrix/{region}/{all}-{region}.bam'
    output:
        'dat/matrix/{region}/base/raw/{all}-{region}.pre.tsv'
    group:
        'bam2hic'
    log:
        'logs/reformatPre/{region}/{all}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        '(samtools view -@ {threads} {input} '
        '| awk -f {SCRIPTS}/bam2pre.awk > {output}) 2> {log} '


def getChromSizes(wc):
    """ Retrieve chromSizes file associated with group or sample. """
    return f'dat/genome/chrom_sizes/{HiC.sample2Cell()[wc.all]}.chrom.sizes',


rule juicerPre:
    input:
        tsv = 'dat/matrix/{region}/base/raw/{all}-{region}.pre.tsv',
        chrom_sizes = getChromSizes
    output:
        'UCSCcompatible/{region}/{all}-{region}.hic'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        resolutions = ','.join([str(bin) for bin in config['resolution']['bins']])
    group:
        'bam2hic'
    log:
        'logs/juicerPre/{region}/{all}.log'
    conda:
        f'{ENVS}/openjdk.yaml'
    shell:
        'java -jar {SCRIPTS}/juicer_tools_1.14.08.jar pre '
        '-c {params.chr} -r {params.resolutions} '
        '{input.tsv} {output} {input.chrom_sizes} &> {log}'


# Reform for HiCcompare input
rule straw:
    input:
        rules.juicerPre.output
    output:
        'dat/matrix/{region}/{bin}/{all}-{region}-{bin}-sutm.txt'
    params:
        # Strip 'chr' as juicer removes by default
        chr = lambda wc: re.sub('chr', '', str(REGIONS['chr'][wc.region])),
        start = lambda wc: REGIONS['start'][wc.region],
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'HiCcompare'
    log:
        'logs/straw/{region}/{bin}/{all}.log'
    conda:
        f'{ENVS}/hic-straw.yaml'
    shell:
        '{SCRIPTS}/run-straw.py NONE {input} '
        '{params.chr}:{params.start}:{params.end} '
        '{params.chr}:{params.start}:{params.end} '
        'BP {wildcards.bin} {output} &> {log}'


rule HiCcompare:
    input:
        'dat/matrix/{region}/{bin}/{group1}-{region}-{bin}-sutm.txt',
        'dat/matrix/{region}/{bin}/{group2}-{region}-{bin}-sutm.txt'
    output:
        all = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}.homer',
        sig = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-sig.homer',
        links = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}.links',
    params:
        dir = lambda wc: f'dat/HiCcompare/{wc.region}/{wc.bin}',
        qcdir = lambda wc: f'qc/HiCcompare/{wc.region}/{wc.bin}',
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        fdr = config['HiCcompare']['fdr']
    group:
        'HiCcompare'
    log:
        'logs/HiCcompare/{region}/{bin}/{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/HiCcompare.yaml'
    shell:
        '{SCRIPTS}/HiCcompare.R {params.dir} {params.qcdir} {params.chr} '
        '{params.start} {params.end} {wildcards.bin} {params.fdr} {input} '
        '&> {log}'


rule multiHiCcompare:
    input:
        group1 = lambda wc: expand(
            'dat/matrix/{{region}}/{{bin}}/{group1}-{rep}-{{region}}-{{bin}}-sutm.txt',
            group1=wc.group1, rep=HiC.groups()[wc.group1]),
        group2 = lambda wc: expand(
            'dat/matrix/{{region}}/{{bin}}/{group2}-{rep}-{{region}}-{{bin}}-sutm.txt',
            group2=wc.group2, rep=HiC.groups()[wc.group2])
    output:
        all = 'dat/multiHiCcompare/{region}/{bin}/{group1}-vs-{group2}.homer',
        sig = 'dat/multiHiCcompare/{region}/{bin}/{group1}-vs-{group2}-sig.homer',
        fdr = 'dat/multiHiCcompare/{region}/{bin}/{group1}-vs-{group2}-fdr.homer',
        links = 'dat/multiHiCcompare/{region}/{bin}/{group1}-vs-{group2}.links'
    params:
        dir = lambda wc: f'dat/multiHiCcompare/{wc.region}/{wc.bin}',
        qcdir = directory('qc/multiHiCcompare'),
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        fdr = config['HiCcompare']['fdr']
    group:
        'HiCcompare'
    log:
        'logs/multiHiCcompare/{region}/{bin}/{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/multiHiCcompare.yaml'
    shell:
        '{SCRIPTS}/multiHiCcompare.R {params.dir} {params.qcdir} '
        '{params.chr} {params.start} {params.end} '
        '{wildcards.bin} {params.fdr} '
        '{input.group1} {input.group2} &> {log}'


rule applyMedianFilter:
    input:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}.homer'
    output:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-logFC.homer'
    params:
        size = config['compareMatrices']['size']
    group:
        'HiCcompare'
    log:
        'logs/applyMedianFilter/{compare}/{region}/{bin}/{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/smoothHiC.py {input} --size {params.size} '
        '> {output} 2> {log}'


rule hicCompareBedgraph:
    input:
        rules.applyMedianFilter.output
    output:
        up = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-up.bedgraph',
        down = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-down.bedgraph'
    params:
        maxDistance = config['compareMatrices']['maxDistance']
    group:
        'HiCcompare'
    log:
        'logs/hicCompareBedgraph/{compare}/{region}/{bin}/{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/hicCompareBedgraph.py --binSize {wildcards.bin} '
        '--maxDistance {params.maxDistance} --upOut {output.up} '
        '--downOut {output.down} {input} &> {log}'


rule homerToH5:
    input:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}.homer'
    output:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}.h5'
    group:
        'HiCcompare'
    log:
        'logs/homerToH5/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '(hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 || touch {output})  &> {log}'


rule filterHiCcompare:
    input:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}.links'
    output:
        up = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-up.links',
        down = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-down.links'
    params:
        p_value = config['HiCcompare']['fdr'],
        log_fc = config['HiCcompare']['logFC'],
    group:
        'HiCcompare'
    log:
        'logs/filterHiCcompare/{compare}/{region}/{bin}/{group1}-vs-{group2}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filterHiCcompare.py '
        '--up {output.up} --down {output.down} '
        '--p_value {params.p_value} --log_fc {params.log_fc} {input} &> {log}'


rule createCompareConfig:
    input:
        mat = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}.h5',
        upBed = rules.hicCompareBedgraph.output.up,
        downBed = rules.hicCompareBedgraph.output.down
    output:
        'plots/{region}/{bin}/HiCcompare/configs/{group1}-vs-{group2}-{compare}-{set}.ini',
    params:
        depth = lambda wc: int(REGIONS['length'][wc.region]),
        colourmap = 'bwr',
        tracks = getTracks,
        vMin = lambda wc: -1 if wc.set == 'fdr' else config['compareMatrices']['vMin'],
        vMax = lambda wc: 1 if wc.set == 'fdr' else config['compareMatrices']['vMax'],
    group:
        'HiCcompare'
    log:
        'logs/createCompareConfig/{compare}/{region}/{bin}/{group1}-{group2}-{set}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.mat} --compare '
        '--sumLogFC {input.upBed},{input.downBed} '
        '{params.tracks} --depth {params.depth} --colourmap {params.colourmap} '
        '--vMin {params.vMin} --vMax {params.vMax} > {output} 2> {log}'


def round_down(wc):
    start = REGIONS['start'][wc.region]
    bin = int(wc.bin)
    return start - (start%bin)


def round_up(wc):
    end = REGIONS['end'][wc.region]
    bin = int(wc.bin)
    return end - (end%bin) + bin


def title(wc):
    title = f'"{wc.group1} vs {wc.group2} - {wc.region} at {wc.bin} bin size - '
    if wc.set == 'sig':
        threshold = config['HiCcompare']['fdr']
        title += f'adj. logFC (FDR <= {threshold})"'
    elif wc.set == 'logFC':
        title += 'adj. logFC"'
    else:
        title += 'FDR"'
    return title

rule plotCompare:
    input:
        rules.createCompareConfig.output
    output:
        'plots/{region}/{bin}/{compare}/{set}/{group1}-vs-{group2}-{region}-{coord}-{bin}-{set}.png'
    params:
        title = title,
        region = setRegion,
        dpi = 600
    group:
        'HiCcompare'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    log:
        'logs/plotAnalysis/{compare}/{region}/{bin}/{group1}-vs-{group2}-{coord}-{set}.log'
    threads:
        THREADS
    shell:
        'export NUMEXPR_MAX_THREADS=1; pyGenomeTracks --tracks {input} '
        '--region {params.region} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


if not ALLELE_SPECIFIC:

    rule sortBam:
        input:
            rules.fixmateBam.output
        output:
            pipe('dat/mapped/{preSample}.sorted.bam')
        params:
            mem = '1G'
        group:
            'deduplicate'
        log:
            'logs/sortBam/{preSample}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            max(math.ceil(THREADS * 0.5), 1)
        shell:
            'samtools sort -@ {threads} -O bam,level=0 '
            '-m {params.mem} {input} > {output} 2> {log}'


    rule deduplicate:
        input:
            rules.sortBam.output
        output:
            bam = temp('dat/mapped/{preSample}.dedup.bam'),
            qc = 'qc/deduplicate/{preSample}.txt'
        group:
            'deduplicate'
        log:
            'logs/deduplicate/{preSample}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            max(math.floor(THREADS * 0.5), 1)
        shell:
            'samtools markdup -@ {threads} '
            '-rsf {output.qc} {input} {output.bam} &> {log}'


    rule mergeBamByCellType:
        input:
            lambda wc: expand('dat/mapped/{preSample}.dedup.bam',
                preSample = HiC.cellTypes()[wc.cellType])
        output:
            pipe('dat/mapped/mergeByCell/{cellType}.merged.bam')
        group:
            'mergeCellType'
        log:
            'logs/mergeBamByCellType/{cellType}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            max(math.ceil(THREADS * 0.5), 1)
        shell:
            'samtools merge -@ {threads} - {input} > {output} 2> {log}'


    rule addReadGroup:
        input:
            rules.mergeBamByCellType.output
        output:
            'dat/mapped/mergeByCell/{cellType}.fixed-RG.bam'
        group:
            'mergeCellType'
        log:
            'logs/addReadGroup/{cellType}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            max(math.floor(THREADS * 0.5), 1)
        shell:
            'samtools addreplacerg -@ {threads} '
            '-r "ID:1\tPL:.\tPU:.\tLB:.\tSM:{wildcards.cellType}" '
            '-O BAM {input} > {output} 2> {log}'


    # GATK PHASE MODE #

    rule createSequenceDictionary:
        input:
            rules.bgzipGenome.output
        output:
            'dat/genome/{cellType}.dict'
        params:
            tmp = config['tmpdir']
        log:
            'logs/gatk/createSequenceDictionary/{cellType}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard CreateSequenceDictionary R={input} O={output} '
            'TMP_DIR={params.tmp} &> {log}'

    # Optionally downsample BAM for base recalibration
    rule downSampleBam:
        input:
            rules.addReadGroup.output
        output:
            temp('dat/mapped/mergeByCell/{cellType}.downSample.bam')
        params:
            fraction = config['gatk']['downSample']
        log:
            'logs/downSampleBam/{cellType}.log'
        threads:
            THREADS
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools view -@ {threads} -s {params.fraction} {input} '
            '> {output} 2> {log}'


    rule indexMergedBam:
        input:
            rules.addReadGroup.output
        output:
            f'{rules.addReadGroup.output}.bai'
        threads:
            THREADS
        log:
            'logs/indexMergedBam/{cellType}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools index -@ {threads} {input} &> {log}'


    def known_sites(input_known):
        input_known_string = ""
        if input_known is not None:
            for known in input_known:
                input_known_string += f' --known-sites {known}'
        return input_known_string


    def baseRecalibratorInput(wc):
        if config['gatk']['downSample'] is not None:
            return rules.downSampleBam.output
        else:
            return rules.addReadGroup.output


    rule baseRecalibrator:
        input:
            bam = baseRecalibratorInput,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'dat/gatk/baseRecalibrator/{cellType}.recal.table'
        params:
            known = known_sites(config['gatk']['all_known']),
            tmp = config['tmpdir'],
            extra = ''
        log:
            'logs/gatk/baseRecalibrator/{cellType}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk BaseRecalibrator {params.known} '
             '--input {input.bam} --reference {input.ref} '
             '--output {output} '
             '--sequence-dictionary {input.ref_dict} '
             '--tmp-dir {params.tmp} {params.extra} &> {log}'


    rule splitIntervals:
        input:
            rules.bgzipGenome.output
        output:
            expand('dat/gatk/splitIntervals/{{cellType}}/{rep}-scattered.interval_list',
                rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])])
        params:
            regions = config['regions'],
            scatterCount = config['gatk']['scatterCount']
        log:
            'logs/gatk/splitIntervals/{cellType}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk SplitIntervals -R {input} -L {params.regions} '
            '--scatter-count {params.scatterCount} '
            '--output dat/gatk/splitIntervals/{wildcards.cellType} &> {log} '


    rule applyBQSR:
        input:
            bam = rules.addReadGroup.output,
            bam_index = rules.indexMergedBam.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output,
            recal_table = rules.baseRecalibrator.output,
            interval = 'dat/gatk/splitIntervals/{cellType}/{rep}-scattered.interval_list'
        output:
            bam = 'dat/mapped/mergeByCell/{cellType}-{rep}.recalibrated.bam',
            index = 'dat/mapped/mergeByCell/{cellType}-{rep}.recalibrated.bai'
        params:
            tmp = config['tmpdir'],
            extra = ''
        group:
            'applyBQSR'
        log:
            'logs/gatk/applyBQSR/{cellType}-{rep}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk ApplyBQSR '
            '--input {input.bam} --reference {input.ref} '
            '--bqsr-recal-file {input.recal_table} --output {output.bam} '
            '--intervals {input.interval} --interval-padding 100 '
            '--tmp-dir {params.tmp} {params.extra} &> {log}'


    rule aggregateApplyBQSR:
        input:
            expand(
                'dat/mapped/mergeByCell/{{cellType}}-{rep}.recalibrated.bam',
                rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])])
        output:
            touch('dat/gatk/.tmp.{cellType}-applyBQSR')
        group:
            'applyBQSR' if config['groupJobs'] else 'aggregateTarget'


    rule haplotypeCaller:
        input:
            bam = rules.applyBQSR.output.bam,
            bam_index = rules.applyBQSR.output.index,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output,
            interval = 'dat/gatk/splitIntervals/{cellType}/{rep}-scattered.interval_list'
        output:
            'dat/gatk/split/{cellType}-{rep}-g.vcf.gz'
        params:
            java_opts = '-Xmx6G',
            min_prune = 2, # Increase to speed up
            downsample = 50, # Decrease to speed up
            tmp = config['tmpdir'],
            extra = ''
        log:
            'logs/gatk/haplotypeCaller/{cellType}-{rep}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        threads:
            min(4, THREADS)
        shell:
            'gatk --java-options {params.java_opts} HaplotypeCaller '
            '{params.extra} --input {input.bam} --output {output} '
            '--max-reads-per-alignment-start {params.downsample} '
            '--min-pruning {params.min_prune} --reference {input.ref} '
            '--intervals {input.interval} --interval-padding 100 '
            '--tmp-dir {params.tmp} -ERC GVCF &> {log}'


    def gatherVCFsInput(wc):
        input = ''
        gvcfs = expand('dat/gatk/split/{cellType}-{rep}-g.vcf.gz',
            rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])],
            cellType=wc.cellType)
        for gvcf in gvcfs:
            input += f' -I {gvcf}'
        return input


    rule gatherVCFs:
        input:
            expand('dat/gatk/split/{{cellType}}-{rep}-g.vcf.gz',
                rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])])
        output:
            'dat/gatk/merged/{cellType}-g.vcf.gz'
        params:
            gvcfs = gatherVCFsInput,
            java_opts = '-Xmx4G'
        group:
            'GATK'
        log:
            'logs/gatk/gatherGVCFs/{cellType}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk --java-options {params.java_opts} GatherVcfs '
             '{params.gvcfs} -O {output} &> {log}'


    rule sortGVCF:
        input:
            rules.gatherVCFs.output
        output:
            'dat/gatk/merged/{cellType}-sorted-g.vcf.gz'
        params:
            tmp = config['tmpdir']
        group:
            'GATK'
        log:
            'logs/picard/sortGVCF/{cellType}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard SortVcf INPUT={input} TMP_DIR={params.tmp} '
            'OUTPUT={output} &> {log}'


    rule indexFeatureFile:
        input:
            rules.sortGVCF.output
        output:
            f'{rules.sortGVCF.output}.tbi'
        params:
            java_opts = '-Xmx4G',
        group:
            'GATK'
        log:
            'logs/gatk/indexFeatureFile/{cellType}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk --java-options {params.java_opts} IndexFeatureFile '
             '-I {input} &> {log}'


    rule genotypeGVCFs:
        input:
            gvcf = rules.sortGVCF.output,
            gvcf_idex = rules.indexFeatureFile.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'dat/gatk/{cellType}.vcf.gz'
        params:
            tmp = config['tmpdir'],
            java_opts = '-Xmx4G'
        group:
            'GATK'
        log:
            'logs/gatk/genotypeGVCFs/{cellType}.log'
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
            'dat/gatk/{cellType}-{mode}.vcf.gz'
        params:
            tmp = config['tmpdir']
        group:
            'GATK'
        log:
            'logs/gatk/selectVariants/{cellType}-{mode}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk SelectVariants --reference {input.ref} '
            '--variant {input.vcf} --select-type-to-include {wildcards.mode} '
            '--tmp-dir {params.tmp} --output {output} &> {log}'


    rule variantRecalibratorSNPs:
        input:
            vcf = 'dat/gatk/{cellType}-SNP.vcf.gz',
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            recal = 'dat/gatk/{cellType}-SNP.vcf.recal',
            tranches = 'dat/gatk/{cellType}-SNP.vcf.tranches'
        params:
            hapmap = f'--resource:hapmap,known=false,training=true,truth=true,'
            f'prior=15.0 {config["gatk"]["hapmap"]}' if config["gatk"]["hapmap"]  else '',
            omni = f'--resource:omni,known=false,training=true,truth=true,'
            f'prior=12.0 {config["gatk"]["omni"]}' if config["gatk"]["omni"]  else '',
            G1K = f'--resource:G1K,known=false,training=true,truth=false,'
            f'prior=10.0 {config["gatk"]["G1K"]}' if config["gatk"]["G1K"]  else '',
            dbsnp = f'--resource:dbsnp,known=true,training=false,truth=false,'
            f'prior=7.0 {config["gatk"]["dbsnp"]}' if config["gatk"]["dbsnp"]  else '',
            trustPoly = '--trust-all-polymorphic' if config['gatk']['trustPoly'] else '',
            max_gaussians = 6,
            tmp = config['tmpdir'],
            java_opts = '-Xmx4G',
            extra = '',  # optional
        group:
            'GATK'
        log:
            'logs/gatk/variantRecalibrator/{cellType}-SNP.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk --java-options {params.java_opts} VariantRecalibrator '
            '--reference {input.ref} --variant {input.vcf} --mode SNP '
            '-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '
            '-an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR -an MQ '
            '--output {output.recal} --tranches-file {output.tranches} '
            '--max-gaussians {params.max_gaussians} '
            '{params.dbsnp} {params.hapmap} {params.omni} '
            '{params.G1K} {params.trustPoly} {params.extra} '
            '--tmp-dir {params.tmp} &> {log}'


    rule variantRecalibratorINDELS:
        input:
            vcf = 'dat/gatk/{cellType}-INDEL.vcf.gz',
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            recal = 'dat/gatk/{cellType}-INDEL.vcf.recal',
            tranches = 'dat/gatk/{cellType}-INDEL.vcf.tranches'
        params:
            mills = f'--resource:mills,known=false,training=true,truth=true,'
            f'prior=12.0 {config["gatk"]["mills"]}' if config["gatk"]["mills"]  else '',
            dbsnp = f'--resource:dbsnp,known=true,training=false,truth=false,'
            f'prior=2.0 {config["gatk"]["dbsnp"]}' if config["gatk"]["dbsnp"]  else '',
            trustPoly = '--trust-all-polymorphic' if config['gatk']['trustPoly'] else '',
            max_gaussians = 4,
            java_opts = '-Xmx4G',
            tmp = config['tmpdir'],
            extra = '',  # optional
        group:
            'GATK'
        log:
            'logs/gatk/variantRecalibrator/{cellType}-INDEL.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk --java-options {params.java_opts} VariantRecalibrator '
            '--reference {input.ref} --variant {input.vcf} --mode INDEL '
            '-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '
            '-an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR '
            '--output {output.recal} --tranches-file {output.tranches} '
            '--max-gaussians {params.max_gaussians} '
            '{params.dbsnp} {params.mills} {params.trustPoly} '
            '{params.extra} --tmp-dir {params.tmp} &> {log}'


    rule applyVQSR:
        input:
            vcf = rules.selectVariants.output,
            tranches = 'dat/gatk/{cellType}-{mode}.vcf.tranches',
            recal = 'dat/gatk/{cellType}-{mode}.vcf.recal',
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output
        output:
            'dat/gatk/{cellType}-{mode}.filt.vcf.gz'
        params:
            sensitivity = lambda wc: 99.5 if wc.mode == 'SNP' else 99.0,
            tmp = config['tmpdir']
        group:
            'GATK'
        log:
            'logs/gatk/applyVQSR/{cellType}-{mode}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk ApplyVQSR --reference {input.ref} --variant {input.vcf} '
            '--tranches-file {input.tranches} --mode {wildcards.mode} '
            '--exclude-filtered --recal-file {input.recal} '
            '--truth-sensitivity-filter-level {params.sensitivity} '
            '--tmp-dir {params.tmp} --output {output} &> {log}'


    rule mergeVCFs:
        input:
            SNP = 'dat/gatk/{cellType}-SNP.filt.vcf.gz',
            INDEL = 'dat/gatk/{cellType}-INDEL.filt.vcf.gz',
        output:
            'dat/gatk/{cellType}-all.filt.vcf.gz'
        params:
            tmp = config['tmpdir']
        group:
            'GATK'
        log:
            'logs/picard/mergeVCFs/{cellType}.log'
        conda:
            f'{ENVS}/picard.yaml'
        shell:
            'picard MergeVcfs INPUT={input.SNP} INPUT={input.INDEL} '
            'TMP_DIR={params.tmp} OUTPUT={output} &> {log}'


    rule splitVCFS:
        input:
            rules.mergeVCFs.output
        output:
            'dat/gatk/{cellType}-all-{region}.filt.vcf'
        params:
            chr = lambda wc: REGIONS['chr'][wc.region],
            start = lambda wc: REGIONS['start'][wc.region] + 1,
            end = lambda wc: REGIONS['end'][wc.region]
        group:
            'GATK'
        log:
            'logs/splitVCFS/{cellType}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view --regions {params.chr}:{params.start}-{params.end} '
            '{input} > {output} 2> {log}'


    # BCFTOOLS PHASE MODE #

    rule mpileup:
        input:
            bam = rules.addReadGroup.output,
            bam_index = rules.indexMergedBam.output,
            genome = rules.bgzipGenome.output,
        output:
            pipe('dat/bcftools/{region}/{cellType}-{region}-mpileup.bcf')
        params:
            chr = lambda wc: REGIONS['chr'][wc.region],
            start = lambda wc: REGIONS['start'][wc.region] + 1,
            end = lambda wc: REGIONS['end'][wc.region]
        group:
            'bcftoolsVariants'
        log:
            'logs/mpileup/{region}/{cellType}.log'
        threads:
            (THREADS - 4) * 0.5
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools mpileup -q 15 --ignore-RG --count-orphans '
            '--regions {params.chr}:{params.start}-{params.end} '
            '--max-depth 100000 --output-type u -f {input.genome} '
            '--threads {threads} {input.bam} > {output} 2> {log} '


    rule callVariants:
        input:
            rules.mpileup.output
        output:
            pipe('dat/bcftools/{region}/{cellType}-{region}-calls.bcf')
        group:
            'bcftoolsVariants'
        log:
            'logs/callVariants/{region}/{cellType}.log'
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
            'dat/bcftools/{region}/{cellType}-{region}.filt.vcf'
        group:
            'bcftoolsVariants'
        log:
            'logs/filterVariants/{region}/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -i "%QUAL>=20" --output-type v '
            '{input} > {output} 2> {log}'


    def hapCut2vcf(wc):
        if PHASE_MODE == 'GATK':
            return rules.splitVCFS.output
        else:
            return rules.filterVariants.output


    rule extractHAIRS:
        input:
            vcf = hapCut2vcf,
            bam = rules.addReadGroup.output
        output:
            'dat/phasing/{region}/{cellType}-{region}.fragments'
        params:
            chr = lambda wc: REGIONS['chr'][wc.region],
            start = lambda wc: REGIONS['start'][wc.region],
            end = lambda wc: REGIONS['end'][wc.region]
        group:
            'hapcut2'
        log:
            'logs/extractHAIRS/{region}/{cellType}.log'
        conda:
            f'{ENVS}/hapcut2.yaml'
        shell:
            'extractHAIRS --hic 1 --bam {input.bam} '
            '--regions {params.chr}:{params.start}-{params.end} '
            '--VCF {input.vcf} --out {output} &> {log}'


    rule hapCut2:
        input:
            fragments = rules.extractHAIRS.output,
            vcf = hapCut2vcf
        output:
            block = 'dat/phasing/{region}/{cellType}-{region}',
            vcf = 'dat/phasing/{region}/{cellType}-{region}.phased.VCF'
        group:
            'hapcut2'
        log:
            'logs/hapCut2/{region}/{cellType}.log'
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
        group:
            'hapcut2'
        log:
            'logs/bgzip_phased/{region}/{cellType}.log'
        conda:
            f'{ENVS}/tabix.yaml'
        shell:
            'bgzip -c {input} > {output} 2> {log}'


    rule indexPhased:
        input:
            rules.bgzipPhased.output
        output:
            f'{rules.bgzipPhased.output}.tbi'
        group:
            'hapcut2'
        log:
            'logs/index_phased/{region}/{cellType}.log'
        conda:
            f'{ENVS}/tabix.yaml'
        shell:
            'tabix {input} &> {log}'


    rule extractBestBlock:
        input:
            rules.hapCut2.output.block
        output:
            'dat/phasing/{region}/{cellType}-{region}.tsv'
        group:
            'hapcut2'
        log:
            'logs/extractBestPhase/{region}/{cellType}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/extractBestHapcut2.py {input} > {output} 2> {log}'


    rule extractVCF:
        input:
            block = rules.extractBestBlock.output,
            vcf = rules.bgzipPhased.output,
            vcf_index = rules.indexPhased.output
        output:
            'dat/phasing/{region}/{cellType}-{region}-best.vcf'
        group:
            'hapcut2'
        log:
            'logs/extractVCF/{region}/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -R {input.block} {input.vcf} '
            '> {output} 2> {log} || touch {output}'


    rule bgzipVCF:
        input:
             rules.extractVCF.output
        output:
            f'{rules.extractVCF.output}.gz'
        group:
            'hapcut2'
        log:
            'logs/bgzipVCF/{region}/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -O z {input} > {output} 2> {log} || touch {output}'


    rule indexVCF:
        input:
            rules.bgzipVCF.output
        output:
            f'{rules.bgzipVCF.output}.csi'
        group:
            'hapcut2'
        log:
            'logs/indexVCF/{region}/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools index -f {input} 2> {log} || touch {output}'


    rule mergeVCFsbyRegion:
        input:
            expand('dat/phasing/{region}/{{cellType}}-{region}-best.vcf.{ext}',
                region=REGIONS.index, ext=['gz', 'gz.csi']),
        output:
            'phasedVCFs/{cellType}-phased.vcf'
        params:
            nonEmpty = nonEmpty
        group:
            'hapcut2' if config['groupJobs'] else 'mergeVCFsbyRegion'
        log:
            'logs/mergeVCFsbyRegion/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            '(bcftools concat --allow-overlaps {params.nonEmpty} > {output} '
            '|| bcftools view {params.nonEmpty} > {output}) 2> {log} '


    rule bcftoolsStats:
        input:
            rules.filterVariants.output
        output:
            'qc/bcftools/{region}/{cellType}-{region}-bcftoolsStats.txt'
        group:
            'bcftoolsVariants'
        log:
            'logs/bcftoolsStats/{region}/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools stats {input} > {output} 2> {log}'


rule sampleReads:
    input:
        'dat/mapped/{preSample}.fixed.bam'
    output:
        'dat/mapped/subsampled/{preSample}-subsample.sam'
    group:
        'filterQC'
    params:
        seed = '42',
        frac = '20'
    threads:
        2 if THREADS > 2 else THREADS
    log:
        'logs/sampleReads/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools view -@ {threads} -s {params.seed}.{params.frac} {input} '
        '> {output} 2> {log}'


rule processHiC:
    input:
        reads = rules.sampleReads.output,
        digest = getRestSites
    output:
        'dat/mapped/subsampled/{preSample}-processed.txt'
    group:
        'filterQC'
    log:
        'logs/process/{preSample}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/processHiC.py {input.digest[0]} {input.reads} '
        '> {output} 2> {log}'


rule plotQC:
    input:
        expand('dat/mapped/subsampled/{sample}-processed.txt',
            sample=HiC.originalSamples())
    output:
        expand('qc/filterQC/{fig}',
               fig=['trans_stats.csv', 'insert_size_frequency.png',
                    'ditag_length.png'])
    params:
        outdir='qc/filterQC'
    group:
        'filterQC' if config['groupJobs'] else 'plotQC'
    log:
        'logs/plotQC/plot_subsample.log'
    conda:
        f'{ENVS}/ggplot2.yaml'
    shell:
        '{SCRIPTS}/plotQC.R {params.outdir} {input} 2> {log}'


rule mergeHicupQC:
    input:
        rules.hicupTruncate.output.summary,
    output:
        'qc/hicup/HiCUP_summary_report-{preSample}.txt'
    log:
        'logs/mergeHicupQC/{preSample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        '{SCRIPTS}/mergeHicupSummary.py --truncater {input} > {output} 2> {log}'

def multiQCconfig():
    if config['multiQCconfig']:
        return f'--config {config["multiQCconfig"]}'
    else:
        return ''

rule multiqc:
    input:
        [expand('qc/fastqc/{sample}-{read}.{mode}_fastqc.zip',
            sample=HiC.originalSamples(),
            read=['R1', 'R2'], mode=['raw', 'trim']),
         expand('qc/cutadapt/{sample}.cutadapt.txt',
            sample=HiC.originalSamples()),
         expand('qc/hicup/HiCUP_summary_report-{sample}.txt',
            sample=HiC.originalSamples()),
         expand('qc/bowtie2/{sample}-{read}.bowtie2.txt',
            sample=HiC.originalSamples(), read=['R1', 'R2']),
        [expand('qc/hicexplorer/{sample}-{region}.{bin}_QC', region=region,
            sample=HiC.samples(), bin=BASE_BIN) for region in regionBin],
         expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
            sample=HiC.originalSamples(), read=['R1', 'R2']) if config['fastq_screen'] else [],
         expand('qc/bcftools/{region}/{cellType}-{region}-bcftoolsStats.txt',
            region=REGIONS.index, cellType=HiC.cellTypes()) if PHASE_MODE=='BCFTOOLS' else []]
    output:
        directory('qc/multiqc')
    params:
        config = multiQCconfig()
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc --outdir {output} --force {params.config} {input} &> {log}'
