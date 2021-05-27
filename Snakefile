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
    'build':             None        ,
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
         'nofill':               False,
         'threads':              4    ,
         'multiplicativeValue':  10000,},
    'compareMatrices':
        {'minZ':         2            ,
         'vMin'       : -4            ,
         'vMax'       :  4            ,
         'size'       :  1            ,
         'maxDistance':  1000000      ,
         'tads'       :  None         ,},
    'gatk':
        {'hapmap'      : None         ,
         'omni'        : None         ,
         'G1K'         : None         ,
         'dbsnp'       : None         ,
         'mills'       : None         ,
         'all_known'   : None         ,
         'maxGaussians': 8            ,
         'trustPoly'   : False        ,
         'downSample'  : None         ,
         'scatterCount': 100          ,},
    'resolution':
        {'base':          5000        ,
         'bins':         [5000, 10000],},
    'plotParams':
        {'distanceNorm'  : False    ,
         'plain'         : True     ,
         'colourmap'     : 'Purples',
         'coordinates'   : None     ,
         'viewpoints'    : None     ,
         'viewpointRange': 500000   ,
         'plotRep'       : True     ,},
    'bigWig'           : {}           ,
    'bed'              : {}           ,
    'localAlignment':    False,
    'fastq_screen':      None,
    'runQC':             True,
    'runPCA':            True,
    'phase':             True,
    'createValidBam':    False,
    'runHiCRep':         True,
    'multiQCconfig':     None,
    'rescalePKL':        False,
    'groupJobs':         False,
    'microC':            False,
    'ASHIC':             True,
}

config = set_config(config, default_config)

workdir: config['workdir']
THREADS = config['threads']
BASE_BIN = config['resolution']['base']
ALLELE_SPECIFIC = True if config['phased_vcf'] else False

HiC = HiCSamples(config['data'], config['restrictionSeqs'], ALLELE_SPECIFIC)
REGIONS = load_regions(config['regions'], adjust=BASE_BIN)

# Remove region-binSize combinations with too few bins
regionBin, binRegion = filterRegions(REGIONS, config['resolution']['bins'], nbins=config['HiCParams']['minBins'])

if ALLELE_SPECIFIC:
    # Turn of phasing if allele specific mode is running
    config['phase'] = False
    # Turn off plotRep for allele specific mode
    if config['ASHIC']:
        if config['plotParams']['plotRep']:
            print('Per replicate plots not availale in allele specific mode '
                  'because replicates are merged for ASHIC. Setting plotRep '
                  'to False.', file=sys.stderr)
            config['plotParams']['plotRep'] = False
        if config['createValidBam']:
            print('createValid bam not available in ASHIC mode. Setting '
                  'createValidBam to False.', file=sys.stderr)
            config['createValidBam'] = False
        if config['runHiCRep']:
            print('HiCRep not available in ASHIC mode. Setting '
                  'runHiCRep to False.', file=sys.stderr)
            config['runHiCRep'] = False

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
    chr = rf'{"|".join(list(REGIONS["chr"].unique().astype(str)))}',
    region = rf'{"|".join(REGIONS.index)}',
    allele = r'[12]',
    rep = r'\d+',
    read = r'R[12]',
    bin = r'\d+',
    mode = r'SNP|INDEL',
    norm = r'raw|KR',
    pm = r'ASHIC|SNPsplit' if ALLELE_SPECIFIC else r'full',
    set = r'logFC|adjIF1|adjIF2',
    adjIF = r'adjIF1|adjIF2',
    compare = r'HiCcompare',
    group = rf'{"|".join(HiC.groups())}',
    group1 = rf'{"|".join(HiC.groups())}',
    group2 = rf'{"|".join(HiC.groups())}',
    sample = rf'{"|".join(HiC.samples(all=True))}',
    all = rf'{"|".join(HiC.samples() + list(HiC.groups()))}',
    combo = r'alt_alt|alt_both-ref|both-ref_both-ref|ref_alt|ref_both-ref|ref_ref'

# Generate dictionary of plot coordinates, may be multple per region
COORDS = load_coords(REGIONS, config['plotParams']['coordinates'], adjust=BASE_BIN)

# Generate dictionary of plot viewpoints
VIEWPOINTS =  load_coords(REGIONS, config['plotParams']['viewpoints'], includeRegions=False)

# Set whether to print a plain HiC map in addiion to a custom
vis = ['plain', 'custom'] if config['plotParams']['plain'] else ['custom']
# Set whether to print a raw HiC map in addiion to a KR
norm = ['raw', 'KR'] if config['plotParams']['plain'] else ['KR']
# Set plot suffix for allele specific mode
if ALLELE_SPECIFIC:
    phaseMode = 'ASHIC' if config['ASHIC'] else 'SNPsplit'
else:
    phaseMode = 'full'

HiC_mode = ([
    [expand('plots/{region}/{bin}/HiCcompare/logFC/{compare}-{region}-{coords}-{bin}-logFC-{pm}.png',
        region=region, coords=COORDS[region], pm=phaseMode,
        compare=HiC.groupCompares(), bin=regionBin[region]) for region in regionBin],
    [expand('plots/{region}/{bin}/viewpoints/HiCcompare/{compare}-{region}-{coords}-{bin}-viewpoint-{pm}.png',
        region=region, coords=VIEWPOINTS[region], pm=phaseMode,
        compare=HiC.groupCompares(), bin=regionBin[region]) for region in regionBin],
    [expand('plots/{region}/{bin}/pyGenomeTracks/{norm}/{group}-{region}-{coords}-{bin}-{vis}-{pm}.png',
        region=region, coords=COORDS[region], norm=norm, pm=phaseMode,
        vis=vis, group=HiC.groups(),
        bin=regionBin[region]) for region in regionBin],
    [expand('plots/{region}/{bin}/viewpoints/{norm}/{preGroup}-{region}-{coords}-{bin}-viewpoint-{pm}.png',
        region=region, coords=VIEWPOINTS[region], norm=norm, pm=phaseMode,
        preGroup=HiC.groups(), bin=regionBin[region]) for region in regionBin],
    [expand('plots/{region}/{bin}/obs_exp/{norm}/{all}-{region}-{bin}-{pm}.png',
        all=(HiC.all() if config['plotParams']['plotRep'] else list(HiC.groups())),
        region=region, bin=regionBin[region], pm=phaseMode,
        norm=norm) for region in regionBin],
     expand('qc/matrixCoverage/{region}/{all}-coverage-{pm}.png',
        all=(HiC.all() if config['plotParams']['plotRep'] else list(HiC.groups())),
        region=regionBin.keys(), pm=phaseMode),
    'qc/hicup/.tmp.aggregatehicupTruncate' if not config['microC'] else []])
# Exclude PCA if not set
if config['runPCA']:
    methods = ['PCA', 'TADinsulation', 'TADboundaries', 'TADdomains']
else:
    methods = ['TADinsulation', 'TADboundaries', 'TADdomains']
rescalePKL = ([
    expand('intervals/{compare}-{dir}-HiCcompare{mode}-{bin}-{pm}.pkl',
            dir=['up', 'down', 'all'],  mode=['', '-count'],
            pm=phaseMode, compare=HiC.groupCompares(), bin=binRegion.keys()),
    expand('intervals/{all}-{bin}-{method}-{pm}.pkl',
            all=(HiC.all() if config['plotParams']['plotRep'] else list(HiC.groups())),
            bin=binRegion.keys(), pm=phaseMode, method=methods)])


def getChromSizes(wc):
    """ Retrieve chromSizes file associated with group or sample. """
    try:
        cellType = HiC.sample2Cell()[wc.all]
    except AttributeError:
        try:
            cellType = HiC.sample2Cell()[wc.preGroup]
        except AttributeError:
            cellType = HiC.sample2Cell()[wc.group1]
            cellType2 = HiC.sample2Cell()[wc.group2]
            if cellType != cellType2:
                sys.stderr.write(
                    f'{wc.group1} and {wc.group2} correspond to different cell '
                    'type. Ensure the chromosome sizes are equal for valid '
                    'bedgraph rescaling of HiCcompare output.\n')
    return f'dat/genome/chrom_sizes/{cellType}.chrom.sizes',


rule all:
    input:
        HiC_mode,
        (expand('phasedVCFs/{cellType}-phased.vcf', cellType=HiC.cellTypes())
         if config['phase'] else []),
        (['qc/multiqc', 'qc/filterQC/ditagLength.png',
          'qc/fastqc/.tmp.aggregateFastqc'] if config['runQC'] else []),
        ([expand('qc/hicrep/{region}-{bin}-hicrep-{pm}.png', region=region,
            bin=regionBin[region], pm=phaseMode) for region in regionBin]
         if config['runHiCRep'] else []),
        (expand('dat/mapped/{sample}-validHiC-{pm}.bam',
            sample=HiC.samples(), pm=phaseMode)
         if (config['createValidBam'] and regionBin) else []),
        rescalePKL if config['rescalePKL'] else []


if ALLELE_SPECIFIC:

    rule filterHomozygous:
        input:
            lambda wc: config['phased_vcf'][wc.cellType]
        output:
            'dat/genome/{cellType}-phasedHet.vcf'
        group:
            'prepareGenome'
        log:
            'logs/filterHomozygous/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -H -m 2 -M 2 -v snps -i \'GT="het"\' --phased '
            '{input} > {output} 2> {log}'


    rule maskPhased:
        input:
            genome = lambda wc: config['genome'][wc.cellType],
            vcf = rules.filterHomozygous.output
        output:
            'dat/genome/{cellType}-masked.fa'
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
            rules.filterHomozygous.output
        output:
            'snpsplit/{cellType}-snpsplit.txt'
        group:
            'prepareGenome'
        log:
            'logs/vcf2SNPsplit/{cellType}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            'python {SCRIPTS}/reformatSNPsplit.py {input} > {output} 2> {log}'


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


rule subsetRestSites:
    input:
        rules.findRestSites.output
    output:
        'dat/genome/{region}/{cellType}-{re}-restSites.bed'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
    group:
        'prepareGenome'
    log:
        'logs/subsetRestSites/{cellType}-{re}-{region}.log'
    shell:
        'python {SCRIPTS}/subsetBED.py {input} '
        '--region {params.chr}:{params.start}-{params.end} '
        '> {output} 2> {log}'


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
        'python {SCRIPTS}/modifyFastQC.py {input} {output} '
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
        'python {SCRIPTS}/modifyCutadapt.py {wildcards.preSample} {input} '
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
        'python {SCRIPTS}/hicup/hicupTruncate.py {params.fill} '
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


def fastqInput(wc):
    if config['microC'] or config['localAlignment']:
        return 'dat/fastq/trimmed/{preSample}-{read}.trim.fastq.gz'
    else:
        return 'dat/fastq/truncated/{preSample}-{read}.trunc.fastq.gz'


rule bowtie2:
    input:
        fastq = fastqInput,
        bt2_index = bowtie2Index
    output:
        sam = pipe('dat/mapped/{preSample}-{read}.sam'),
        qc = 'qc/bowtie2/{preSample}-{read}.bowtie2.txt'
    params:
        index = bowtie2Basename,
        cellType = lambda wc: HiC.sample2Cell()[wc.preSample],
        sensitivity = 'sensitive',
        local = '--local' if config['localAlignment'] else ''
    group:
        'bowtie2'
    log:
        'logs/bowtie2/{preSample}-{read}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        THREADS - 2 if THREADS > 1 else 1
    shell:
        'bowtie2 -x {params.index} -U {input.fastq} {params.local} '
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
        'awk -f {SCRIPTS}/addReadFlag.awk -v flag={params.flag} {input} '
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
        'samtools collate -Ou {input} {params.tmpPrefix} > {output} 2> {log}'


rule fixmateBam:
    input:
        rules.collateBam.output
    output:
        'dat/mapped/{preSample}.fixed.bam'
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
        'dat/mapped/{preSample}.hic.bam'
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
        expand('dat/snpsplit/{{preSample}}.hic.{ext}',
            ext = ['SNPsplit_report.txt', 'SNPsplit_sort.txt']),
        bam = expand('dat/snpsplit/{{preSample}}.hic.{ext}',
            ext = ['G1_G1.bam', 'G1_UA.bam',
                   'G2_G2.bam', 'G2_UA.bam',
                   'G1_G2.bam', 'UA_UA.bam'])
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
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} -n {output} {input} &> {log}'


rule mergeALLSNPsplit:
    input:
        rules.SNPsplit.output.bam
    output:
        temp('dat/snpsplit/merged/{preSample}.ASHIC.bam')
    group:
        'mergeALLSNPsplit'
    log:
        'logs/mergeSNPsplit/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} -n {output} {input} &> {log}'


def splitInput(wc):
    if ALLELE_SPECIFIC:
        if config['ASHIC']:
            return 'dat/snpsplit/merged/{sample}.ASHIC.bam'
        else:
            return 'dat/snpsplit/merged/{sample}.hic.bam'
    else:
        return 'dat/mapped/{sample}.fixed.bam'


rule splitPairedReads:
    input:
        splitInput
    output:
        temp('dat/mapped/split/{sample}-{read}-{pm}.bam')
    params:
        flag = lambda wc: '0x40' if wc.read == 'R1' else '0x80'
    group:
        'prepareBAM'
    log:
        'logs/splitPairedReads/{sample}-{read}-{pm}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -@ {threads} -f {params.flag} -b {input} '
        '> {output} 2> {log}'


def getRestSites(wc):
    """ Retrieve restSite files associated with sample wildcard """
    try:
        sample = wc.sample
    except AttributeError:
        sample = wc.preSample
    try:
        region=f'{wc.region}/'
    except AttributeError:
        region=''
    return expand('dat/genome/{region}{cellType}-{re}-restSites.bed',
        region=region,
        cellType=HiC.sample2Cell()[sample],
        re=HiC.restrictionSeqs()[sample])


def getRestrictionSeqs(wc):
    sequences = HiC.restrictionSeqs(removeCut=True)[wc.sample].values()
    return ' '.join(sequences)


def getDanglingSequences(wc):
    sequences = HiC.restrictionSeqs(dangling=True)[wc.sample].values()
    return ' '.join(sequences)

if config['microC']:
    rule buildBaseMatrix:
        input:
            bams = expand('dat/mapped/split/{{sample}}-{read}-{{pm}}.bam', read=['R1', 'R2']),
        output:
            hic = f'dat/matrix/{{region}}/base/raw/{{sample}}-{{region}}.{BASE_BIN}-{{pm}}.h5',
            bam = 'dat/matrix/{region}/{sample}-{region}-{pm}.bam',
            qc = directory(f'qc/hicexplorer/{{sample}}-{{region}}.{BASE_BIN}-{{pm}}_QC')
        params:
            bin = BASE_BIN,
            chr = lambda wc: REGIONS['chr'][wc.region],
            start = lambda wc: REGIONS['start'][wc.region] + 1,
            end = lambda wc: REGIONS['end'][wc.region],
            inputBufferSize = 400000,
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
            'logs/buildBaseMatrix/{sample}-{region}-{pm}.log'
        threads:
            max(2, config['HiCParams']['threads'])
        conda:
            f'{ENVS}/hicexplorer3.4.yaml'
        shell:
            'hicBuildMatrix --samFiles {input.bams} '
            '--region {params.chr}:{params.start}-{params.end} '
            '--maxLibraryInsertSize {params.maxLibraryInsertSize} '
            '--minDistance {params.minDistance} '
            '--minMappingQuality {params.minMappingQuality} '
            '--removeSelfLigation {params.removeSelfLigation} '
            '--inputBufferSize {params.inputBufferSize} '
            '{params.keepSelfCircles} '
            '{params.skipDuplicationCheck} --binSize {params.bin} '
            '--outFileName {output.hic} --outBam {output.bam} '
            '--QCfolder {output.qc} --threads {threads} '
            '&> {log} || mkdir -p {output.qc}; touch {output.hic} {output.bam}'
else:
    rule buildBaseMatrix:
        input:
            bams = expand('dat/mapped/split/{{sample}}-{read}-{{pm}}.bam', read=['R1', 'R2']),
            restSites = getRestSites
        output:
            hic = f'dat/matrix/{{region}}/base/raw/{{sample}}-{{region}}.{BASE_BIN}-{{pm}}.h5',
            bam = 'dat/matrix/{region}/{sample}-{region}-{pm}.bam',
            qc = directory(f'qc/hicexplorer/{{sample}}-{{region}}.{BASE_BIN}-{{pm}}_QC')
        params:
            bin = BASE_BIN,
            chr = lambda wc: REGIONS['chr'][wc.region],
            start = lambda wc: REGIONS['start'][wc.region] + 1,
            end = lambda wc: REGIONS['end'][wc.region],
            reSeqs = getRestrictionSeqs,
            inputBufferSize = 400000,
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
            'logs/buildBaseMatrix/{sample}-{region}-{pm}.log'
        threads:
            max(2, config['HiCParams']['threads'])
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
            '--inputBufferSize {params.inputBufferSize} '
            '{params.keepSelfCircles} '
            '{params.skipDuplicationCheck} --binSize {params.bin} '
            '--outFileName {output.hic} --outBam {output.bam} '
            '--QCfolder {output.qc} --threads {threads} '
            '&> {log} || mkdir -p {output.qc}; touch {output.hic} {output.bam}'


rule mergeValidHiC:
    input:
        expand('dat/matrix/{region}/{{sample}}-{region}-{{pm}}.bam',
            region=regionBin.keys())
    output:
        'dat/mapped/{sample}-validHiC-{pm}.bam'
    log:
        'logs/mergeValidHiC/{sample}-{pm}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {input} '
        '2> {log} || touch {output}'


rule reformatASHIC:
    input:
        'dat/mapped/{preSample}-validHiC-ASHIC.bam'
    output:
        readPairs = expand('dat/ashic/readPairs/{{preSample}}/{chr}_{combo}',
            chr=list(REGIONS['chr'].unique()),
            combo=['alt_alt', 'alt_both-ref', 'both-ref_both-ref',
                   'ref_alt', 'ref_both-ref', 'ref_ref']),
        dir = directory('dat/ashic/readPairs/{preSample}/')
    params:
        prefix = lambda wc: f'dat/ashic/readPairs/{wc.preSample}/'
    log:
        'logs/reformatASHIC/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        '(touch {output.readPairs} && awk -f {SCRIPTS}/processSNPsplit.awk '
        '-v prefix={params.prefix} '
        '<(samtools cat {input} | samtools view -@ {threads})) &> {log} '
        '|| touch {output.readPairs}; mkdir -p {output.dir}'


rule mergeASHIC:
    input:
        lambda wc: expand(
            'dat/ashic/readPairs/{preGroup}-{rep}/{chr}_{combo}',
            preGroup=wc.preGroup, rep=HiC.originalGroups()[wc.preGroup],
            chr=wc.chr, combo=wc.combo)
    output:
        'dat/ashic/readPairs/{preGroup}_{chr}_{combo}'
    log:
        'logs/mergeAshic/{preGroup}-{chr}-{combo}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'cat {input} > {output} 2> {log}'


rule ASHICbin:
    input:
        lambda wc: expand(
            'dat/ashic/readPairs/{{preGroup}}_{chr}_{combo}',
            chr=REGIONS['chr'][wc.region],
            combo=['alt_alt', 'alt_both-ref', 'both-ref_both-ref',
                   'ref_alt', 'ref_both-ref', 'ref_ref']),
        chromSizes = getChromSizes
    output:
        directory(f'dat/ashic/binned/{{preGroup}}-{{region}}-{BASE_BIN}-binned/')
    params:
        bin = BASE_BIN,
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region],
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'ASHIC'
    log:
        'logs/ASHICbin/{preGroup}-{region}.log'
    conda:
        f'{ENVS}/ashic.yaml'
    shell:
        'ashic-data binning --res {params.bin} --chrom {params.chr} '
        '--c1 1 --p1 2 --a1 3 --c2 4 --p2 5 --a2 6 '
        '--genome {input.chromSizes} '
        '--start {params.start} --end {params.end} '
        'dat/ashic/readPairs/{wildcards.preGroup} {output} &> {log}'


rule ASHICpack:
    input:
        rules.ASHICbin.output
    output:
        directory(f'dat/ashic/packed/{{preGroup}}-{{region}}-{BASE_BIN}-packed/')
    params:
        diag = 0,
        perc = 2
    group:
        'ASHIC'
    log:
        'logs/ASHICpack/{preGroup}-{region}.log'
    conda:
        f'{ENVS}/ashic.yaml'
    shell:
        'ashic-data pack --diag {params.diag} --perc {params.perc} {input} '
        '{output} &> {log}'


rule ASHIC:
    input:
        rules.ASHICpack.output
    output:
        directory(f'dat/ashic/imputed/{{preGroup}}-{{region}}-{BASE_BIN}-imputed/')
    params:
        model = 'ASHIC-ZIPM',
        diag = 0,
        maxIter = 100,
        seed = 0
    group:
        'ASHIC'
    log:
        'logs/ASHIC/{preGroup}-{region}.log'
    conda:
        f'{ENVS}/ashic.yaml'
    shell:
        'ashic -i {input}/*.pickle -o {output} --model {params.model} '
        '--diag {params.diag} --seed {params.seed} '
        '--max-iter {params.maxIter} &> {log}'


rule ASHIC2homer:
    input:
        rules.ASHIC.output
    output:
        f'dat/matrix/{{region}}/base/raw/{{preGroup}}_a{{allele}}-{{region}}-ASHIC.{BASE_BIN}.homer'
    params:
        bin = BASE_BIN,
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region],
        mat = lambda wc: 't_mm.txt' if wc.allele == '1' else 't_pp.txt'
    group:
        'ASHIC'
    log:
        'logs/ASHIC2homer/{preGroup}-a{allele}-{region}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/ashic2homer.py {input}/matrices/{params.mat} '
        '--chrom {params.chr} --start {params.start} --binSize {params.bin} '
        '> {output} 2> {log} '


rule ASHIChomerToH5:
    input:
        rules.ASHIC2homer.output
    output:
        f'dat/matrix/{{region}}/base/raw/{{preGroup}}_a{{allele}}-{{region}}.{BASE_BIN}-ASHIC.h5'
    group:
        'ASHIC'
    log:
        'logs/ASHIChomerToH5/{preGroup}-a{allele}-{region}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 &> {log}'


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
            'dat/matrix/{{region}}/base/raw/{{group}}-{rep}-{{region}}.{bin}-{{pm}}.h5',
            rep=HiC.groups()[wc.group], bin=BASE_BIN),
    output:
        f'dat/matrix/{{region}}/base/raw/{{group}}-{{region}}.{BASE_BIN}-{{pm}}.h5'
    params:
        nonEmpty = nonEmpty
    log:
        'logs/sumReplicates/{group}-{region}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicSumMatrices --matrices {params.nonEmpty} --outFileName {output} '
        '&> {log} || touch {output}'


rule normCounts01:
    input:
        f'dat/matrix/{{region}}/base/raw/{{all}}-{{region}}.{BASE_BIN}-{{pm}}.h5'
    output:
        f'dat/matrix/{{region}}/base/raw/{{all}}-{{region}}.{BASE_BIN}-{{pm}}-norm01.h5'
    log:
        'logs/normCounts01/{all}-{region}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicNormalize --matrices {input} --normalize norm_range '
        '--outFileName {output} &> {log} || touch {output}'


rule normCountsConstant:
    input:
        rules.normCounts01.output
    output:
        f'dat/matrix/{{region}}/base/raw/{{all}}-{{region}}.{BASE_BIN}-{{pm}}-normConstant.h5'
    params:
        multiplicativeValue = config['HiCParams']['multiplicativeValue']
    log:
        'logs/normCountsConstant/{all}-{region}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicNormalize --matrices {input} --normalize multiplicative '
        '--multiplicativeValue {params.multiplicativeValue} '
        '--outFileName {output} &> {log} || touch {output}'


rule mergeBins:
    input:
        f'dat/matrix/{{region}}/base/raw/{{all}}-{{region}}.{BASE_BIN}-{{pm}}-normConstant.h5'
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}.h5'
    params:
        nbins = lambda wc: int(int(wc.bin) / BASE_BIN)
    log:
        'logs/mergeBins/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log} || touch {output}'


rule correctMatrix:
    input:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}.h5'
    output:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.h5'
    group:
        'processHiC'
    log:
        'logs/correctMatrix/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicCorrectMatrix correct --matrix {input} --correctionMethod KR '
        '--outFileName {output} &> {log} || touch {output}'


rule TADinsulation:
    input:
        rules.correctMatrix.output
    output:
        expand(
            'dat/matrix/{{region}}/{{bin}}/tads/{{all}}-{{region}}-{{bin}}-{{pm}}{ext}',
            ext = ['_boundaries.bed', '_boundaries.gff', '_domains.bed',
                   '_score.bedgraph', '_zscore_matrix.h5']),
        score = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-{pm}_tad_score.bm'
    params:
        method = 'fdr',
        bin = lambda wc: wc.bin,
        region = lambda wc: wc.region,
        all = lambda wc: wc.all,
        min_depth = lambda wc: int(wc.bin) * 3,
        max_depth = lambda wc: int(wc.bin) * 10,
        prefix = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-{pm}'
    group:
        'processHiC'
    threads:
        THREADS
    log:
        'logs/TADinsulation/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicFindTADs --matrix {input} '
        '--minDepth {params.min_depth} --maxDepth {params.max_depth} '
        '--step {wildcards.bin} --outPrefix {params.prefix} '
        '--correctForMultipleTesting {params.method} '
        '--numberOfProcessors {threads} &> {log} || touch {output}'


rule rescaleTADinsulation:
    input:
        bedgraphs = lambda wc: expand(
            'dat/matrix/{region}/{{bin}}/tads/{{all}}-{region}-{{bin}}-{{pm}}_score.bedgraph',
            region=binRegion[wc.bin]),
        chromSizes = getChromSizes
    output:
        'intervals/{all}-{bin}-TADinsulation-{pm}.pkl'
    params:
        regions = config['regions'],
        name = lambda wc: f'{wc.all}-{wc.bin}-TADinsulation',
    group:
        'rescaleBedgraphs'
    log:
        'logs/rescaleTADinsulations/{all}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/rescaleBedgraph.py sum --name {params.name} '
        '--out {output} --binSize {wildcards.bin} --regions {params.regions} '
        '--filetype bedgraph {input.chromSizes} <(cat {input.bedgraphs}) '
        '&> {log}'


rule rescaleTADdomains:
    input:
        bedgraphs = lambda wc: expand(
            'dat/matrix/{region}/{{bin}}/tads/{{all}}-{region}-{{bin}}-{{pm}}-ontad_domains.bed',
            region=binRegion[wc.bin]),
        chromSizes = getChromSizes
    output:
        'intervals/{all}-{bin}-TADdomains-{pm}.pkl'
    params:
        name = lambda wc: f'{wc.all}-{wc.bin}-TADdomains',
        regions = config['regions']
    group:
        'rescaleBedgraphs'
    log:
        'logs/rescaleTADdomains/{all}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/rescaleBedgraph.py count --name {params.name} '
        '--out {output} --binSize {wildcards.bin} --regions {params.regions} '
        '{input.chromSizes} <(cat {input.bedgraphs}) &> {log}'


rule rescaleTADboundaries:
    input:
        bedgraphs = lambda wc: expand(
            'dat/matrix/{region}/{{bin}}/tads/{{all}}-{region}-{{bin}}-{{pm}}-ontad_boundaries.bed',
            region=binRegion[wc.bin]),
        chromSizes = getChromSizes
    output:
        'intervals/{all}-{bin}-TADboundaries-{pm}.pkl'
    params:
        name = lambda wc: f'{wc.all}-{wc.bin}-TADboundaries',
        regions = config['regions']
    group:
        'rescaleBedgraphs'
    log:
        'logs/rescaleTADboundaries/{all}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/rescaleBedgraph.py count --name {params.name} '
        '--out {output} --binSize {wildcards.bin} --regions {params.regions} '
        '{input.chromSizes} <(cat {input.bedgraphs}) &> {log}'


rule detectLoops:
    input:
        rules.correctMatrix.output
    output:
        'dat/matrix/{region}/{bin}/loops/unmod/{all}-{region}-{bin}-{pm}.bedgraph'
    params:
        peakWidth = 6,
        windowSize = 10,
        pValuePre = 0.1,
        peakInteractionsThreshold = 10,
        obsExpThreshold = 1.5,
        pValue = 0.1,
        maxLoopDistance = 2000000
    group:
        'processHiC'
    log:
        'logs/detectLoops/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    threads:
        THREADS
    shell:
        'hicDetectLoops --matrix {input} --outFileName {output} '
        '--maxLoopDistance {params.maxLoopDistance} '
        '--windowSize {params.windowSize} '
        '--peakWidth {params.peakWidth} '
        '--pValuePreselection {params.pValuePre} '
        '--obsExpThreshold {params.obsExpThreshold} '
        '--pValue {params.pValue} '
        '--peakInteractionsThreshold {params.peakInteractionsThreshold} '
        '--threads 1 --threadsPerChromosome {threads} '
        '&> {log} || touch {output} && touch {output} '

# Hacky fix to correct loop intervals that extend too far
rule clampLoops:
    input:
        rules.detectLoops.output
    output:
        'dat/matrix/{region}/{bin}/loops/{all}-{region}-{bin}-{pm}.bedgraph'
    params:
        minPos = lambda wc: REGIONS['start'][wc.region],
        maxPos = lambda wc: REGIONS['end'][wc.region]
    group:
        'processHiC'
    log:
        'logs/clampLoops/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/clampLoops.py '
        '{params.minPos} {params.maxPos} {input} > {output} 2> {log}'


rule hicPCA:
    input:
        rules.correctMatrix.output
    output:
        'dat/matrix/{region}/{bin}/PCA/{all}-{region}-{bin}-{pm}.bedgraph'
    params:
        method = "dist_norm",
        format = 'bedgraph'
    group:
        'processHiC'
    log:
        'logs/hicPCA/{all}-{region}-{bin}-{pm}.log'
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
        'dat/matrix/{region}/{bin}/PCA/{all}-{region}-{bin}-fix-{pm}.bedgraph'
    params:
        pos = lambda wc: REGIONS['end'][wc.region]
    group:
        'processHiC'
    log:
        'logs/fixBedgraph/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/fixBedgraph.py {input} --pos {params.pos} '
        '> {output} 2> {log}'


rule rescalePCA:
    input:
        bedgraphs = lambda wc: expand(
            'dat/matrix/{region}/{{bin}}/PCA/{{all}}-{region}-{{bin}}-fix-{{pm}}.bedgraph',
            region=binRegion[wc.bin]),
        chromSizes = getChromSizes
    output:
        'intervals/{all}-{bin}-PCA-{pm}.pkl'
    params:
        regions = config['regions'],
        name = lambda wc: f'{wc.all}-{wc.bin}-PCA',
    group:
        'rescaleBedgraphs'
    log:
        'logs/rescalePCA/{all}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/rescaleBedgraph.py sum --name {params.name} '
        '--out {output} --binSize {wildcards.bin} --regions {params.regions} '
        '--filetype bedgraph {input.chromSizes} <(cat {input.bedgraphs}) '
        '&> {log}'


rule reformatHomer:
    input:
        'dat/matrix/{region}/{bin}/{method}/{all}-{region}-{bin}-{pm}.h5'
    output:
        'dat/matrix/{region}/{bin}/{method}/{all}-{region}-{bin}-{pm}.gz'
    group:
        'processHiC'
    log:
        'logs/reformatHomer/{all}-{region}-{bin}-{method}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat h5 --outputFormat homer &> {log} || touch {output}'


rule reformatNxN:
    input:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.gz'
    output:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.nxn.tsv'
    group:
        'processHiC'
    log:
        'logs/reformatNxN/{region}/{bin}/{all}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/reformatNxN.py <(zcat {input}) '
        '> {output} 2> {log} || touch {output}'


rule plotCoverage:
    input:
        lambda wc: expand('dat/matrix/{{region}}/{bin}/raw/{{all}}-{{region}}-{bin}-{{pm}}.gz',
            bin=regionBin[wc.region])
    output:
        'qc/matrixCoverage/{region}/{all}-coverage-{pm}.png'
    params:
        dpi = 300,
        nBins = 10000,
        fontSize = 12,
        nonEmpty = nonEmpty
    log:
        'logs/plotCoverage/{region}/{all}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/contactCoverage.py {params.nonEmpty} '
        '--out {output} --dpi {params.dpi} --nBins {params.nBins} '
        '--fontSize {params.fontSize} &> {log}'


rule OnTAD:
    input:
        rules.reformatNxN.output
    output:
        bed = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad-{pm}.bed',
        tad = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad-{pm}.tad'
    params:
        chr = lambda wc: re.sub('chr', '', str(REGIONS['chr'][wc.region])),
        length = lambda wc: REGIONS['length'][wc.region],
        outprefix = 'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-ontad-{pm}'
    group:
        'processHiC'
    log:
        'logs/OnTAD/{region}/{bin}/{all}-{pm}.log'
    shell:
        '{SCRIPTS}/OnTAD {input} -o {params.outprefix} -bedout chr{params.chr} '
        '{params.length} {wildcards.bin} &> {log} || touch {output}'


def setTrim(wc):
    if str(REGIONS['chr'][wc.region]).startswith('chr'):
        trim = ''
    else:
        trim = '--trimChr'
    return trim


rule reformatDomains:
    input:
        rules.OnTAD.output.bed
    output:
        'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-{pm}-ontad_domains.bed'
    params:
        scale = lambda wc: REGIONS['start'][wc.region] - int(wc.bin),
        trimChr = setTrim
    group:
        'processHiC'
    log:
        'logs/reformatDomains/{region}/{bin}/{all}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/reformatDomains.py {input} --scale {params.scale} '
        '{params.trimChr} > {output} 2> {log}'


rule domain2boundaries:
    input:
        rules.reformatDomains.output
    output:
        'dat/matrix/{region}/{bin}/tads/{all}-{region}-{bin}-{pm}-ontad_boundaries.bed'
    group:
        'processHiC'
    log:
        'logs/domain2boundaries/{region}/{bin}/{all}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/TADdomains2boundaries.py {wildcards.bin} {input} '
        '> {output} 2> {log}'


rule computeStripeScore:
    input:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.gz'
    output:
        forward = 'dat/matrix/{region}/{bin}/stripes/{all}-{region}-{bin}-forwardStripe-{pm}.bedgraph',
        rev = 'dat/matrix/{region}/{bin}/stripes/{all}-{region}-{bin}-reverseStripe-{pm}.bedgraph'
    group:
        'processHiC'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/computeStripeScore/{all}-{region}-{bin}-{pm}.log'
    shell:
        'python {SCRIPTS}/computeStripeScore.py {input} '
        '{output.forward} {output.rev} &> {log} '


rule distanceNormalise:
    input:
        'dat/matrix/{region}/{bin}/{norm}/{all}-{region}-{bin}-{pm}.h5'
    output:
        'dat/matrix/{region}/{bin}/{norm}/obs_exp/{all}-{region}-{bin}-{pm}.h5'
    params:
        method = 'obs_exp'
    group:
        'processHiC'
    log:
        'logs/distanceNormalise/{all}-{region}-{bin}-{norm}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicTransform -m {input} --method {params.method} -o {output} '
        '&> {log} || touch {output}'


def getTracks(wc):
    """ Build track command for generate config """
    command = ''
    for title, track in config['bigWig'].items():
        command += f'--bigWig {title},{track} '
    for title, track in config['bed'].items():
        command += f'--bed {title},{track} '
    return command


def getMatrix(wc):
    """ Return either normal or obs_exp matrix """
    if config['plotParams']['distanceNorm']:
        return 'dat/matrix/{region}/{bin}/{norm}/obs_exp/{group}-{region}-{bin}-{pm}.h5'
    else:
        return 'dat/matrix/{region}/{bin}/{norm}/{group}-{region}-{bin}-{pm}.h5'

def getPCAinput(wc):
    if config['runPCA']:
        return 'dat/matrix/{region}/{bin}/PCA/{group}-{region}-{bin}-fix-{pm}.bedgraph'
    else:
        return []

def getPCAparams(wc):
    if config['runPCA']:
        pca = f'dat/matrix/{wc.region}/{wc.bin}/PCA/{wc.group}-{wc.region}-{wc.bin}-fix-{wc.pm}.bedgraph'
        return f'--bigWig PCA1,{pca}'
    else:
        return ''

rule createConfig:
    input:
        matrix = getMatrix,
        loops = 'dat/matrix/{region}/{bin}/loops/{group}-{region}-{bin}-{pm}.bedgraph',
        insulations = 'dat/matrix/{region}/{bin}/tads/{group}-{region}-{bin}-{pm}_tad_score.bm',
        tads = 'dat/matrix/{region}/{bin}/tads/{group}-{region}-{bin}-{pm}-ontad_domains.bed',
        pca = getPCAinput,
        #forwardStripe = 'dat/matrix/{region}/{bin}/stripes/{group}-{region}-{bin}-forwardStripe.bedgraph',
        #reverseStripe = 'dat/matrix/{region}/{bin}/stripes/{group}-{region}-{bin}-reverseStripe.bedgraph'
    output:
        'plots/{region}/{bin}/pyGenomeTracks/{norm}/configs/{group}-{region}-{bin}-{vis}-{pm}.ini'
    params:
        tracks = getTracks,
        depth = lambda wc: int(REGIONS['length'][wc.region]),
        colourmap = config['plotParams']['colourmap'],
        vMin = '--vMin 0' if config['plotParams']['distanceNorm'] else '',
        vMax = '--vMax 2' if config['plotParams']['distanceNorm'] else '',
        log = '' if config['plotParams']['distanceNorm'] else '--log',
        plain = lambda wc: '--plain' if wc.vis == 'plain' else '',
        pca = getPCAparams
    group:
        'processHiC'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/createConfig/{group}-{region}-{bin}-{norm}-{vis}-{pm}.log'
    shell:
        'python {SCRIPTS}/generate_config.py --matrix {input.matrix} '
        '{params.log} --colourmap {params.colourmap} {params.tracks} '
        '--depth {params.depth} {params.vMin} {params.vMax} {params.plain} '
        '--insulations {input.insulations} --loops {input.loops} '
        '{params.pca} --tads {input.tads} > {output} 2> {log}'


def setRegion(wc):
    """ Replace underscores with : and - for valid --region argument. """
    region = list(wc.coord)
    # Find indices of all underscores
    inds = [i for i,c in enumerate(region) if c == '_']
    # Replace penultimate and last underscore with : and -
    region[inds[-1]] = '-'
    region[inds[-2]] = ':'
    return ''.join(region)


def setMatrixTitle(wc):
    if config['build'] is not None:
        build = config['build'].replace('"', '') # Double quotes disallowed
        build = f' ({build})'
    else:
        build = ''
    try:
        name = wc.group
    except AttributeError:
        name = wc.all
    title = f'"{name} : {wc.region}{build} at {wc.bin} bin size ({wc.norm} - {wc.pm})"',
    return title


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        'plots/{region}/{bin}/pyGenomeTracks/{norm}/{group}-{region}-{coord}-{bin}-{vis}-{pm}.png'
    params:
        region = setRegion,
        title = setMatrixTitle,
        dpi = 600
    group:
        'processHiC'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    log:
        'logs/plotHiC/{group}-{coord}-{region}-{bin}-{norm}-{vis}-{pm}.log'
    threads:
        THREADS
    shell:
        'export NUMEXPR_MAX_THREADS=1; pyGenomeTracks --tracks {input} '
        '--region {params.region} --outFileName {output} '
        '--title {params.title} --dpi {params.dpi} &> {log}'


def makeViewRegion(wc):
    """ Pad either side of viewpoint to define view range then
        replace underscores with : and - for valid --region argument. """
    region = list(wc.coord)
    # Find indices of all underscores
    inds = [i for i,c in enumerate(region) if c == '_']
    end = int(''.join(region[inds[-1]+1:]))
    start = int(''.join(region[inds[-2]+1:inds[-1]]))
    chrom = ''.join(region[:inds[-2]])
    size = config['plotParams']['viewpointRange']
    # Ensure padding does not exceed region length
    start = max(start - size, REGIONS.loc[wc.region, 'start'] + 1)
    end = min(end + size, REGIONS.loc[wc.region, 'end'])
    return f'{chrom}:{start}-{end}'


rule runViewpoint:
    input:
        'dat/matrix/{region}/{bin}/{norm}/{group}-{region}-{bin}-{pm}.h5'
    output:
        bedgraph = 'dat/viewpoints/{region}/{bin}/{norm}/{group}-{region}-{coord}-{bin}-{pm}.bedgraph',
        plot = temp('dat/viewpoints/{region}/{bin}/{norm}/{group}-{region}-{coord}-{bin}-{pm}.png')
    params:
        referencePoint = setRegion,
        region = makeViewRegion,
    group:
        'processHiC'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    log:
        'logs/runViewpoint/{group}-{coord}-{region}-{bin}-{norm}-{pm}.log'
    shell:
        '(hicPlotViewpoint --matrix {input} --region {params.region} '
        '--outFileName {output.plot} --referencePoint {params.referencePoint} '
        '--interactionOutFileName {output.bedgraph} '
        '&& mv {output.bedgraph}_*.bedgraph {output.bedgraph}) &> {log}'


rule plotViewpoint:
    input:
        rules.runViewpoint.output.bedgraph
    output:
        'plots/{region}/{bin}/viewpoints/{norm}/{group}-{region}-{coord}-{bin}-viewpoint-{pm}.png'
    params:
        dpi = 600,
        build = f'--build {config["build"]}' if config['build'] else ''
    group:
        'processHiC'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/plotViewpoint/{group}-{coord}-{region}-{bin}-{norm}-{pm}.log'
    shell:
        'python {SCRIPTS}/plotViewpoint.py {input} --out {output} '
        '--dpi {params.dpi} {params.build} &> {log}'


rule plotMatrix:
    input:
        rules.distanceNormalise.output
    output:
        'plots/{region}/{bin}/obs_exp/{norm}/{all}-{region}-{bin}-{pm}.png'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        title = setMatrixTitle,
        dpi = 600,
        colour = 'YlGn'
    log:
        'logs/plotMatrix/{all}-{region}-{bin}-{norm}-{pm}.log'
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
        'dat/matrix/{region}/{bin}/raw/{sample}-{region}-{bin}-{pm}.gz'
    output:
        'dat/matrix/{region}/{bin}/raw/{sample}-{region}-{bin}-{pm}.nxnp3.tsv'
    group:
        'HiCRep'
    log:
        'logs/reformatNxN3p/{sample}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/reformatNxN3p.py {wildcards.bin} {wildcards.region} '
        '<(zcat {input}) > {output} 2> {log}'


rule HiCRep:
    input:
        'dat/matrix/{region}/{bin}/raw/{sample1}-{region}-{bin}-{pm}.nxnp3.tsv',
        'dat/matrix/{region}/{bin}/raw/{sample2}-{region}-{bin}-{pm}.nxnp3.tsv',
    output:
        'qc/hicrep/data/{sample1}-vs-{sample2}-{region}-{bin}-hicrep-{pm}.csv'
    params:
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'HiCRep'
    log:
        'logs/HiCRep/{sample1}-vs-{sample2}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicrep.yaml'
    shell:
        'Rscript {SCRIPTS}/runHiCRep.R {output} {wildcards.bin} '
        '{params.start} {params.end} {input} &> {log}'


rule plotHiCRep:
    input:
        expand('qc/hicrep/data/{compare}-{{region}}-{{bin}}-hicrep-{{pm}}.csv',
            compare=HiC.sampleCompares())
    output:
        'qc/hicrep/{region}-{bin}-hicrep-{pm}.png'
    params:
        dpi = 600,
    group:
        'HiCRep'
    log:
        'logs/plotHiCRep/{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/plotHiCRep.py --out {output} --dpi {params.dpi} '
        '{input} &> {log}'


rule mergeBamByReplicate:
    input:
        lambda wc: expand(
            'dat/matrix/{{region}}/{{group}}-{rep}-{{region}}.bam',
            rep = HiC.groups()[wc.group]),
    output:
        'dat/matrix/{region}/{group}-{region}.bam'
    params:
        nonEmpty = nonEmpty
    log:
        'logs/mergeBamByReplicate/{region}/{group}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {params.nonEmpty} 2> {log}'


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


rule reformatSUTM:
    input:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}.gz'
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}-sutm.txt'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region],
    group:
        'HiCcompare'
    log:
        'logs/reformatSUTM/{region}/{bin}/{all}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/homer2sutm.py {input} --chrom {params.chr} '
        '--start {params.start} --binSize {wildcards.bin} '
        '> {output} 2> {log}'


rule HiCcompare:
    input:
        'dat/matrix/{region}/{bin}/raw/{group1}-{region}-{bin}-{pm}-sutm.txt',
        'dat/matrix/{region}/{bin}/raw/{group2}-{region}-{bin}-{pm}-sutm.txt'
    output:
        all = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{pm}.homer',
        links = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{pm}.links',
        adjIF1 = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-adjIF1-{pm}.homer',
        adjIF2 = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-adjIF2-{pm}.homer'
    params:
        dir = lambda wc: f'dat/HiCcompare/{wc.region}/{wc.bin}',
        qcdir = lambda wc: f'qc/HiCcompare/{wc.region}/{wc.bin}',
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        suffix = lambda wc: wc.pm,
    group:
        'HiCcompare'
    log:
        'logs/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{pm}.log'
    conda:
        f'{ENVS}/HiCcompare.yaml'
    shell:
        'Rscript {SCRIPTS}/HiCcompare.R {params.dir} {params.qcdir} '
        '{params.chr} {params.start} {params.end} '
        '{wildcards.bin} {params.suffix} {input} &> {log}'


rule applyMedianFilter:
    input:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{pm}.homer'
    output:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-logFC-{pm}.homer'
    params:
        size = config['compareMatrices']['size']
    group:
        'HiCcompare'
    log:
        'logs/applyMedianFilter/{compare}/{region}/{bin}/{group1}-vs-{group2}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/smoothHiC.py {input} --size {params.size} '
        '> {output} 2> {log}'


rule hicCompareBedgraph:
    input:
        rules.applyMedianFilter.output
    output:
        all = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-all-{pm}.bedgraph',
        up = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-up-{pm}.bedgraph',
        down = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-down-{pm}.bedgraph'
    params:
        maxDistance = config['compareMatrices']['maxDistance']
    group:
        'HiCcompare'
    log:
        'logs/hicCompareBedgraph/{compare}/{region}/{bin}/{group1}-vs-{group2}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/hicCompareBedgraph.py --Z '
        '--maxDistance {params.maxDistance} --allOut {output.all} '
        '--upOut {output.up} --downOut {output.down} {input} &> {log}'


rule rescaleHiCcompare:
    input:
        bedgraphs = lambda wc: expand(
            'dat/{{compare}}/{region}/{{bin}}/{{group1}}-vs-{{group2}}-{{dir}}-{{pm}}.bedgraph',
            region=binRegion[wc.bin]),
        chromSizes = getChromSizes
    output:
        'intervals/{group1}-vs-{group2}-{dir}-{compare}-{bin}-{pm}.pkl'
    params:
        regions = config['regions'],
        name = lambda wc: f'{wc.group1}-vs-{wc.group2}-{wc.dir}-{wc.compare}-{wc.bin}',
    group:
        'rescaleBedgraphs'
    log:
        'logs/rescaleHiCcompare/{group1}-vs-{group2}-{dir}-{compare}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/rescaleBedgraph.py sum --name {params.name} '
        '--out {output} --binSize {wildcards.bin} --regions {params.regions} '
        '--filetype bedgraph {input.chromSizes} <(cat {input.bedgraphs}) '
        '&> {log}'


rule rescaleHiCcompareCount:
    input:
        bedgraphs = lambda wc: expand(
            'dat/{{compare}}/{region}/{{bin}}/{{group1}}-vs-{{group2}}-{{dir}}-{{pm}}.bedgraph',
            region=binRegion[wc.bin]),
        chromSizes = getChromSizes
    output:
        'intervals/{group1}-vs-{group2}-{dir}-{compare}-count-{bin}-{pm}.pkl'
    params:
        threshold = config['compareMatrices']['minZ'],
        regions = config['regions'],
        name = lambda wc: f'{wc.group1}-vs-{wc.group2}-{wc.dir}-{wc.compare}-{wc.bin}-peak',
    group:
        'rescaleBedgraphs'
    log:
        'logs/rescaleHiCcompareBinary/{group1}-vs-{group2}-{dir}-{compare}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/rescaleBedgraph.py count --name {params.name} '
        '--out {output} --binSize {wildcards.bin} --regions {params.regions} '
        '--filetype bedgraph --threshold {params.threshold} '
        '{input.chromSizes} <(cat {input.bedgraphs}) &> {log}'


rule homerToH5:
    input:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}-{pm}.homer'
    output:
        'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}-{pm}.h5'
    group:
        'HiCcompare'
    log:
        'logs/homerToH5/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '(hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 || touch {output})  &> {log}'


if config['compareMatrices']['tads'] is not None:
    rule filterTADs:
        input:
            config['compareMatrices']['tads']
        output:
            'dat/tads/{region}/{region}-referenceTADs.bed'
        params:
            chr = lambda wc: REGIONS['chr'][wc.region],
            start = lambda wc: REGIONS['start'][wc.region],
            end = lambda wc: REGIONS['end'][wc.region]
        group:
            'HiCcompare'
        log:
            'logs/filterTADs/{region}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            'python {SCRIPTS}/filterTADs.py --chrom {params.chr} '
            '--start {params.start} --end {params.end} '
            '{input} > {output} 2> {log}'


def setControl(wc):
    """ Set control matrix as the other IF1 relative to the target """
    if wc.adjIF == 'adjIF1':
        adj = 'adjIF2'
    else:
        adj = 'adjIF1'
    return f'dat/HiCcompare/{{region}}/{{bin}}/{{group1}}-vs-{{group2}}-{adj}-{{pm}}.h5'


def setDomains(wc):
    """ Set TADdomains as same as target matrix e.g. IF1 = group1 """
    if config['compareMatrices']['tads'] is not None:
        return f'dat/tads/{wc.region}/{wc.region}-referenceTADs.bed'
    elif wc.adjIF == 'adjIF1':
        group = wc.group1
    else:
        group = wc.group2
    return f'dat/matrix/{{region}}/{{bin}}/tads/{group}-{{region}}-{{bin}}-{{pm}}-ontad_domains.bed'


rule differentialTAD:
    input:
        target = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{adjIF}-{pm}.h5',
        control = setControl,
        tadDomains = setDomains
    output:
        accepted = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}_accepted.diff_tad',
        rejected = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}_rejected.diff_tad'
    params:
        pValue = 0.05,
        mode = 'all',
        modeReject = 'one',
        chr = lambda wc: REGIONS['chr'][wc.region],
        prefix = lambda wc: f'dat/tads/{wc.region}/{wc.bin}/{wc.group1}-vs-{wc.group2}-{wc.region}-{wc.bin}-{wc.adjIF}-{wc.pm}'
    group:
        'HiCcompare'
    log:
        'logs/differentialTAD/{region}/{bin}/{group1}-vs-{group2}-{adjIF}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicDifferentialTAD --targetMatrix {input.target} '
        '--controlMatrix {input.control} '
        '--tadDomains {input.tadDomains} '
        '--pValue {params.pValue} --mode {params.mode} '
        '--modeReject {params.modeReject} --outFileNamePrefix {params.prefix} '
        ' &> {log} || touch {output} '


rule reformatDifferentialTAD:
    input:
        'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}_{result}.diff_tad'
    output:
        'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}_{result}_domains.bed'
    group:
        'HiCcompare'
    log:
        'logs/reformatDifferentialTAD/{region}/{bin}/{group1}-vs-{group2}-{adjIF}-{pm}-{result}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/reformatDomains.py {input} --scale 0 '
        '> {output} 2> {log}'


rule createCompareConfig:
    input:
        mat = 'dat/{compare}/{region}/{bin}/{group1}-vs-{group2}-{set}-{pm}.h5',
        upBed = rules.hicCompareBedgraph.output.up,
        downBed = rules.hicCompareBedgraph.output.down,
        tads1 = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF1-{pm}_rejected_domains.bed',
        tads2 = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF2-{pm}_rejected_domains.bed'
    output:
        'plots/{region}/{bin}/HiCcompare/configs/{group1}-vs-{group2}-{compare}-{set}-{pm}.ini',
    params:
        depth = lambda wc: int(REGIONS['length'][wc.region]),
        colourmap = 'bwr',
        tracks = getTracks,
        threshold = config['compareMatrices']['minZ'],
        sumLogFC_title = f'"sum(logFC) ({config["compareMatrices"]["maxDistance"]:.1e}bp) threshold = {config["compareMatrices"]["minZ"]}"',
        vMin = config['compareMatrices']['vMin'],
        vMax = config['compareMatrices']['vMax'],
    group:
        'HiCcompare'
    log:
        'logs/createCompareConfig/{compare}/{region}/{bin}/{group1}-{group2}-{set}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/generate_config.py --matrix {input.mat} --compare '
        '--sumLogFC {input.upBed} {input.downBed} '
        '--sumLogFC_title {params.sumLogFC_title} '
        '--sumLogFC_hline {params.threshold} '
        '--tads {input.tads1} {input.tads2} '
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


def setCompareTitle(wc):
    if config['build'] is not None:
        build = config['build'].replace('"', '') # Double quotes disallowed
        build = f' ({build})'
    else:
        build = ''
    title = (f'"{wc.group1} vs {wc.group2} - {wc.region}{build} at '
             f'{wc.bin} bin size - adj. logFC - {wc.pm}"')
    return title


rule plotCompare:
    input:
        rules.createCompareConfig.output
    output:
        'plots/{region}/{bin}/{compare}/{set}/{group1}-vs-{group2}-{region}-{coord}-{bin}-{set}-{pm}.png'
    params:
        title = setCompareTitle,
        region = setRegion,
        dpi = 600
    group:
        'HiCcompare'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    log:
        'logs/plotAnalysis/{compare}/{region}/{bin}/{group1}-vs-{group2}-{coord}-{set}-{pm}.log'
    threads:
        THREADS
    shell:
        'export NUMEXPR_MAX_THREADS=1; pyGenomeTracks --tracks {input} '
        '--region {params.region} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


rule runCompareViewpoint:
    input:
        'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{set}-{pm}.h5'
    output:
        plot = 'dat/HiCcompare/{region}/{bin}/viewpoint/{group1}-vs-{group2}-{set}-{region}-{coord}-{bin}-{pm}.png',
        bedgraph = 'dat/HiCcompare/{region}/{bin}/viewpoint/{group1}-vs-{group2}-{set}-{region}-{coord}-{bin}-{pm}.bedgraph',
    params:
        referencePoint = setRegion,
        region = makeViewRegion,
        dpi = 600
    group:
        'HiCcompare'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    log:
        'logs/runCompareViewpoint/{group1}-vs-{group2}-{set}-{region}-{coord}-{bin}-{pm}.log'
    shell:
        '(hicPlotViewpoint --matrix {input} --region {params.region} '
        '--outFileName {output.plot} --referencePoint {params.referencePoint} '
        '--interactionOutFileName {output.bedgraph} --dpi {params.dpi} '
        '&& mv {output.bedgraph}_*.bedgraph {output.bedgraph}) &> {log}'


rule plotCompareViewpoint:
    input:
        IF1 = 'dat/HiCcompare/{region}/{bin}/viewpoint/{group1}-vs-{group2}-adjIF1-{region}-{coord}-{bin}-{pm}.bedgraph',
        IF2 = 'dat/HiCcompare/{region}/{bin}/viewpoint/{group1}-vs-{group2}-adjIF2-{region}-{coord}-{bin}-{pm}.bedgraph',
    output:
        'plots/{region}/{bin}/viewpoints/HiCcompare/{group1}-vs-{group2}-{region}-{coord}-{bin}-viewpoint-{pm}.png',
    params:
        dpi = 600,
        build = f'--build {config["build"]}' if config['build'] else ''
    group:
        'HiCcompare'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/plotCompareViewpoint/{group1}-vs-{group2}-{region}-{coord}-{bin}-{pm}.log'
    shell:
        'python {SCRIPTS}/plotViewpoint.py {input} --out {output} '
        '--dpi {params.dpi} {params.build} &> {log}'


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
            'samtools view -b -@ {threads} -s {params.fraction} {input} '
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
            maxGaussians = config['gatk']['maxGaussians'],
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
            '--max-gaussians {params.maxGaussians} '
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
            maxGaussians = config['gatk']['maxGaussians'],
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
            '--max-gaussians {params.maxGaussians} '
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
            end = lambda wc: REGIONS['end'][wc.region],
            maxfragments = 20000000
        group:
            'hapcut2'
        log:
            'logs/extractHAIRS/{region}/{cellType}.log'
        conda:
            f'{ENVS}/hapcut2.yaml'
        shell:
            'extractHAIRS --hic 1 --bam {input.bam} '
            '--maxfragments {params.maxfragments} '
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
            'python {SCRIPTS}/extractBestHapcut2.py {input} > {output} 2> {log}'


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
            'mergeVCFsbyRegion'
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
        temp('dat/mapped/subsampled/{preSample}.subsampled.sam')
    group:
        'filterQC'
    params:
        nLines = 1000000 * 2
    log:
        'logs/sampleReads/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'head -n {params.nLines} <(samtools view {input}) > {output} 2> {log}'


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
        'python {SCRIPTS}/processHiC.py {input.digest[0]} {input.reads} '
        '> {output} 2> {log}'


rule plotQC:
    input:
        expand('dat/mapped/subsampled/{sample}-processed.txt',
            sample=HiC.originalSamples())
    output:
        ditagOut = 'qc/filterQC/ditagLength.png',
        insertOut = 'qc/filterQC/insertSizeFrequency.png'
    group:
        'filterQC' if config['groupJobs'] else 'plotQC'
    log:
        'logs/plotQC.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/plotQC.py {input} --insertOut {output.insertOut} '
        '--ditagOut {output.ditagOut} &> {log}'


rule mergeHicupQC:
    input:
        rules.hicupTruncate.output.summary,
    output:
        'qc/hicup/HiCUP_summary_report-{preSample}.txt'
    log:
        'logs/mergeHicupQC/{preSample}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/mergeHicupSummary.py --truncater {input} '
        '> {output} 2> {log}'

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
            sample=HiC.originalSamples()) if not config['microC'] else [],
         expand('qc/bowtie2/{sample}-{read}.bowtie2.txt',
            sample=HiC.originalSamples(), read=['R1', 'R2']),
         expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
            sample=HiC.originalSamples(), read=['R1', 'R2']) if config['fastq_screen'] else [],
         expand('qc/bcftools/{region}/{cellType}-{region}-bcftoolsStats.txt',
            region=REGIONS.index, cellType=HiC.cellTypes()) if PHASE_MODE=='BCFTOOLS' else []],
        [expand('qc/hicexplorer/{sample}-{region}.{bin}-{pm}_QC', region=region,
            sample=HiC.samples(), bin=BASE_BIN, pm=phaseMode) for region in regionBin] if not config['ASHIC'] else [],
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
