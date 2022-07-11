#!/usr/bin/env python3

import os
import re
import sys
import math
import random
import tempfile
import itertools
import numpy as np
from snake_setup import set_config, load_regions, load_coords, filterRegions, HiCSamples, getValidBins

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
    'data':              ''          ,
    'phasedVCF':         None        ,
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
        {'minBins':              50    ,
         'minDistance':          300   ,
         'maxLibraryInsertSize': 1000  ,
         'minMappingQuality':    15    ,
         'keepSelfLigation':     False ,
         'keepSelfCircles':      False ,
         'nofill':               False ,
         'makeBam':              False ,
         'compartmentScore':     False ,
         'threads':              4     ,
         'multiplicativeValue':  10000 ,
         'microC':               False ,},
    'loops':
        {'detectLoops':               False  ,
         'peakWidth':                 6      ,
         'windowSize':                10     ,
         'pValuePreselection':        0.1    ,
         'peakInteractionsThreshold': 10     ,
         'obsExpThreshold':           1.5    ,
         'pValue':                    0.025  ,
         'maxLoopDistance':           2000000,
         'plot':                      False  ,
         'expected':                 'mean'  ,},
    'compareMatrices':
        {'colourmap'    : 'bwr'        ,
         'zThreshold'   : 2            ,
         'vMax'         : 1            ,
         'tads'         : None         ,
         'loops'        : None         ,
         'allPairs'     : False        ,
         'minSum'       : 25           ,},
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
        {'base':          10000        ,
         'bins':         [10000, 20000],},
    'plotParams':
        {'distanceNorm'  : False    ,
         'colourmap'     : 'Purples',
         'coordinates'   : None     ,
         'viewpoints'    : None     ,
         'viewpointRange': 500000   ,
         'plotRep'       : True     ,
         'vLines'        : []       ,
         'includeRegions': True     ,
         'filetype'      : 'svg'    ,
         'geneLabelFontSize': 12    ,
         'maxLabels'        : 60    ,},
    'Genes':
        {'gff3'        :  []        ,
         'typeKey'     : 'gene_type',
         'label'       : 'gene_name',
         'geneID'      : 'gene_id'  ,},
    'bigWig'           : {}        ,
    'bed'              : {}        ,
    'QC':
        {'runQC'        : True    ,
         'flipSNP'      : False   ,
         'QCsample'     : 100000  ,
         'fastqScreen'  : None    ,
         'multiQCconfig': None    ,
         'runHiCRep'    : True    ,
         'HiCRep_bin'   : 150000  ,},
    'phase':            False,

}

config = set_config(config, default_config)

workdir: config['workdir']
THREADS = workflow.cores

BASE_BIN = config['resolution']['base']
ALLELE_SPECIFIC = True if config['phasedVCF'] else False

HiC = HiCSamples(
    config['data'], config['restrictionSeqs'], ALLELE_SPECIFIC,
    allPairs=config['compareMatrices']['allPairs'])
adjust = 1 if not config['resolution']['bins'] else max(config['resolution']['bins'])
REGIONS = load_regions(config['regions'], adjust=adjust)

# Remove region-binSize combinations with too few bins
regionBin, binRegion = filterRegions(
    REGIONS, config['resolution']['bins'],
    nbins=config['HiCParams']['minBins'])

hicrepRegionsBin, _ = filterRegions(
    REGIONS, [config['QC']['HiCRep_bin']],
    nbins=config['HiCParams']['minBins'])


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
    pm = r'SNPsplit' if ALLELE_SPECIFIC else r'full',
    set = r'logFC|adjIF1|adjIF2',
    adjIF = r'adjIF1|adjIF2',
    filter = r'noFilter|medianFilter',
    compare = r'HiCcompare',
    type = rf'{config["plotParams"]["filetype"]}',
    group = rf'{"|".join(HiC.groups())}',
    group1 = rf'{"|".join(HiC.groups())}',
    group2 = rf'{"|".join(HiC.groups())}',
    sample = rf'{"|".join(HiC.samples(all=True))}',
    all = rf'{"|".join(HiC.samples() + list(HiC.groups()))}',
    combo = r'alt_alt|alt_both-ref|both-ref_both-ref|ref_alt|ref_both-ref|ref_ref',

# Generate dictionary of plot coordinates, may be multple per region
COORDS = load_coords(
    REGIONS, config['plotParams']['coordinates'],
    adjust=adjust,
    includeRegions=config['plotParams']['includeRegions'])

# Generate dictionary of plot viewpoints
VIEWPOINTS = load_coords(
    REGIONS, config['plotParams']['viewpoints'], includeRegions=False)

if config['phase'] and not ALLELE_SPECIFIC:
    if config['gatk']['all_known']:
        PHASE_MODE = 'GATK'
    else:
        PHASE_MODE = 'BCFTOOLS'
else:
    PHASE_MODE = None

# Set plot suffix for allele specific mode
pm = 'SNPsplit' if ALLELE_SPECIFIC else 'full'


HiC_mode = ([
    [expand('plots/{region}/{bin}/HiCsubtract/{filter}/{compare}-{region}-{coords}-{bin}-LOESSdiff-{filter}-{pm}.{type}',
        region=region, coords=COORDS[region], pm=pm, compare=HiC.groupCompares(), filter=['medianFilter', 'noFilter'],
        bin=regionBin[region], type=config['plotParams']['filetype']) for region in regionBin],
    [expand('plots/{region}/{bin}/pyGenomeTracks/{group}-{region}-{coords}-{bin}-{pm}.{type}',
        region=region, coords=COORDS[region], pm=pm, group=HiC.groups(),
        bin=regionBin[region], type=config['plotParams']['filetype']) for region in regionBin],
    [expand('plots/{region}/{bin}/HiCsubtract/configs/{compare}-{coords}-LOESSdiff-{filter}-{pm}.ini',
        region=region, coords=COORDS[region], pm=pm, compare=HiC.groupCompares(), filter=['medianFilter', 'noFilter'],
        bin=regionBin[region], type=config['plotParams']['filetype']) for region in regionBin],
    [expand('plots/{region}/{bin}/pyGenomeTracks/configs/{group}-{region}-{coords}-{bin}-{pm}.ini',
        region=region, coords=COORDS[region], pm=pm, group=HiC.groups(),
        bin=regionBin[region], type=config['plotParams']['filetype']) for region in regionBin],
    [expand('plots/{region}/{bin}/viewpoints/HiCcompare/{compare}-{region}-{coords}-{bin}-viewpoint-{pm}.{type}',
        region=region, coords=VIEWPOINTS[region], pm=pm, compare=HiC.groupCompares(),
        bin=regionBin[region], type=config['plotParams']['filetype']) for region in regionBin],
    [expand('plots/{region}/{bin}/viewpoints/{preGroup}-{region}-{coords}-{bin}-viewpoint-{pm}.{type}',
        region=region, coords=VIEWPOINTS[region], pm=pm, preGroup=HiC.groups(),
        bin=regionBin[region], type=config['plotParams']['filetype']) for region in regionBin],
    [expand('plots/{region}/{bin}/obs_exp/{all}-{region}-{bin}-{pm}.{type}',
        all=(HiC.all() if config['plotParams']['plotRep'] else list(HiC.groups())),
        region=region, bin=regionBin[region], pm=pm,
        type=config['plotParams']['filetype']) for region in regionBin],
     expand('referenceTADs/{all}-{bin}-ontad-{pm}.bed',
        all=(HiC.all() if config['plotParams']['plotRep'] else list(HiC.groups())),
        bin=getValidBins(regionBin), pm=pm),
     expand('qc/matrixCoverage/{region}/{all}-coverage-{pm}.{type}',
        all=(HiC.all() if config['plotParams']['plotRep'] else list(HiC.groups())),
        region=regionBin.keys(), pm=pm, type=config['plotParams']['filetype'])])


methods = ['TADinsulation', 'TADboundaries', 'TADdomains']


def getChromSizes(wc):
    """ Retrieve chromSizes file associated with group or sample. """
    try:
        cellType = HiC.sample2Cell()[wc.all]
    except AttributeError:
        try:
            cellType = HiC.sample2Cell()[wc.preGroup]
        except AttributeError:
            try:
                cellType = HiC.sample2Cell()[wc.sample]
            except AttributeError:
                cellType = HiC.sample2Cell()[wc.group1]
                cellType2 = HiC.sample2Cell()[wc.group2]
                if cellType != cellType2:
                    sys.stderr.write(
                        f'{wc.group1} and {wc.group2} correspond to different '
                        'cell type. Ensure the chromosome sizes are equal for '
                        'valid bedgraph rescaling of HiCcompare output.\n')
    return f'dat/genome/chrom_sizes/{cellType}.chrom.sizes',

if config['HiCParams']['compartmentScore']:
    include: 'CscoreTool.snake'
    compartmentOutput = (
        [expand('dat/Cscore/{region}/{bin}/{group}-{region}-{bin}-{pm}-Cscore{ext}',
            region=region, bin=regionBin[region], group=list(HiC.groups()),
            pm=pm, ext=['bias.txt', '_hh.txt', '_cscore.txt', '_cscore.bedgraph'])
        for region in regionBin])
else:
    compartmentOutput = []

if len(HiC.sampleCompares()) == 0:
    config['QC']['runHiCRep'] = False

rule all:
    input:
        HiC_mode,
        compartmentOutput,
        (expand('phasedVCFs/{cellType}-phased.vcf', cellType=HiC.cellTypes())
            if PHASE_MODE is not None else []),
        ([f'qc/filterQC/insertSizeFrequency.{config["plotParams"]["filetype"]}',
        'qc/multiqc'] if config['QC']['runQC'] else []),
        ([expand('qc/hicrep/{region}-{bin}-hicrep-{pm}.{type}',
            bin=hicrepRegionsBin[region], pm=pm, region=region,
            type=config['plotParams']['filetype'])
            for region in hicrepRegionsBin] if config['QC']['runHiCRep'] else [])


if ALLELE_SPECIFIC:

    rule filterHomozygous:
        input:
            lambda wc: config['phasedVCF'][wc.cellType]
        output:
            'dat/genome/{cellType}-phasedHet.vcf'
        group:
            'prepareGenome'
        log:
            'logs/filterHomozygous/{cellType}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools view -H -m 2 -M 2 -v snps -g het --phased '
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
        params:
            flip = lambda wc: '--flipSNP' if config['QC']['flipSNP'] else ''
        group:
            'prepareGenome'
        log:
            'logs/vcf2SNPsplit/{cellType}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            'python {SCRIPTS}/reformatSNPsplit.py {params.flip} {input} '
            '> {output} 2> {log}'


    rule SNPcoverage:
        input:
            rules.filterHomozygous.output
        output:
            'dat/genome/SNPcoverage/{cellType}-{bin}-phasedHet.bed'
        group:
            'prepareGenome'
        log:
            'logs/SNPcoverage/{cellType}-{bin}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            'python {SCRIPTS}/SNPcoverage.py {input} --binSize {wildcards.bin} '
            '> {output} 2> {log}'


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
        'logs/indexGenome/{cellType}.log'
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


rule gff3ToGenePred:
    input:
        config['Genes']['gff3']
    output:
        genePred = f'annotation/{config["build"]}.unsorted.genePred',
        attributes = f'annotation/{config["build"]}.attrs'
    log:
        'logs/gff3ToGenePred.log'
    conda:
        f'{ENVS}/UCSCtools.yaml'
    shell:
        'gff3ToGenePred -attrsOut={output.attributes} {input} '
        '{output.genePred} &> {log}'


rule sortGenePred:
    input:
        f'annotation/{config["build"]}.unsorted.genePred'
    output:
        f'annotation/{config["build"]}.genePred'
    log:
        'logs/sortGenePred.log'
    conda:
        f'{ENVS}/UCSCtools.yaml'
    shell:
        'sort -k2,2 -k4n,4n {input} > {output} 2> {log}'


rule genePredToBed:
    input:
        rules.sortGenePred.output
    output:
        f'annotation/{config["build"]}.bed12'
    log:
        'logs/genePredToBed.log'
    conda:
        f'{ENVS}/UCSCtools.yaml'
    shell:
        'genePredToBed {input} {output} &> {log}'


rule processBED12:
    input:
        bed12 = rules.genePredToBed.output,
        attrs = rules.gff3ToGenePred.output.attributes
    output:
        f'annotation/{config["build"]}-coloured.bed12'
    params:
        label = config['Genes']['label'],
        typeKey = config['Genes']['typeKey'],
        geneID = config['Genes']['geneID']
    group:
        'prepareGenome'
    log:
        'logs/processBED12.log'
    shell:
        'python {SCRIPTS}/processBED12.py --label {params.label} '
        '--typeKey {params.typeKey} --geneID {params.geneID} '
        '{input.bed12} {input.attrs} > {output} 2> {log}'


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
        'logs/fastQC/{preSample}-{read}.log'
    conda:
        f'{ENVS}/fastqc.yaml'
    shell:
        'python3 {SCRIPTS}/fastqc.py {input} --htmlOut {output.html} '
        '--dataOut {output.zip} &> {log}'


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
        'dat/fastq/{preSample}-{read}.trimmed.fastq.gz'
    output:
        html = 'qc/fastqc/{preSample}-{read}.trim_fastqc.html',
        zip = 'qc/fastqc/{preSample}-{read}.trim_fastqc.zip'
    group:
        'fastqc'
    log:
        'logs/fastQCTrimmed/{preSample}-{read}.log'
    conda:
        f'{ENVS}/fastqc.yaml'
    shell:
        'python3 {SCRIPTS}/fastqc.py {input} --htmlOut {output.html} '
        '--dataOut {output.zip} &> {log}'


if config['QC']['fastqScreen'] is not None:

    rule fastQScreen:
        input:
            'dat/fastq/{preSample}-{read}.trimmed.fastq.gz'
        output:
            txt = 'qc/fastqScreen/{preSample}-{read}_screen.txt',
            png = 'qc/fastqScreen/{preSample}-{read}.fastqScreen.png'
        params:
            config = config['QC']['fastqScreen'],
            subset = 100000,
        log:
            'logs/fastqScreen/{preSample}-{read}.log'
        conda:
            f'{ENVS}/fastqScreen.yaml'
        threads:
            THREADS
        shell:
            'python {SCRIPTS}/fastqScreen.py {input} {params.config} '
            '--subset {params.subset} --threads {threads} '
            '--plotOut {output.png} --dataOut {output.txt} &> {log}'


rule cutadapt:
    input:
        lambda wc: HiC.path(wc.preSample, ['R1', 'R2'])
    output:
        trimmed = [temp('dat/fastq/{preSample}-R1.trimmed.fastq.gz'),
                   temp('dat/fastq/{preSample}-R2.trimmed.fastq.gz')],
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
        truncated = [temp('dat/fastq/{preSample}-R1.truncated.fastq.gz'),
                     temp('dat/fastq/{preSample}-R2.truncated.fastq.gz')],
        summary = 'qc/hicup/HiCUP_summary_report-{preSample}.txt'
    params:
        re1 = lambda wc: list(HiC.restrictionSeqs()[wc.preSample].values())[0],
        fill = '--nofill' if config['HiCParams']['nofill'] else '',
        tmpdir = config['tmpdir']
    group:
        'hicupTruncate'
    threads:
        2 if THREADS > 2 else THREADS
    log:
        'logs/hicupTruncate/{preSample}.log'
    conda:
        f'{ENVS}/hicup.yaml'
    shell:
        'python {SCRIPTS}/hicupTruncate.py {params.fill} --re1 {params.re1} '
        '--hicup {SCRIPTS}/hicup_truncater --output {output.truncated} '
        '--threads {threads} --tmpdir {params.tmpdir} '
        '{input} > {output.summary} 2> {log}'


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
    mode = 'trimmed' if config['HiCParams']['microC'] else 'truncated'
    return f'dat/fastq/{wc.preSample}-{wc.read}.{mode}.fastq.gz'


rule bowtie2:
    input:
        fastq = fastqInput,
        bt2_index = bowtie2Index
    output:
        sam = pipe('dat/mapped/{preSample}-{read}.sam'),
        qc = 'qc/bowtie2/{preSample}-{read}.bowtie2.txt'
    params:
        index = bowtie2Basename,
        sensitivity = 'sensitive',
        cellType = lambda wc: HiC.sample2Cell()[wc.preSample]
    group:
        'bowtie2'
    log:
        'logs/bowtie2/{preSample}-{read}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        THREADS - 2 if THREADS > 2 else 1
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
        'awk -f {SCRIPTS}/addReadFlag.awk -v flag={params.flag} '
        '{input} > {output} 2> {log}'


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


rule mergeBam:
    input:
        'dat/mapped/{preSample}-R1-addFlag.bam',
        'dat/mapped/{preSample}-R2-addFlag.bam'
    output:
        temp('dat/mapped/{preSample}-merged.bam')
    group:
        'prepareBAM'
    log:
        'logs/mergeBam/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} -n {output} {input} &> {log}'


rule fixmateBam:
    input:
        rules.mergeBam.output
    output:
        pipe('dat/mapped/{preSample}.fixed.bam')
    group:
        'prepareBAM'
    log:
        'logs/fixmateBam/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools fixmate -p -O SAM {input} {output} &> {log}'


rule samblaster:
    input:
        rules.fixmateBam.output
    output:
        sam = pipe('dat/mapped/{preSample}-dedup.sam'),
        qc = 'qc/deduplicate/{preSample}.txt'
    group:
        'prepareBAM'
    log:
        'logs/samblaster/{preSample}.log'
    conda:
        f'{ENVS}/samblaster.yaml'
    shell:
        'samblaster --removeDups --input {input} > {output.sam} '
        '2> {log} && cp {log} {output.qc}'


rule sam2bam2:
    input:
        rules.samblaster.output.sam
    output:
        temp('dat/mapped/{preSample}-dedup.bam')
    group:
        'prepareBAM'
    log:
        'logs/sam2bam2/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools view -b {input} > {output} 2> {log}'


rule removeUnmapped:
    input:
        rules.sam2bam2.output
    output:
        'dat/mapped/{preSample}.hic.bam'
    log:
        'logs/removeUnmapped/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools view -b -@ {threads} -F 2316 {input} > {output} 2> {log}'


rule samtoolsStats:
    input:
        rules.sam2bam2.output
    output:
        'qc/samtools/stats/{preSample}.stats.txt'
    group:
        'samQC'
    log:
        'logs/samtoolsStats/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'


def SNPsplitInput(wc):
    """ Retrieve cell type associated with sample. """
    cellType = HiC.sample2Cell()[wc.preSample]
    return f'snpsplit/{cellType}-snpsplit.txt'


rule SNPsplit:
    input:
        bam = rules.removeUnmapped.output,
        snps = SNPsplitInput
    output:
        report = 'dat/snpsplit/{preSample}.hic.SNPsplit_report.txt',
        bam = 'dat/snpsplit/{preSample}.hic.allele_flagged.bam'
    params:
        outdir = 'dat/snpsplit/'
    group:
        'SNPsplit'
    log:
        'logs/SNPsplit/{preSample}.log'
    conda:
        f'{ENVS}/snpsplit.yaml'
    shell:
        'SNPsplit {input.bam} --snp_file {input.snps} --hic '
        '-skip_tag2sort --output_dir {params.outdir} &> {log}'


rule tag2sort:
    input:
        rules.SNPsplit.output.bam
    output:
        report = 'dat/snpsplit/{preSample}.hic.SNPsplit_sort.txt',
        bam = expand('dat/snpsplit/{{preSample}}.hic.{ext}',
            ext = ['G1_G1.bam', 'G1_UA.bam', 'G2_G2.bam',
                   'G2_UA.bam', 'G1_G2.bam', 'UA_UA.bam'])
    params:
        outdir = 'dat/snpsplit/',
        bam = lambda wc: f'{wc.preSample}.hic.allele_flagged.bam'
    group:
        'SNPsplit'
    log:
        'logs/tag2sort/{preSample}.log'
    conda:
        f'{ENVS}/snpsplit.yaml'
    shell:
        'tag2sort {params.bam} --hic --output_dir {params.outdir} &> {log}'


rule mergeSNPsplit:
    input:
        'dat/snpsplit/{preGroup}-{rep}.hic.G{allele}_G{allele}.bam',
        'dat/snpsplit/{preGroup}-{rep}.hic.G{allele}_UA.bam'
    output:
        'dat/snpsplit/merged/{preGroup}_a{allele}-{rep}.hic.bam'
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


def splitInput(wc):
    if ALLELE_SPECIFIC:
        return 'dat/snpsplit/merged/{sample}.hic.bam'
    else:
        return 'dat/mapped/{sample}.hic.bam'


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
    if config['HiCParams']['microC']:
        return ['dat/genome/restSites-empty.bed']
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


def outBam():
    if config['HiCParams']['makeBam']:
        return 'dat/matrix/{sample}-{pm}.bam'
    else:
        return []


rule buildBaseMatrix:
    input:
        bams = ancient(expand('dat/mapped/split/{{sample}}-{read}-{{pm}}.bam', read=['R1', 'R2'])),
        restSites = getRestSites,
        chromSizes = getChromSizes
    output:
        hic = f'dat/matrix/base/raw/{{sample}}.{BASE_BIN}-{{pm}}.h5',
        bam = outBam(),
        qc = directory(f'qc/hicexplorer/{{sample}}.{BASE_BIN}-{{pm}}_QC'),
        qc_table = f'qc/hicexplorer/{{sample}}.{BASE_BIN}-{{pm}}_QC/QC_table.txt'
    params:
        bin = BASE_BIN,
        reSeqs = getRestrictionSeqs,
        inputBufferSize = 400000,
        danglingSequences = getDanglingSequences,
        maxLibraryInsertSize = config['HiCParams']['maxLibraryInsertSize'],
        minMappingQuality = config['HiCParams']['minMappingQuality'],
        minDistance = config['HiCParams']['minDistance'],
        keepSelfLigation = (
            '--keepSelfLigation' if config['HiCParams']['keepSelfLigation'] else ''),
        keepSelfCircles = (
            '--keepSelfCircles' if config['HiCParams']['keepSelfCircles'] else ''),
    group:
        'buildBaseMatrix'
    log:
        'logs/buildBaseMatrix/{sample}-{pm}.log'
    threads:
        max(2, config['HiCParams']['threads'])
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicBuildMatrix --samFiles {input.bams} '
        '--restrictionCutFile {input.restSites} '
        '--restrictionSequence {params.reSeqs} '
        '--maxLibraryInsertSize {params.maxLibraryInsertSize} '
        '--minDistance {params.minDistance} '
        '--minMappingQuality {params.minMappingQuality} '
        '--danglingSequence {params.danglingSequences} '
        '--inputBufferSize {params.inputBufferSize} '
        '--chromosomeSizes {input.chromSizes} '
        '{params.keepSelfCircles} {params.keepSelfLigation} '
        '--binSize {params.bin} '
        '--outFileName {output.hic} --skipDuplicationCheck '
        '--QCfolder {output.qc} --threads {threads} '
        f'{"--outBam {output.bam} " if config["HiCParams"]["makeBam"] else ""}'
        '&> {log}'


rule adjustMatrix:
    input:
        rules.buildBaseMatrix.output.hic
    output:
        f'dat/matrix/{{region}}/base/raw/{{sample}}-{{region}}.{BASE_BIN}-{{pm}}.h5'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'buildBaseMatrix'
    log:
        'logs/adjustMatrix/{sample}-{region}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicAdjustMatrix --matrix {input} --outFileName {output} '
        '--regions <(echo -e \'{params.chr}\t{params.start}\t{params.end}\') '
        '&> {log}'


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
        '&> {log}'


rule mergeBins:
    input:
        f'dat/matrix/{{region}}/base/raw/{{all}}-{{region}}.{BASE_BIN}-{{pm}}.h5'
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}-raw.h5'
    params:
        nbins = lambda wc: int(int(wc.bin) / BASE_BIN)
    log:
        'logs/mergeBins/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log}'


rule normCounts01:
    input:
        rules.mergeBins.output
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}-norm01.h5'
    log:
        'logs/normCounts01/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicNormalize --matrices {input} --normalize norm_range '
        '--outFileName {output} &> {log}'


rule normCountsConstant:
    input:
        rules.normCounts01.output
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}.h5'
    params:
        multiplicativeValue = config['HiCParams']['multiplicativeValue']
    log:
        'logs/normCountsConstant/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicNormalize --matrices {input} --normalize multiplicative '
        '--multiplicativeValue {params.multiplicativeValue} '
        '--outFileName {output} &> {log}'


rule correctMatrix:
    input:
        rules.normCountsConstant.output
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
        '--outFileName {output} &> {log}'


rule TADinsulation:
    input:
        rules.correctMatrix.output
    output:
        expand(
            'dat/tads/{{region}}/{{bin}}/{{all}}-{{region}}-{{bin}}-{{pm}}{ext}',
            ext = ['_boundaries.bed', '_boundaries.gff', '_domains.bed',
                   '_score.bedgraph', '_zscore_matrix.h5']),
        score = 'dat/tads/{region}/{bin}/{all}-{region}-{bin}-{pm}_tad_score.bm'
    params:
        method = 'fdr',
        bin = lambda wc: wc.bin,
        region = lambda wc: wc.region,
        all = lambda wc: wc.all,
        min_depth = lambda wc: int(wc.bin) * 3,
        max_depth = lambda wc: int(wc.bin) * 10,
        prefix = 'dat/tads/{region}/{bin}/{all}-{region}-{bin}-{pm}'
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
        '--numberOfProcessors {threads} &> {log}'


rule detectLoops:
    input:
        rules.correctMatrix.output
    output:
        'dat/loops/{region}/{bin}/unmod/{all}-{region}-{bin}-{pm}.bedgraph'
    params:
        peakWidth = config['loops']['peakWidth'],
        windowSize = config['loops']['windowSize'],
        pValuePreselection = config['loops']['pValuePreselection'],
        peakInteractionsThreshold = config['loops']['peakInteractionsThreshold'],
        obsExpThreshold = config['loops']['obsExpThreshold'],
        pValue = config['loops']['pValue'],
        maxLoopDistance = config['loops']['maxLoopDistance'],
        expected = config['loops']['expected']
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
        '--pValuePreselection {params.pValuePreselection} '
        '--obsExpThreshold {params.obsExpThreshold} '
        '--pValue {params.pValue} '
        '--peakInteractionsThreshold {params.peakInteractionsThreshold} '
        '--expected {params.expected} '
        '--threads 1 --threadsPerChromosome {threads} '
        '&> {log} || touch {output}'

# Fix to correct loop intervals that extend too far
rule clampLoops:
    input:
        rules.detectLoops.output
    output:
        'dat/loops/{region}/{bin}/{all}-{region}-{bin}-{pm}.bedgraph'
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
        'python {SCRIPTS}/clampLoops.py {params.minPos} {params.maxPos} '
        '{input} > {output} 2> {log}'


rule H5_to_NxN:
    input:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.h5'
    output:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.nxn.tsv'
    group:
        'processHiC'
    log:
        'logs/H5_to_NxN/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/H5_to_NxN.py {input} > {output} 2> {log}'


rule plotCoverage:
    input:
        lambda wc: expand('dat/matrix/{{region}}/{bin}/raw/{{all}}-{{region}}-{bin}-{{pm}}.h5',
            bin=regionBin[wc.region])
    output:
        'qc/matrixCoverage/{region}/{all}-coverage-{pm}.{type}'
    params:
        dpi = 600,
        nBins = 10000,
        fontSize = 12,
        nonEmpty = nonEmpty
    log:
        'logs/plotCoverage/{all}-{region}-{pm}-{type}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/contactCoverage.py {params.nonEmpty} '
        '--out {output} --dpi {params.dpi} --nBins {params.nBins} '
        '--fontSize {params.fontSize} &> {log}'


rule OnTAD:
    input:
        rules.H5_to_NxN.output
    output:
        'dat/tads/{region}/{bin}/{all}-{region}-{bin}-ontad-{pm}.tad'
    params:
        chr = lambda wc: re.sub('chr', '', str(REGIONS['chr'][wc.region])),
        length = lambda wc: REGIONS['length'][wc.region],
        outprefix = 'dat/tads/{region}/{bin}/{all}-{region}-{bin}-ontad-{pm}'
    group:
        'processHiC'
    log:
        'logs/OnTAD/{all}-{region}-{bin}-{pm}.log'
    shell:
        '{SCRIPTS}/OnTAD {input} -o {params.outprefix} -bedout chr{params.chr} '
        '{params.length} {wildcards.bin} &> {log}'


rule reformatOnTAD:
    input:
        rules.OnTAD.output
    output:
        'dat/tads/{region}/{bin}/{all}-{region}-{bin}-ontad-{pm}.bed'
    params:
        chrom = lambda wc: REGIONS['chr'][wc.region],
        scale = lambda wc: REGIONS['start'][wc.region],
        maxPos = lambda wc: REGIONS['end'][wc.region]
    group:
        'processHiC'
    log:
        'logs/reformatOnTAD/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/reformatOnTAD.py {input} --scale {params.scale} '
        '--binSize {wildcards.bin} --chrom {params.chrom} '
        '--maxPos {params.maxPos} > {output} 2> {log}'


rule mergeOnTAD:
    input:
        expand('dat/tads/{region}/{{bin}}/{{all}}-{region}-{{bin}}-ontad-{{pm}}.bed',
            region=REGIONS.index)
    output:
        'referenceTADs/{all}-{bin}-ontad-{pm}.bed'
    log:
        'logs/mergeOnTAD/{all}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools sort -i <(cat {input}) > {output} 2> {log}'


rule distanceNormalise:
    input:
        'dat/matrix/{region}/{bin}/KR/{all}-{region}-{bin}-{pm}.h5'
    output:
        'dat/matrix/{region}/{bin}/KR/obs_exp/{all}-{region}-{bin}-{pm}.h5'
    params:
        method = 'obs_exp'
    group:
        'processHiC'
    log:
        'logs/distanceNormalise/{all}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicTransform -m {input} --method {params.method} -o {output} &> {log}'


def getTracks(wc):
    """ Build track command for generate config """
    command = ''
    if config['Genes']['gff3']:
        command += f'--genes Genes,{rules.processBED12.output},3 '
    if isinstance(config['bigWig'], dict):
        for title, track in config['bigWig'].items():
            command += f'--bigWig {title},{track},3 '
    if isinstance(config['bed'], dict):
        for title, track in config['bed'].items():
            command += f'--genes {title},{track},1.5 '
    if config['plotParams']['vLines']:
        command += f'--vLines {config["plotParams"]["vLines"]} '
    return command


def getMatrix(wc):
    """ Return either normal or obs_exp matrix """
    if config['plotParams']['distanceNorm']:
        return 'dat/matrix/{region}/{bin}/KR/obs_exp/{group}-{region}-{bin}-{pm}.h5'
    else:
        return 'dat/matrix/{region}/{bin}/KR/{group}-{region}-{bin}-{pm}.h5'

def getCscoreInput(wc):
    if config['HiCParams']['compartmentScore']:
        return f'dat/Cscore/{wc.region}/{wc.bin}/{wc.group}-{wc.region}-{wc.bin}-{wc.pm}-Cscore_cscore.bed'
    else:
        return []

def getCscoreParams(wc):
    if config['HiCParams']['compartmentScore']:
        cscore = f'dat/Cscore/{wc.region}/{wc.bin}/{wc.group}-{wc.region}-{wc.bin}-{wc.pm}-Cscore_cscore.bed'
        return f'--CScore {cscore}'
    else:
        return ''

def getGenesInput(wc):
    if config['Genes']['gff3']:
        return f'{rules.processBED12.output}'
    else:
        return []

def getDepth(wc):
    chrom, start, end = wc.coord.split('_')
    return int(end) - int(start)

def getLoops(wc):
    if config['loops']['detectLoops']:
        return f'dat/loops/{wc.region}/{wc.bin}/{wc.group}-{wc.region}-{wc.bin}-{wc.pm}.bedgraph'
    else:
        return []

def getLoopParams(wc):
    if config['loops']['plot'] and config['loops']['detectLoops']:
        return f'--loops {getLoops(wc)}'
    else:
        return ''

rule createConfig:
    input:
        matrix = getMatrix,
        loops = getLoops,
        insulations = 'dat/tads/{region}/{bin}/{group}-{region}-{bin}-{pm}_tad_score.bm',
        tads = 'dat/tads/{region}/{bin}/{group}-{region}-{bin}-ontad-{pm}.bed',
        cscore = getCscoreInput,
        vLines = config['plotParams']['vLines'],
        genes = getGenesInput
    output:
        'plots/{region}/{bin}/pyGenomeTracks/configs/{group}-{region}-{coord}-{bin}-{pm}.ini'
    params:
        tracks = getTracks,
        depth = getDepth,
        loops = getLoopParams,
        colourmap = config['plotParams']['colourmap'],
        geneLabelFontSize = config['plotParams']['geneLabelFontSize'],
        maxLabels = config['plotParams']['maxLabels'],
        vMin = '--vMin 0' if config['plotParams']['distanceNorm'] else '',
        vMax = '--vMax 2' if config['plotParams']['distanceNorm'] else '',
        log = '' if config['plotParams']['distanceNorm'] else '--log',
        cscore = getCscoreParams
    group:
        'processHiC'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/createConfig/{group}-{region}-{coord}-{bin}-{pm}.log'
    shell:
        'python {SCRIPTS}/generate_config.py --matrix {input.matrix} '
        '{params.log} --colourmap {params.colourmap} {params.tracks} '
        '--depth {params.depth} {params.vMin} {params.vMax} '
        '--geneLabelFontSize {params.geneLabelFontSize} '
        '--maxLabels {params.maxLabels} '
        '--insulation {input.insulations} {params.loops} '
        '{params.cscore} --tads {input.tads} | uniq > {output} 2> {log}'


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
    title = f'"{name} : {wc.region}{build} at {wc.bin} bin size (KR - {wc.pm})"',
    return title


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        'plots/{region}/{bin}/pyGenomeTracks/{group}-{region}-{coord}-{bin}-{pm}.{type}'
    params:
        region = setRegion,
        title = setMatrixTitle,
        dpi = 600
    group:
        'processHiC'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    log:
        'logs/plotHiC/{group}-{coord}-{region}-{bin}-{pm}-{type}.log'
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
        'dat/matrix/{region}/{bin}/KR/{group}-{region}-{bin}-{pm}.h5'
    output:
        bedgraph = 'dat/viewpoints/{region}/{bin}/KR/{group}-{region}-{coord}-{bin}-{pm}.bedgraph',
        plot = temp('dat/viewpoints/{region}/{bin}/KR/{group}-{region}-{coord}-{bin}-{pm}.png')
    params:
        referencePoint = setRegion,
        region = makeViewRegion,
    group:
        'processHiC'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    log:
        'logs/runViewpoint/{group}-{coord}-{region}-{bin}-{pm}.log'
    shell:
        '(hicPlotViewpoint --matrix {input} --region {params.region} '
        '--outFileName {output.plot} --referencePoint {params.referencePoint} '
        '--interactionOutFileName {output.bedgraph} '
        '&& mv {output.bedgraph}_*.bedgraph {output.bedgraph}) &> {log}'


rule plotViewpoint:
    input:
        rules.runViewpoint.output.bedgraph
    output:
        'plots/{region}/{bin}/viewpoints/{group}-{region}-{coord}-{bin}-viewpoint-{pm}.{type}'
    params:
        dpi = 600,
        build = f'--build {config["build"]}' if config['build'] else ''
    group:
        'processHiC'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/plotViewpoint/{group}-{coord}-{region}-{bin}-{pm}-{type}.log'
    shell:
        'python {SCRIPTS}/plotViewpoint.py {input} --out {output} '
        '--dpi {params.dpi} {params.build} &> {log}'


rule plotMatrix:
    input:
        rules.distanceNormalise.output
    output:
        'plots/{region}/{bin}/obs_exp/{all}-{region}-{bin}-{pm}.{type}'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        title = setMatrixTitle,
        dpi = 600,
        colour = 'YlGn'
    log:
        'logs/plotMatrix/{all}-{region}-{bin}-{pm}-{type}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicPlotMatrix --matrix {input} --outFileName {output} '
        '--region {params.chr}:{params.start}-{params.end} '
        '--colorMap {params.colour} --title {params.title} '
        '--vMin 0 --vMax 2 --dpi {params.dpi} &> {log}'


rule H5_to_NxN3p:
    input:
        'dat/matrix/{region}/{bin}/raw/{sample}-{region}-{bin}-{pm}-raw.h5'
    output:
        'dat/matrix/{region}/{bin}/raw/{sample}-{region}-{bin}-{pm}.nxn3p.tsv'
    group:
        'HiCRep'
    log:
        'logs/H5_to_NxN3p/{sample}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/H5_to_NxN3p.py {input} > {output} 2> {log}'


rule HiCRep:
    input:
        'dat/matrix/{region}/{bin}/raw/{sample1}-{region}-{bin}-{pm}.nxn3p.tsv',
        'dat/matrix/{region}/{bin}/raw/{sample2}-{region}-{bin}-{pm}.nxn3p.tsv',
    output:
        'qc/hicrep/data/{sample1}-vs-{sample2}-{region}-{bin}-hicrep-{pm}.csv'
    params:
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'HiCRep'
    log:
        'logs/HiCRep/{sample1}-vs-{sample2}-{region}-{bin}-{pm}.log'
    params:
        THREADS
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
        'qc/hicrep/{region}-{bin}-hicrep-{pm}.{type}'
    params:
        dpi = 600,
    group:
        'HiCRep'
    log:
        'logs/plotHiCRep/{region}-{bin}-{pm}-{type}.log'
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
        'logs/mergeBamByReplicate/{group}-{region}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        THREADS
    shell:
        'samtools merge -@ {threads} {output} {params.nonEmpty} 2> {log}'


rule H5_to_SUTM:
    input:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}-raw.h5'
    output:
        'dat/matrix/{region}/{bin}/raw/{all}-{region}-{bin}-{pm}-sutm.txt'
    group:
        'HiCcompare'
    log:
        'logs/H5_to_SUTM/{all}-{region}-{bin}-{pm}.log'
    conda:
         f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/H5_to_SUTM.py {input} > {output} 2> {log}'


rule HiCcompare:
    input:
        'dat/matrix/{region}/{bin}/raw/{group1}-{region}-{bin}-{pm}-sutm.txt',
        'dat/matrix/{region}/{bin}/raw/{group2}-{region}-{bin}-{pm}-sutm.txt'
    output:
        adjIF1 = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF1-{pm}.2d.txt',
        adjIF2 = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF2-{pm}.2d.txt',
        adjIF1sutm = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF1-{pm}.sutm',
        adjIF2sutm = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF2-{pm}.sutm'
    params:
        dir = lambda wc: f'dat/HiCcompare/{wc.region}/{wc.bin}',
        qcdir = lambda wc: f'qc/HiCcompare/{wc.region}/{wc.bin}',
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region],
        region = lambda wc: wc.region,
        suffix = lambda wc: wc.pm
    group:
        'HiCcompare'
    log:
        'logs/HiCcompare/{group1}-vs-{group2}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/HiCcompare.yaml'
    shell:
        'Rscript {SCRIPTS}/HiCcompare.R {params.dir} {params.qcdir} '
        '{params.chr} {params.start} {params.end} {params.region} '
        '{wildcards.bin} {params.suffix} {input} &> {log}'


rule Text2DToH5:
    input:
        chromSizes = getChromSizes,
        matrix = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{set}-{pm}.2d.txt'
    output:
        'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{set}-{pm}-allChrom.h5'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region]
    group:
        'HiCcompare'
    log:
        'logs/Text2DToH5/{group1}-vs-{group2}-{region}-{bin}-{set}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        '(hicConvertFormat --matrices {input.matrix} --outFileName {output} '
        '--inputFormat 2D-text --outputFormat h5 --resolutions {wildcards.bin} '
        '--chromosomeSizes <(grep -w {params.chr} {input.chromSizes})) &> {log}'


rule adjustCompareMatrix:
    input:
        rules.Text2DToH5.output
    output:
        'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{set}-{pm}.h5'
    params:
        chr = lambda wc: REGIONS['chr'][wc.region],
        start = lambda wc: REGIONS['start'][wc.region] + 1,
        end = lambda wc: REGIONS['end'][wc.region]
    group:
        'HiCcompare'
    log:
        'logs/adjustCompareMatrix/{group1}-vs-{group2}-{region}-{bin}-{set}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicAdjustMatrix --matrix {input} --outFileName {output} '
        '--regions <(echo -e \'{params.chr}\t{params.start}\t{params.end}\') '
        '&> {log}'


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
    return f'dat/HiCcompare/{{region}}/{{bin}}/{{group1}}-vs-{{group2}}-{{region}}-{{bin}}-{adj}-{{pm}}.h5'


def setDomains(wc):
    """ Set TADdomains as same as target matrix e.g. IF1 = group1 """
    if config['compareMatrices']['tads'] is not None:
        return f'dat/tads/{wc.region}/{wc.region}-referenceTADs.bed'
    elif wc.adjIF == 'adjIF1':
        group = wc.group1
    else:
        group = wc.group2
    return f'dat/tads/{{region}}/{{bin}}/{group}-{{region}}-{{bin}}-ontad-{{pm}}.bed'


rule processDiffTAD:
    input:
        allTADs = setDomains,
        matrix = 'dat/HiCsubtract/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-LOESSdiff-medianFilter-{pm}.h5',
        raw1 = 'dat/matrix/{region}/{bin}/raw/{group1}-{region}-{bin}-{pm}-raw.h5',
        raw2 = 'dat/matrix/{region}/{bin}/raw/{group2}-{region}-{bin}-{pm}-raw.h5'
    output:
        outPkl = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}.pkl',
        outDiff = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}-diffTAD.bed',
        out = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}-allTAD.bed'
    params:
        threshold = config['compareMatrices']['zThreshold'],
        name = 'ASTAD' if ALLELE_SPECIFIC else 'diffTAD'
    group:
        'HiCcompare'
    log:
        'logs/processDiffTAD/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/processDiffTAD.py {input.matrix} '
        '{input.raw1} {input.raw2} {input.allTADs} '
        '--threshold {params.threshold} --name {params.name} '
        '--outDiff {output.outDiff} --outPickle {output.outPkl} '
        '> {output.out} 2> {log}'


def getLoopsInput(wc):
    if config['compareMatrices']['loops'] is None:
        loops = ([
            'dat/loops/{region}/{bin}/{group1}-{region}-{bin}-{pm}.bedgraph',
            'dat/loops/{region}/{bin}/{group2}-{region}-{bin}-{pm}.bedgraph'
        ])
    else:
        loops = config['compareMatrices']['loops']
    return loops


rule scoreLoopDiff:
    input:
        loops = getLoopsInput,
        matrix = 'dat/HiCsubtract/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-LOESSdiff-medianFilter-{pm}.h5'
    output:
        interactOut = 'dat/loops/diff/{group1}-vs-{group2}-{region}-LOESSdiff-{bin}-{pm}.interact',
        linksUp = 'dat/loops/diff/{group1}-vs-{group2}-{region}-LOESSdiff-{bin}-{pm}-linksUp.links',
        linksDown = 'dat/loops/diff/{group1}-vs-{group2}-{region}-LOESSdiff-{bin}-{pm}-linksDown.links'
    params:
        nBins = 20,
        maxLineWidth = 3
    log:
        'logs/scoreLoopDiff/{group1}-vs-{group2}-LOESSdiff-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/scoreLoopDiff.py {input.loops} '
        '--maxLineWidth {params.maxLineWidth} --nBins {params.nBins} '
        '--matrix {input.matrix} --interactOut {output.interactOut} '
        '--linksUp {output.linksUp} --linksDown {output.linksDown} &> {log}'


rule runCompareViewpoint:
    input:
        'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{set}-{pm}.h5'
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
        'plots/{region}/{bin}/viewpoints/HiCcompare/{group1}-vs-{group2}-{region}-{coord}-{bin}-viewpoint-{pm}.{type}',
    params:
        dpi = 600,
        build = f'--build {config["build"]}' if config['build'] else ''
    group:
        'HiCcompare'
    conda:
        f'{ENVS}/python3.yaml'
    log:
        'logs/plotCompareViewpoint/{group1}-vs-{group2}-{region}-{coord}-{bin}-{pm}-{type}.log'
    shell:
        'python {SCRIPTS}/plotViewpoint.py {input} --out {output} '
        '--dpi {params.dpi} {params.build} &> {log}'


rule distanceNormaliseAdjIF:
    input:
        'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-{pm}.h5'
    output:
        'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-{adjIF}-obsExp-{pm}.h5'
    params:
        method = 'obs_exp'
    group:
        'HiCcompare'
    log:
        'logs/distanceNormaliseAdjIF/{group1}-vs-{group2}-{adjIF}-{region}-{bin}-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicTransform -m {input} --method {params.method} -o {output} &> {log}'


rule HiCsubtract:
    input:
        mat1 = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF1-obsExp-{pm}.h5',
        mat2 = 'dat/HiCcompare/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF2-obsExp-{pm}.h5',
        raw1 = 'dat/matrix/{region}/{bin}/raw/{group1}-{region}-{bin}-{pm}-raw.h5',
        raw2 = 'dat/matrix/{region}/{bin}/raw/{group2}-{region}-{bin}-{pm}-raw.h5'
    output:
        out = 'dat/HiCsubtract/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-LOESSdiff-noFilter-{pm}.h5',
        outFilt = 'dat/HiCsubtract/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-LOESSdiff-medianFilter-{pm}.h5'
    params:
        minSum = config['compareMatrices']['minSum']
    log:
        'logs/HiCsubtract/{group1}-vs-{group2}-{bin}-{region}-LOESSdiff-{pm}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'python {SCRIPTS}/compareHiC.py {input.mat1} {input.mat2} '
        '--outMatrix {output.out} '
        '--outMatrixFilter {output.outFilt} --minSum {params.minSum} '
        '--raw {input.raw1} {input.raw2} &> {log}'


def getSNPcoverage(wc):
    cellType1 = HiC.sample2Cell()[wc.group1]
    cellType2 = HiC.sample2Cell()[wc.group2]
    bin = wc.bin
    if ALLELE_SPECIFIC and (cellType1 == cellType2):
        return f'dat/genome/SNPcoverage/{cellType1}-{bin}-phasedHet.bed'
    return []

def getSNPcommand(wc):
    cellType1 = HiC.sample2Cell()[wc.group1]
    cellType2 = HiC.sample2Cell()[wc.group2]
    bin = wc.bin
    command = ''
    if ALLELE_SPECIFIC and (cellType1 == cellType2):
        track = f'dat/genome/SNPcoverage/{cellType1}-{bin}-phasedHet.bed'
        command += f'--SNPdensity {track} '
    return command

def getVmax(wc):
    vMax = config['compareMatrices']['vMax']
    # Reduce vMin slightly to account for median filter
    if wc.filter == 'medianFilter':
        vMax *= 0.75
    return vMax

def getVmin(wc):
    return -getVmax(wc)

def getCscoreInputSubtact(wc):
    if config['HiCParams']['compartmentScore']:
        return ([
            f'dat/Cscore/{wc.region}/{wc.bin}/{wc.group1}-{wc.region}-{wc.bin}-{wc.pm}-Cscore_cscore.bed',
            f'dat/Cscore/{wc.region}/{wc.bin}/{wc.group2}-{wc.region}-{wc.bin}-{wc.pm}-Cscore_cscore.bed'])
    else:
        return []

def getCscoreParams(wc):
    if config['HiCParams']['compartmentScore']:
        cscore = getCscoreInputSubtact(wc)
        return f'--CScore {cscore[0]} {cscore[1]}'
    else:
        return ''

def getLinksInput(wc):
    if ((config['compareMatrices']['loops'] is None)
            and (not config['loops']['detectLoops'])):
        return []
    else:
        return [f'dat/loops/diff/{wc.group1}-vs-{wc.group2}-{wc.region}-LOESSdiff-{wc.bin}-{wc.pm}-linksUp.links',
                f'dat/loops/diff/{wc.group1}-vs-{wc.group2}-{wc.region}-LOESSdiff-{wc.bin}-{wc.pm}-linksDown.links']

def getLinksParams(wc):
    links = getLinksInput(wc)
    if (links == []):
        return ''
    else:
        return f'--links {links[0]} {links[1]}'

rule createSubtractConfig:
    input:
        mat = 'dat/HiCsubtract/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-LOESSdiff-{filter}-{pm}.h5',
        tads = 'dat/tads/{region}/{bin}/{group1}-vs-{group2}-{region}-{bin}-adjIF1-{pm}-allTAD.bed',
        cscore = getCscoreInputSubtact,
        linksInput = getLinksInput,
        vLines = config['plotParams']['vLines'],
        SNPcoverage = getSNPcoverage,
        genes = getGenesInput
    output:
        ini = 'plots/{region}/{bin}/HiCsubtract/configs/{group1}-vs-{group2}-{coord}-LOESSdiff-{filter}-{pm}.ini',
        tmpLinks = temp('plots/{region}/{bin}/HiCsubtract/configs/{group1}-vs-{group2}-{coord}-LOESSdiff-{filter}-{pm}.tmp.links')
    params:
        vMin = getVmin,
        vMax = getVmax,
        depth = getDepth,
        tracks = getTracks,
        links = getLinksParams,
        cscore = getCscoreParams,
        SNPcoverage = getSNPcommand,
        maxLabels = config['plotParams']['maxLabels'],
        colourmap = config['compareMatrices']['colourmap'],
        geneLabelFontSize = config['plotParams']['geneLabelFontSize']
    group:
        'plotHiCsubtract'
    log:
        'logs/createSubtractConfig/{group1}-vs-{group2}-{bin}-{region}-{coord}-LOESSdiff-{filter}-{pm}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/generate_config.py --compare '
        '--matrix {input.mat} --vMin {params.vMin} --vMax {params.vMax} '
        '--tads {input.tads} {params.SNPcoverage} '
        '{params.links} --tmpLinks {output.tmpLinks} '
        '--geneLabelFontSize {params.geneLabelFontSize} '
        '--maxLabels {params.maxLabels} '
        '{params.cscore} '
        '--depth {params.depth} --colourmap {params.colourmap} '
        '{params.tracks} | uniq > {output.ini} 2> {log}'


def setSubtractTitle(wc):
    if config['build'] is not None:
        build = config['build'].replace('"', '') # Double quotes disallowed
        build = f' ({build})'
    else:
        build = ''
    title = (f'"{wc.group1} vs {wc.group2} - {wc.region}{build} at '
             f'{wc.bin} bin size - LOESSdiff - {wc.pm}"')
    return title


rule plotSubtract:
    input:
        ini = rules.createSubtractConfig.output.ini,
        tmpLinks = rules.createSubtractConfig.output.tmpLinks
    output:
        'plots/{region}/{bin}/HiCsubtract/{filter}/{group1}-vs-{group2}-{region}-{coord}-{bin}-LOESSdiff-{filter}-{pm}.{type}'
    params:
        title = setSubtractTitle,
        region = setRegion,
        dpi = 600
    group:
        'plotHiCsubtract'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    log:
        'logs/plotSubtract/{group1}-vs-{group2}-{bin}-{region}-{coord}-LOESSdiff-{filter}-{pm}-{type}.log'
    threads:
        THREADS
    shell:
        'export NUMEXPR_MAX_THREADS=1; pyGenomeTracks --tracks {input.ini} '
        '--region {params.region} --outFileName {output} '
        '--title {params.title} --dpi {params.dpi} &> {log}'


rule reformatPre:
    input:
        'dat/matrix/{region}/{all}-{region}.bam'
    output:
        'dat/matrix/{region}/base/raw/{all}-{region}.pre.tsv'
    group:
        'bam2hic'
    log:
        'logs/reformatPre/{all}-{region}.log'
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
        'logs/juicerPre/{all}-{region}.log'
    conda:
        f'{ENVS}/openjdk.yaml'
    shell:
        'java -jar {SCRIPTS}/juicer_tools_1.14.08.jar pre '
        '-c {params.chr} -r {params.resolutions} '
        '{input.tsv} {output} {input.chrom_sizes} &> {log}'


if not ALLELE_SPECIFIC:

    rule sortBam:
        input:
            'dat/mapped/{preSample}-dedup.bam'
        output:
            'dat/mapped/{preSample}-sort.bam'
        params:
            mem = '1G'
        log:
            'logs/sortBam/{preSample}.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            THREADS
        shell:
            'samtools sort -@ {threads} -m {params.mem} {input} '
            '> {output} 2> {log}'


    rule mergeBamByCellType:
        input:
            lambda wc: expand('dat/mapped/{preSample}-sort.bam',
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
            'samtools merge -u -@ {threads} - {input} > {output} 2> {log}'


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
            'logs/createSequenceDictionary/{cellType}.log'
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


    rule splitIntervals:
        input:
            genome = rules.bgzipGenome.output,
            index = rules.indexGenome.output,
            seqDict = rules.createSequenceDictionary.output
        output:
            expand('dat/gatk/splitIntervals/{{cellType}}/{rep}-scattered.interval_list',
                rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])])
        params:
            regions = config['regions'],
            scatterCount = config['gatk']['scatterCount']
        log:
            'logs/splitIntervals/{cellType}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
            'gatk SplitIntervals -R {input.genome} -L {params.regions} '
            '--scatter-count {params.scatterCount} '
            '--output dat/gatk/splitIntervals/{wildcards.cellType} &> {log} '


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
            ref_dict = rules.createSequenceDictionary.output,
            interval = 'dat/gatk/splitIntervals/{cellType}/{rep}-scattered.interval_list'
        output:
            'dat/gatk/baseRecalibrator/{cellType}-{rep}.recal.table'
        params:
            known = known_sites(config['gatk']['all_known']),
            tmp = config['tmpdir'],
            extra = ''
        log:
            'logs/baseRecalibrator/{cellType}-{rep}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk BaseRecalibrator {params.known} '
             '--input {input.bam} --reference {input.ref} '
             '--output {output} --intervals {input.interval} '
             '--sequence-dictionary {input.ref_dict} '
             '--tmp-dir {params.tmp} {params.extra} &> {log}'


    def gatherRecal(wc):
        splitRecal = expand(
            'dat/gatk/baseRecalibrator/{cellType}-{rep}.recal.table',
            rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])],
            cellType=wc.cellType
        )
        inputGather = ""
        for i in splitRecal:
            inputGather += f' --input {i}'
        return inputGather


    rule GatherBQSRReports:
        input:
            expand('dat/gatk/baseRecalibrator/{{cellType}}-{rep}.recal.table',
                rep=[str(i).zfill(4) for i in range(config['gatk']['scatterCount'])])
        output:
            'dat/gatk/baseRecalibrator/{cellType}.recal.table'
        params:
            input = gatherRecal
        log:
            'logs/GatherBQSRReports/{cellType}.log'
        conda:
            f'{ENVS}/gatk.yaml'
        shell:
             'gatk GatherBQSRReports --input {params.input} '
             '--output {output} &> {log}'


    rule applyBQSR:
        input:
            bam = rules.addReadGroup.output,
            bam_index = rules.indexMergedBam.output,
            ref = rules.bgzipGenome.output,
            ref_index = rules.indexGenome.output,
            ref_dict = rules.createSequenceDictionary.output,
            recal_table = rules.GatherBQSRReports.output,
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
            'logs/applyBQSR/{cellType}-{rep}.log'
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
            'logs/haplotypeCaller/{cellType}-{rep}.log'
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
            'logs/gatherGVCFs/{cellType}.log'
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
            'logs/sortGVCF/{cellType}.log'
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
            'logs/indexFeatureFile/{cellType}.log'
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
            'logs/genotypeGVCFs/{cellType}.log'
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
            'logs/selectVariants/{cellType}-{mode}.log'
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
            'logs/variantRecalibratorSNPs/{cellType}.log'
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
            maxGaussians = 4, # config['gatk']['maxGaussians'],
            java_opts = '-Xmx4G',
            tmp = config['tmpdir'],
            extra = '',  # optional
        group:
            'GATK'
        log:
            'logs/variantRecalibratorINDELS/{cellType}.log'
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
            'logs/applyVQSR/{cellType}-{mode}.log'
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
            'logs/mergeVCFs/{cellType}.log'
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
            'logs/mpileup/{cellType}-{region}.log'
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
            'logs/callVariants/{cellType}-{region}.log'
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
            'logs/filterVariants/{cellType}-{region}.log'
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
            'logs/extractHAIRS/{cellType}-{region}.log'
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
            'logs/hapCut2/{cellType}-{region}.log'
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
            'logs/bgzipPhased/{cellType}-{region}.log'
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
            'logs/indexPhased/{cellType}-{region}.log'
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
            'logs/extractBestPhase/{cellType}-{region}.log'
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
            'logs/extractVCF/{cellType}-{region}.log'
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
            'logs/bgzipVCF/{cellType}-{region}.log'
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
            'logs/indexVCF/{cellType}-{region}.log'
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
            'logs/bcftoolsStats/{cellType}-{region}.log'
        conda:
            f'{ENVS}/bcftools.yaml'
        shell:
            'bcftools stats {input} > {output} 2> {log}'


rule sampleReads:
    input:
        'dat/mapped/{preSample}.hic.bam'
    output:
        temp('dat/mapped/{preSample}.subsampled.sam')
    group:
        'filterQC'
    params:
        nLines = config['QC']['QCsample'] * 2
    log:
        'logs/sampleReads/{preSample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'cat <(samtools view -H {input}) '
        '<(samtools view {input} | head -n {params.nLines}) '
        '> {output} 2> {log}'


rule processHiC:
    input:
        rules.sampleReads.output,
    output:
        'dat/mapped/{preSample}-processed.pkl'
    group:
        'filterQC'
    log:
        'logs/process/{preSample}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/processHiC.py {wildcards.preSample} {output} '
        '{input} &> {log}'


rule plotQC:
    input:
        expand('dat/mapped/{sample}-processed.pkl',
            sample=HiC.originalSamples())
    output:
        'qc/filterQC/insertSizeFrequency.{type}'
    params:
        dpi = 600
    group:
        'filterQC'
    log:
        'logs/plotQC-{type}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/plotQC.py {output} {input} '
        '--dpi {params.dpi} &> {log}'


def multiQCconfig():
    if config['QC']['multiQCconfig']:
        return f'--config {config["QC"]["multiQCconfig"]}'
    else:
        return ''

rule multiqc:
    input:
        [expand('qc/fastqc/{sample}-{read}.{mode}_fastqc.zip',
            sample=HiC.originalSamples(),
            read=['R1', 'R2'], mode=['raw', 'trim']),
         expand('qc/cutadapt/{sample}.cutadapt.txt',
            sample=HiC.originalSamples()),
         expand('qc/fastqScreen/{sample}-{read}_screen.txt',
            sample=HiC.originalSamples(), read=['R1', 'R2']) if config['QC']['fastqScreen'] else [],
         expand('qc/bowtie2/{sample}-{read}.bowtie2.txt',
            sample=HiC.originalSamples(), read=['R1', 'R2']),
         expand('qc/hicup/HiCUP_summary_report-{sample}.txt',
            sample=HiC.originalSamples()) if not config['HiCParams']['microC'] else [],
         expand('qc/bcftools/{region}/{cellType}-{region}-bcftoolsStats.txt',
            region=REGIONS.index, cellType=HiC.cellTypes()) if PHASE_MODE=='BCFTOOLS' else []],
         expand('qc/deduplicate/{sample}.txt', sample=HiC.originalSamples()),
         expand('qc/samtools/{tool}/{sample}.{tool}.txt',
            sample=HiC.originalSamples(), tool=['stats']),
         expand('qc/hicexplorer/{sample}.{bin}-{pm}_QC',
            sample=HiC.samples(), bin=BASE_BIN, pm=pm),
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
