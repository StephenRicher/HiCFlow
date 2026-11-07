nextflow.enable.dsl = 2

include { FASTQC_RAW } from './modules/local/fastqc'
include { CUTADAPT } from './modules/local/cutadapt'
include { FASTQC_TRIMMED } from './modules/local/fastqc_trimmed'
include { HICUP_TRUNCATE } from './modules/local/hicup_truncate'
include { BOWTIE2_BUILD } from './modules/local/bowtie2_build'
include { BOWTIE2_ALIGN } from './modules/local/bowtie2_align'
include { PREPARE_BAM } from './modules/local/prepare_bam'
include { SAMTOOLS_STATS } from './modules/local/samtools_stats'

def createSamplesChannel(samplesheetPath) {
    def samplesheet = file(samplesheetPath, checkIfExists: true)
    def samplesheetDir = samplesheet.parent

    Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def missingValues = ['experiment', 'cell_type', 'group', 'rep', 'R1', 'R2']
                .findAll { column -> !row[column] }

            if (missingValues) {
                error "Missing required values in ${samplesheet}: ${missingValues.join(', ')}"
            }

            def sampleId = "${row.group}-${row.rep}"
            def read1 = samplesheetDir.resolve(row.R1).normalize()
            def read2 = samplesheetDir.resolve(row.R2).normalize()

            if (!read1.exists()) {
                error "R1 FASTQ does not exist for ${sampleId}: ${read1}"
            }

            if (!read2.exists()) {
                error "R2 FASTQ does not exist for ${sampleId}: ${read2}"
            }

            def meta = [
                id         : sampleId,
                experiment : row.experiment,
                cell_type  : row.cell_type,
                group      : row.group,
                replicate  : row.rep as Integer
            ]

            tuple(meta, read1, read2)
        }
}

workflow {
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }

    if (!params.genome) {
        error 'Missing required genome configuration'
    }

    log.info """
    HiCFlow Nextflow
    =================
    Sample sheet: ${params.samplesheet}
    Output directory: ${params.outdir}
    """

    samples = createSamplesChannel(params.samplesheet)

    genomes = Channel.fromList(
        params.genome.collect { cellType, referencePath ->
            tuple(
                cellType,
                file(referencePath, checkIfExists: true)
            )
        }
    )

    BOWTIE2_BUILD(genomes)

    FASTQC_RAW(samples)

    CUTADAPT(
        samples,
        params.cutadapt
    )

    FASTQC_TRIMMED(CUTADAPT.out.reads)

    HICUP_TRUNCATE(
        CUTADAPT.out.reads,
        params.hicup
    )

    readsByCellType = HICUP_TRUNCATE.out.reads.map { meta, read1, read2 ->
        tuple(meta.cell_type, meta, read1, read2)
    }

    alignmentInputs = readsByCellType
        .combine(BOWTIE2_BUILD.out.index, by: 0)
        .map { cellType, meta, read1, read2, indexFiles ->
            tuple(cellType, meta, read1, read2, indexFiles)
        }

    BOWTIE2_ALIGN(alignmentInputs)

    PREPARE_BAM(BOWTIE2_ALIGN.out.alignments)

    SAMTOOLS_STATS(PREPARE_BAM.out.bam)
}
