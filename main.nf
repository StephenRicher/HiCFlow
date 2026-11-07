nextflow.enable.dsl = 2

include { FASTQC_RAW } from './modules/local/fastqc'

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

    log.info """
    HiCFlow Nextflow
    =================
    Sample sheet: ${params.samplesheet}
    Output directory: ${params.outdir}
    """

    samples = createSamplesChannel(params.samplesheet)

    FASTQC_RAW(samples)
}
