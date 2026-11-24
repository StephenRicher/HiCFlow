nextflow.enable.dsl = 2

include { DOWNLOAD_FASTQ } from './modules/download_fastq'
include { DOWNLOAD_REFERENCE } from './modules/download_reference'
include { PREPARE_REFERENCE } from './modules/prepare_reference'
include { ALIGN_AND_FILTER } from './modules/align_and_filter'
include { MAKE_BIGWIG } from './modules/make_bigwig'
include { MERGE_BIGWIGS } from './modules/merge_bigwigs'
include { DOWNLOAD_CORIGAMI_RESOURCES } from './modules/download_corigami_resources'
include { CORIGAMI_PREDICT } from './modules/corigami_predict'
include { CORIGAMI_EDIT } from './modules/corigami_edit'
include { COMPARE_CORIGAMI_PREDICTIONS } from './modules/compare_corigami_predictions'

workflow {
    if (!params.samples) {
        error 'Missing required parameter: --samples'
    }

    if (!params.outdir) {
        error 'Missing required parameter: --outdir'
    }

    samples = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def missing = ['gsm', 'cell_line', 'assay', 'replicate']
                .findAll { column -> !row[column] }

            if (missing) {
                error "Missing sample values: ${missing.join(', ')}"
            }

            def meta = [
                id       : "${row.cell_line}_${row.assay}_rep${row.replicate}",
                cell_line: row.cell_line,
                assay    : row.assay,
                replicate: row.replicate as Integer
            ]

            tuple(meta, row.gsm)
        }

    DOWNLOAD_FASTQ(samples)

    DOWNLOAD_REFERENCE(
        Channel.value(params.referenceUrl)
    )

    PREPARE_REFERENCE(
        DOWNLOAD_REFERENCE.out.reference
    )

    alignmentInputs = DOWNLOAD_FASTQ.out.reads
        .combine(PREPARE_REFERENCE.out.reference)
        .map {
            meta,
            read1,
            read2,
            reference,
            index ->

            tuple(
                meta,
                read1,
                read2,
                reference,
                index
            )
        }

    ALIGN_AND_FILTER(alignmentInputs)

    MAKE_BIGWIG(
        ALIGN_AND_FILTER.out.bam
    )

    mergedBigWigInputs = MAKE_BIGWIG.out.bigwigs
        .map { meta, bigWig ->
            tuple(
                meta.cell_line,
                meta.assay,
                bigWig
            )
        }
        .groupTuple(by: [0, 1])

    MERGE_BIGWIGS(mergedBigWigInputs)

    if (params.runCorigami) {
        if (!params.deletions) {
            error 'Missing required parameter: --deletions'
        }

        DOWNLOAD_CORIGAMI_RESOURCES(
            Channel.value(params.corigamiResourcesUrl)
        )

        finalSignals = MERGE_BIGWIGS.out.bigwigs
            .map { cellLine, assay, bigWig ->
                tuple(
                    "${cellLine}|${assay}",
                    cellLine,
                    assay,
                    bigWig
                )
            }

        atacSignals = finalSignals
            .filter { key, cellLine, assay, bigWig ->
                assay == 'ATAC'
            }
            .map { key, cellLine, assay, bigWig ->
                tuple(cellLine, bigWig)
            }

        ctcfSignals = finalSignals
            .filter { key, cellLine, assay, bigWig ->
                assay == 'CTCF'
            }
            .map { key, cellLine, assay, bigWig ->
                tuple(cellLine, bigWig)
            }

        cellSignals = ctcfSignals
            .join(atacSignals)
            .map { cellLine, ctcf, atac ->
                tuple(cellLine, ctcf, atac)
            }

        deletions = Channel
            .fromPath(params.deletions, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def deletionStart = row.start as Integer
                def deletionEnd = row.end as Integer

                tuple(
                    row.cell_line,
                    row.clone,
                    row.allele as Integer,
                    row.chromosome,
                    params.predictionWindowStart as Integer,
                    deletionStart,
                    deletionEnd
                )
            }

        deletionsWithSignals = deletions
            .combine(cellSignals, by: 0)
            .map {
                cellLine,
                clone,
                allele,
                chromosome,
                windowStart,
                deletionStart,
                deletionEnd,
                ctcf,
                atac ->

                tuple(
                    cellLine,
                    clone,
                    allele,
                    chromosome,
                    windowStart,
                    deletionStart,
                    deletionEnd,
                    ctcf,
                    atac
                )
            }

        deletionInputs = deletionsWithSignals
            .combine(DOWNLOAD_CORIGAMI_RESOURCES.out.resources)
            .map {
                cellLine,
                clone,
                allele,
                chromosome,
                windowStart,
                deletionStart,
                deletionEnd,
                ctcf,
                atac,
                resources ->

                tuple(
                    cellLine,
                    clone,
                    allele,
                    chromosome,
                    windowStart,
                    deletionStart,
                    deletionEnd,
                    ctcf,
                    atac,
                    resources
                )
            }

        CORIGAMI_EDIT(
            deletionInputs,
            params.padding
        )

        wildTypeInputs = deletionInputs
            .map {
                cellLine,
                clone,
                allele,
                chromosome,
                windowStart,
                deletionStart,
                deletionEnd,
                ctcf,
                atac,
                resources ->

                tuple(
                    "${cellLine}|${chromosome}|${windowStart}",
                    cellLine,
                    chromosome,
                    windowStart,
                    ctcf,
                    atac,
                    resources
                )
            }
            .unique { key, cellLine, chromosome, windowStart, ctcf, atac, resources ->
                key
            }
            .map { key, cellLine, chromosome, windowStart, ctcf, atac, resources ->
                tuple(
                    cellLine,
                    chromosome,
                    windowStart,
                    ctcf,
                    atac,
                    resources
                )
            }

        CORIGAMI_PREDICT(wildTypeInputs)

        wildTypePredictions = CORIGAMI_PREDICT.out.predictions
            .map {
                cellLine,
                chromosome,
                windowStart,
                prediction ->

                tuple(
                    "${cellLine}|${chromosome}|${windowStart}",
                    prediction
                )
            }

        comparisonInputs = CORIGAMI_EDIT.out.predictions
            .map {
                cellLine,
                clone,
                allele,
                chromosome,
                windowStart,
                deletionStart,
                deletionEnd,
                prediction ->

                tuple(
                    "${cellLine}|${chromosome}|${windowStart}",
                    cellLine,
                    clone,
                    allele,
                    chromosome,
                    windowStart,
                    deletionStart,
                    deletionEnd,
                    prediction
                )
            }
            .combine(wildTypePredictions, by: 0)
            .map {
                key,
                cellLine,
                clone,
                allele,
                chromosome,
                windowStart,
                deletionStart,
                deletionEnd,
                deletionPrediction,
                wildTypePrediction ->

                tuple(
                    cellLine,
                    clone,
                    allele,
                    chromosome,
                    windowStart,
                    deletionStart,
                    deletionEnd,
                    wildTypePrediction,
                    deletionPrediction
                )
            }

        comparisonSettings = [
            windowSize          : 2_097_152,
            predictionPercentile: 99.5,
            differencePercentile: 99.0,
            differenceLimit     : params.corigamiDifferenceLimit,
            dpi                 : params.corigamiPlotDpi
        ]

        COMPARE_CORIGAMI_PREDICTIONS(
            comparisonInputs,
            comparisonSettings
        )
    }
}
