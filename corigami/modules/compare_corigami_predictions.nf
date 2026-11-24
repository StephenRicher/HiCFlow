process COMPARE_CORIGAMI_PREDICTIONS {
    tag "${cellLine}:${clone}:allele${allele}"

    cpus 1
    memory '2 GB'
    time '2h'

    publishDir "${params.outdir}/corigami/comparisons/${cellLine}/${clone}/allele${allele}",
        mode: 'copy'

    input:
    tuple val(cellLine),
        val(clone),
        val(allele),
        val(chromosome),
        val(windowStart),
        val(deletionStart),
        val(deletionEnd),
        path(wildTypePrediction),
        path(deletionPrediction)

    val settings

    output:
    tuple val(cellLine),
        val(clone),
        val(allele),
        val(chromosome),
        val(windowStart),
        path("${cellLine}-${clone}-allele${allele}-wild-type.npy"),
        path("${cellLine}-${clone}-allele${allele}-deletion.npy"),
        path("${cellLine}-${clone}-allele${allele}-deletion-minus-WT.npy"),
        path("${cellLine}-${clone}-allele${allele}-wild-type.svg"),
        path("${cellLine}-${clone}-allele${allele}-deletion.svg"),
        path("${cellLine}-${clone}-allele${allele}-deletion-minus-WT.svg"),
        path("${cellLine}-${clone}-allele${allele}-viewpoint.tsv"),
        path("${cellLine}-${clone}-allele${allele}-viewpoint.svg"),
        path("${cellLine}-${clone}-allele${allele}-metadata.json"),
        emit: comparisons

    script:
    def prefix = "${cellLine}-${clone}-allele${allele}"

    def differenceLimitArgument = settings.differenceLimit != null
        ? "--difference-limit ${settings.differenceLimit}"
        : ''

    def deletionLabel = clone in ['004', '005']
        ? 'DIRAS3 deletion'
        : 'B1 deletion'
    """
    set -euo pipefail

    /opt/conda/envs/corigami/bin/python \
        /workflow/scripts/compare_corigami_predictions.py \
        --wild-type ${wildTypePrediction} \
        --deletion ${deletionPrediction} \
        --cell-line '${cellLine}' \
        --clone '${clone}' \
        --allele '${allele}' \
        --chromosome '${chromosome}' \
        --window-start ${windowStart} \
        --window-size ${settings.windowSize} \
        --deletion-start ${deletionStart} \
        --deletion-end ${deletionEnd} \
        --output-prefix ${prefix} \
        --deletion-label '${deletionLabel}' \
        --prediction-percentile ${settings.predictionPercentile} \
        --difference-percentile ${settings.differencePercentile} \
        --dpi ${settings.dpi} \
        ${differenceLimitArgument}
    """
}
