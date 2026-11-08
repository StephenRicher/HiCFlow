process PLOT_SUBTRACTION {
    tag "${group2.id}-minus-${group1.id}:${coordinate.id}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/subtraction/medianFilter",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/subtraction/configs",
        mode: 'copy',
        pattern: '*.ini'

    publishDir "${params.outdir}/logs/plot_subtraction",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path(subtractionMatrix),
        path(tads),
        val(coordinate)

    path genes
    path bigWigs
    val settings

    output:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        val(coordinate),
        path("${group2.id}-minus-${group1.id}-${region.name}-${coordinate.id}-${resolution}-medianFilter.${settings.filetype}"),
        emit: plots

    path "${group2.id}-minus-${group1.id}-${region.name}-${coordinate.id}-${resolution}-medianFilter.ini",
        emit: configs

    path "${group2.id}-minus-${group1.id}-${region.name}-${coordinate.id}-${resolution}.plot_subtraction.log",
        emit: logs

    script:
    def prefix = "${group2.id}-minus-${group1.id}-${region.name}-${coordinate.id}-${resolution}"
    def depth = coordinate.end - coordinate.start
    def vMax = settings.vMax * 0.75
    def vMin = -vMax

    def bigWigArguments = settings.bigWigTitles
        .withIndex()
        .collect { bigWigTitle, index ->
            "--bigWig '${bigWigTitle},${bigWigs[index]},3'"
        }
        .join(' ')

    def vLinesArgument = settings.vlines
        ? "--vLines '${settings.vlines}'"
        : ''

    def title = "${group2.id} - ${group1.id} : ${region.name} (${settings.genomeBuild}) at ${resolution} bp; " +
        "positive = higher in ${group2.id}, negative = higher in ${group1.id}"

    """
    set -euo pipefail

    export NUMEXPR_MAX_THREADS=1

    python /opt/hicflow/workflow/scripts/generate_config.py \
        --compare \
        --matrix ${subtractionMatrix} \
        --vMin ${vMin} \
        --vMax ${vMax} \
        --tads ${tads} \
        --depth ${depth} \
        --colourmap '${settings.colourmap}' \
        --geneLabelFontSize ${settings.geneLabelFontSize} \
        --maxLabels ${settings.maxLabels} \
        --genes 'Genes,${genes},3' \
        ${bigWigArguments} \
        ${vLinesArgument} \
        > ${prefix}-medianFilter.ini \
        2> ${prefix}.plot_subtraction.log

    pyGenomeTracks \
        --tracks ${prefix}-medianFilter.ini \
        --region '${coordinate.chromosome}:${coordinate.start + 1}-${coordinate.end}' \
        --outFileName ${prefix}-medianFilter.${settings.filetype} \
        --title '${title}' \
        --dpi ${settings.dpi} \
        >> ${prefix}.plot_subtraction.log \
        2>&1
    """
}
