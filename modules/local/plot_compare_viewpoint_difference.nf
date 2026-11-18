process PLOT_COMPARE_VIEWPOINT_DIFFERENCE {
    tag "${comparison.group2}-minus-${comparison.group1}:${viewpoint.id}:${resolution}"

    label 'process_low'

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/viewpoints/comparisons",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/viewpoints/${region.name}/${resolution}/comparisons",
        mode: 'copy',
        pattern: '*.tsv'

    publishDir "${params.outdir}/logs/plot_compare_viewpoint_difference",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(comparison),
        val(region),
        val(resolution),
        val(viewpoint),
        path(group1Bedgraph),
        path(group2Bedgraph)

    val settings

    output:
    tuple val(comparison),
        val(region),
        val(resolution),
        val(viewpoint),
        path("${comparison.group2}-minus-${comparison.group1}-${region.name}-${viewpoint.id}-${resolution}-viewpoint.${settings.filetype}"),
        emit: plots

    tuple val(comparison),
        val(region),
        val(resolution),
        val(viewpoint),
        path("${comparison.group2}-minus-${comparison.group1}-${region.name}-${viewpoint.id}-${resolution}-viewpoint.tsv"),
        emit: tables

    path "${comparison.group2}-minus-${comparison.group1}-${region.name}-${viewpoint.id}-${resolution}.plot_compare_viewpoint_difference.log",
        emit: logs

    script:
    def prefix = "${comparison.group2}-minus-${comparison.group1}-${region.name}-${viewpoint.id}-${resolution}"
    def buildArgument = settings.genomeBuild
        ? "--build '${settings.genomeBuild}'"
        : ''

    def title = "${comparison.group2} - ${comparison.group1}; " +
        "positive = higher in ${comparison.group2}, " +
        "negative = higher in ${comparison.group1}"

    def referencePoint = "${viewpoint.chromosome}:${viewpoint.start + 1}-${viewpoint.end}"
    """
    set -euo pipefail

    python /opt/hicflow/workflow/scripts/plotCompareViewpointDifference.py \
        ${group1Bedgraph} \
        ${group2Bedgraph} \
        --group1 '${comparison.group1}' \
        --group2 '${comparison.group2}' \
        --viewpoint '${referencePoint}' \
        --table ${prefix}-viewpoint.tsv \
        --out ${prefix}-viewpoint.${settings.filetype} \
        --title '${title}' \
        --dpi ${settings.dpi} \
        ${buildArgument} \
        > ${prefix}.plot_compare_viewpoint_difference.log \
        2>&1
    """
}
