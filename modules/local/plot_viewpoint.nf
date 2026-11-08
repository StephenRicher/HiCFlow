process PLOT_VIEWPOINT {
    tag "${meta.id}:${viewpoint.id}:${resolution}"

    label 'process_low'

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/viewpoints",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/viewpoints/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*.bedgraph'

    publishDir "${params.outdir}/logs/plot_viewpoint",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix),
        val(viewpoint)

    val settings

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        val(viewpoint),
        path("${meta.id}-${region.name}-${viewpoint.id}-${resolution}-viewpoint.${settings.filetype}"),
        emit: plots

    tuple val(meta),
        val(region),
        val(resolution),
        val(viewpoint),
        path("${meta.id}-${region.name}-${viewpoint.id}-${resolution}-viewpoint.bedgraph"),
        emit: bedgraphs

    path "${meta.id}-${region.name}-${viewpoint.id}-${resolution}.plot_viewpoint.log",
        emit: logs

    script:
    def prefix = "${meta.id}-${region.name}-${viewpoint.id}-${resolution}"
    def viewStart = Math.max(viewpoint.start - settings.viewpointRange, region.start + 1)
    def viewEnd = Math.min(viewpoint.end + settings.viewpointRange, region.end)
    def referencePoint = "${viewpoint.chromosome}:${viewpoint.start + 1}-${viewpoint.end}"
    def viewRegion = "${viewpoint.chromosome}:${viewStart}-${viewEnd}"
    def buildArgument = settings.genomeBuild ? "--build '${settings.genomeBuild}'" : ''

    """
    set -euo pipefail

    hicPlotViewpoint \
        --matrix ${matrix} \
        --region '${viewRegion}' \
        --referencePoint '${referencePoint}' \
        --outFileName ${prefix}.temporary.png \
        --interactionOutFileName ${prefix}-viewpoint \
        > ${prefix}.plot_viewpoint.log \
        2>&1

    generated_bedgraph=\$(find . \
        -maxdepth 1 \
        -name '${prefix}-viewpoint*.bedgraph' \
        -type f \
        | head -n 1)

    if [ -z "\${generated_bedgraph}" ]; then
        echo "hicPlotViewpoint did not produce a bedGraph" >&2
        exit 1
    fi

    mv \
        "\${generated_bedgraph}" \
        ${prefix}-viewpoint.bedgraph

    python /opt/hicflow/workflow/scripts/plotViewpoint.py \
        ${prefix}-viewpoint.bedgraph \
        --viewpoint '${referencePoint}' \
        --out ${prefix}-viewpoint.${settings.filetype} \
        --dpi ${settings.dpi} \
        ${buildArgument} \
        >> ${prefix}.plot_viewpoint.log \
        2>&1
    """
}
