process PLOT_OBS_EXP_MATRIX {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_low'

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/obs_exp",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/logs/plot_obs_exp_matrix",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix)

    val settings

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-full-obs-exp.${settings.filetype}"),
        emit: plots

    path "${meta.id}-${region.name}-${resolution}.plot_obs_exp_matrix.log",
        emit: logs

    script:
    def start = region.start + 1
    def title = "${meta.id} : ${region.name} (${settings.genomeBuild}) at ${resolution} bin size (KR - full)"

    """
    hicPlotMatrix \
        --matrix ${matrix} \
        --outFileName ${meta.id}-${region.name}-${resolution}-full-obs-exp.${settings.filetype} \
        --region '${region.chromosome}:${start}-${region.end}' \
        --colorMap YlGn \
        --title '${title}' \
        --vMin 0 \
        --vMax 2 \
        --dpi ${settings.dpi} \
        > ${meta.id}-${region.name}-${resolution}.plot_obs_exp_matrix.log \
        2>&1
    """
}
