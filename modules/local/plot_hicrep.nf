process PLOT_HICREP {
    tag "${region.name}:${resolution}"

    label 'process_low'

    publishDir "${params.outdir}/qc/hicrep",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/logs/plot_hicrep",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(region),
        val(resolution),
        path(scoreFiles)

    val settings

    output:
    tuple val(region),
        val(resolution),
        path("${region.name}-${resolution}-hicrep-full.${settings.filetype}"),
        emit: plots

    path "${region.name}-${resolution}.plot_hicrep.log",
        emit: logs

    script:
    def scoreArguments = scoreFiles.join(' ')

    """
    python /opt/hicflow/workflow/scripts/plotHiCRep.py \
        --out ${region.name}-${resolution}-hicrep-full.${settings.filetype} \
        --dpi ${settings.dpi} \
        ${scoreArguments} \
        > ${region.name}-${resolution}.plot_hicrep.log \
        2>&1
    """
}
