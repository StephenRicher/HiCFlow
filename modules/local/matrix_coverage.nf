process MATRIX_COVERAGE {
    tag "${meta.id}:${region.name}"

    label 'process_low'

    publishDir "${params.outdir}/qc/matrix_coverage/${region.name}",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/logs/matrix_coverage",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolutions),
        path(matrices)

    val settings

    output:
    tuple val(meta),
        val(region),
        path("${meta.id}-${region.name}-coverage.${settings.filetype}"),
        emit: plots

    path "${meta.id}-${region.name}.matrix_coverage.log",
        emit: logs

    script:
    def matrixArguments = matrices.join(' ')
    def resolutionArguments = resolutions.join(' ')

    """
    python /opt/hicflow/workflow/scripts/contactCoverage.py \
        ${matrixArguments} \
        --resolutions ${resolutionArguments} \
        --sample '${meta.id}' \
        --out ${meta.id}-${region.name}-coverage.${settings.filetype} \
        --dpi ${settings.dpi} \
        --nBins ${settings.nBins} \
        --fontSize ${settings.fontSize} \
        > ${meta.id}-${region.name}.matrix_coverage.log \
        2>&1
    """
}
