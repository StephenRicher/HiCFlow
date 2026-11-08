process ADJUST_MATRIX {
    tag "${meta.id}:${region.name}"

    label 'process_low'

    publishDir "${params.outdir}/matrix/regions/base/raw", mode: 'copy', pattern: '*.h5'
    publishDir "${params.outdir}/logs/adjust_matrix", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(matrix), val(region)

    output:
    tuple val(meta),
        val(region),
        path("${meta.id}-${region.name}.${params.resolution.base}-full.h5"),
        emit: matrices

    path "${meta.id}-${region.name}.adjust_matrix.log", emit: logs

    script:
    """
    printf '%s\\t%s\\t%s\\n' \
        '${region.chromosome}' \
        '${region.start + 1}' \
        '${region.end}' \
        > region.bed

    hicAdjustMatrix \
        --matrix ${matrix} \
        --outFileName ${meta.id}-${region.name}.${params.resolution.base}-full.h5 \
        --regions region.bed \
        > ${meta.id}-${region.name}.adjust_matrix.log \
        2>&1
    """
}
