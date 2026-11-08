process DISTANCE_NORMALISE {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/matrix/regions/${resolution}/KR/obs_exp",
        mode: 'copy',
        pattern: '*-obs-exp.h5'

    publishDir "${params.outdir}/logs/distance_normalise",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix)

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-full-obs-exp.h5"),
        emit: matrices

    path "${meta.id}-${region.name}-${resolution}.distance_normalise.log",
        emit: logs

    script:
    """
    hicTransform \
        --matrix ${matrix} \
        --method obs_exp \
        --outFileName ${meta.id}-${region.name}-${resolution}-full-obs-exp.h5 \
        > ${meta.id}-${region.name}-${resolution}.distance_normalise.log \
        2>&1
    """
}
