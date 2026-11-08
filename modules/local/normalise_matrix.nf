process NORMALISE_MATRIX {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/matrix/regions/${resolution}/raw",
        mode: 'copy',
        pattern: '*-normalised.h5'

    publishDir "${params.outdir}/matrix/regions/${resolution}/KR",
        mode: 'copy',
        pattern: '*-KR.h5'

    publishDir "${params.outdir}/logs/normalise_matrix",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix)

    val multiplicativeValue

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-full-normalised.h5"),
        emit: normalised

    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-full-KR.h5"),
        emit: corrected

    path "${meta.id}-${region.name}-${resolution}.normalise_matrix.log",
        emit: logs

    script:
    """
    set -euo pipefail

    exec > ${meta.id}-${region.name}-${resolution}.normalise_matrix.log 2>&1

    hicNormalize \
        --matrices ${matrix} \
        --normalize norm_range \
        --outFileName ${meta.id}-${region.name}-${resolution}-norm01.h5

    hicNormalize \
        --matrices ${meta.id}-${region.name}-${resolution}-norm01.h5 \
        --normalize multiplicative \
        --multiplicativeValue ${multiplicativeValue} \
        --outFileName ${meta.id}-${region.name}-${resolution}-full-normalised.h5

    hicCorrectMatrix correct \
        --matrix ${meta.id}-${region.name}-${resolution}-full-normalised.h5 \
        --correctionMethod KR \
        --outFileName ${meta.id}-${region.name}-${resolution}-full-KR.h5
    """
}
