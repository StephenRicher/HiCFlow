process SUM_REPLICATES {
    tag "${groupMeta.id}:${region.name}"

    label 'process_low'

    publishDir "${params.outdir}/matrix/regions/base/raw", mode: 'copy', pattern: '*.h5'
    publishDir "${params.outdir}/logs/sum_replicates", mode: 'copy', pattern: '*.log'

    input:
    tuple val(groupMeta), val(region), path(matrices)

    output:
    tuple val(groupMeta),
        val(region),
        path("${groupMeta.id}-${region.name}.${params.resolution.base}-full.h5"),
        emit: matrices

    path "${groupMeta.id}-${region.name}.sum_replicates.log", emit: logs

    script:
    def matrixArguments = matrices.join(' ')

    """
    set -euo pipefail

    hicSumMatrices \
        --matrices ${matrixArguments} \
        --outFileName ${groupMeta.id}-${region.name}.${params.resolution.base}-full.h5 \
        > ${groupMeta.id}-${region.name}.sum_replicates.log \
        2>&1
    """
}
