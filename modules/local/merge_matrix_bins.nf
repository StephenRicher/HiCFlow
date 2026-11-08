process MERGE_MATRIX_BINS {
    tag "${meta.id}:${region.name}:${targetResolution}"

    label 'process_low'

    publishDir "${params.outdir}/matrix/regions/${targetResolution}/raw",
        mode: 'copy',
        pattern: '*.h5'

    publishDir "${params.outdir}/logs/merge_matrix_bins",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(baseResolution),
        val(targetResolution),
        path(matrix)

    output:
    tuple val(meta),
        val(region),
        val(targetResolution),
        path("${meta.id}-${region.name}-${targetResolution}-full-raw.h5"),
        emit: matrices

    path "${meta.id}-${region.name}-${targetResolution}.merge_matrix_bins.log",
        emit: logs

    script:
    def numberOfBins = targetResolution.intdiv(baseResolution)

    if (targetResolution % baseResolution != 0) {
        error "Resolution ${targetResolution} is not divisible by base resolution ${baseResolution}"
    }

    """
    hicMergeMatrixBins \
        --matrix ${matrix} \
        --numBins ${numberOfBins} \
        --outFileName ${meta.id}-${region.name}-${targetResolution}-full-raw.h5 \
        > ${meta.id}-${region.name}-${targetResolution}.merge_matrix_bins.log \
        2>&1
    """
}
