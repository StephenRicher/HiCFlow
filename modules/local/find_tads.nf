process FIND_TADS {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/tads/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*_boundaries.bed|*_boundaries.gff|*_domains.bed|*_score.bedgraph|*_zscore_matrix.h5'

    publishDir "${params.outdir}/logs/find_tads",
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
        path("${meta.id}-${region.name}-${resolution}-full_boundaries.bed"),
        path("${meta.id}-${region.name}-${resolution}-full_boundaries.gff"),
        path("${meta.id}-${region.name}-${resolution}-full_domains.bed"),
        path("${meta.id}-${region.name}-${resolution}-full_score.bedgraph"),
        path("${meta.id}-${region.name}-${resolution}-full_zscore_matrix.h5"),
        emit: results

    path "${meta.id}-${region.name}-${resolution}.find_tads.log",
        emit: logs

    script:
    def prefix = "${meta.id}-${region.name}-${resolution}-full"
    def minDepth = resolution * 3
    def maxDepth = resolution * 10

    """
    hicFindTADs \
        --matrix ${matrix} \
        --minDepth ${minDepth} \
        --maxDepth ${maxDepth} \
        --step ${resolution} \
        --outPrefix ${prefix} \
        --correctForMultipleTesting fdr \
        --numberOfProcessors ${task.cpus} \
        > ${meta.id}-${region.name}-${resolution}.find_tads.log \
        2>&1
    """
}
