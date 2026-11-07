process SAMTOOLS_STATS {
    tag "${meta.id}"

    label 'process_low'

    publishDir "${params.outdir}/qc/samtools/stats", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
        path("${meta.id}.stats.txt"),
        emit: reports

    path "${meta.id}.samtools_stats.log", emit: logs

    script:
    """
    samtools stats \
        ${bam} \
        > ${meta.id}.stats.txt \
        2> ${meta.id}.samtools_stats.log
    """
}
