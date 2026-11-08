process SPLIT_PAIRED_READS {
    tag "${meta.id}"

    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
        path("${meta.id}-R1.bam"),
        path("${meta.id}-R2.bam"),
        emit: reads

    path "${meta.id}.split_paired_reads.log", emit: logs

    script:
    """
    set -euo pipefail

    samtools view \
        -@ ${task.cpus} \
        -f 0x40 \
        -b \
        ${bam} \
        > ${meta.id}-R1.bam \
        2> ${meta.id}.split_paired_reads.log

    samtools view \
        -@ ${task.cpus} \
        -f 0x80 \
        -b \
        ${bam} \
        > ${meta.id}-R2.bam \
        2>> ${meta.id}.split_paired_reads.log
    """
}
