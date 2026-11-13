process PREPARE_BAM {
    tag "${meta.id}"

    label 'process_medium'

    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.hic.bam'
    publishDir "${params.outdir}/qc/deduplicate", mode: 'copy', pattern: '*.deduplicate.txt'
    publishDir "${params.outdir}/logs/prepare_bam", mode: 'copy', pattern: '*.prepare_bam.log'

    input:
    tuple val(meta), path(read1Sam), path(read2Sam)

    output:
    tuple val(meta),
        path("${meta.id}.hic.bam"),
        emit: bam

    tuple val(meta),
        path("${meta.id}.deduplicate.txt"),
        emit: reports

    path "${meta.id}.prepare_bam.log", emit: logs

    script:
    """
    set -euo pipefail

    exec 2> ${meta.id}.prepare_bam.log

    gawk \
        -f /opt/hicflow/workflow/scripts/addReadFlag.awk \
        -v flag=0x41 \
        ${read1Sam} \
        > ${meta.id}-R1-addFlag.sam

    samtools view \
        -u \
        ${meta.id}-R1-addFlag.sam \
        > ${meta.id}-R1-addFlag.bam

    rm -f ${meta.id}-R1-addFlag.sam

    gawk \
        -f /opt/hicflow/workflow/scripts/addReadFlag.awk \
        -v flag=0x81 \
        ${read2Sam} \
        > ${meta.id}-R2-addFlag.sam

    samtools view \
        -u \
        ${meta.id}-R2-addFlag.sam \
        > ${meta.id}-R2-addFlag.bam

    rm -f ${meta.id}-R2-addFlag.sam

    samtools merge \
        -@ ${task.cpus} \
        -n \
        ${meta.id}-merged.bam \
        ${meta.id}-R1-addFlag.bam \
        ${meta.id}-R2-addFlag.bam

    rm -f \
        ${meta.id}-R1-addFlag.bam \
        ${meta.id}-R2-addFlag.bam

    samtools fixmate \
        -p \
        -O SAM \
        ${meta.id}-merged.bam \
        ${meta.id}.fixed.sam

    rm -f ${meta.id}-merged.bam

    samblaster \
        --removeDups \
        --input ${meta.id}.fixed.sam \
        > ${meta.id}-dedup.sam \
        2> ${meta.id}.deduplicate.txt

    rm -f ${meta.id}.fixed.sam

    samtools view \
        -b \
        -@ ${task.cpus} \
        -F 2316 \
        ${meta.id}-dedup.sam \
        > ${meta.id}.hic.bam

    rm -f ${meta.id}-dedup.sam
    """
}
