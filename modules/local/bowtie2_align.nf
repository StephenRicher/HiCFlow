process BOWTIE2_ALIGN {
    tag "${meta.id}"

    label 'process_medium'

    publishDir "${params.outdir}/qc/bowtie2", mode: 'copy', pattern: '*.bowtie2.txt'
    publishDir "${params.outdir}/logs/bowtie2", mode: 'copy', pattern: '*.bowtie2.log'

    input:
    tuple val(cellType), val(meta), path(read1), path(read2), path(indexFiles)

    output:
    tuple val(meta),
        path("${meta.id}-R1.sam"),
        path("${meta.id}-R2.sam"),
        emit: alignments

    tuple val(meta),
        path("${meta.id}-R1.bowtie2.txt"),
        path("${meta.id}-R2.bowtie2.txt"),
        emit: reports

    path "${meta.id}-R1.bowtie2.log", emit: logsR1
    path "${meta.id}-R2.bowtie2.log", emit: logsR2

    script:
    """
    bowtie2 \
        -x ${cellType} \
        -U ${read1} \
        --reorder \
        --threads ${task.cpus} \
        --sensitive \
        > ${meta.id}-R1.sam \
        2> ${meta.id}-R1.bowtie2.log

    cp \
        ${meta.id}-R1.bowtie2.log \
        ${meta.id}-R1.bowtie2.txt

    bowtie2 \
        -x ${cellType} \
        -U ${read2} \
        --reorder \
        --threads ${task.cpus} \
        --sensitive \
        > ${meta.id}-R2.sam \
        2> ${meta.id}-R2.bowtie2.log

    cp \
        ${meta.id}-R2.bowtie2.log \
        ${meta.id}-R2.bowtie2.txt
    """
}
