process CUTADAPT {
    tag "${meta.id}"

    label 'process_low'

    publishDir "${params.outdir}/fastq/trimmed", mode: 'copy', pattern: '*.fastq.gz'
    publishDir "${params.outdir}/qc/cutadapt", mode: 'copy', pattern: '*.cutadapt.txt'

    input:
    tuple val(meta), path(read1), path(read2)
    val settings

    output:
    tuple val(meta),
        path("${meta.id}-R1.trimmed.fastq.gz"),
        path("${meta.id}-R2.trimmed.fastq.gz"),
        emit: reads

    tuple val(meta),
        path("${meta.id}.cutadapt.txt"),
        emit: reports

    script:
    """
    cutadapt \
        --cores ${task.cpus} \
        -a '${settings.forwardAdapter}' \
        -A '${settings.reverseAdapter}' \
        --overlap ${settings.overlap} \
        --error-rate ${settings.errorRate} \
        --minimum-length ${settings.minimumLength} \
        --quality-cutoff ${settings.qualityCutoff} \
        --output ${meta.id}-R1.trimmed.fastq.gz\
        --paired-output ${meta.id}-R2.trimmed.fastq.gz \
        ${read1} \
        ${read2} \
        > ${meta.id}.cutadapt.txt
    """
}
