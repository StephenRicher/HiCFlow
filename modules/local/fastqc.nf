process FASTQC_RAW {
    tag "${meta.id}"

    label 'process_low'

    publishDir "${params.outdir}/qc/fastqc/raw", mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta),
        path("*_fastqc.html"),
        path("*_fastqc.zip"),
        emit: reports

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${read1} \
        ${read2}
    """
}
