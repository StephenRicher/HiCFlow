process HICUP_TRUNCATE {
    tag "${meta.id}"

    label 'process_low'

    publishDir "${params.outdir}/qc/hicup", mode: 'copy', pattern: 'HiCUP_summary_report-*.txt'
    publishDir "${params.outdir}/logs/hicup_truncate", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(read1), path(read2)
    val settings

    output:
    tuple val(meta),
        path("${meta.id}-R1.truncated.fastq.gz"),
        path("${meta.id}-R2.truncated.fastq.gz"),
        emit: reads

    tuple val(meta),
        path("HiCUP_summary_report-${meta.id}.txt"),
        emit: reports

    path "${meta.id}.hicup_truncate.log", emit: logs

    script:
    def nofill = settings.nofill ? '--nofill' : ''

    """
    python /opt/hicflow/workflow/scripts/hicupTruncate.py \
        ${nofill} \
        --re1 '${settings.restrictionSequence}' \
        --hicup /opt/hicflow/workflow/scripts/hicup_truncater \
        --output \
            ${meta.id}-R1.truncated.fastq.gz \
            ${meta.id}-R2.truncated.fastq.gz \
        --threads ${task.cpus} \
        --tmpdir "\$PWD" \
        ${read1} \
        ${read2} \
        > HiCUP_summary_report-${meta.id}.txt \
        2> ${meta.id}.hicup_truncate.log
    """
}
