process PROCESS_HIC_QC {
    tag "${meta.id}"

    label 'process_low'

    publishDir "${params.outdir}/logs/process_hic_qc",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta), path(bam)
    val samplePairs

    output:
    tuple val(meta),
        path("${meta.id}.processed.pkl"),
        emit: processed

    path "${meta.id}.process_hic_qc.log",
        emit: logs

    script:
    def alignmentLines = samplePairs * 2

    """
    set -euo pipefail

    samtools view -H ${bam} \
        > ${meta.id}.subsampled.sam

    set +o pipefail

    samtools view ${bam} \
        | head -n ${alignmentLines} \
        >> ${meta.id}.subsampled.sam

    samtools_status=\${PIPESTATUS[0]}

    set -o pipefail

    if [ "\${samtools_status}" -ne 0 ] && [ "\${samtools_status}" -ne 141 ]; then
        echo "samtools view failed with exit status \${samtools_status}" >&2
        exit "\${samtools_status}"
    fi

    python /opt/hicflow/workflow/scripts/processHiC.py \
        '${meta.id}' \
        ${meta.id}.processed.pkl \
        ${meta.id}.subsampled.sam \
        > ${meta.id}.process_hic_qc.log \
        2>&1
    """
}
