process MULTIQC {
    tag 'multiqc'

    label 'process_low'

    publishDir "${params.outdir}/qc/multiqc",
        mode: 'copy'

    publishDir "${params.outdir}/logs/multiqc",
        mode: 'copy',
        pattern: '*.log'

    input:
    path qcInputs
    path multiqcConfig

    output:
    path 'multiqc_report.html',
        emit: report

    path 'multiqc_report_data',
        emit: data

    path 'multiqc.log',
        emit: log

    script:
    """
    multiqc \
        . \
        --config ${multiqcConfig} \
        --filename multiqc_report.html \
        --outdir . \
        --force \
        > multiqc.log \
        2>&1
    """
}
