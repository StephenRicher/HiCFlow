process PLOT_HIC_QC {
    tag 'insert-size-frequency'

    label 'process_low'

    publishDir "${params.outdir}/qc/filter_qc",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/logs/plot_hic_qc",
        mode: 'copy',
        pattern: '*.log'

    input:
    path processedFiles
    val settings

    output:
    path "insertSizeFrequency.${settings.filetype}",
        emit: plot

    path 'insertSizeFrequency.plot_hic_qc.log',
        emit: log

    script:
    def inputArguments = processedFiles.join(' ')

    """
    python /opt/hicflow/workflow/scripts/plotQC.py \
        insertSizeFrequency.${settings.filetype} \
        ${inputArguments} \
        --dpi ${settings.dpi} \
        > insertSizeFrequency.plot_hic_qc.log \
        2>&1
    """
}
