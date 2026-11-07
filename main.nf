nextflow.enable.dsl = 2

workflow {
    log.info """
    HiCFlow Nextflow
    =================
    Output directory: ${params.outdir}
    """

    Channel
        .of('HiCFlow Nextflow initialisation complete')
        .view()
}
