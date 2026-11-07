nextflow.enable.dsl = 2

params.outdir = 'results'

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
