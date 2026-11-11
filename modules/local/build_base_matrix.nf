process BUILD_BASE_MATRIX {
    tag "${meta.id}"

    label 'process_medium'

    publishDir "${params.outdir}/qc/hicexplorer", mode: 'copy', pattern: '*_QC'
    publishDir "${params.outdir}/logs/build_base_matrix", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta),
        path(read1Bam),
        path(read2Bam),
        path(chromSizes),
        path(restrictionSites)

    val settings

    output:
    tuple val(meta),
        path("${meta.id}.${settings.baseResolution}-full.h5"),
        emit: matrices

    tuple val(meta),
        path("${meta.id}.${settings.baseResolution}-full_QC"),
        emit: qc

    path "${meta.id}.build_base_matrix.log", emit: logs

    script:
    def keepSelfLigation = settings.keepSelfLigation ? '--keepSelfLigation' : ''
    def keepSelfCircles = settings.keepSelfCircles ? '--keepSelfCircles' : ''

    """
    hicBuildMatrix \
        --samFiles \
            ${read1Bam} \
            ${read2Bam} \
        --restrictionCutFile ${restrictionSites} \
        --restrictionSequence '${settings.restrictionSequence}' \
        --danglingSequence '${settings.danglingSequence}' \
        --maxLibraryInsertSize ${settings.maxLibraryInsertSize} \
        --minDistance ${settings.minDistance} \
        --minMappingQuality ${settings.minMappingQuality} \
        --inputBufferSize ${settings.inputBufferSize} \
        --chromosomeSizes ${chromSizes} \
        ${keepSelfLigation} \
        ${keepSelfCircles} \
        --binSize ${settings.baseResolution} \
        --outFileName ${meta.id}.${settings.baseResolution}-full.h5 \
        --skipDuplicationCheck \
        --QCfolder ${meta.id}.${settings.baseResolution}-full_QC \
        --threads ${task.cpus} \
        > ${meta.id}.build_base_matrix.log \
        2>&1
    """
}
