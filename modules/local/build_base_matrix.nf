process BUILD_BASE_MATRIX {
    tag "${meta.id}:${region.name}"

    label 'process_medium'

    publishDir "${params.outdir}/qc/hicexplorer",
        mode: 'copy',
        pattern: '*_QC'

    publishDir "${params.outdir}/logs/build_base_matrix",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        path(read1Bam),
        path(read2Bam),
        path(chromSizes),
        path(restrictionSites)

    val settings

    output:
    tuple val(meta),
        val(region),
        path("${meta.id}-${region.name}.${settings.baseResolution}-full.h5"),
        emit: matrices

    tuple val(meta),
        val(region),
        path("${meta.id}-${region.name}.${settings.baseResolution}-full_QC"),
        emit: qc

    path "${meta.id}-${region.name}.build_base_matrix.log",
        emit: logs

    script:
    def keepSelfLigation = settings.keepSelfLigation ? '--keepSelfLigation' : ''
    def keepSelfCircles = settings.keepSelfCircles ? '--keepSelfCircles' : ''
    def matrixRegion = "${region.chromosome}:${region.start + 1}-${region.end}"
    def prefix = "${meta.id}-${region.name}.${settings.baseResolution}-full"

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
        --region '${matrixRegion}' \
        ${keepSelfLigation} \
        ${keepSelfCircles} \
        --binSize ${settings.baseResolution} \
        --outFileName ${prefix}.h5 \
        --skipDuplicationCheck \
        --QCfolder ${prefix}_QC \
        --threads ${task.cpus} \
        > ${meta.id}-${region.name}.build_base_matrix.log \
        2>&1
    """
}
