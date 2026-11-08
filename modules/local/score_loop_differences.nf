process SCORE_LOOP_DIFFERENCES {
    tag "${group2.id}-minus-${group1.id}:${region.name}:${resolution}"

    label 'process_low'

    publishDir "${params.outdir}/loops/${region.name}/${resolution}/differential",
        mode: 'copy',
        pattern: '*.interact'

    publishDir "${params.outdir}/loops/${region.name}/${resolution}/differential",
        mode: 'copy',
        pattern: '*.links'

    publishDir "${params.outdir}/logs/score_loop_differences",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path(subtractionMatrix),
        path(group1Loops),
        path(group2Loops)

    val settings

    output:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-loops.interact"),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-higher-in-${group2.id}.links"),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-higher-in-${group1.id}.links"),
        emit: differences

    path "${group2.id}-minus-${group1.id}-${region.name}-${resolution}.score_loop_differences.log",
        emit: logs

    script:
    def prefix = "${group2.id}-minus-${group1.id}-${region.name}-${resolution}"

    """
    python /opt/hicflow/workflow/scripts/scoreLoopDiff.py \
        ${group1Loops} \
        ${group2Loops} \
        --matrix ${subtractionMatrix} \
        --maxLineWidth ${settings.maximumLineWidth} \
        --nBins ${settings.differenceBins} \
        --interactOut ${prefix}-loops.interact \
        --linksUp ${prefix}-higher-in-${group2.id}.links \
        --linksDown ${prefix}-higher-in-${group1.id}.links \
        > ${prefix}.score_loop_differences.log \
        2>&1
    """
}
