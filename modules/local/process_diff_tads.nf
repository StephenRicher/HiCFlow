process PROCESS_DIFF_TADS {
    tag "${group2.id}-minus-${group1.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/tads/${region.name}/${resolution}/differential",
        mode: 'copy',
        pattern: '*.bed'

    publishDir "${params.outdir}/tads/${region.name}/${resolution}/differential",
        mode: 'copy',
        pattern: '*.pkl'

    publishDir "${params.outdir}/logs/process_diff_tads",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path(subtractionMatrix),
        path(rawMatrix1),
        path(rawMatrix2),
        path(tads1),
        path(tads2)

    val threshold

    output:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-${group1.id}-allTAD.bed"),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-${group2.id}-allTAD.bed"),
        emit: allTads

    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-${group1.id}-diffTAD.bed"),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-${group2.id}-diffTAD.bed"),
        emit: differentialTads

    path "${group2.id}-minus-${group1.id}-${region.name}-${resolution}-${group1.id}.pkl",
        emit: group1Data

    path "${group2.id}-minus-${group1.id}-${region.name}-${resolution}-${group2.id}.pkl",
        emit: group2Data

    path "${group2.id}-minus-${group1.id}-${region.name}-${resolution}.process_diff_tads.log",
        emit: logs

    script:
    def subtractionPrefix = "${group2.id}-minus-${group1.id}-${region.name}-${resolution}"

    """
    set -euo pipefail

    exec > ${subtractionPrefix}.process_diff_tads.log 2>&1

    python /opt/hicflow/workflow/scripts/processDiffTAD.py \
        ${subtractionMatrix} \
        ${rawMatrix1} \
        ${rawMatrix2} \
        ${tads1} \
        --threshold ${threshold} \
        --name diffTAD \
        --outDiff ${subtractionPrefix}-${group1.id}-diffTAD.bed \
        --outPickle ${subtractionPrefix}-${group1.id}.pkl \
        > ${subtractionPrefix}-${group1.id}-allTAD.bed

    python /opt/hicflow/workflow/scripts/processDiffTAD.py \
        ${subtractionMatrix} \
        ${rawMatrix1} \
        ${rawMatrix2} \
        ${tads2} \
        --threshold ${threshold} \
        --name diffTAD \
        --outDiff ${subtractionPrefix}-${group2.id}-diffTAD.bed \
        --outPickle ${subtractionPrefix}-${group2.id}.pkl \
        > ${subtractionPrefix}-${group2.id}-allTAD.bed
    """
}
