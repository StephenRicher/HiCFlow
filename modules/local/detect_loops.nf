process DETECT_LOOPS {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/loops/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*.bedgraph'

    publishDir "${params.outdir}/logs/detect_loops",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix)

    val settings

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-loops.bedgraph"),
        emit: loops

    path "${meta.id}-${region.name}-${resolution}.detect_loops.log",
        emit: logs

    script:
    def prefix = "${meta.id}-${region.name}-${resolution}"
    def expectedArgument = settings.expected
        ? "--expected '${settings.expected}'"
        : ''

    """
    set -uo pipefail

    hicDetectLoops \
        --matrix ${matrix} \
        --outFileName ${prefix}-unclamped.bedgraph \
        --maxLoopDistance ${settings.maxLoopDistance} \
        --windowSize ${settings.windowSize} \
        --peakWidth ${settings.peakWidth} \
        --pValuePreselection ${settings.pValuePreselection} \
        --obsExpThreshold ${settings.obsExpThreshold} \
        --pValue ${settings.pValue} \
        --peakInteractionsThreshold ${settings.peakInteractionsThreshold} \
        ${expectedArgument} \
        --threads 1 \
        --threadsPerChromosome ${task.cpus} \
        > ${prefix}.detect_loops.log \
        2>&1

    status=\$?

    if [ "\${status}" -ne 0 ]; then
        echo "hicDetectLoops exited with status \${status}; publishing an empty loop file" \
            >> ${prefix}.detect_loops.log

        : > ${prefix}-unclamped.bedgraph
    fi

    python /opt/hicflow/workflow/scripts/clampLoops.py \
        ${region.start} \
        ${region.end} \
        ${prefix}-unclamped.bedgraph \
        > ${prefix}-loops.bedgraph \
        2>> ${prefix}.detect_loops.log
    """
}
