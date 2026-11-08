process CSCORE {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/cscore/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*-Cscore_cscore.bed'

    publishDir "${params.outdir}/cscore/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*-switch-input.bedgraph'

    publishDir "${params.outdir}/logs/cscore",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(baseMatrix),
        path(reference),
        path(referenceIndex)

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-Cscore_cscore.bed"),
        emit: scores

    path "${meta.id}-${region.name}-${resolution}-Cscore-reorient.bedgraph",
        emit: bedgraphs

    path "${meta.id}-${region.name}-${resolution}.cscore.log",
        emit: logs

    script:
    def prefix = "${meta.id}-${region.name}-${resolution}-Cscore"
    def chromosome = region.chromosome.replaceFirst('^chr', '')
    def cscoreChromosome = "chr${chromosome}"

    """
    set -euo pipefail

    exec > ${meta.id}-${region.name}-${resolution}.cscore.log 2>&1

    printf '%s\\t%s\\t%s\\t%s\\n' \
        '${region.chromosome}' \
        '${region.start}' \
        '${region.end}' \
        '${region.name}' \
        > region.bed

    python /opt/hicflow/workflow/scripts/generateEqualLengthBed.py \
        region.bed \
        ${resolution} \
        > windows.bed

    samtools faidx \
        ${reference} \
        '${region.chromosome}:${region.start}-${region.end}' \
        > region.fa

    python /opt/hicflow/workflow/scripts/makeHiCsummary.py \
        ${baseMatrix} \
        | gzip \
        > matrix.HiCsummary.gz

    /opt/hicflow/workflow/scripts/CscoreTool1.1 \
        windows.bed \
        <(zcat matrix.HiCsummary.gz) \
        ${prefix} \
        ${task.cpus} \
        1000000 \
        ${cscoreChromosome}

    python /opt/hicflow/workflow/scripts/orientCscore.py \
        ${prefix}_cscore.bedgraph \
        region.fa \
        > ${prefix}-reorient.bedgraph

    python /opt/hicflow/workflow/scripts/Bedgraph_to_Bed.py \
        ${prefix}-reorient.bedgraph \
        > ${prefix}_cscore.bed
    """
}
