process ONTAD {
    tag "${meta.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/tads/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*-ontad-full.tad'

    publishDir "${params.outdir}/tads/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*-ontad-full.bed'

    publishDir "${params.outdir}/logs/ontad",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix)

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        path("${meta.id}-${region.name}-${resolution}-ontad-full.bed"),
        emit: domains

    path "${meta.id}-${region.name}-${resolution}-ontad-full.tad",
        emit: raw

    path "${meta.id}-${region.name}-${resolution}.ontad.log",
        emit: logs

    script:
    def prefix = "${meta.id}-${region.name}-${resolution}-ontad-full"
    def regionLength = region.end - region.start
    def onTadChromosome = "chr${region.chromosome.replaceFirst('^chr', '')}"

    """
    set -euo pipefail

    exec > ${meta.id}-${region.name}-${resolution}.ontad.log 2>&1

    python /opt/hicflow/workflow/scripts/H5_to_NxN.py \
        ${matrix} \
        > ${meta.id}-${region.name}-${resolution}.nxn.tsv

    /opt/hicflow/workflow/scripts/OnTAD \
        ${meta.id}-${region.name}-${resolution}.nxn.tsv \
        -o ${prefix} \
        -bedout ${onTadChromosome} \
        ${regionLength} \
        ${resolution}

    python /opt/hicflow/workflow/scripts/reformatOnTAD.py \
        ${prefix}.tad \
        --scale ${region.start} \
        --binSize ${resolution} \
        --chrom '${region.chromosome}' \
        --maxPos ${region.end} \
        > ${prefix}.bed
    """
}
