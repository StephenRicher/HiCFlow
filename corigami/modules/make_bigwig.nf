process MAKE_BIGWIG {
    tag "${meta.id}"

    publishDir "${params.outdir}/signal/replicates",
        mode: 'copy',
        pattern: '*.CPM.hg38.bw'

    input:
    tuple val(meta),
        path(bam),
        path(bai)

    output:
    tuple val(meta),
        path("${meta.id}.CPM.hg38.bw"),
        emit: bigwigs

    script:
    """
    set -euo pipefail

    bamCoverage \
        --bam ${bam} \
        --outFileName ${meta.id}.CPM.hg38.bw \
        --outFileFormat bigwig \
        --binSize 10 \
        --normalizeUsing CPM \
        --minMappingQuality 30 \
        --numberOfProcessors ${task.cpus}
    """
}
