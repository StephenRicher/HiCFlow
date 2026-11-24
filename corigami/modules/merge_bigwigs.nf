process MERGE_BIGWIGS {
    tag "${cellLine}:${assay}"

    publishDir "${params.outdir}/signal/final",
        mode: 'copy',
        pattern: '*.CPM.hg38.bw'

    input:
    tuple val(cellLine),
        val(assay),
        path(replicateBigWigs)

    output:
    tuple val(cellLine),
        val(assay),
        path("${cellLine}_${assay}_signal.CPM.hg38.bw"),
        emit: bigwigs

    script:
    def orderedBigWigs = replicateBigWigs.sort { first, second ->
        first.name <=> second.name
    }

    if (orderedBigWigs.size() != 2) {
        error "Expected two ${cellLine} ${assay} BigWigs, found ${orderedBigWigs.size()}"
    }

    """
    set -euo pipefail

    bigwigCompare \
        --bigwig1 ${orderedBigWigs[0]} \
        --bigwig2 ${orderedBigWigs[1]} \
        --operation mean \
        --binSize 10 \
        --outFileName ${cellLine}_${assay}_signal.CPM.hg38.bw

    test -s ${cellLine}_${assay}_signal.CPM.hg38.bw
    """
}
