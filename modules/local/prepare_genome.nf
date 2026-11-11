process PREPARE_GENOME {
    tag "${cellType}"

    label 'process_medium'

    publishDir "${params.outdir}/genome",
        mode: 'copy',
        pattern: '*.chrom.sizes'

    publishDir "${params.outdir}/genome",
        mode: 'copy',
        pattern: '*-restSites.bed'

    publishDir "${params.outdir}/logs/prepare_genome",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(cellType), path(reference)
    val restrictionSequence

    output:
    tuple val(cellType),
        path("${cellType}.fa.gz"),
        path("${cellType}.fa.gz.fai"),
        path("${cellType}.chrom.sizes"),
        path("${cellType}-restSites.bed"),
        emit: genome

    path "${cellType}.prepare_genome.log", emit: logs

    script:
    def searchPattern = restrictionSequence.replace('^', '')

    """
    set -euo pipefail

    {
        zcat -f ${reference} \
            | bgzip \
            > ${cellType}.fa.gz

        samtools faidx ${cellType}.fa.gz

        cut \
            -f 1,2 \
            ${cellType}.fa.gz.fai \
            > ${cellType}.chrom.sizes

        hicFindRestSite \
            --fasta ${cellType}.fa.gz \
            --searchPattern '${searchPattern}' \
            --outFile ${cellType}-restSites.bed
    } > ${cellType}.prepare_genome.log 2>&1
    """
}
