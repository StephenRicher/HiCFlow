process PREPARE_GENES {
    tag "${settings.genomeBuild}"

    label 'process_medium'

    publishDir "${params.outdir}/annotation",
        mode: 'copy',
        pattern: '*-coloured.bed12'

    publishDir "${params.outdir}/logs/prepare_genes",
        mode: 'copy',
        pattern: '*.log'

    input:
    path gff3
    val settings

    output:
    path "${settings.genomeBuild}-coloured.bed12",
        emit: annotation

    path "${settings.genomeBuild}.prepare_genes.log",
        emit: logs

    script:
    """
    set -euo pipefail

    exec > ${settings.genomeBuild}.prepare_genes.log 2>&1

    gff3ToGenePred \
        -attrsOut=${settings.genomeBuild}.attrs \
        ${gff3} \
        ${settings.genomeBuild}.unsorted.genePred

    sort \
        -k2,2 \
        -k4n,4n \
        ${settings.genomeBuild}.unsorted.genePred \
        > ${settings.genomeBuild}.genePred

    genePredToBed \
        ${settings.genomeBuild}.genePred \
        ${settings.genomeBuild}.bed12

    python /opt/hicflow/workflow/scripts/processBED12.py \
        --label '${settings.label}' \
        --typeKey '${settings.typeKey}' \
        --geneID '${settings.geneID}' \
        ${settings.genomeBuild}.bed12 \
        ${settings.genomeBuild}.attrs \
        > ${settings.genomeBuild}-coloured.bed12
    """
}
