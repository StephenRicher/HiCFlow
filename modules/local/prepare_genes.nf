process PREPARE_GENES {
    tag "${settings.genomeBuild}"

    label 'process_low'

    publishDir "${params.outdir}/annotation",
        mode: 'copy',
        pattern: '*-coloured.bed12'

    publishDir "${params.outdir}/logs/prepare_genes",
        mode: 'copy',
        pattern: '*.log'

    input:
    path gff3
    path regions
    val settings

    output:
    path "${settings.genomeId}-coloured.bed12",
        emit: annotation

    path "${settings.genomeId}.prepare_genes.log",
        emit: logs

    script:
    """
    set -euo pipefail

    python /opt/hicflow/workflow/scripts/gff3_to_bed12.py \
        ${gff3} \
        ${regions} \
        ${settings.genomeId}-coloured.bed12 \
        > ${settings.genomeId}.prepare_genes.log \
        2>&1
    """
}
