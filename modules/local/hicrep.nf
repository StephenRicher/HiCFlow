process HICREP {
    tag "${meta1.id}-vs-${meta2.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/qc/hicrep/data",
        mode: 'copy',
        pattern: '*.csv'

    publishDir "${params.outdir}/logs/hicrep",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta1),
        val(meta2),
        val(region),
        val(resolution),
        path(matrix1),
        path(matrix2)

    output:
    tuple val(meta1),
        val(meta2),
        val(region),
        val(resolution),
        path("${meta1.id}-vs-${meta2.id}-${region.name}-${resolution}-hicrep-full.csv"),
        emit: scores

    path "${meta1.id}-vs-${meta2.id}-${region.name}-${resolution}.hicrep.log",
        emit: logs

    script:
    def prefix = "${meta1.id}-vs-${meta2.id}-${region.name}-${resolution}"
    def start = region.start + 1

    """
    set -euo pipefail

    python /opt/hicflow/workflow/scripts/H5_to_NxN3p.py \
        ${matrix1} \
        > ${meta1.id}.nxn3p.tsv \
        2> ${prefix}.hicrep.log

    python /opt/hicflow/workflow/scripts/H5_to_NxN3p.py \
        ${matrix2} \
        > ${meta2.id}.nxn3p.tsv \
        2>> ${prefix}.hicrep.log

    Rscript /opt/hicflow/workflow/scripts/runHiCRep.R \
        ${prefix}-hicrep-full.csv \
        ${resolution} \
        ${start} \
        ${region.end} \
        ${meta1.id}.nxn3p.tsv \
        ${meta2.id}.nxn3p.tsv \
        >> ${prefix}.hicrep.log \
        2>&1
    """
}
