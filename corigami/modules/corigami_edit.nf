process CORIGAMI_EDIT {
    tag "${cellLine}:${clone}:allele${allele}"

    cpus 2
    memory '7 GB'
    time '24h'
    maxForks 1

    publishDir(
        "${params.outdir}/corigami/predictions/deletions/${cellLine}/${clone}/allele${allele}",
        mode: 'copy',
        pattern: '*-deletion.npy'
    )

    input:
    tuple val(cellLine),
        val(clone),
        val(allele),
        val(chromosome),
        val(windowStart),
        val(deletionStart),
        val(deletionEnd),
        path(ctcf),
        path(atac),
        path(resources)

    val padding

    output:
    tuple val(cellLine),
        val(clone),
        val(allele),
        val(chromosome),
        val(windowStart),
        val(deletionStart),
        val(deletionEnd),
        path("${cellLine}-${clone}-allele${allele}-${chromosome}-${windowStart}-deletion.npy"),
        emit: predictions

    script:
    def deletionWidth = deletionEnd - deletionStart

    if (deletionWidth <= 0) {
        error "Invalid deletion: ${chromosome}:${deletionStart}-${deletionEnd}"
    }

    """
    set -euo pipefail

    model="${resources}/model_weights/corigami_base.ckpt"
    sequence_dir="${resources}/data/hg38/dna_sequence"

    if [ ! -s "\${model}" ]; then
        echo "Missing C.Origami model: \${model}" >&2
        exit 1
    fi

    if [ ! -s "\${sequence_dir}/chr1.fa.gz" ]; then
        echo "Missing C.Origami hg38 sequence files: \${sequence_dir}" >&2
        exit 1
    fi

    /opt/conda/envs/corigami/bin/python -c \
        "import pyBigWig; pyBigWig.open('${ctcf}').close(); pyBigWig.open('${atac}').close()"

    mkdir -p editing-output

    /opt/conda/envs/corigami/bin/corigami-edit \
        --out editing-output \
        --celltype '${cellLine}_${clone}_allele${allele}' \
        --chr '${chromosome}' \
        --start ${windowStart} \
        --model "\${model}" \
        --seq "\${sequence_dir}" \
        --ctcf ${ctcf} \
        --atac ${atac} \
        --del-start ${deletionStart} \
        --del-width ${deletionWidth} \
        --padding '${padding}'

    prediction=\$(find editing-output \
        -type f \
        -path '*/deletion/npy/*.npy' \
        -print \
        -quit)

    if [ -z "\${prediction}" ]; then
        echo "C.Origami did not produce a deletion prediction matrix" >&2
        exit 1
    fi

    cp \
        "\${prediction}" \
        ${cellLine}-${clone}-allele${allele}-${chromosome}-${windowStart}-deletion.npy
    """
}
