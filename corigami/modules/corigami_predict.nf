process CORIGAMI_PREDICT {
    tag "${cellLine}:${chromosome}:${windowStart}"

    cpus 2
    memory '12 GB'
    time '24h'

    publishDir "${params.outdir}/corigami/predictions/wild-type",
        mode: 'copy'

    publishDir "${params.outdir}/corigami/predictions/wild-type",
        mode: 'copy',
        pattern: '*-wild-type.npy'

    input:
    tuple val(cellLine),
        val(chromosome),
        val(windowStart),
        path(ctcf),
        path(atac),
        path(resources)

    output:
    tuple val(cellLine),
        val(chromosome),
        val(windowStart),
        path("${cellLine}-${chromosome}-${windowStart}-wild-type.npy"),
        emit: predictions


    script:
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

    mkdir -p prediction-output

    /opt/conda/envs/corigami/bin/python -c \
        "import pyBigWig; pyBigWig.open('${ctcf}').close(); pyBigWig.open('${atac}').close()"

    /opt/conda/envs/corigami/bin/corigami-predict \
        --out prediction-output \
        --celltype '${cellLine}' \
        --chr '${chromosome}' \
        --start ${windowStart} \
        --model "\${model}" \
        --seq "\${sequence_dir}" \
        --ctcf ${ctcf} \
        --atac ${atac}

    prediction=\$(find prediction-output \
        -type f \
        -path '*/prediction/npy/*.npy' \
        -print \
        -quit)

    if [ -z "\${prediction}" ]; then
        echo "C.Origami did not produce a WT prediction matrix" >&2
        exit 1
    fi

    cp \
        "\${prediction}" \
        ${cellLine}-${chromosome}-${windowStart}-wild-type.npy
    """
}
