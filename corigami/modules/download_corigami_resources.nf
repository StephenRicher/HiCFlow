process DOWNLOAD_CORIGAMI_RESOURCES {
    tag 'C.Origami resources'

    publishDir "${params.outdir}/corigami/resources",
        mode: 'copy',
        pattern: 'resource-manifest.tsv'

    input:
    val resourcesUrl

    output:
    path 'corigami_data',
        emit: resources

    path 'resource-manifest.tsv',
        emit: manifest

    script:
    """
    set -euo pipefail

    curl \
        --fail \
        --location \
        --retry 10 \
        --retry-delay 10 \
        --retry-all-errors \
        --output corigami_data.tar.gz.partial \
        '${resourcesUrl}'

    mv \
        corigami_data.tar.gz.partial \
        corigami_data.tar.gz

    tar \
        --extract \
        --gzip \
        --file corigami_data.tar.gz

    resource_root=\$(find . \
        -maxdepth 2 \
        -type d \
        -name 'corigami_data' \
        | head -n 1)

    if [ -z "\${resource_root}" ]; then
        echo "Could not locate extracted corigami_data directory" >&2
        exit 1
    fi

    if [ "\${resource_root}" != './corigami_data' ]; then
        mv "\${resource_root}" ./corigami_data
    fi

    model=\$(find corigami_data \
        -type f \
        \\( -name '*.ckpt' -o -name '*.pt' -o -name '*.pth' \\) \
        | head -n 1)

    if [ -z "\${model}" ]; then
        echo "Could not locate the pretrained C.Origami model" >&2
        exit 1
    fi

    printf 'resource\\tpath\\n' \
        > resource-manifest.tsv

    printf 'model\\t%s\\n' "\${model}" \
        >> resource-manifest.tsv
    """
}
