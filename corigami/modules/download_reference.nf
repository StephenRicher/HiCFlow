process DOWNLOAD_REFERENCE {
    tag 'hg38'

    input:
    val referenceUrl

    output:
    path 'hg38.fa.gz',
        emit: reference

    script:
    """
    set -euo pipefail

    curl \
        --fail \
        --location \
        --retry 10 \
        --retry-delay 10 \
        --retry-all-errors \
        --output hg38.fa.gz.partial \
        '${referenceUrl}'

    mv hg38.fa.gz.partial hg38.fa.gz

    gzip -t hg38.fa.gz
    """
}
