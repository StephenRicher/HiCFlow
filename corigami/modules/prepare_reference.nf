process PREPARE_REFERENCE {
    tag 'hg38'

    input:
    path reference

    output:
    tuple path('hg38.fa'),
        path('hg38.*.bt2*'),
        emit: reference

    script:
    """
    set -euo pipefail

    gzip -dc ${reference} \
        > hg38.fa

    bowtie2-build \
        --threads ${task.cpus} \
        hg38.fa \
        hg38
    """
}
