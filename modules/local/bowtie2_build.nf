process BOWTIE2_BUILD {
    tag "${cellType}"

    label 'process_medium'

    input:
    tuple val(cellType), path(reference)

    output:
    tuple val(cellType),
        path("${cellType}*.bt2"),
        emit: index

    path "${cellType}.bowtie2-build.log", emit: logs

    script:
    """
    bowtie2-build \
        --threads ${task.cpus} \
        ${reference} \
        ${cellType} \
        > ${cellType}.bowtie2-build.log \
        2>&1
    """
}
