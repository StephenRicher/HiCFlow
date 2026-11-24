process ALIGN_AND_FILTER {
    tag "${meta.id}"

    publishDir "${params.outdir}/qc/fastp",
        mode: 'copy',
        pattern: '*.fastp.*'

    publishDir "${params.outdir}/qc/bowtie2",
        mode: 'copy',
        pattern: '*.bowtie2.log'

    publishDir "${params.outdir}/bam",
        mode: 'copy',
        pattern: '*.filtered.bam*'

    input:
    tuple val(meta),
        path(read1),
        path(read2),
        path(reference),
        path(index)

    output:
    tuple val(meta),
        path("${meta.id}.filtered.bam"),
        path("${meta.id}.filtered.bam.bai"),
        emit: bam

    path "${meta.id}.fastp.html",
        emit: fastpHtml

    path "${meta.id}.fastp.json",
        emit: fastpJson

    path "${meta.id}.bowtie2.log",
        emit: bowtieLog

    script:
    """
    set -euo pipefail

    if gzip -l ${read2} \
        | awk 'NR == 2 {exit !(\$2 > 0)}'
    then
        fastp \
            --thread ${task.cpus} \
            --in1 ${read1} \
            --in2 ${read2} \
            --out1 ${meta.id}-R1.trimmed.fastq.gz \
            --out2 ${meta.id}-R2.trimmed.fastq.gz \
            --html ${meta.id}.fastp.html \
            --json ${meta.id}.fastp.json

        bowtie2 \
            --threads ${task.cpus} \
            --sensitive \
            -x hg38 \
            -1 ${meta.id}-R1.trimmed.fastq.gz \
            -2 ${meta.id}-R2.trimmed.fastq.gz \
            2> ${meta.id}.bowtie2.log \
            | samtools view \
                -@ ${task.cpus} \
                -b \
                -q 30 \
                -F 2316 \
                - \
            | samtools sort \
                -@ ${task.cpus} \
                -n \
                -o ${meta.id}.name-sorted.bam \
                -

        samtools fixmate \
            -m \
            ${meta.id}.name-sorted.bam \
            ${meta.id}.fixmate.bam

        samtools sort \
            -@ ${task.cpus} \
            -o ${meta.id}.coordinate-sorted.bam \
            ${meta.id}.fixmate.bam

        samtools markdup \
            -@ ${task.cpus} \
            -r \
            ${meta.id}.coordinate-sorted.bam \
            ${meta.id}.filtered.bam

        rm -f \
            ${meta.id}-R1.trimmed.fastq.gz \
            ${meta.id}-R2.trimmed.fastq.gz \
            ${meta.id}.name-sorted.bam \
            ${meta.id}.fixmate.bam \
            ${meta.id}.coordinate-sorted.bam
    else
        fastp \
            --thread ${task.cpus} \
            --in1 ${read1} \
            --out1 ${meta.id}-R1.trimmed.fastq.gz \
            --html ${meta.id}.fastp.html \
            --json ${meta.id}.fastp.json

        bowtie2 \
            --threads ${task.cpus} \
            --sensitive \
            -x hg38 \
            -U ${meta.id}-R1.trimmed.fastq.gz \
            2> ${meta.id}.bowtie2.log \
            | samtools view \
                -@ ${task.cpus} \
                -b \
                -q 30 \
                -F 2308 \
                - \
            | samtools sort \
                -@ ${task.cpus} \
                -o ${meta.id}.coordinate-sorted.bam \
                -

        samtools markdup \
            -@ ${task.cpus} \
            -r \
            ${meta.id}.coordinate-sorted.bam \
            ${meta.id}.filtered.bam

        rm -f \
            ${meta.id}-R1.trimmed.fastq.gz \
            ${meta.id}.coordinate-sorted.bam
    fi

    samtools index \
        -@ ${task.cpus} \
        ${meta.id}.filtered.bam
    """
}
