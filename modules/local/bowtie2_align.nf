process BOWTIE2_ALIGN {
    tag "${meta.id}"

    label 'process_medium'

    publishDir "${params.outdir}/bam",
        mode: 'copy',
        pattern: '*.hic.bam'

    publishDir "${params.outdir}/qc/bowtie2",
        mode: 'copy',
        pattern: '*.bowtie2.txt'

    publishDir "${params.outdir}/qc/deduplicate",
        mode: 'copy',
        pattern: '*.deduplicate.txt'

    publishDir "${params.outdir}/logs/bowtie2",
        mode: 'copy',
        pattern: '*.bowtie2.log'

    publishDir "${params.outdir}/logs/prepare_bam",
        mode: 'copy',
        pattern: '*.prepare_bam.log'

    input:
    tuple val(cellType), val(meta), path(read1), path(read2), path(indexFiles)

    output:
    tuple val(meta),
        path("${meta.id}.hic.bam"),
        emit: bam

    tuple val(meta),
        path("${meta.id}-R1.bowtie2.txt"),
        path("${meta.id}-R2.bowtie2.txt"),
        emit: bowtieReports

    tuple val(meta),
        path("${meta.id}.deduplicate.txt"),
        emit: deduplicationReports

    path "${meta.id}-R1.bowtie2.log",
        emit: logsR1

    path "${meta.id}-R2.bowtie2.log",
        emit: logsR2

    path "${meta.id}.prepare_bam.log",
        emit: prepareBamLogs

    script:
    def alignContig = params.align_contig ?: ''

    """
    set -euo pipefail

    : > ${meta.id}.prepare_bam.log

    bowtie2 \
        -x ${cellType} \
        -U ${read1} \
        --reorder \
        --threads ${task.cpus} \
        --sensitive \
        2> ${meta.id}-R1.bowtie2.log \
        | gawk \
            -f /opt/hicflow/workflow/scripts/addReadFlag.awk \
            -v flag=0x41 \
        | gawk \
            -v target='${alignContig}' \
            'BEGIN { FS = OFS = "\\t" } /^@/ || target == "" || \$3 == target' \
        | samtools view \
            -1 \
            -o ${meta.id}-R1-addFlag.bam \
            - \
            2>> ${meta.id}.prepare_bam.log

    cp \
        ${meta.id}-R1.bowtie2.log \
        ${meta.id}-R1.bowtie2.txt

    bowtie2 \
        -x ${cellType} \
        -U ${read2} \
        --reorder \
        --threads ${task.cpus} \
        --sensitive \
        2> ${meta.id}-R2.bowtie2.log \
        | gawk \
            -f /opt/hicflow/workflow/scripts/addReadFlag.awk \
            -v flag=0x81 \
        | gawk \
            -v target='${alignContig}' \
            'BEGIN { FS = OFS = "\\t" } /^@/ || target == "" || \$3 == target' \
        | samtools view \
            -1 \
            -o ${meta.id}-R2-addFlag.bam \
            - \
            2>> ${meta.id}.prepare_bam.log

    cp \
        ${meta.id}-R2.bowtie2.log \
        ${meta.id}-R2.bowtie2.txt

    samtools merge \
        -@ ${task.cpus} \
        -n \
        ${meta.id}-merged.bam \
        ${meta.id}-R1-addFlag.bam \
        ${meta.id}-R2-addFlag.bam \
        2>> ${meta.id}.prepare_bam.log

    rm -f \
        ${meta.id}-R1-addFlag.bam \
        ${meta.id}-R2-addFlag.bam

    samtools view \
        -h \
        ${meta.id}-merged.bam \
        2>> ${meta.id}.prepare_bam.log \
        | gawk '
            BEGIN {
                FS = OFS = "\t"
            }

            /^@/ {
                print
                next
            }

            function flush_group(    i) {
                if (has_r1 && has_r2) {
                    for (i = 1; i <= count; i++) {
                        print records[i]
                    }
                }

                delete records
                count = 0
                has_r1 = 0
                has_r2 = 0
            }

            current_name != "" && \$1 != current_name {
                flush_group()
            }

            {
                current_name = \$1
                records[++count] = \$0

                if (and(\$2, 64)) {
                    has_r1 = 1
                }

                if (and(\$2, 128)) {
                    has_r2 = 1
                }
            }

            END {
                flush_group()
            }
        ' \
        | samtools fixmate \
            -p \
            -O SAM \
            - \
            - \
            2>> ${meta.id}.prepare_bam.log \
        | samblaster \
            --removeDups \
            --input /dev/stdin \
            2> ${meta.id}.deduplicate.txt \
        | samtools view \
            -b \
            -@ ${task.cpus} \
            -F 2316 \
            -o ${meta.id}.hic.bam \
            - \
            2>> ${meta.id}.prepare_bam.log

    rm -f ${meta.id}-merged.bam
    """
}
