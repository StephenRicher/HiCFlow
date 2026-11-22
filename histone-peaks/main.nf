nextflow.enable.dsl = 2

params.samples = '/workflow/samples.tsv'
params.referenceUrl = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz'
params.captureRegion = 'chr1:66299999-70699997'
params.outdir = '/data/results'

process DOWNLOAD_FASTQ {
    tag "${meta.id}"

    cpus 2
    memory '4 GB'
    time '24h'
    maxForks 2

    publishDir "${params.outdir}/metadata",
        mode: 'copy',
        pattern: '*.sra.tsv'

    input:
    tuple val(meta), val(gsm)

    output:
    tuple val(meta),
        path("${meta.id}-R1.fastq.gz"),
        path("${meta.id}-R2.fastq.gz"),
        emit: reads

    path "${meta.id}.sra.tsv",
        emit: metadata

    script:
    """
    set -euo pipefail

    pysradb gsm-to-srr ${gsm} \
        > ${meta.id}.sra.tsv

    mapfile -t runs < <(
        awk -F '\\t' '
            NR > 1 && \$2 ~ /^SRR/ {
                print \$2
            }
        ' ${meta.id}.sra.tsv \
            | sort -u
    )

    if [ "\${#runs[@]}" -eq 0 ]; then
        echo "No SRR accessions found for ${gsm}" >&2
        exit 1
    fi

    mkdir fastq

    for run in "\${runs[@]}"; do
        prefetch \
            --max-size 100G \
            "\${run}"

        mkdir -p fastq/tmp

        fasterq-dump \
            --split-files \
            --threads ${task.cpus} \
            --temp fastq/tmp \
            --disk-limit 500G \
            --disk-limit-tmp 500G \
            --outdir fastq \
            "\${run}"

        rm -rf fastq/tmp
    done

    mapfile -t read1_files < <(
        find fastq \
            -maxdepth 1 \
            -type f \
            -name '*_1.fastq' \
            | sort
    )

    mapfile -t read2_files < <(
        find fastq \
            -maxdepth 1 \
            -type f \
            -name '*_2.fastq' \
            -size +0c \
            | sort
    )

    mapfile -t unpaired_files < <(
        find fastq \
            -maxdepth 1 \
            -type f \
            -name '*.fastq' \
            ! -name '*_1.fastq' \
            ! -name '*_2.fastq' \
            | sort
    )

    if [ "\${#read2_files[@]}" -gt 0 ]; then
        if [ "\${#read1_files[@]}" -ne "\${#read2_files[@]}" ]; then
            echo \
                "Mixed or incomplete paired-end FASTQs for ${gsm}: " \
                "\${#read1_files[@]} R1 and \${#read2_files[@]} R2 files" \
                >&2
            exit 1
        fi

        if [ "\${#unpaired_files[@]}" -gt 0 ]; then
            echo "Mixed paired-end and single-end runs for ${gsm}" >&2
            exit 1
        fi

        cat "\${read1_files[@]}" \
            | gzip -c \
            > ${meta.id}-R1.fastq.gz

        cat "\${read2_files[@]}" \
            | gzip -c \
            > ${meta.id}-R2.fastq.gz
    else
        single_files=(
            "\${read1_files[@]}"
            "\${unpaired_files[@]}"
        )

        if [ "\${#single_files[@]}" -eq 0 ]; then
            echo "No FASTQ files produced for ${gsm}" >&2
            exit 1
        fi

        cat "\${single_files[@]}" \
            | gzip -c \
            > ${meta.id}-R1.fastq.gz

        gzip -c </dev/null \
            > ${meta.id}-R2.fastq.gz
    fi

    rm -rf fastq SRR*
    """
}

process DOWNLOAD_REFERENCE {
    tag 'hg19'

    cpus 1
    memory '1 GB'
    time '4h'

    input:
    val referenceUrl

    output:
    path 'hg19.fa.gz',
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
        --output hg19.fa.gz.partial \
        '${referenceUrl}'

    mv hg19.fa.gz.partial hg19.fa.gz

    gzip -t hg19.fa.gz
    """
}

process PREPARE_REFERENCE {
    tag 'hg19'

    cpus 4
    memory '6 GB'
    time '24h'

    input:
    path reference

    output:
    tuple path('hg19.fa'),
        path('hg19.*.bt2*'),
        emit: reference

    script:
    """
    set -euo pipefail

    gzip -dc ${reference} \
        > hg19.fa

    bowtie2-build \
        --threads ${task.cpus} \
        hg19.fa \
        hg19
    """
}

process ALIGN_AND_FILTER {
    tag "${meta.id}"

    cpus 4
    memory '6 GB'
    time '24h'

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
            -x hg19 \
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
            -x hg19 \
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

process MAKE_BIGWIG {
    tag "${meta.id}"

    cpus 2
    memory '4 GB'
    time '8h'

    publishDir "${params.outdir}/signal/replicates",
        mode: 'copy',
        pattern: '*.CPM.bw'

    input:
    tuple val(meta),
        path(bam),
        path(bai)

    val captureRegion

    output:
    tuple val(meta),
        path("${meta.id}.CPM.bw"),
        emit: bigwigs

    script:
    def matcher = captureRegion =~ /^([^:]+):(\d+)-(\d+)$/

    if (!matcher.matches()) {
        error "Invalid capture region: ${captureRegion}"
    }

    def chromosome = matcher[0][1]
    def start = matcher[0][2]
    def end = matcher[0][3]
    def deepToolsRegion = "${chromosome}:${start}:${end}"

    """
    set -euo pipefail

    bamCoverage \
        --bam ${bam} \
        --outFileName ${meta.id}.CPM.bw \
        --outFileFormat bigwig \
        --binSize 25 \
        --normalizeUsing CPM \
        --minMappingQuality 30 \
        --region '${deepToolsRegion}' \
        --numberOfProcessors ${task.cpus}
    """
}

process MERGE_BIGWIGS {
    tag "${cellLine}:${mark}"

    cpus 1
    memory '2 GB'
    time '4h'

    publishDir "${params.outdir}/signal/final",
        mode: 'copy',
        pattern: '*.CPM.bw'

    publishDir "${params.outdir}/signal/final",
        mode: 'copy',
        pattern: '*.tsv'

    input:
    tuple val(cellLine),
        val(mark),
        path(replicateBigWigs)

    output:
    tuple val(cellLine),
        val(mark),
        path("${cellLine}_${mark}_signal.CPM.bw"),
        path("${cellLine}_${mark}_signal_provenance.tsv"),
        emit: bigwigs

    script:
    def orderedBigWigs = replicateBigWigs.sort { first, second ->
        first.name <=> second.name
    }

    def replicateCount = orderedBigWigs.size()

    if (replicateCount < 1 || replicateCount > 2) {
        error "Expected one or two BigWigs, found ${replicateCount}"
    }

    def outputPrefix = cellLine + '_' + mark

    if (replicateCount == 2) {
        def firstBigWig = orderedBigWigs[0]
        def secondBigWig = orderedBigWigs[1]

        """
        set -euo pipefail

        bigwigCompare \
            --bigwig1 ${firstBigWig} \
            --bigwig2 ${secondBigWig} \
            --operation mean \
            --binSize 25 \
            --outFileName ${outputPrefix}_signal.CPM.bw

        printf 'cell_line\\tmark\\treplicates\\ttrack_type\\n' \
            > ${outputPrefix}_signal_provenance.tsv

        printf '%s\\t%s\\t%s\\t%s\\n' \
            '${cellLine}' \
            '${mark}' \
            '${replicateCount}' \
            'mean of two CPM-normalised replicates' \
            >> ${outputPrefix}_signal_provenance.tsv

        test -s ${outputPrefix}_signal.CPM.bw
        """
    } else {
        def singleBigWig = orderedBigWigs[0]

        """
        set -euo pipefail

        cp \
            ${singleBigWig} \
            ${outputPrefix}_signal.CPM.bw

        printf 'cell_line\\tmark\\treplicates\\ttrack_type\\n' \
            > ${outputPrefix}_signal_provenance.tsv

        printf '%s\\t%s\\t%s\\t%s\\n' \
            '${cellLine}' \
            '${mark}' \
            '${replicateCount}' \
            'single CPM-normalised replicate' \
            >> ${outputPrefix}_signal_provenance.tsv

        test -s ${outputPrefix}_signal.CPM.bw
        """
    }
}

process CALL_PEAKS {
    tag "${meta.id}"

    cpus 1
    memory '4 GB'
    time '8h'

    publishDir "${params.outdir}/peaks/replicates",
        mode: 'copy',
        pattern: '*.peaks.bed'

    publishDir "${params.outdir}/logs/macs2",
        mode: 'copy',
        pattern: '*.macs2.log'

    input:
    tuple val(meta),
        path(bam),
        path(bai)

    output:
    tuple val(meta),
        path("${meta.id}.peaks.bed"),
        emit: peaks

    path "${meta.id}.macs2.log",
        emit: logs

    script:
    def peakArguments = meta.mark == 'CTCF'
        ? '--qvalue 0.01'
        : '--broad --broad-cutoff 0.1 --qvalue 0.05'

    def peakFile = meta.mark == 'CTCF'
        ? "${meta.id}_peaks.narrowPeak"
        : "${meta.id}_peaks.broadPeak"

    """
    set -euo pipefail

    if [ "\$(samtools view -c -f 1 ${bam})" -gt 0 ]; then
        format=BAMPE
    else
        format=BAM
    fi

    echo "Calling ${meta.mark} peaks using MACS2 format \${format}" \
        > ${meta.id}.macs2.log

    macs2 callpeak \
        --treatment ${bam} \
        --format "\${format}" \
        --gsize hs \
        --name ${meta.id} \
        --nomodel \
        --shift 0 \
        --extsize 200 \
        ${peakArguments} \
        >> ${meta.id}.macs2.log \
        2>&1

    test -s ${peakFile}

    cut -f1-3 ${peakFile} \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        > ${meta.id}.peaks.bed
    """
}

process CONSENSUS_PEAKS {
    tag "${cellLine}:${mark}"

    cpus 1
    memory '2 GB'
    time '2h'

    publishDir "${params.outdir}/peaks/final",
        mode: 'copy',
        pattern: '*.{bed,tsv}'

    input:
    tuple val(cellLine),
        val(mark),
        path(replicatePeaks)

    val captureRegion

    output:
    tuple val(cellLine),
        val(mark),
        path("${cellLine}_${mark}_peaks.hg19.bed"),
        path("${cellLine}_${mark}_chr1_capture.hg19.bed"),
        path("${cellLine}_${mark}_provenance.tsv"),
        emit: peaks

    script:
    def orderedPeaks = replicatePeaks.sort { first, second ->
        first.name <=> second.name
    }

    def replicateCount = orderedPeaks.size()

    if (replicateCount < 1 || replicateCount > 2) {
        error "Expected one or two replicates, found ${replicateCount}"
    }

    def captureMatch = captureRegion =~ /^([^:]+):(\d+)-(\d+)$/

    if (!captureMatch.matches()) {
        error "Invalid capture region: ${captureRegion}"
    }

    def chromosome = captureMatch[0][1]
    def start = captureMatch[0][2]
    def end = captureMatch[0][3]
    def outputPrefix = [cellLine, mark].join('_')
    def firstPeak = orderedPeaks[0]
    def secondPeak = replicateCount == 2 ? orderedPeaks[1] : orderedPeaks[0]
    def trackType = replicateCount == 2
        ? 'replicate-overlap consensus'
        : 'single replicate'

    """
    set -euo pipefail

    if [ ${replicateCount} -eq 2 ]; then
        bedtools intersect \
            -a ${firstPeak} \
            -b ${secondPeak} \
            -u \
            | cut -f1-3 \
            | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
            | bedtools merge \
                -i - \
            > ${outputPrefix}_peaks.hg19.bed
    else
        cut -f1-3 ${firstPeak} \
            | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
            | bedtools merge \
                -i - \
            > ${outputPrefix}_peaks.hg19.bed
    fi

    printf '%s\\t%s\\t%s\\n' \
        '${chromosome}' \
        '${start}' \
        '${end}' \
        > capture-region.bed

    bedtools intersect \
        -a ${outputPrefix}_peaks.hg19.bed \
        -b capture-region.bed \
        -u \
        > ${outputPrefix}_chr1_capture.hg19.bed

    printf 'cell_line\\tmark\\treplicates\\ttrack_type\\n' \
        > ${outputPrefix}_provenance.tsv

    printf '%s\\t%s\\t%s\\t%s\\n' \
        '${cellLine}' \
        '${mark}' \
        '${replicateCount}' \
        '${trackType}' \
        >> ${outputPrefix}_provenance.tsv

    test -s ${outputPrefix}_peaks.hg19.bed
    test -s ${outputPrefix}_chr1_capture.hg19.bed
    """
}

workflow {
    samples = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id       : "${row.cell_line}_${row.mark}_rep${row.replicate}",
                cell_line: row.cell_line,
                mark     : row.mark,
                replicate: row.replicate as Integer
            ]

            tuple(meta, row.gsm)
        }

    DOWNLOAD_FASTQ(samples)

    DOWNLOAD_REFERENCE(
        Channel.value(params.referenceUrl)
    )

    PREPARE_REFERENCE(
        DOWNLOAD_REFERENCE.out.reference
    )

    alignmentInputs = DOWNLOAD_FASTQ.out.reads
        .combine(PREPARE_REFERENCE.out.reference)
        .map {
            meta,
            read1,
            read2,
            reference,
            index ->

            tuple(
                meta,
                read1,
                read2,
                reference,
                index
            )
        }

    ALIGN_AND_FILTER(alignmentInputs)

    CALL_PEAKS(ALIGN_AND_FILTER.out.bam)

    MAKE_BIGWIG(
        ALIGN_AND_FILTER.out.bam,
        params.captureRegion
    )

    bigWigInputs = MAKE_BIGWIG.out.bigwigs
        .map { meta, bigWig ->
            tuple(
                meta.cell_line,
                meta.mark,
                bigWig
            )
        }
        .groupTuple(by: [0, 1])

    MERGE_BIGWIGS(bigWigInputs)

    consensusInputs = CALL_PEAKS.out.peaks
        .map { meta, peaks ->
            tuple(
                meta.cell_line,
                meta.mark,
                peaks
            )
        }
        .groupTuple(by: [0, 1])

    CONSENSUS_PEAKS(
        consensusInputs,
        params.captureRegion
    )
}
