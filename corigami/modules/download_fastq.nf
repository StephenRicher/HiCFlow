process DOWNLOAD_FASTQ {
    tag "${meta.id}"

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

    mkdir -p fastq tmp

    for run in "\${runs[@]}"; do
        prefetch \
            --max-size 100G \
            --output-directory sra \
            "\${run}"

        fasterq-dump \
            --split-files \
            --threads ${task.cpus} \
            --temp tmp \
            --outdir fastq \
            "sra/\${run}"
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
                "Incomplete paired-end FASTQs for ${gsm}: " \
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

    rm -rf fastq tmp sra
    """
}
