#!/usr/bin/env bash

set -euo pipefail

results_dir="${1:-results}"

require_file() {
    local path="$1"

    if [[ ! -s "${path}" ]]; then
        echo "Missing or empty file: ${path}" >&2
        exit 1
    fi
}

require_count() {
    local expected="$1"
    local directory="$2"
    local pattern="$3"
    local actual

    if [[ ! -d "${directory}" ]]; then
        echo "Missing directory: ${directory}" >&2
        exit 1
    fi

    actual="$(
        find "${directory}" \
            -maxdepth 1 \
            -type f \
            -name "${pattern}" \
            | wc -l \
            | tr -d ' '
    )"

    if [[ "${actual}" != "${expected}" ]]; then
        echo \
            "Expected ${expected} files matching ${pattern} in ${directory}, found ${actual}" \
            >&2
        exit 1
    fi
}

require_count 4 "${results_dir}/bam" '*.hic.bam'

require_count 6 \
    "${results_dir}/matrix/regions/3000/raw" \
    '*-full-raw.h5'

require_count 6 \
    "${results_dir}/matrix/regions/3000/KR" \
    '*-full-KR.h5'

require_count 6 \
    "${results_dir}/tads/chr3L/3000" \
    '*_score.bedgraph'

require_count 6 \
    "${results_dir}/tads/chr3L/3000" \
    '*-ontad-full.bed'

require_count 12 \
    "${results_dir}/plots/chr3L/3000/pyGenomeTracks" \
    '*.svg'

require_count 2 \
    "${results_dir}/plots/chr3L/3000/subtraction/medianFilter" \
    '*.svg'

require_count 12 \
    "${results_dir}/plots/chr3L/3000/viewpoints" \
    '*.svg'

require_count 12 \
    "${results_dir}/viewpoints/chr3L/3000" \
    '*.bedgraph'

require_count 4 \
    "${results_dir}/plots/chr3L/3000/viewpoints/hiccompare" \
    '*.svg'

require_count 2 \
    "${results_dir}/plots/chr3L/3000/viewpoints/comparisons" \
    '*.svg'

require_count 2 \
    "${results_dir}/viewpoints/chr3L/3000/comparisons" \
    '*.tsv'

require_count 2 \
    "${results_dir}/qc/hicrep/data" \
    '*.csv'

require_count 6 \
    "${results_dir}/qc/matrix_coverage/chr3L" \
    '*.svg'

require_file "${results_dir}/annotation/BDGP6-coloured.bed12"
require_file "${results_dir}/qc/filter_qc/insertSizeFrequency.svg"
require_file "${results_dir}/qc/hicrep/chr3L-3000-hicrep-full.svg"
require_file "${results_dir}/qc/multiqc/multiqc_report.html"

subtraction_metadata="${results_dir}/matrix/compare/chr3L/3000/G1S-minus-AS-chr3L-3000-subtraction.tsv"

require_file "${subtraction_metadata}"

grep -Fqx $'operation\tG1S - AS' "${subtraction_metadata}"
grep -Fqx $'positive_values\tHigher in G1S' "${subtraction_metadata}"
grep -Fqx $'negative_values\tHigher in AS' "${subtraction_metadata}"

while IFS= read -r table; do
    require_file "${table}"

    header="$(head -n 1 "${table}")"

    if [[ "${header}" != $'chromosome\tstart\tend\tAS\tG1S\tdifference' ]]; then
        echo "Unexpected viewpoint comparison header in ${table}: ${header}" >&2
        exit 1
    fi

    awk -F $'\t' '
        NR == 1 {
            next
        }

        NF != 6 {
            exit 1
        }

        {
            expected = $5 - $4
            difference = $6 - expected

            if (difference < -1e-8 || difference > 1e-8) {
                exit 1
            }

            rows++
        }

        END {
            if (rows == 0) {
                exit 1
            }
        }
    ' "${table}"
done < <(
    find "${results_dir}/viewpoints/chr3L/3000/comparisons" \
        -maxdepth 1 \
        -type f \
        -name 'G1S-minus-AS-*.tsv' \
        | sort
)

echo "HiCFlow example validation passed"
