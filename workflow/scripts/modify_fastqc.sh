#!/usr/bin/env bash

main() {
    local fastqc_in="${1}"
    local fastqc_out="${2}"
    local sample="${3}"

    # Retrieve fastqc data path to modify
    local path=$(unzip -l "${fastqc_in}" | grep 'fastqc_data.txt' | awk '{print $NF}')

    cp "${fastqc_in}" "${fastqc_out}"

    if [ ! -z "${sample}" ]; then
      TEMP="$(mktemp)"
      unzip -o "${fastqc_in}" "${path}"
      awk -v OFS='\t' -v sample="${sample}" '
        $0 ~ /^Filename/ {$2 = sample} {print $0}' "${path}" > "${TEMP}"
      mv "${TEMP}" "${path}"
      zip -mq "${fastqc_out}" "${path}"
      rmdir $(dirname "${path}")
    fi
}

main "${@}"
