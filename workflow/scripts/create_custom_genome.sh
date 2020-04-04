#!/usr/bin/env bash

main() {
    local capture_regions="${1}"
    local genome="${2}"

    while IFS=$'\t' read -r chr start end region; do

        samtools faidx "${genome}" "${chr}":$((start+1))-"${end}" \
        | sed "1 s/^.*$/>${region}/"

    done <"${capture_regions}"

}

main "${@}"
