#!/usr/bin/env bash

main() {
    threads=1

    while getopts ':r:c:s:e:q:@:' flag; do
        case "${flag}" in
            r) readonly region="${OPTARG}" ;;
            c) readonly chr="${OPTARG}" ;;
            s) readonly start="${OPTARG}" ;;
            e) readonly end="${OPTARG}" ;;
            q) readonly qc="${OPTARG}" ;;
            @) threads="${OPTARG}" ;;
            *) usage ;;
        esac
    done
    readonly threads
    shift "$((OPTIND-1))"

    local input="${1}"

    local region_length=$((end - start))
    samtools view -h -@ "${threads}" "${input}" \
      | grep -v "@SQ" \
      | sed "2i @SQ\tSN:${region}\tLN:${region_length}" \
      | awk -v OFS='\t' -v chr="${chr}" \
            -v start="${start}" -v region="${region}" '
          !/^@/ {$3=region; $4=$4-start+1; $8=$8-start+1} {print}' \
      | samtools sort -@ "${threads}"

    local count=$(($(samtools view -c -@ "${threads}" "${input}") / 2))
    local pairs_per_kb=$(("${count}" / ("${region_length}" / 1000)))

    printf '%s\t%s\t%s\t%s\t%s\t\n%s\t%s\t%s\t%s\t%s\t\n' \
            "sample" "capture_region" "valid_hic_pairs" \
            "region_length" "hic_pairs_per_kb" \
            "${input}" "${region}" "${count}" \
            "${region_length}" "${pairs_per_kb}" > "${qc}"
}

main "${@}"
