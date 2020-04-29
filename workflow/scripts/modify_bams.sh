#!/usr/bin/env bash

main() {

    while getopts ':r:c:s:e:q:@:' flag; do
        case "${flag}" in
            r) readonly region="${OPTARG}" ;;
            c) readonly chr="${OPTARG}" ;;
            s) readonly start="${OPTARG}" ;;
            e) readonly end="${OPTARG}" ;;
            *) usage ;;
        esac
    done
    shift "$((OPTIND-1))"

    local input="${1}"

    local region_length=$((end - start))
    samtools view -h "${input}" \
      | grep -v "@SQ" \
      | sed "2i @SQ\tSN:${region}\tLN:${region_length}" \
      | awk -v OFS='\t' -v chr="${chr}" \
            -v start="${start}" -v region="${region}" '
          !/^@/ {$3=region; $4=$4-start+1; $8=$8-start+1} {print}' \
      | samtools view -b

}

main "${@}"
