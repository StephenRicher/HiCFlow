#!/usr/bin/env bash

main() {
    iternum=1000
    upper_threshold=5
    while getopts ':p:o:u:i:' flag; do
        case "${flag}" in
            p) readonly diagnostic_plot="${OPTARG}" ;;
            o) readonly corrected_matrix="${OPTARG}" ;;
            u) upper_threshold="${OPTARG}" ;;
            i) iternum="${OPTARG}" ;;
            *) usage ;;
        esac
    done
    readonly iternum
    readonly upper_threshold
    shift "$((OPTIND-1))"

    local input="${1}"

    local lower_threshold=$(hicCorrectMatrix diagnostic_plot \
      --matrix "${input}" --plotName "${diagnostic_plot}" 2>&1 \
    | grep "mad threshold" | awk '{print $NF}')

    if [ -z "${lower_threshold}" ]; then
      >&2 echo "No valid lower threshold detected."
      exit 1
    fi

    hicCorrectMatrix correct \
        --matrix "${input}" --correctionMethod ICE --iterNum "${iternum}" \
        --filterThreshold "${lower_threshold}" "${upper_threshold}" \
        --outFileName "${corrected_matrix}"
}

main "${@}"
