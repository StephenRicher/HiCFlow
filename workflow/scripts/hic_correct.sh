#!/usr/bin/env bash

main() {
    threads=1
    iternum=1000
    while getopts ':p:o:i:@:' flag; do
        case "${flag}" in
            p) readonly diagnostic_plot="${OPTARG}" ;;
            o) readonly corrected_matrix="${OPTARG}" ;;
            i) iternum="${OPTARG}" ;;
            @) threads="${OPTARG}" ;;
            *) usage ;;
        esac
    done
    readonly threads
    readonly iternum
    shift "$((OPTIND-1))"

    local input="${1}"

    local lower_threshold=$(hicCorrectMatrix diagnostic_plot \
      --matrix "${input}" --plotName "${diagnostic_plot}" 2>&1 \
    | grep "mad threshold" | awk '{print $NF}')

    if [ -z "${lower_threshold}" ]; then
      fail "No valid lower threshold detected."
    fi

    hicCorrectMatrix correct \
        --matrix "${input}" --correctionMethod ICE \
        --iterNum "${iternum}" --filterThreshold "${lower_threshold}" 5 \
        --outFileName "${corrected_matrix}"
}


fail() {
    local red='\033[0;31m'
    local no_colour='\033[0m'

    tput setaf 1
    >&2 echo "Error in "${0}": "${FUNCNAME[1]}"."
    all_empty "${@}" || >&2 echo "${1}"
    tput sgr0
    exit "${2-1}"
}

main "${@}"
