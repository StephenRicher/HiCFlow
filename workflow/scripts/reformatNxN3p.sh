#!/usr/bin/env bash

main() {
  while getopts ':b:r:' flag; do
      case "${flag}" in
          b) local binsize="${OPTARG}" ;;
          r) local region="${OPTARG}" ;;
          *) usage ;;
      esac
  done
  shift "$((OPTIND-1))"

  local input="${1}"

  zcat -f "${input}" \
    | tail -n +2 \
    | cut -f 3- \
    | awk -v OFS='\t' -v bin="${binsize}" -v chr="${region}" '
        {start=(NR-1)*bin; end=start+bin; print chr, start, end, $0}'
}

main "${@}"
