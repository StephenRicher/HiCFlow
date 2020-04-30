#!/usr/bin/env bash

main() {
  local start="${1}"
  local bin="${2}"
  local input="${3}"

  # Ontad BED coordinates are the END of the bin region
  # Therefore we rescale the intervals to start = n - bin, end = n
  tail -n+2 "${input}" \
      | awk -v start="${start}" -v bin="${bin}" -v OFS='\t' '
          {$2=$2+start; $3=$3+start}
          {print $1, $2-bin, $2, $1, $3-bin, $3, 0}'

}

main "${@}"
