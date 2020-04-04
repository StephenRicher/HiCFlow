#!/usr/bin/env bash

main() {
  local start="${1}"
  local input="${2}"

  tail -n+2 "${input}" \
      | awk -v start="${start}" -v OFS='\t' '
          {$2=$2+start; $3=$3+start; $7=$7+start; $8=$8+start}
          {print $1, $2, $3, $1, $7, $8, 0}'

}

main "${@}"
