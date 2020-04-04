#!/usr/bin/env bash

main() {
  local input="${1}"

  local top_score_line=$(grep -n BLOCK "${input}" \
                          | sort -k 6 -nr \
                          | cut -f 1 -d ':' \
                          | head -n 1)

  awk -v n="${top_score_line}" '
    NR<n {next} NR==n {print;next} /^BLOCK/ {exit} {print}' "${input}"
}

main "${@}"
