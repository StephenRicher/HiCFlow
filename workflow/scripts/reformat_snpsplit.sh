#!/usr/bin/env bash

main() {
  local input="${1}"

  awk -v OFS=$'\t' '
    substr($10,1,3)=="0|1" {print $3, $1, $2, 1, $4"/"$5}
    substr($10,1,3)=="1|0" {print $3, $1, $2, 1, $5"/"$4}' "${input}"
}

main "${@}"
