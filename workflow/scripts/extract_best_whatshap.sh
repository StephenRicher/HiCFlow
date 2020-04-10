#!/usr/bin/env bash

main() {
  local vcf="${1}"
  local log="${2}"

  block=$(grep 'Largest component' "${log}" | awk '{print $(NF-2)}')

  grep -e '^#' -e "${block}"$ "${vcf}"

}

main "${@}"
