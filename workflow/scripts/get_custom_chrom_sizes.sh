#!/usr/bin/env bash

main() {
  sam="${1}"

  samtools view -H "${1}" \
    | awk '$1 ~ /^@SQ/' \
    | sed 's/:/ /g' \
    | awk -v OFS='\t' '{print $3, $5}'
}

main "${@}"
