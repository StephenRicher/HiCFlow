#!/usr/bin/env bash

main() {
  local fastq_r1="${1}"
  local fastq_r2="${2}"

  paste <(zcat -f "${fastq_r1}" | tr '\t' ' ') \
        <(zcat -f "${fastq_r2}" | tr '\t' ' ') \
  | paste - - - - \
  | awk -v OFS="\n" -v FS="\t" '
      {print($1, $3, $5, $7, $2, $4, $6, $8)}' \
  | gzip --stdout
}

main "${@}"
