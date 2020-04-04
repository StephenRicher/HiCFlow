#!/usr/bin/env bash

# Interleave SAM file - keep only @HQ and @PG and remove duplicates from header
main() {
  read1="${1}"
  read2="${2}"

  paste -d '\n' "${read1}" "${read2}" \
    | awk '/^@HQ|^@SQ/ {if(!seen[$0]++) {print $0}}
           /^@/ {next}
           {print $0}'
}

main "${@}"
