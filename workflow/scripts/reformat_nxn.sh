#!/usr/bin/env bash

main() {
  local input="${1}"

  zcat "${input}" \
    | cut -f 3- \
    | tail -n +2
}

main "${@}"
