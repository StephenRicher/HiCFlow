#!/usr/bin/env bash

main() {
  local configs=("${@}")

  # Output everything before 'End Sample Specific' in each file.
  for config in "${configs[@]}"; do
    sed '/End Sample Specific/q' "${config}"
  done
  # Print eveything after 'End Sample Specific in one of the files'
  sed -n '/End Sample Specific/,$p' "${configs[0]}"
}

main "${@}"
