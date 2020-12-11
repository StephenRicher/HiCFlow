#!/usr/bin/env bash

# Randomly sample N pairs of reads from SAM file.
# Assumes SAM compliant format (no whitespace) and contains
# primary read pairs only, grouped by name.

sam="${1}"
nSamples=${2:-1000000}
grep -v '^#' "${1}" | paste -sd " \n" | shuf -n "${nSamples}" | tr ' ' '\n'
