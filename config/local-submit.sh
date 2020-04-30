#!/usr/bin/env bash

config=/media/stephen/Data/HiC-subsample/config/local-config.yaml
snakemake --use-conda -kp --notemp --configfile "${config}" \
--conda-prefix /media/stephen/Data/HiC-subsample/analysis/.snakemake "${@}"
