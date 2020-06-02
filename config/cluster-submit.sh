#!/usr/bin/env bash

config=/home/u/sr467/scratch/projects/HiC-new-pipeline/config/config.yaml
snakemake --profile slurm --use-conda -kp --notemp --configfile "${config}" \
--conda-prefix /home/u/sr467/scratch/projects/HiC-pipeline/.snakemake \
--cluster-config /home/u/sr467/scratch/projects/HiC-new-pipeline/config/cluster-config.yaml \
--wrapper-prefix file:/home/u/sr467/scratch/snakemake-wrappers/ "${@}"
