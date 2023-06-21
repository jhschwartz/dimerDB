#!/bin/bash

. ./set_vars.sh

$conda_env/bin/python -c "from datetime import datetime; print(datetime.utcnow())" > build_start.tmp;

$conda_env/bin/python -m snakemake --profile .smk_profile_slurm_expanse;

rm build_start.tmp;
