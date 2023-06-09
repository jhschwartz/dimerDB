#!/bin/bash

. ./set_vars.sh

$conda_env/bin/python -m snakemake --profile .smk_profile_slurm_expanse 
