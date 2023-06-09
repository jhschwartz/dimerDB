#!/bin/bash



$conda_env/bin/python -m snakemake -s snakefile_after_download_only --profile $smk_profile

