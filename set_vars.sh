#!/bin/bash

export basepath=$(pwd)
export conda_env="$basepath/env"
export lib_path="$basepath/lib"

sed -e "s|\!BASEPATH\!|$basepath|" mod/config.yaml.mod > config.yaml;

#export smk_profile=$(realpath .smk_profile_slurm_expanse);
export smk_profile=$(realpath .smk_profile_slurm_lh);

#export use_tmpdir="/scratch/job_{JOBID}";

