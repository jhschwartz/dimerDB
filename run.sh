#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH -p sigbio
#SBATCH -t 24:00:00
#SBATCH --mem=5G

/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python -m snakemake -c1
