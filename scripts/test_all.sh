#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH --mem=5G
#SBATCH -t 1:00:00
#SBATCH -p sigbio
#SBATCH --ntasks=8

cd test;
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python -m unittest -v;
