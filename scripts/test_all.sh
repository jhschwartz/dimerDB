#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH --mem=5G
#SBATCH -t 10:00
#SBATCH -p sigbio

cd test;
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python -m unittest -v -ff;
