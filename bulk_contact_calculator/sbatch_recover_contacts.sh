#!/bin/sh

#SBATCH -A mia174
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t "24:00:00"
#SBATCH --mem=5G
#SBATCH -o recover_contacts.log

env/bin/python recover_contacts.py
