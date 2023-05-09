#!/bin/sh

#SBATCH -A sigbio_project7
#SBATCH -p sigbio
#SBATCH --mem=5G
#SBATCH -t "4:00:00"
#SBATCH -c 8
#SBATCH -o pgm.log

/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python extract_incontact_dimers.py -u UPI000000105F -o PGMcontacts.tmp;
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_cluster.py -c 8 -l 0.02 -r 0.9 -i PGMcontacts.tmp;

