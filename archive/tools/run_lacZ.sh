#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH -p sigbio
#SBATCH --mem="5G"
#SBATCH -t "48:00:00"
#SBATCH -c 1 
#SBATCH -o lacZ2.log

#/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python extract_incontact_dimers.py -u UPI0000049419 -o LACZcontacts.txt;
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_cluster.py -c 8 -l 0.02 -r 0.9 -i LACZcontacts.txt;

