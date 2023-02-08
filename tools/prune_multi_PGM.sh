#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH -p sigbio_largemem
#SBATCH --exclusive
#SBATCH -t "48:00:00"
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=8
#SBATCH -o multiprune_PGM.log

/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.1 > PGM_thresh_10 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.2 > PGM_thresh_20 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.3 > PGM_thresh_30 & 
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.4 > PGM_thresh_40 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.5 > PGM_thresh_50 & 
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.6 > PGM_thresh_60 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.7 > PGM_thresh_70 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i PGMcontacts.tmp -c 8 -t 0.8 > PGM_thresh_80 &


wait
