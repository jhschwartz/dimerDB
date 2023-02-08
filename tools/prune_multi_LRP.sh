#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH -p sigbio_largemem
#SBATCH --exclusive
#SBATCH -t "48:00:00"
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=8
#SBATCH -o multiprune.log

/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i lrpdimers.txt -c 8 -t 0.2 > LRP_thresh_20 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i lrpdimers.txt -c 8 -t 0.3 > LRP_thresh_30 & 
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i lrpdimers.txt -c 8 -t 0.4 > LRP_thresh_40 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i lrpdimers.txt -c 8 -t 0.5 > LRP_thresh_50 & 
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i lrpdimers.txt -c 8 -t 0.6 > LRP_thresh_60 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i lrpdimers.txt -c 8 -t 0.7 > LRP_thresh_70 &


wait
