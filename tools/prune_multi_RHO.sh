#!/bin/sh

#SBATCH -A sigbio_project1
#SBATCH -p sigbio_largemem
#SBATCH --exclusive
#SBATCH -t "48:00:00"
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=8
#SBATCH -o multiprune_RHO.log

/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.1 > RHO_thresh_10 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.2 > RHO_thresh_20 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.3 > RHO_thresh_30 & 
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.4 > RHO_thresh_40 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.5 > RHO_thresh_50 & 
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.6 > RHO_thresh_60 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.7 > RHO_thresh_70 &
/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python manual_prune.py -i rhodimers.txt -c 8 -t 0.8 > RHO_thresh_80 &


wait
