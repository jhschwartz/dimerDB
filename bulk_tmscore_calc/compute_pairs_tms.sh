#!/bin/sh

pairsfile=$1;

# await chengxin's new version of usalign...
../bin/USalign/USalign -fast -mm 1 -ter 1 -outfmt 2 -dir ../lib/rcsb/ $pairsfile
