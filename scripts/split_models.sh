#!/bin/sh

pdb_assembly_file=$1;
pdb_name=$2;
split_chain_exe=$3;


awk '
    BEGIN{ base_fn = $pdb_name "-model1.pdb"; n = 1 }
    {
       print > fn
       if ($0 ~ "ENDMDL") {
           close (fn)
           n++
           fn = $pdb_name "-model" n ".pdb"
       }
    }
'
