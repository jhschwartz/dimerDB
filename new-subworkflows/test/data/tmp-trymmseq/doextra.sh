#!/bin/sh


#for i in 30 40 50 60 70 80;
#for i in 30 70; 
#do
#    mmseqs easy-cluster shorternonredseqs.fasta seqs_$i tmpdir --min-seq-id 0.$i --cov-mode 0;
#done

i=80
mmseqs easy-cluster shorternonredseqs.fasta seqs_$i tmpdir --min-seq-id 0.$i --cov-mode 0;
