#!/bin/sh
#echo $1
#a=$(echo $1 | awk '{print $1}')
#b=$(echo $1 | awk '{print $2}')
#mkdir -p $a
#touch $a/members_chains.txt
#echo "$b" >> $a/members_chains.txt

f=$(ls -d 6*)
for d in $f
do
#    sort $d/members_chains.txt > $d/tmp
#    rm $d/members_chains.txt
#    mv $d/tmp $d/members_chains.txt
#     mems="$d/members_chains.txt";
#     cat $mems | xargs -I{} grep {} ../../../../out_all_dimers_chains/expected_dimers.txt > $d/members_dimers.txt;
    #cat $d/members_chains.txt > $d/members_dimers.txt
    #echo "cat $d/members_dimers.txt > $d/representative.tsv"
    awk  '{print $2}' $d/representation.tsv | sort -u > $d/representatives.txt
    #mv $d/representative.tsv $d/representation.tsv
done

#echo $f | xargs -I{} bash -e "sort {}/members_chains.txt > {}/tmp"

