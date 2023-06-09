#!/bin/sh
# file to create a test environment that relies on a tiny subsample of the PDB
# instead of the entire libraru output of subworkflow0

workdir=$1;
testpath=$(pwd);

set -e;



# create after-download snakefile
grep -v "out_sub0" $testpath/../snakefile > $workdir/snakefile_after_download_only;

# copy config
cp $testpath/../config.yaml $workdir;

cd $workdir;

# edit paths
sed -i "s|basepath:.*|basepath: $workdir|" config.yaml;
sed -i "s|intermediates_dir:.*|intermediates_dir: $workdir/intermediates|" config.yaml;
sed -i "s|lib:.*|lib: $workdir/lib|" config.yaml;
sed -i "s|tmscore_db:.*|tmscore_db: $workdir/lib/tmdb|" config.yaml;
sed -i "s|out_dir:.*|out_dir: $workdir/homodimerDB|" config.yaml;



# copy subworkflows
mkdir workflow;
cp $testpath/../workflow/*.smk workflow/;

# copy scripts
mkdir scripts;
cp $testpath/../scripts/*.py scripts;
cp $testpath/../scripts/*.sh scripts;


# copy bin
cp $testpath/../bin bin/ -r; 


cd $testpath;

