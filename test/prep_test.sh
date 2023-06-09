#!/bin/sh
# file to create a test environment that relies on a tiny subsample of the PDB
# instead of the entire libraru output of subworkflow0


set -e;

testpath=$(pwd);
mkdir tmp;


# create after-download snakefile
grep -v "out_sub0" ../snakefile > tmp/snakefile_after_download_only;

# copy config
cp ../config.yaml tmp/;

cd tmp;

# edit config to work with tmp
tmppath=$(pwd);

# edit paths
sed -i "s|basepath:.*|basepath: $tmppath|" config.yaml;
sed -i "s|intermediates_dir:.*|intermediates_dir: $tmppath/intermediates|" config.yaml;
sed -i "s|lib:.*|lib: $testpath/lib|" config.yaml;
sed -i "s|tmscore_db:.*|tmscore_db: $testpath/lib/tmdb|" config.yaml;
sed -i "s|out_dir:.*|out_dir: $tmppath/homodimerDB|" config.yaml;



# copy subworkflows
mkdir workflow;
cp ../../workflow/*.smk workflow/;

# copy scripts
mkdir scripts;
cp ../../scripts/*.py scripts;
cp ../../scripts/*.sh scripts;


# copy bin
cp ../../bin bin/ -r; 


cd $testpath;

