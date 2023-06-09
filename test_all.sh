#!/bin/bash 

set -e;

. ./set_vars.sh;

base=$(pwd);
python=$(realpath env/bin/python);

echo "Beginning tests of scripts....";
echo "---------"
cd scripts/test && $python -m unittest -vff;
echo "---------"
echo "Passed tests of scripts.";
echo "";

echo "Beginning tests of subworkflows....";
echo "---------"
cd $base/workflow/test && $python -m unittest -vff;
echo "---------"
echo "Passed tests of subworkflows.";
echo "";

echo "Beginning end-to-end test....";
echo "---------"
cd $base/test && $python -m unittest -vff;
echo "---------"
echo "Passed end-to-end test.";
echo "";

echo "Passed all tests!";



