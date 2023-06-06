#!/bin/sh

module load cpu;
conda activate ./env;
export OMP_NUM_THREADS=2;

