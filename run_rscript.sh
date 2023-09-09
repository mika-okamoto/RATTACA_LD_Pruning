#!/bin/bash 
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16
#PBS -N sample_size3
#PBS -j oe
#PBS -q hotel
#PBS -o /projects/ps-palmer/bbjohnson/rattaca/ld_prune/sample_size3.out
#PBS -e /projects/ps-palmer/bbjohnson/rattaca/ld_prune/sample_size3.err

source activate r_env
cd /projects/ps-palmer/bbjohnson/rattaca/ld_prune/scripts
Rscript var_sample_size.r
