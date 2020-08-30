#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 8    # CPUs
#SBATCH --mem-per-cpu=4G    # memory
#SBATCH -p pierce    # partition
#SBATCH -t 0-06:00    # runtime
#SBATCH -o snakemake.out    # output file
#SBATCH -e snakemake.err    # error file

snakemake -j 8 --use-conda flagstat_CM_CN_TP
