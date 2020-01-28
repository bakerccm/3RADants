#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 4    # CPUs
#SBATCH --mem-per-cpu=16G    # memory
#SBATCH -p pierce    # partition
#SBATCH -t 0-06:00    # runtime
#SBATCH -o snakemake.out    # output file
#SBATCH -e snakemake.err    # error file

module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01

snakemake -j 4
