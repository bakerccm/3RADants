#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 5    # CPUs
#SBATCH --mem-per-cpu=32G    # memory
#SBATCH -p shared    # partition
#SBATCH -t 0-06:00    # runtime
#SBATCH -o snakemake.out    # output file
#SBATCH -e snakemake.err    # error file
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bakerccm@gmail.com

snakemake -j 5 dereplicate_all
