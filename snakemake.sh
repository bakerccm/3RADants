#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 40    # CPUs
#SBATCH --mem=10G    # memory
#SBATCH -p shared    # partition
#SBATCH -t 0-09:00    # runtime
#SBATCH -o snakemake.out    # output file
#SBATCH -e snakemake.err    # error file
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bakerccm@gmail.com

snakemake -j 40 --use-conda all_mapped_sample_stats
