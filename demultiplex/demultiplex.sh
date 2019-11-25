#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 1    # CPUs
#SBATCH --mem-per-cpu=16G    # memory
#SBATCH -p pierce    # partition
#SBATCH -t 0-12:00    # runtime
#SBATCH -o demultiplex.out    # output file
#SBATCH -e demultiplex.err    # error file
#SBATCH --mail-type=BEGIN    # notifications: BEGIN,END,FAIL,ALL
#SBATCH --mail-type=END    # notifications: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bakerccm@gmail.com

# scripts adapted from those supplied by Jack Boyle
#------------------------------------------------------------------------------#
# This is step 1 in the stacks pipeline. It demultiplexes the gbs data.
#------------------------------------------------------------------------------#

module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01

process_radtags -P -p ../data/Rawdata/pool3_CGAGACTA_CGAGACTA  -b ../barcodes/sample_tags_plate3.tsv -o ./pool3_CGAGACTA_CGAGACTA/ \
    -c -q -r --inline_inline --renz_1 nheI --renz_2 ecoRI
   
process_radtags -P -p ../data/Rawdata/pool6_ACTCCATC_ACTCCATC  -b ../barcodes/sample_tags_plate6.tsv -o ./pool6_ACTCCATC_ACTCCATC/ \
    -c -q -r --inline_inline --renz_1 nheI --renz_2 ecoRI
