# snakefile for processing 3RAD data from Brendan and PJ's ant and plant samples

# Chris Baker
# bakerccm@gmail.com
# 29 January 2020

# scripts adapted from those supplied by Jack Boyle

################
# preamble

# loads the stacks module prior to any shell commands
# see https://snakemake.readthedocs.io/en/stable/project_info/faq.html?highlight=shell.prefix#i-want-to-configure-the-behavior-of-my-shell-for-all-rules-how-can-that-be-achieved-with-snakemake
shell.prefix("module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01;")

# get sample names from barcode files
import pandas as pd
SAMPLES={}
PLATES=['plate' + str(i) for i in range(1,9)]
for p in PLATES:
    SAMPLES[p] = pd.read_table("barcodes/sample_tags_" + p + ".tsv", header = None, names = ['row_index','col_index','sample'])['sample'].tolist()

################
rule all:
    input:
        "demultiplexed/plate1",
        "demultiplexed/plate2",
        "demultiplexed/plate3",
        "demultiplexed/plate4",
        "demultiplexed/plate5",
        "demultiplexed/plate6",
        "demultiplexed/plate7",
        "demultiplexed/plate8"

rule dereplicated_files:
    input:
        # file lists can be python lists
        # supply multiple lists separated by , or concatenate with +
        ["dereplicated/" + plate + "/" + sample + ".1.1.fq.gz" for plate in SAMPLES.keys() for sample in SAMPLES[plate]],
        ["dereplicated/" + plate + "/" + sample + ".2.2.fq.gz" for plate in SAMPLES.keys() for sample in SAMPLES[plate]]

################
# step 2 in the stacks pipeline
# remove clones using UMIs, which are in the fastq headers

rule remove_clones:
    input:
        # note that we're just ignoring the remainder files for now
        file1="demultiplexed/{plate}/{sample}.1.fq.gz",
        file2="demultiplexed/{plate}/{sample}.2.fq.gz"
    output:
        # unclear why it adds the extra numbers, but it does
        file1="dereplicated/{plate}/{sample}.1.1.fq.gz",
        file2="dereplicated/{plate}/{sample}.2.2.fq.gz"
    shell:
        """
        clone_filter \
        -1 {input.file1} -2 {input.file2} \
        -o dereplicated/{wildcards.plate}/{wildcards.sample} \
        --null_index -i gzfastq --oligo_len_2 8
        """

        #"clone_filter -1 pair_1 -2 pair_2] -o out_dir [-i type] [-y type] [-D] [-h]

        #"clone_filter [-f in_file | -p in_dir [-P] [-I] | -1 pair_1 -2 pair_2] -o out_dir [-i type] [-y type] [-D] [-h]

        #"clone_filter -P -p ./raw/ -i gzfastq -o ./filtered/ --null_index --oligo_len_2 8"

################
# step 1 in the stacks pipeline
# demultiplex each pool

# assumes raw data files have already been renamed using data/create_links.sh
# to conform to naming convention expected by process_radtags

# this should probably be changed to find files rather than folders

# use --retain_header to retain i5 index for use as UMI in clone_filter

rule demultiplex_plate:
    input:
        sequences="data/links/{plate}",
        barcodes="barcodes/sample_tags_{plate}.tsv"
    output:
        directory("demultiplexed/{plate}")
    shell:
        """
        rm -rf "demultiplexed/{wildcards.plate}"
        mkdir "demultiplexed/{wildcards.plate}"
        process_radtags -P -p {input.sequences} -b {input.barcodes} -o {output} -c -q -r --inline_inline --renz_1 nheI --renz_2 ecoRI --retain_header
        """
################
# step 0
# create links from input files
# to conform to naming convention expected by process_radtags

# this should probably be changed to find files rather than folders

#rule make_fastq_links:
#    input:
#        "data/Rawdata/ ?? fastq files"
#    output:
#        "data/links/{plate}",
#    shell:


################
