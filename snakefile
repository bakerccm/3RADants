# snakefile for processing 3RAD data from Brendan and PJ's ant and plant samples

# Chris Baker
# bakerccm@gmail.com
# 29 January 2020

# scripts adapted from those supplied by Jack Boyle

################
# preamble

# forces the stacks module to be loaded prior to running any shell commands
# see https://snakemake.readthedocs.io/en/stable/project_info/faq.html?highlight=shell.prefix#i-want-to-configure-the-behavior-of-my-shell-for-all-rules-how-can-that-be-achieved-with-snakemake
shell.prefix("module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01;")

# get sample names from barcode files
import pandas as pd
SAMPLES = {}
PLATES = ['plate' + str(plate_no) for plate_no in range(1,9)]  # plate1 through plate8
for plate in PLATES:
    SAMPLES[plate] = pd.read_table("barcodes/sample_tags_" + plate + ".tsv", header = None,
        names = ['row_index','col_index','sample'])['sample'].tolist()

################
# rules to run parts of the pipeline

# dereplicate rules

rule dereplicate_Brendan:
    input:
        # file lists can be python lists
        # supply multiple lists separated by , or concatenate with +
        ["dereplicated/" + plate + "/" + sample + ".1.1.fq.gz" for plate in ['plate4','plate5','plate6', 'plate7','plate8'] for sample in SAMPLES[plate]],
        ["dereplicated/" + plate + "/" + sample + ".2.2.fq.gz" for plate in ['plate4','plate5','plate6', 'plate7','plate8'] for sample in SAMPLES[plate]]

rule dereplicate_PJ:
    input:
        ["dereplicated/" + plate + "/" + sample + ".1.1.fq.gz" for plate in ['plate1','plate2','plate3'] for sample in SAMPLES[plate]],
        ["dereplicated/" + plate + "/" + sample + ".2.2.fq.gz" for plate in ['plate1','plate2','plate3'] for sample in SAMPLES[plate]]

rule dereplicate_all:
    input:
        ["dereplicated/" + plate + "/" + sample + ".1.1.fq.gz" for plate in PLATES for sample in SAMPLES[plate]],
        ["dereplicated/" + plate + "/" + sample + ".2.2.fq.gz" for plate in PLATES for sample in SAMPLES[plate]]

# demultiplex rules

rule demultiplex_Brendan:
    input:
        expand("demultiplexed/{plate}", plate = ['plate4','plate5','plate6', 'plate7','plate8'])

rule demultiplex_PJ:
    input:
        expand("demultiplexed/{plate}", plate = ['plate1','plate2','plate3'])

rule demultiplex_all:
    input:
        expand("demultiplexed/{plate}", plate = PLATES)

################
# step 2 in the stacks pipeline
# remove clones using UMIs, which are in the fastq headers

rule dereplicate_sample:
    input:
        # note that we're just ignoring the remainder files for now ... is this what we want?
        flag="demultiplexed/{plate}/demultiplex.done",
        file1="demultiplexed/{plate}/{sample}.1.fq.gz",
        file2="demultiplexed/{plate}/{sample}.2.fq.gz"
    output:
        # unclear why it adds the extra numbers, but it does
        file1="dereplicated/{plate}/{sample}.1.1.fq.gz",
        file2="dereplicated/{plate}/{sample}.2.2.fq.gz"
    shell:
        "clone_filter -1 {input.file1} -2 {input.file2} -o dereplicated/{wildcards.plate} --null_index -i gzfastq --oligo_len_2 8"

################
# step 1 in the stacks pipeline
# demultiplex each pool

# assumes raw data files have already been renamed using data/create_links.sh
# to conform to naming convention expected by process_radtags

# use --retain_header to retain i5 index for use as UMI in clone_filter

rule demultiplex_plate:
    input:
        sequences="data/links/{plate}",
        barcodes="barcodes/sample_tags_{plate}.tsv"
    output:
        touch("demultiplexed/{plate}/demultiplex.done")
    shell:
        "process_radtags -P -p {input.sequences} -b {input.barcodes} -o {output} -c -q -r --inline_inline --renz_1 nheI --renz_2 ecoRI --retain_header"
################
# step 0
# create links from input files
# to conform to naming convention expected by process_radtags

#rule make_fastq_links:
#    input:
#        "data/Rawdata/ ?? fastq files"
#    output:
#        "data/links/{plate}",
#    shell:


################
