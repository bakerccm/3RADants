# scripts adapted from those supplied by Jack Boyle

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

################
# step 2 in the stacks pipeline
# remove clones using UMIs, which are in the fastq headers

rule remove_clones:
    input:
        file1="demultiplexed/plate3/COR051.1.fq.gz",
        file2="demultiplexed/plate3/COR051.2.fq.gz"
    output:
        dereplicated/plate3/
    shell:
        "clone_filter -1 pair_1 -2 pair_2] -o out_dir [-i type] [-y type] [-D] [-h]
        #"clone_filter [-f in_file | -p in_dir [-P] [-I] | -1 pair_1 -2 pair_2] -o out_dir [-i type] [-y type] [-D] [-h]
        #"clone_filter -P -p ./raw/ -i gzfastq -o ./filtered/ --null_index --oligo_len_2 8"


        clone_filter \
        -1 demultiplexed/plate3/COR051.1.fq.gz \
        -2 demultiplexed/plate3/COR051.2.fq.gz \
        -o demultiplexed/plate3 --null_index -i gzfastq \
        --oligo_len_2 8



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
