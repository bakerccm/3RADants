# scripts adapted from those supplied by Jack Boyle

rule all:
    input:
        "demultiplex/plate1_demultiplexed",
        "demultiplex/plate2_demultiplexed",
        "demultiplex/plate3_demultiplexed",
        "demultiplex/plate4_demultiplexed",
        "demultiplex/plate5_demultiplexed",
        "demultiplex/plate6_demultiplexed",
        "demultiplex/plate7_demultiplexed",
        "demultiplex/plate8_demultiplexed"

################
# step 1 in the stacks pipeline
# demultiplex each pool

# assumes raw data files have already been renamed using data/create_links.sh
# to conform to naming convention expected by process_radtags

rule demultiplex_plate:
    input:
        sequences="data/links/{plate}",
        barcodes="barcodes/sample_tags_{plate}.tsv"
    output:
        directory("demultiplex/{plate}_demultiplexed")
    shell:
        """
        rm -rf "demultiplex/{wildcards.plate}_demultiplexed"
        mkdir "demultiplex/{wildcards.plate}_demultiplexed"
        process_radtags -P -p {input.sequences} -b {input.barcodes} -o {output} -c -q -r --inline_inline --renz_1 nheI --renz_2 ecoRI
        """

################
