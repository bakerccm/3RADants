# snakefile for processing 3RAD data from Brendan and PJ's ant and plant samples

# Chris Baker
# bakerccm@gmail.com
# 29 January 2020

# scripts adapted from those supplied by Jack Boyle

################
# preamble

# forces the stacks module to be loaded prior to running any shell commands
# see https://snakemake.readthedocs.io/en/stable/project_info/faq.html?highlight=shell.prefix#i-want-to-configure-the-behavior-of-my-shell-for-all-rules-how-can-that-be-achieved-with-snakemake
#    shell.prefix("module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01;")
# actually don't do that because it slows down steps that don't need it, e.g. creating symlinks

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

# make links to dereplicated data by ant species

# unclear that grouping by ant sp is the best way to do this ... consider grouping later e.g. when we create stacks
rule link_Brendan:
    input:
        ["dereplicated_byant/" + sample[0:2] + "/" + sample + ".1.1.fq.gz" for plate in ['plate4','plate5','plate6', 'plate7','plate8'] for sample in SAMPLES[plate]],
        ["dereplicated_byant/" + sample[0:2] + "/" + sample + ".2.2.fq.gz" for plate in ['plate4','plate5','plate6', 'plate7','plate8'] for sample in SAMPLES[plate]]

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
# step 1
# demultiplex each plate

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
        """
        module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01
        process_radtags -P -p {input.sequences} -b {input.barcodes} -o {output} -c -q -r --inline_inline --renz_1 nheI --renz_2 ecoRI --retain_header
        """

################
# step 2
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
        """
        module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01
        clone_filter -1 {input.file1} -2 {input.file2} -o dereplicated/{wildcards.plate} --null_index -i gzfastq --oligo_len_2 8
        """

################
# step 3
# prepare ant genomes from Richard in preparation for read mapping

rule index_ant_genomes:
    input:
        "genomes/CN.done", "genomes/TP.done"
        #"genomes/CM.done", "genomes/CN.done", "genomes/TP.done"

# names of original genome files from Richard
genome_names={'CM': 'Cmimosae_FINAL_VER2.1', 'CN': 'Cnig.1', 'TP': 'Tpen_r3.1'}

# temporarily decompress ant genome fasta.gz file since bowtie2-build requires fasta (not fasta.gz) format
rule decompress_ant_genome:
    input:
        "genomes/{orig_genome}.fasta.gz"
    output:
        temp("genomes/{orig_genome}.fasta")
    shell:
        "zcat {input} > {output}"

# this runs in ~10 min for CN, 17 min for TP with one thread
# 2GB memory should be sufficient
rule index_ant_genome:
    input:
        lambda wildcards: "genomes/" + genome_names[wildcards.genome] + ".fasta"
    output:
        touch("genomes/{genome}.done")
    threads:
        1
    benchmark:
        "genomes/{genome}.benchmark.txt"
    shell:
        """
        module load bowtie2/2.3.2-fasrc02
        bowtie2-build --threads {threads} {input} genomes/{wildcards.genome}
        """

################
# step 4
# map reads to indexed ant genomes using bowtie, then quality-filter the matches,
# sort them, and output them as a .bam file

# note that samples need to be mapped to the right ant genome for the species
# organize reads by ant species for this purpose (is this a good idea??)
rule organize_reads_by_antsp:
    input:
        # gets plate label based on sample number
        file1 = lambda wildcards: "dereplicated/" + [plate for plate, samples in SAMPLES.items() if wildcards.sample in samples][0] + "/" + wildcards.sample + ".1.1.fq.gz",
        file2 = lambda wildcards: "dereplicated/" + [plate for plate, samples in SAMPLES.items() if wildcards.sample in samples][0] + "/" + wildcards.sample + ".2.2.fq.gz"
    output:
        # unclear why it adds the extra numbers, but it does
        file1="dereplicated_byant/{antsp}/{sample}.1.1.fq.gz",
        file2="dereplicated_byant/{antsp}/{sample}.2.2.fq.gz"
    shell:
        # note that file paths for the link are relative to link location, not working directory
        """
        ln -s ../../{input.file1} {output.file1}
        ln -s ../../{input.file2} {output.file2}
        """

# conda environment samtools should already be created:
#     conda create --name samtools1.10
#     conda activate samtools1.10
#     # set conda channels
#     conda config --add channels defaults
#     conda config --add channels bioconda
#     conda config --add channels conda-forge
#     conda install samtools==1.10

# genomes should already be indexed at genomes/CM.*, genomes/CN.* and genomes/TP.*

rule map_sort_partial:
    input:
        "mapped/CN/CN.NMW.D7.post.1.bam",
        "mapped/CN/CN.NMW.D17.post.4.bam",
        "mapped/CN/CN.NMW.D131.post.3.bam",
        "mapped/CN/CN.NMW.D131.pre.1.bam",
        "mapped/CN/CN.SMW.885.pre.4.bam",
        "mapped/CN/CN.NMW.D131.pre.1.bam",
        "mapped/CN/CN.NMW.D17.post.3.bam",
        "mapped/CN/CN.NMW.D23.pre.4.bam",
        "mapped/CN/CN.SMW.840.pre.2.bam",
        "mapped/CN/CN.SMW.860.post.3.bam"

# maps reads to indexed ant genome
# ~2 min per sample with four cores, ~10 min per sample with 1 core
# 2GB memory is plenty (probably uses more like 600MB)
rule map_to_genome:
    input:
        genome="genomes/{antsp}.done",
        fastq1="dereplicated_byant/{antsp}/{sample}.1.1.fq.gz",
        fastq2="dereplicated_byant/{antsp}/{sample}.2.2.fq.gz"
    output:
        temp("mapped/{antsp}/{sample}.sam") # make this temporary since the bam file after sorting is much smaller
    params:
        # I think we should expect fragments from about 330bp to 540 bp
        min_length=300,
        max_length=600
    threads: 4
    benchmark:
        "mapped/{antsp}/{sample}.map.benchmark.txt"
    shell:
        # Command from Jack specified --end-to-end and --very-sensitive-local but these seem mutually exclusive.
        # Instead try --end-to-end and --very-sensitive, per Jack's suggestion by email 5 Feb 2020.
        # can't use module load bowtie2/2.2.6-fasrc01 as I think --threads was only odded to bowtie-build in 2.2.7
        """
        module load bowtie2/2.3.2-fasrc02
        bowtie2 --end-to-end --very-sensitive -p {threads} -I {params.min_length} -X {params.max_length} \
        -x genomes/{wildcards.antsp} -1 {input.fastq1} -2 {input.fastq2} -S {output}
        """

# sorts mapped reads (.sam file from map_to_genome is discarded after this completes)
# note: need to run snakemake --use-conda
rule sort_mapped_reads:
    input:
        "mapped/{antsp}/{sample}.sam"
    output:
        "mapped/{antsp}/{sample}.bam"
    benchmark:
        "mapped/{antsp}/{sample}.sort.benchmark.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -hq 5 {input} | samtools sort - -o {output}"

################
