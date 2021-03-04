# snakefile for processing 3RAD data from Brendan and PJ's ant and plant samples

# Chris Baker
# bakerccm@gmail.com
# 30 August 2020

# scripts adapted from those supplied by Jack Boyle

################
## get config information

# config file
# config['plates'] gives the plate numbers to run (i.e. [4,5,6,7,8] to run Brendan's samples)
configfile: 'config/config.yaml'

# get metadata
import pandas as pd
# mappings to original raw data files
RAWDATA = pd.read_csv("config/rawdata_filenames.csv", header = 0, index_col = 'link')
# sample metadata
SAMPLES = pd.read_csv("config/sample_tags.csv", header = 0, index_col = 'sample')

#list(SAMPLES[SAMPLES['plate'] == 3].index) # e.g. get sample names for plate 3

#list(SAMPLES[SAMPLES['plate'].isin([1,2,3])].index)

################
# reformat sample metadata files for use with stacks: one for each plate of samples

# - double braces in shell commands are necessary to prevent bash braces from being interpreted as indicating snakemake wildcards
# - note addition of G and T at end of tags (between tag and restriction site)
# - note also replacement of '/' with '-' since only letters, numbers, '.', '-' and '_' are allowed by process_radtags

rule all_reformat_metadata:
    input:
        expand("out/barcodes/sample_tags_plate{plate}.tsv", plate = config['plates'])

rule reformat_metadata:
    input:
        "config/sample_tags.csv"
    output:
        "out/barcodes/sample_tags_plate{plate}.tsv"
    shell:
        '''
        grep "^{wildcards.plate}," {input} | awk -v FS=, -v OFS="\t" '{{print $3"G",$4"T",$5}}' | tr '/' '-' > {output}
        '''

################
# step 0
# create links from input files
# to conform to naming convention expected by process_radtags

rule all_fastq_links:
    input:
        expand("out/data/{link}", link = list(RAWDATA.index))

rule fastq_link:
   input:
       lambda wildcards: 'data/' + RAWDATA.loc[wildcards.link]['original']
   output:
       "out/data/{link}"
   shell:
        # note use of -r to get relative link is not available in all versions of ln
        "ln -sr {input} {output}"

################
# step 1
# demultiplex each plate

# assumes raw data files have already been renamed using data/create_links.sh
# to conform to naming convention expected by process_radtags

# use --retain_header to retain i5 index for use as UMI in clone_filter

# this remove adaptors and trims everything to the same length
# R1 reads with shorter barcodes get a few extra bp trimmed at the 3' expand
# so that all the R1 reads end up being 140bp (since the longest barcode is 9+1 bp).

rule demultiplex_all:
    input:
        expand("out/demultiplexed/{sample}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index), read = [1,2])

# generates a demultiplex rule for each plate
# - list of output files varies by plate but this can't be left as a wildcard in the output expand() if a single rule is used for all plates
# - alternatively, if you write a rule that just specifies a single sample's outputs rather than all the output files for a plate, the rule runs once per sample rather than once per plate
# - see discussion at https://stackoverflow.com/questions/41135801/snakemake-best-practice-for-demultiplexing
# - note: requires snakemake 5.31.0 or later for 'name' keyword to work
for p in config['plates']:
    rule:
        name:
            'demultiplex_plate' + str(p)
        input:
            ["out/data/plate" + str(p) + "/plate" + str(p) + "_" + r + "_" + s + ".fastq.gz" for r in ["R1", "R2"] for s in ["001", "002", "003"]],
            barcodes = "out/barcodes/sample_tags_plate" + str(p) + ".tsv"
        output:
            # ignores remainder files
            expand("out/demultiplexed/{sample}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES['plate'] == p].index), read = [1,2])
        params:
            sequences = "out/data/plate" + str(p), # folder containing soft links to original data
            renz_1 = config['renz_1'],
            renz_2 = config['renz_2']
        shell:
            '''
            module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01
            process_radtags -P -p {params.sequences} -b {input.barcodes} -o out/demultiplexed -c -q -r --inline_inline --renz_1 {params.renz_1} --renz_2 {params.renz_2} --retain_header
            '''

################
# step 2
# remove clones using UMIs, which are in the fastq headers

rule dereplicate_all:
    input:
        expand("out/dereplicated/{sample}.{read}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index), read = [1,2])

rule dereplicate_sample:
    input:
        # note that we're just ignoring the remainder files for now ... is this what we want?
        file1="out/demultiplexed/{sample}.1.fq.gz",
        file2="out/demultiplexed/{sample}.2.fq.gz"
    output:
        # unclear why it adds the extra numbers, but it does
        file1="out/dereplicated/{sample}.1.1.fq.gz",
        file2="out/dereplicated/{sample}.2.2.fq.gz"
    shell:
        """
        module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01
        clone_filter -1 {input.file1} -2 {input.file2} -o out/dereplicated --null_index -i gzfastq --oligo_len_2 8
        """

################
# step 3
# prepare ant genomes from Richard in preparation for read mapping

# suffixes for filenames output by bowtie2-build
BOWTIE_SUFFIXES = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

rule ant_genomes:
    input:
        expand("out/genomes/{genome}.{suffix}", genome = ["CM", "CN", "TP"], suffix = BOWTIE_SUFFIXES)

# temporarily decompress ant genome fasta.gz file since bowtie2-build requires fasta (not fasta.gz) format
# note use of config['genome_data'] to get original file names
rule decompress_genome:
    input:
        lambda wildcards: config['genome_data'][wildcards.genome]
    output:
        temp("out/genomes/{genome}.fasta")
    shell:
        "zcat {input} > {output}"

# this runs in ~10 min for CN, 17 min for TP with one thread
# 2GB memory should be sufficient
rule index_genome:
    input:
        "out/genomes/{genome}.fasta"
    output:
        expand("out/genomes/{{genome}}.{suffix}", suffix = BOWTIE_SUFFIXES)
    threads:
        1
    shell:
        """
        module load bowtie2/2.3.2-fasrc02
        bowtie2-build --threads {threads} {input} out/genomes/{wildcards.genome}
        """

################
# step 4
# map reads to indexed ant genomes using bowtie, then quality-filter the matches,
# sort them, and output them as a .bam file

# conda environment samtools should already be created:
#     conda create --name samtools1.10
#     conda activate samtools1.10
#     # set conda channels
#     conda config --add channels defaults
#     conda config --add channels bioconda
#     conda config --add channels conda-forge
#     conda install samtools==1.10

# need to get these samples names automatically
rule all_mapped_samples:
    input:
        expand("out/mapped/{sample}.bam", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index))

# maps reads to indexed ant genome
# ~2 min per sample with four cores, ~10 min per sample with 1 core
# 2GB memory is plenty (probably uses more like 600MB)
rule map_to_genome:
    input:
        # first two characters of sample name (i.e. wildcards.sample[0:2]) denote ant species
        lambda wildcards: ["out/genomes/" + config['genome_mapping'][wildcards.sample[0:2]] + "." + suffix for suffix in BOWTIE_SUFFIXES],
        fastq1="out/dereplicated/{sample}.1.1.fq.gz",
        fastq2="out/dereplicated/{sample}.2.2.fq.gz"
    output:
        temp("out/mapped/{sample}.sam") # make this temporary since the bam file after sorting is much smaller
    params:
        genome=lambda wildcards: "out/genomes/" + config['genome_mapping'][wildcards.sample[0:2]], # stem for genome files to supply to bowtie2 command
        min_length=config['bowtie2']['min_length'],
        max_length=config['bowtie2']['max_length']
    threads: 4
    #benchmark:
    #    "out/mapped/{antsp}/{sample}.map.benchmark.txt"
    log:
        "out/mapped/{sample}.log"
    shell:
        # Command from Jack specified --end-to-end and --very-sensitive-local but these seem mutually exclusive.
        # Instead try --end-to-end and --very-sensitive, per Jack's suggestion by email 5 Feb 2020.
        # can't use module load bowtie2/2.2.6-fasrc01 as I think --threads was only odded to bowtie-build in 2.2.7
        """
        module load bowtie2/2.3.2-fasrc02
        bowtie2 --end-to-end --very-sensitive -p {threads} -I {params.min_length} -X {params.max_length} \
        -x {params.genome} -1 {input.fastq1} -2 {input.fastq2} -S {output} 2>{log}
        """

# sorts mapped reads (.sam file from map_to_genome is discarded after this completes)
# note: need to run snakemake --use-conda
rule sort_mapped_reads:
    input:
        "out/mapped/{sample}.sam"
    output:
        "out/mapped/{sample}.bam"
    #benchmark:
    #    "out/mapped/{antsp}/{sample}.sort.benchmark.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -hq 5 {input} | samtools sort - -o {output}"

################

# examine stats for mappings to genomes

rule all_mapped_sample_stats:
    input:
        expand("out/mapped/{sample}.flagstat", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index))

rule mapped_sample_stats:
    input:
        "out/mapped/{antsp}/{sample}.bam"
    output:
        flagstat="out/mapped/{antsp}/{sample}.flagstat",
        idxstat="out/mapped/{antsp}/{sample}.idxstat",
        stat="out/mapped/{antsp}/{sample}.stat"
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        samtools flagstat -O tsv {input} >{output.flagstat}
        samtools idxstat {input} >{output.idxstat}
        samtools stat {input} >{output.stat}
        '''

################
