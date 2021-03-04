# snakefile for processing 3RAD data from Brendan Kenyan ant samples

# Chris Baker
# bakerccm@gmail.com
# 4 March 2021

# scripts adapted from those supplied by Jack Boyle

################
## get config information

# config file
# config['plates'] gives the plate numbers to run (i.e. [4,5,6,7,8] to run Brendan's samples)
configfile: 'config/config.yaml'

import pandas as pd

# get mappings from original raw data files to soft links used by process_radtags
RAWDATA = pd.read_csv("config/rawdata_filenames.csv", header = 0, index_col = 'link')

# sample metadata
# list(SAMPLES[SAMPLES['plate'] == 4].index) # gives sample names for plate 4
# list(SAMPLES[SAMPLES['plate'].isin([4,5])].index) # gives sample names for plates 4 and 5
# list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index) # gives sample names for plates specified in config file
SAMPLES = pd.read_csv("config/sample_tags.csv", header = 0, index_col = 'sample')

################
# default rule
# currently the same as rule all_mapped_sample_stats as this should request all files in the pipeline so far
rule all:
    input:
        expand("out/mapped/{sample}.flagstat", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index))

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
# create links from original raw data files to conform to naming convention expected by process_radtags

rule all_fastq_links:
    input:
        list(RAWDATA[RAWDATA['plate'].isin(config['plates'])].index)

rule fastq_link:
   input:
       lambda wildcards: RAWDATA.loc["out/data/plate" + wildcards.plate + "/" + wildcards.file + ".fastq.gz"]['original']
   output:
       "out/data/plate{plate}/{file}.fastq.gz"
   shell:
        # note use of -r to get relative link is not available in all versions of ln
        "ln -sr {input} {output}"

################
# pipeline step 1: demultiplex each plate

# this remove adaptors and trims everything to the same length
# R1 reads with shorter barcodes get a few extra bp trimmed at the 3' expand
# so that all the R1 reads end up being 140bp (since the longest barcode is 9+1 bp).

# use --retain_header to retain i5 index for use as UMI in clone_filter

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
# pipeline step 2: remove clones using UMIs, which are in the fastq headers following process_radtags

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
# pipeline step 3: prepare ant genomes from Richard ahead of read mapping

# suffixes for filenames output by bowtie2-build
BOWTIE_SUFFIXES = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

rule prepare_genomes:
    input:
        expand("out/genomes/{genome}.{suffix}", genome = ["CM", "CN", "TP"], suffix = BOWTIE_SUFFIXES)

# temporarily decompress ant genome fasta.gz file since bowtie2-build requires fasta (not fasta.gz) format
# note use of config['genome_data'] to get original genome file names
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
# pipeline step 4a: map reads to indexed ant genomes using bowtie, then quality-filter the matches, sort them, and output them as a .bam file

rule map_all_samples:
    input:
        expand("out/mapped/{sample}.bam", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index))

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
# note: need to use --use-conda
rule sort_mapped_reads:
    input:
        "out/mapped/{sample}.sam"
    output:
        "out/mapped/{sample}.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -hq 5 {input} | samtools sort - -o {output}"

################
# pipeline step 4b: examine stats for mappings to genomes

rule all_mapped_sample_stats:
    input:
        expand("out/mapped/{sample}.flagstat", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index))

rule mapped_sample_stats:
    input:
        "out/mapped/{sample}.bam"
    output:
        flagstat="out/mapped/{sample}.flagstat",
        idxstat="out/mapped/{sample}.idxstat",
        stat="out/mapped/{sample}.stat"
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        samtools flagstat -O tsv {input} >{output.flagstat}
        samtools idxstat {input} >{output.idxstat}
        samtools stat {input} >{output.stat}
        '''

################
