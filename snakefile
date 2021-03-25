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
# list(SAMPLES[SAMPLES.plate == 4].index) # gives sample names for plate 4
# list(SAMPLES[SAMPLES.plate.isin([4,5])].index) # gives sample names for plates 4 and 5
# list(SAMPLES[SAMPLES.plate.isin(config['plates'])].index) # gives sample names for plates specified in config file
SAMPLES = pd.read_csv("config/sample_tags.csv", header = 0, index_col = 'sample')

################
# default rule
rule all:
    input:
        expand("out/mapped/{sample}.{ext}", sample = list(SAMPLES[SAMPLES.plate.isin(config['plates'])].index), ext = ["flagstat", "idxstat", "stats"]),
        ["out/populations_1/" + species + "/populations.log" for species in SAMPLES.species.unique() if not pd.isna(species)]

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
        list(RAWDATA[RAWDATA.plate.isin(config['plates'])].index)

rule fastq_link:
   input:
       lambda wildcards: RAWDATA.loc["out/data/plate" + wildcards.plate + "/" + wildcards.file + ".fastq.gz"].original
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
        expand("out/demultiplexed/{sample}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES.plate.isin(config['plates'])].index), read = [1,2])

# generates a demultiplex rule for each plate
# - list of output files varies by plate but this can't be left as a wildcard in the output expand() if a single rule is used for all plates
# - alternatively, if you write a rule that just specifies a single sample's outputs rather than all the output files for a plate, the rule runs once per sample rather than once per plate
# - see discussion at https://stackoverflow.com/questions/41135801/snakemake-best-practice-for-demultiplexing
# - note: requires snakemake 5.31.0 or later for 'name' keyword to work
# runs in ~1 hour per plate (32GB works, 16GB seems to be not enough)
for p in config['plates']:
    rule:
        name:
            'demultiplex_plate' + str(p)
        input:
            ["out/data/plate" + str(p) + "/plate" + str(p) + "_" + r + "_" + s + ".fastq.gz" for r in ["R1", "R2"] for s in ["001", "002", "003"]],
            barcodes = "out/barcodes/sample_tags_plate" + str(p) + ".tsv"
        output:
            # ignores remainder files
            expand("out/demultiplexed/{sample}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES.plate == p].index), read = [1,2])
        params:
            sequences = "out/data/plate" + str(p), # folder containing soft links to original data
            renz_1 = config['renz_1'],
            renz_2 = config['renz_2']
        envmodules:
            "gcc/7.1.0-fasrc01",
            "stacks/2.4-fasrc01"
        shell:
            "process_radtags -P -p {params.sequences} -b {input.barcodes} -o out/demultiplexed -c -q -r --inline_inline --renz_1 {params.renz_1} --renz_2 {params.renz_2} --retain_header"

################
# pipeline step 2: remove clones using UMIs, which are in the fastq headers following process_radtags

# runs on all samples in 5:10h with 1 core, 24Gb (but probably uses less)
rule dereplicate_all:
    input:
        expand("out/dereplicated/{sample}.{read}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES.plate.isin(config['plates'])].index), read = [1,2])

# runs in ~30-60 s per sample, 1 core, 24Gb per core (but probably uses less)
rule dereplicate_sample:
    input:
        # note that we're just ignoring the remainder files for now ... is this what we want?
        file1="out/demultiplexed/{sample}.1.fq.gz",
        file2="out/demultiplexed/{sample}.2.fq.gz"
    output:
        # unclear why it adds the extra numbers, but it does
        file1="out/dereplicated/{sample}.1.1.fq.gz",
        file2="out/dereplicated/{sample}.2.2.fq.gz"
    envmodules:
        "gcc/7.1.0-fasrc01",
        "stacks/2.4-fasrc01"
    shell:
        "clone_filter -1 {input.file1} -2 {input.file2} -o out/dereplicated --null_index -i gzfastq --oligo_len_2 8"

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

# this runs in ~10-20 min per genome
# 2GB memory should be sufficient
rule index_genome:
    input:
        "out/genomes/{genome}.fasta"
    output:
        expand("out/genomes/{{genome}}.{suffix}", suffix = BOWTIE_SUFFIXES)
    envmodules:
        "bowtie2/2.3.2-fasrc02"
    threads:
        1
    shell:
        "bowtie2-build --threads {threads} {input} out/genomes/{wildcards.genome}"

################
# pipeline step 4a: map reads to indexed ant genomes using bowtie, then quality-filter the matches, sort them, and output them as a .bam file

rule map_all_samples:
    input:
        expand("out/mapped/{sample}.bam", sample = list(SAMPLES[SAMPLES.plate.isin(config['plates'])].index))

# 1.5 hours with 40 cores, total 10Gb memory to map and sort all 5 plates of samples, and calculate stats
# i.e. steps 4a and 4b combined
rule map_to_genome:
    input:
        # get ant species from SAMPLES and then use dict from config.yaml to determine which genome to map to
        lambda wildcards: ["out/genomes/" + config['genome_mapping'][SAMPLES.loc[wildcards.sample].species] + "." + suffix for suffix in BOWTIE_SUFFIXES],
        fastq1="out/dereplicated/{sample}.1.1.fq.gz",
        fastq2="out/dereplicated/{sample}.2.2.fq.gz"
    output:
        temp("out/mapped/{sample}.sam") # make this temporary since the bam file after sorting is much smaller
    params:
        genome=lambda wildcards: "out/genomes/" + config['genome_mapping'][SAMPLES.loc[wildcards.sample].species], # stem for genome files to supply to bowtie2 command
        min_length=config['bowtie2']['min_length'],
        max_length=config['bowtie2']['max_length']
    envmodules:
        "bowtie2/2.3.2-fasrc02"
    threads: 4
    log:
        "out/mapped/{sample}.log"
    shell:
        # Command from Jack specified --end-to-end and --very-sensitive-local but these seem mutually exclusive.
        # Instead try --end-to-end and --very-sensitive, per Jack's suggestion by email 5 Feb 2020.
        # can't use module load bowtie2/2.2.6-fasrc01 as I think --threads was only added to bowtie-build in 2.2.7
        '''
        bowtie2 --end-to-end --very-sensitive -p {threads} -I {params.min_length} -X {params.max_length} \
        -x {params.genome} -1 {input.fastq1} -2 {input.fastq2} -S {output} 2>{log}
        '''

# sorts mapped reads (.sam file from map_to_genome is discarded after this completes)
# note: need to use --use-conda
# ~ 1 min per sample, 1GB memory is sufficient
rule sort_mapped_reads:
    input:
        "out/mapped/{sample}.sam"
    output:
        "out/mapped/{sample}.bam"
    conda:
        "envs/samtools.yaml"
    log:
        "out/mapped/{sample}_sort.log"
    shell:
        "samtools view -hq 5 {input} | samtools sort - -o {output} 2>{log}"

################
# pipeline step 4b: examine stats for mappings to genomes

rule all_mapped_sample_stats:
    input:
        expand("out/mapped/{sample}.{ext}", sample = list(SAMPLES[SAMPLES.plate.isin(config['plates'])].index), ext = ["flagstat", "idxstat", "stats"])

rule mapped_sample_index:
    input:
        "out/mapped/{sample}.bam"
    output:
        "out/mapped/{sample}.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule mapped_sample_flagstat:
    input:
        "out/mapped/{sample}.bam"
    output:
        "out/mapped/{sample}.flagstat"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools flagstat -O tsv {input} >{output}"

rule mapped_sample_idxstat:
    input:
        bam="out/mapped/{sample}.bam",
        bai="out/mapped/{sample}.bam.bai"
    output:
        "out/mapped/{sample}.idxstat"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools idxstat {input.bam} >{output}"

rule mapped_sample_stats:
    input:
        "out/mapped/{sample}.bam"
    output:
        "out/mapped/{sample}.stats"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools stats {input} >{output}"

################
# make population maps for stacks
# (each ant species gets a separate map, but all individuals are labelled as being in the same population here)

rule all_population_maps:
    input:
        ["out/population_maps/" + species + ".tsv" for species in SAMPLES.species.unique() if not pd.isna(species)]

rule make_population_maps:
    input:
        "config/sample_tags.csv"
    output:
        "out/popmaps_1/{species}.tsv"
    run:
        popmap = SAMPLES[SAMPLES.species == wildcards.species].index.to_frame(index=False)
        popmap['population'] = 1
        popmap.to_csv("out/popmaps_1/" + wildcards.species + ".tsv", index = False, sep = "\t", header = False)

################
# run gstacks and populations
# these rules replace a single run of ref_map.pl:
#   ref_map.pl --samples out/mapped --popmap out/population_maps/TP.tsv -o out/gstacks -T 8 -d

# runs for ~20 min with 8 cores, total 32G memory (2 threads per job)
rule all_gstacks:
    input:
        ["out/gstacks/" + species + "/catalog.fa.gz" for species in SAMPLES.species.unique() if not pd.isna(species)],
        ["out/gstacks/" + species + "/catalog.calls" for species in SAMPLES.species.unique() if not pd.isna(species)]

rule gstacks:
    input:
        lambda wildcards: ["out/mapped/" + sample + ".bam" for sample in list(SAMPLES[SAMPLES.species == wildcards.species].index)],
        popmap = "out/popmaps_1/{species}.tsv"
    output:
        "out/gstacks/{species}/catalog.fa.gz",
        "out/gstacks/{species}/catalog.calls"
    envmodules:
        "gcc/7.1.0-fasrc01",
        "stacks/2.4-fasrc01"
    threads: 2
    shell:
        "gstacks -I out/mapped -M {input.popmap} -O out/gstacks/{wildcards.species} -t {threads}"

# runs for ~2 min with 8 cores, total 32G memory (2 threads per job)
rule all_populations_1:
    input:
        ["out/populations_1/" + species + "/populations.log" for species in SAMPLES.species.unique() if not pd.isna(species)]

# produces genepop output (among other things) which works with Jack's genepopper.pl script
rule populations_1:
    input:
        "out/gstacks/{species}/catalog.fa.gz",
        "out/gstacks/{species}/catalog.calls",
        popmap = "out/population_maps/{species}.tsv"
    output:
        "out/populations_1/{species}/populations.log",
        "out/populations_1/{species}/populations.snps.genepop"
    envmodules:
        "gcc/7.1.0-fasrc01",
        "stacks/2.4-fasrc01"
    threads: 2
    shell:
        '''
        populations -P out/gstacks/{wildcards.species} -M {input.popmap} \
        -O out/populations_1/{wildcards.species} -t {threads} --genepop
        '''

rule genepopper_1:
    input:
        "out/populations_1/{species}/populations.snps.genepop"
    output:
        log = "out/populations_1/{species}/genepopper.log",
        popmap = "out/popmaps_2/{species}.tsv"
    params:
        cutoff = 0.1
    shell:
        "code/genepopper.pl -i {input} -o {output.log} -c {params.cutoff} -p {output.popmap}"



################
