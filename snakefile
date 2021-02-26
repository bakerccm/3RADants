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

# sample metadata
import pandas as pd
SAMPLES = pd.read_csv("config/sample_tags.csv", header = 0, index_col = 'sample')

#list(SAMPLES[SAMPLES['plate'] == 3].index) # e.g. get sample names for plate 3

#list(SAMPLES[SAMPLES['plate'].isin([1,2,3])].index)

# get sample names for Brendan's ants -- all and by ant species
# not sure what the best way to organize this is yet ...
#BRENDAN_SAMPLES={}
#BRENDAN_SAMPLES['all'] = [sample for plate, samples in SAMPLES.items() if plate in PLATES_BRENDAN for sample in samples]
#for antsp in ['CM', 'CN', 'CS', 'TP']:
#    BRENDAN_SAMPLES[antsp] =[sample for sample in BRENDAN_SAMPLES['all'] if sample[0:2]==antsp]

################
# reformat sample metadata files for use with stacks: one for each plate of samples

# - double braces in shell commands are necessary to prevent bash braces from being interpreted as indicating snakemake wildcards
# - note addition of G and T at end of tags (between tag and restriction site)
# - note also replacement of '/' with '-' since only letters, numbers, '.', '-' and '_' are allowed by process_radtags

rule all_reformat_metadata:
    input:
        expand("out/barcodes/sample_tags_plate{plate}.tsv", plate = config['plates']['all'])

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

# rules to run parts of the pipeline

# dereplicate rules

# rule dereplicate_Brendan:
#     input:
#         # file lists can be python lists
#         # supply multiple lists separated by , or concatenate with +
#         ["out/dereplicated/" + plate + "/" + sample + ".1.1.fq.gz" for plate in PLATES_BRENDAN for sample in SAMPLES[plate]],
#         ["out/dereplicated/" + plate + "/" + sample + ".2.2.fq.gz" for plate in PLATES_BRENDAN for sample in SAMPLES[plate]]

# demultiplex rules

# rule demultiplex_Brendan:
#     input:
#         expand("out/demultiplexed/{plate}", plate = PLATES_BRENDAN)

# make links to dereplicated data by ant species

# unclear that grouping by ant sp is the best way to do this ... consider grouping later e.g. when we create stacks
# rule link_Brendan:
#     input:
#         ["out/dereplicated_byant/" + sample[0:2] + "/" + sample + ".1.1.fq.gz" for plate in PLATES_BRENDAN for sample in SAMPLES[plate]],
#         ["out/dereplicated_byant/" + sample[0:2] + "/" + sample + ".2.2.fq.gz" for plate in PLATES_BRENDAN for sample in SAMPLES[plate]]

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

# this remove adaptors and trims everything to the same length
# R1 reads with shorter barcodes get a few extra bp trimmed at the 3' expand
# so that all the R1 reads end up being 140bp (since the longest barcode is 9+1 bp).

rule demultiplex_brendan:
    input:
        expand("out/demultiplexed/{sample}.{read}.fq.gz", sample = list(SAMPLES[SAMPLES['plate'].isin(config['plates'])].index), read = [1,2])

rule demultiplex_plate:
    input:
        sequences = lambda wildcards: "data/links/plate" + str(SAMPLES.loc[wildcards.sample,'plate']),
        barcodes = lambda wildcards: "out/barcodes/sample_tags_plate" + str(SAMPLES.loc[wildcards.sample,'plate']) + ".tsv"
    output:
        # ignores remainder files
        "out/demultiplexed/{sample}.1.fq.gz",
        "out/demultiplexed/{sample}.2.fq.gz"
    params:
        renz_1 = config['renz_1'],
        renz_2 = config['renz_2']
    shell:
        """
        module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01
        process_radtags -P -p {input.sequences} -b {input.barcodes} -o out/demultiplexed -c -q -r --inline_inline --renz_1 {renz_1} --renz_2 {renz_2} --retain_header
        """

################
# step 2
# remove clones using UMIs, which are in the fastq headers

# rule dereplicate_sample:
#     input:
#         # note that we're just ignoring the remainder files for now ... is this what we want?
#         flag="out/demultiplexed/{plate}/demultiplex.done",
#         file1="out/demultiplexed/{plate}/{sample}.1.fq.gz",
#         file2="out/demultiplexed/{plate}/{sample}.2.fq.gz"
#     output:
#         # unclear why it adds the extra numbers, but it does
#         file1="out/dereplicated/{plate}/{sample}.1.1.fq.gz",
#         file2="out/dereplicated/{plate}/{sample}.2.2.fq.gz"
#     shell:
#         """
#         module load gcc/7.1.0-fasrc01 stacks/2.4-fasrc01
#         clone_filter -1 {input.file1} -2 {input.file2} -o out/dereplicated/{wildcards.plate} --null_index -i gzfastq --oligo_len_2 8
#         """

################
# step 3
# prepare ant genomes from Richard in preparation for read mapping

# rule index_ant_genomes:
#     input:
#         "genomes/CM.done", "genomes/CN.done", "genomes/TP.done"

# names of original genome files from Richard
# genome_names={'CM': 'Cmimosae_FINAL_VER2.1', 'CN': 'Cnig.1', 'TP': 'Tpen_r3.1'}

# temporarily decompress ant genome fasta.gz file since bowtie2-build requires fasta (not fasta.gz) format
# rule decompress_ant_genome:
#     input:
#         "genomes/{orig_genome}.fasta.gz"
#     output:
#         temp("genomes/{orig_genome}.fasta")
#     shell:
#         "zcat {input} > {output}"

# this runs in ~10 min for CN, 17 min for TP with one thread
# 2GB memory should be sufficient
# rule index_ant_genome:
#     input:
#         lambda wildcards: "genomes/" + genome_names[wildcards.genome] + ".fasta"
#     output:
#         touch("genomes/{genome}.done")
#     threads:
#         1
#     benchmark:
#         "genomes/{genome}.benchmark.txt"
#     shell:
#         """
#         module load bowtie2/2.3.2-fasrc02
#         bowtie2-build --threads {threads} {input} genomes/{wildcards.genome}
#         """

################
# step 4
# map reads to indexed ant genomes using bowtie, then quality-filter the matches,
# sort them, and output them as a .bam file

# note that samples need to be mapped to the right ant genome for the species
# organize reads by ant species for this purpose (is this a good idea??)
# rule organize_reads_by_antsp:
#     input:
#         # gets plate label based on sample number
#         file1 = lambda wildcards: "out/dereplicated/" + [plate for plate, samples in SAMPLES.items() if wildcards.sample in samples][0] + "/" + wildcards.sample + ".1.1.fq.gz",
#         file2 = lambda wildcards: "out/dereplicated/" + [plate for plate, samples in SAMPLES.items() if wildcards.sample in samples][0] + "/" + wildcards.sample + ".2.2.fq.gz"
#     output:
#         # unclear why it adds the extra numbers, but it does
#         file1="out/dereplicated_byant/{antsp}/{sample}.1.1.fq.gz",
#         file2="out/dereplicated_byant/{antsp}/{sample}.2.2.fq.gz"
#     shell:
#         # note that file paths for the link are relative to link location, not working directory
#         """
#         ln -s ../../{input.file1} {output.file1}
#         ln -s ../../{input.file2} {output.file2}
#         """

# conda environment samtools should already be created:
#     conda create --name samtools1.10
#     conda activate samtools1.10
#     # set conda channels
#     conda config --add channels defaults
#     conda config --add channels bioconda
#     conda config --add channels conda-forge
#     conda install samtools==1.10

# genomes should already be indexed at genomes/CM.*, genomes/CN.* and genomes/TP.*


# need to get these samples names automatically
# rule map_sort:
#     input:
#         expand("out/mapped/CM/{sample}.bam", sample = BRENDAN_SAMPLES['CM']),
#         expand("out/mapped/CN/{sample}.bam", sample = BRENDAN_SAMPLES['CN']),
#         #expand("out/mapped/CS/{sample}.bam", sample = BRENDAN_SAMPLES['CS']),
#         expand("out/mapped/TP/{sample}.bam", sample = BRENDAN_SAMPLES['TP'])

# maps reads to indexed ant genome
# ~2 min per sample with four cores, ~10 min per sample with 1 core
# 2GB memory is plenty (probably uses more like 600MB)
# rule map_to_genome:
#     input:
#         genome=lambda wildcards: "genomes/{antsp}.done",
#         fastq1="out/dereplicated_byant/{antsp}/{sample}.1.1.fq.gz",
#         fastq2="out/dereplicated_byant/{antsp}/{sample}.2.2.fq.gz"
#     output:
#         temp("out/mapped/{antsp}/{sample}.sam") # make this temporary since the bam file after sorting is much smaller
#     params:
#         # I think we should expect fragments from about 330bp to 540 bp
#         min_length=300,
#         max_length=600
#     threads: 4
#     #benchmark:
#     #    "out/mapped/{antsp}/{sample}.map.benchmark.txt"
#     log:
#         "out/mapped/{antsp}/{sample}.log"
#     shell:
#         # Command from Jack specified --end-to-end and --very-sensitive-local but these seem mutually exclusive.
#         # Instead try --end-to-end and --very-sensitive, per Jack's suggestion by email 5 Feb 2020.
#         # can't use module load bowtie2/2.2.6-fasrc01 as I think --threads was only odded to bowtie-build in 2.2.7
#         """
#         module load bowtie2/2.3.2-fasrc02
#         bowtie2 --end-to-end --very-sensitive -p {threads} -I {params.min_length} -X {params.max_length} \
#         -x genomes/{wildcards.antsp} -1 {input.fastq1} -2 {input.fastq2} -S {output} 2>{log}
#         """

# sorts mapped reads (.sam file from map_to_genome is discarded after this completes)
# note: need to run snakemake --use-conda
# rule sort_mapped_reads:
#     input:
#         "out/mapped/{antsp}/{sample}.sam"
#     output:
#         "out/mapped/{antsp}/{sample}.bam"
#     #benchmark:
#     #    "out/mapped/{antsp}/{sample}.sort.benchmark.txt"
#     conda:
#         "envs/samtools.yaml"
#     shell:
#         "samtools view -hq 5 {input} | samtools sort - -o {output}"

################

# examine mappings to genomes

# rule flagstat_CM_CN_TP:
#     input:
#         expand("out/mapped/CM/{sample}.flagstat", sample = BRENDAN_SAMPLES['CM']),
#         expand("out/mapped/CN/{sample}.flagstat", sample = BRENDAN_SAMPLES['CN']),
#         #expand("out/mapped/CS/{sample}.flagstat", sample = BRENDAN_SAMPLES['CS']),
#         expand("out/mapped/TP/{sample}.flagstat", sample = BRENDAN_SAMPLES['TP'])

# rule flagstat_mapped_sample:
#     input:
#         "out/mapped/{antsp}/{sample}.bam"
#     output:
#         "out/mapped/{antsp}/{sample}.flagstat"
#     conda:
#         "envs/samtools.yaml"
#     shell:
#         "samtools flagstat -O tsv {input} >{output}"

# rule idxstat_mapped_sample:
#     input:
#         "out/mapped/{antsp}/{sample}.bam"
#     output:
#         "out/mapped/{antsp}/{sample}.idxstat"
#     conda:
#         "envs/samtools.yaml"
#     shell:
#         "samtools idxstat {input} >{output}"

# rule stat_mapped_sample:
#     input:
#         "out/mapped/{antsp}/{sample}.bam"
#     output:
#         "out/mapped/{antsp}/{sample}.stat"
#     conda:
#         "envs/samtools.yaml"
#     shell:
#         "samtools stat {input} >{output}"

################
