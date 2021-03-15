# Information on sample metadata

## Barcodes and sample ID information

`metadata/sample_tags.csv` [md5 84392c84ffe13b9c9f7477ce0d7ed656] is a UTF-8 CSV file giving information on inline barcodes used in this 3RAD sequencing.

Use this to extract information on barcodes for each sample well within each plate.

## Datasets covered by this pipeline

This pipeline implements the steps of a fairly generic 3RAD processing pipeline, but it was designed with a specific dataset in mind: a set of 8 pools (plates), 1-3 comprising PJ's samples, and 4-8 comprising Brendan's samples. **Within this dataset, the pipeline is primarily designed to process Brendan's samples in plates 4-8.** The pipeline would require modification to process PJ's samples properly, or other datasets. It is, however, designed so such modifications are likely to be relatively minor.

- Some of PJ's samples are non-uniquely labelled. There are cases where multiple wells (either on the same plate or on different plates) have the same sample ID. Presumably these are different WGA products or something like that. This pipeline is not currently designed to deal with non-uniquely named samples. Depending on what makes biological sense, the pipeline should probably be modified to either concatenate reads from those samples, or rename those samples with unique names, perhaps using plate and well identifiers in conjunction with sample names.

- The species of each sample is listed in sample_tags.csv for Brendan's samples. This is currently empty for PJ's samples. This information is used by the pipeline (in conjunction with genome_mapping in `config.yaml`) to identify the genome to map reads to, and to determine which samples to group for `gstacks`.

- PJ's samples need to be mapped to different genomes than Brendan's. It is likely that this information can simply be added to genome_mapping in `config.yaml` if the species of each sample can be extracted from the sample name or from `sample_tags.csv`.

- Other parameters such as the restriction enzymes and length filtering parameters are currently specified in ```config.yaml``` and apply across all samples. These values may need to be modified if processing other samples, or potentially made contingent on the sample if the values do not apply across all samples.

## Reformatting metadata

To export a metadata file for each plate of samples as required by stacks:

```bash
for plate in {1..8}; do
grep "^${plate}," metadata/sample_tags.csv | awk -v FS=, -v OFS="\t" '{print $3"G",$4"T",$5}' | tr '/' '-' > barcodes/sample_tags_plate${plate}.tsv
done
```

Note addition of G and T at end of tags (between tag and restriction site).

Note also replacement of '/' with '-' since only letters, numbers, '.', '-' and '_' are allowed by ```process_radtags```.

This conversion is carried out by the rule `all_reformat_metadata` and `reformat_metadata`, but note that shell commands in `snakefile` look different e.g. due to use of snakemake wildcards and double braces to escape braces that need to be passed to `bash`.
