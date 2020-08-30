# Information on sample metadata

## Barcodes and sample ID information

- ```metadata/sample_tags.csv``` [md5 84392c84ffe13b9c9f7477ce0d7ed656] is a UTF-8 CSV file giving information on inline barcodes used in this 3RAD sequencing.

  - use this to extract information on barcodes for each sample well within each plate.


- N.B.!!!! Some of PJ's samples are non-uniquely labelled. There are cases where multiple wells (either on the same plate
or on different plates) have the same sample ID. Presumably these are different WGA products or something like that.
We should think about whether we want to concatenate these as part of the pipeline, and if so make sure that the concatenation
works properly including across plates. Or whether we want to give these fully unique identifiers. (DUR001_b, DUR001_a, etc, or
perhaps with duplicates identified by appending original well numbers).

- For now though: the pipeline appears to demultiplex PJ's samples OK, so just construct the rest of the pipeline for Brendan's
ants and revisit PJs's data later as necessary.

## Reformatting metadata

- To export a single metadata file for PJ for plates 1 through 3:

  ```bash
  grep "^[123]," metadata/sample_tags.csv | awk -v FS=, -v OFS="\t" '{print $1,$3,$4,$5}' > barcodes/sample_tags_PJ.tsv
  ```

  This can be performed by executing the snakemake rule ```reformat_metadata_PJ```.

- To export a metadata file for each plate of samples as required by stacks:

  ```bash
  for plate in {1..8}; do
    grep "^${plate}," metadata/sample_tags.csv | awk -v FS=, -v OFS="\t" '{print $3"G",$4"T",$5}' | tr '/' '-' > barcodes/sample_tags_plate${plate}.tsv
  done
  ```

  Note addition of G and T at end of tags (between tag and restriction site).
  
  Note also replacement of '/' with '-' since only letters, numbers, '.', '-' and '_' are allowed by ```process_radtags```.

  Note that shell commands in ```snakefile``` look different e.g. due to use of snakemake wildcards and double braces to escape braces that need to be passed to ```bash```.
