(1) barcodes and sample ID information

sample_tags.csv [md5 84392c84ffe13b9c9f7477ce0d7ed656] is a UTF-8 CSV file giving information on inline barcodes used in this 3RAD sequencing.
--> use this to extract information on barcodes for each sample well within each plate.

(2) reformat for use with stacks

    # export a file for PJ
        grep "^[123]," sample_tags.csv | awk -v FS=, -v OFS="\t" '{print $1,$3,$4,$5}' > sample_tags_PJ.tsv

    # one file for each plate of samples
    # note addition of G and T at end of tags (between tag and restriction site)
        for plate in {1..8}; do
            grep "^${plate}," sample_tags.csv | awk -v FS=, -v OFS="\t" '{print $3"G",$4"T",$5}' > sample_tags_plate${plate}.tsv
        done
