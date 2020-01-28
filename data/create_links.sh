#!/bin/bash

DATADIR="/n/piercefs/protected/Users/cbaker/3RADants/data"

# these commands create soft links to raw data files to rename them prior to demultiplexing

# plate 1

	mkdir -p "$DATADIR/links/plate1"
	cd "$DATADIR/links/plate1"
	
	for read in {1,2}; do
		ln -s ../../Rawdata/pool1_ACGTTACC_ACGTTACC/pool1_ACGTTACC_ACGTTACC_CKDL190141872-1A-ACGTTACC-ACGTTACC_H3VMVBBXX_L7_${read}.fq.gz plate1_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool1_ATAAGGCG_ATAAGGCG/pool1_ATAAGGCG_ATAAGGCG_CKDL190141872-1A-ATAAGGCG-ATAAGGCG_H3VMVBBXX_L7_${read}.fq.gz plate1_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool1_CAGCGATT_CAGCGATT/pool1_CAGCGATT_CAGCGATT_CKDL190141872-1A-CAGCGATT-CAGCGATT_H3VMVBBXX_L7_${read}.fq.gz plate1_R${read}_003.fastq.gz
	done

# plate 2

	mkdir -p "$DATADIR/links/plate2"
	cd "$DATADIR/links/plate2"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool2_CTGTGTTG_CTGTGTTG/pool2_CTGTGTTG_CTGTGTTG_CKDL190141872-1A-CTGTGTTG-CTGTGTTG_H3VMVBBXX_L7_${read}.fq.gz plate2_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool2_CTTACCTG_CTTACCTG/pool2_CTTACCTG_CTTACCTG_CKDL190141872-1A-CTTACCTG-CTTACCTG_H3VMVBBXX_L7_${read}.fq.gz plate2_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool2_TAGTGACC_TAGTGACC/pool2_TAGTGACC_TAGTGACC_CKDL190141872-1A-TAGTGACC-TAGTGACC_H3VMVBBXX_L7_${read}.fq.gz plate2_R${read}_003.fastq.gz
	done

# plate 3

	mkdir -p "$DATADIR/links/plate3"
	cd "$DATADIR/links/plate3"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool3_CGAGACTA_CGAGACTA/pool3_CGAGACTA_CGAGACTA_CKDL190141872-1A-CGAGACTA-CGAGACTA_H3VMVBBXX_L7_${read}.fq.gz plate3_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool3_CGTTGCAA_CGTTGCAA/pool3_CGTTGCAA_CGTTGCAA_CKDL190141872-1A-CGTTGCAA-CGTTGCAA_H3VMVBBXX_L7_${read}.fq.gz plate3_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool3_TGAGGTGT_TGAGGTGT/pool3_TGAGGTGT_TGAGGTGT_CKDL190141872-1A-TGAGGTGT-TGAGGTGT_H3VMVBBXX_L7_${read}.fq.gz plate3_R${read}_003.fastq.gz
	done

# plate 4

	mkdir -p "$DATADIR/links/plate4"
	cd "$DATADIR/links/plate4"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool4_GACATGGT_GACATGGT/pool4_GACATGGT_GACATGGT_CKDL190141872-1A-GACATGGT-GACATGGT_H3VMVBBXX_L7_${read}.fq.gz plate4_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool4_GATCCATG_GATCCATG/pool4_GATCCATG_GATCCATG_CKDL190141872-1A-GATCCATG-GATCCATG_H3VMVBBXX_L7_${read}.fq.gz plate4_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool4_GATTCAGC_GATTCAGC/pool4_GATTCAGC_GATTCAGC_CKDL190141872-1A-GATTCAGC-GATTCAGC_H3VMVBBXX_L7_${read}.fq.gz plate4_R${read}_003.fastq.gz
	done

# plate 5

	mkdir -p "$DATADIR/links/plate5"
	cd "$DATADIR/links/plate5"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool5_GCATGTCT_GCATGTCT/pool5_GCATGTCT_GCATGTCT_CKDL190141872-1A-GCATGTCT-GCATGTCT_H3VMVBBXX_L7_${read}.fq.gz plate5_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool5_GCCTATCA_GCCTATCA/pool5_GCCTATCA_GCCTATCA_CKDL190141872-1A-GCCTATCA-GCCTATCA_H3VMVBBXX_L7_${read}.fq.gz plate5_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool5_TCACGTTC_TCACGTTC/pool5_TCACGTTC_TCACGTTC_CKDL190141872-1A-TCACGTTC-TCACGTTC_H3VMVBBXX_L7_${read}.fq.gz plate5_R${read}_003.fastq.gz
	done

# plate 6

	mkdir -p "$DATADIR/links/plate6"
	cd "$DATADIR/links/plate6"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool6_AACAACCG_AACAACCG/pool6_AACAACCG_AACAACCG_CKDL190141872-1A-AACAACCG-AACAACCG_H3VMVBBXX_L7_${read}.fq.gz plate6_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool6_ACTCCATC_ACTCCATC/pool6_ACTCCATC_ACTCCATC_CKDL190141872-1A-ACTCCATC-ACTCCATC_H3VMVBBXX_L7_${read}.fq.gz plate6_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool6_TGTGCGTT_TGTGCGTT/pool6_TGTGCGTT_TGTGCGTT_CKDL190141872-1A-TGTGCGTT-TGTGCGTT_H3VMVBBXX_L7_${read}.fq.gz plate6_R${read}_003.fastq.gz
	done

# plate 7

	mkdir -p "$DATADIR/links/plate7"
	cd "$DATADIR/links/plate7"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool7_ACTCGTTG_ACTCGTTG/pool7_ACTCGTTG_ACTCGTTG_CKDL190141872-1A-ACTCGTTG-ACTCGTTG_H3VMVBBXX_L7_${read}.fq.gz plate7_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool7_TAGTTGCG_TAGTTGCG/pool7_TAGTTGCG_TAGTTGCG_CKDL190141872-1A-TAGTTGCG-TAGTTGCG_H3VMVBBXX_L7_${read}.fq.gz plate7_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool7_TGTGACTG_TGTGACTG/pool7_TGTGACTG_TGTGACTG_CKDL190141872-1A-TGTGACTG-TGTGACTG_H3VMVBBXX_L7_${read}.fq.gz plate7_R${read}_003.fastq.gz
	done

# plate 8

	mkdir -p "$DATADIR/links/plate8"
	cd "$DATADIR/links/plate8"

	for read in {1,2}; do
		ln -s ../../Rawdata/pool8_AAGAGCCA_AAGAGCCA/pool8_AAGAGCCA_AAGAGCCA_CKDL190141872-1A-AAGAGCCA-AAGAGCCA_H3VMVBBXX_L7_${read}.fq.gz plate8_R${read}_001.fastq.gz
		ln -s ../../Rawdata/pool8_CCTATGGT_CCTATGGT/pool8_CCTATGGT_CCTATGGT_CKDL190141872-1A-CCTATGGT-CCTATGGT_H3VMVBBXX_L7_${read}.fq.gz plate8_R${read}_002.fastq.gz
		ln -s ../../Rawdata/pool8_CGAAGAAC_CGAAGAAC/pool8_CGAAGAAC_CGAAGAAC_CKDL190141872-1A-CGAAGAAC-CGAAGAAC_H3VMVBBXX_L7_${read}.fq.gz plate8_R${read}_003.fastq.gz
	done
