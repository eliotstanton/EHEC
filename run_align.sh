#!/bin/bash

# Align reads from FRIK804 to prophage sequences:
scripts/AlignReads.pl \
	-o analysis/align/phi9phi10 \
	-s FRIK804 \
	-t 2097886 \
	-u 2187895 \
	FASTA/FRIK804.fasta \
/media/eliot/09E628311F540FB8/FASTQ/FRIK804_CGATGT_L001_R1_001.fastq,/media/eliot/09E628311F540FB8/FASTQ/FRIK804_CGATGT_L001_R2_001.fastq
#	FASTQ/FRIK804.aligned_1.fastq,FASTQ/FRIK804.aligned_2.fastq

rm analysis/align/phi9phi10/*.ebwt
rm analysis/align/phi9phi10/*.sam

# Align reads from FRIK1275 to prophage sequences:
scripts/AlignReads.pl \
	-o analysis/align/phi9phi10 \
	-s FRIK1275 \
	-t 2097886 \
	-u 2187895 \
	FASTA/FRIK804.fasta \
/media/eliot/09E628311F540FB8/FASTQ/FRIK1275_TGACCA_L001_R1_001.fastq,/media/eliot/09E628311F540FB8/FASTQ/FRIK1275_TGACCA_L001_R2_001.fastq
#	FASTQ/FRIK1275.aligned_1.fastq,FASTQ/FRIK1275.aligned_2.fastq

rm analysis/align/phi9phi10/*.ebwt
rm analysis/align/phi9phi10/*.sam

# Align reads from FRIK1625 to prophage sequences:
scripts/AlignReads.pl \
	-o analysis/align/phi9phi10 \
	-s FRIK1625 \
	-t 2097886 \
	-u 2187895 \
	FASTA/FRIK804.fasta \
/media/eliot/09E628311F540FB8/FASTQ/FRIK1625_ACAGTG_L001_R1_001.fastq,/media/eliot/09E628311F540FB8/FASTQ/FRIK1625_ACAGTG_L001_R2_001.fastq
#	FASTQ/FRIK1625.aligned_1.fastq,FASTQ/FRIK1625.aligned_2.fastq

rm analysis/align/phi9phi10/*.ebwt
rm analysis/align/phi9phi10/*.sam
