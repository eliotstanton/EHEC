#!/bin/bash

mkdir lineage2

# Sakai vs. EC4115:
scripts/Synteny.pl \
	-o lineage2/Sakai/EC4115 \
	FASTA/Sakai.fasta,FASTA/EC4115.fasta \
	features/Sakai.features,features/EC4115.features

# Sakai vs. FRIK2455:
scripts/Synteny.pl \
	-o lineage2/Sakai/FRIK2455 \
	FASTA/Sakai.fasta,FASTA/FRIK2455.fasta \
	features/Sakai.features,features/FRIK2455.features

# FRIK2455 vs. EC4115:
scripts/Synteny.pl \
	-o lineage2/FRIK2455/EC4115 \
	FASTA/FRIK2455.fasta,FASTA/EC4115.fasta \
	features/FRIK2455.features,features/EC4115.features
