#!/bin/bash

# Construct map of synteny shared between Sakai and FRIK804
scripts/Synteny.pl \
	-o analysis/synteny/FRIK804-Sakai/chromosome \
	FASTA/FRIK804.fasta,FASTA/Sakai.fasta \
	features/FRIK804.features,features/Sakai.features

# Construct map of synteny shared between Mu-like prophage:
scripts/Synteny.pl \
	-o analysis/synteny/FRIK804-Sakai/Mu \
	FASTA/FRIK804.Φ804-3.fasta,FASTA/Sakai.Sp18.fasta \
	features/FRIK804.Φ804-3.features,features/Sakai.Sp18.features

# Construct map of synteny shared between PLE804-1 and SpLE1:
scripts/Synteny.pl \
	-o analysis/synteny/FRIK804-Sakai/indel-1 \
	FASTA/FRIK804.indel-1.fasta,FASTA/Sakai.SpLE1.fasta \
	features/FRIK804.indel-1.features,features/Sakai.SpLE1.features

# Construct map of synteny shared between stx2 prophage in FRIK804 and Sakai:
scripts/Synteny.pl \
	-o analysis/synteny/FRIK804-Sakai/stx2 \
	FASTA/FRIK804.Φ804-6.fasta,FASTA/Sakai.Sp5.fasta \
	features/FRIK804.Φ804-6.features,features/Sakai.Sp5.features

# Construct map of synteny shared between stx2 prophage in FRIK804, EDL933, and
# Sakai:
scripts/Synteny.pl \
	-o analysis/synteny/FRIK804-EDL933-Sakai/stx2 \
	FASTA/FRIK804.Φ804-6.fasta,FASTA/EDL933.933W.fasta,FASTA/Sakai.Sp5.fasta \
	features/FRIK804.Φ804-6.features,features/EDL933.933W.features,features/Sakai.Sp5.features
