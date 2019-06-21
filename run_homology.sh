#!/bin/bash

# Examine homolohy present between ΦFRIK804-7 and ΦFRIK804-15 in FRIK804:
#scripts/HomologyLite.pl \
#	-m 100 \
#	-o analysis/homology/Φ804-7-Φ804-15 \
#	FASTA/FRIK804.Φ804-7.fasta,FASTA/FRIK804.Φ804-15.fasta

#scripts/WriteORFs.pl \
#	features/FRIK804.Φ804-7.features\
#	analysis/homology/Φ804-7-Φ804-15/Φ804-7.ORFs.svg

#scripts/WriteORFs.pl \
#	features/FRIK804.Φ804-15.features\
#	analysis/homology/Φ804-7-Φ804-15/Φ804-15.ORFs.svg

# Examine homology present between ΦFRIK804-9 and ΦFRIK804-10 in FRIK804:
#scripts/HomologyLite.pl \
#	-m 100 \
#	-o analysis/homology/FRIK804-phi9-FRIK804-phi10 \
#	FASTA/FRIK804.Φ804-9.fasta,FASTA/FRIK804.Φ804-10.fasta

#scripts/WriteORFs.pl \
#	features/FRIK804.Φ804-9.features\
#	analysis/homology/Φ804-9-Φ804-10/Φ804-9.ORFs.svg

#scripts/WriteORFs.pl \
#	features/FRIK804.Φ804-10.features\
#	analysis/homology/Φ804-9-Φ804-10/Φ804-10.ORFs.svg

# Examine homology present within FRIK804 at 100 bp:
scripts/HomologyAnalyzer.pl \
        -m 100 \
        -o analysis/homology/FRIK804.100 \
        FASTA/FRIK804.fasta \
        features/FRIK804.features

# Examine homology present within MG1655 at 100 bp:
scripts/HomologyAnalyzer.pl \
        -m 100 \
        -o analysis/homology/MG1655.100 \
        FASTA/MG1655.fasta \
        features/MG1655.features
#rm analysis/homology/MG1655.100/*.hash

exit

# Examine homology present within FRIK804 at 500 bp:
scripts/HomologyAnalyzer.pl \
	-m 500 \
	-o analysis/homology/FRIK804.500 \
	FASTA/FRIK804.fasta \
	features/FRIK804.features
rm analysis/homology/FRIK804.500/*.hash

# Examine homology present within Sakai at 100 bp:
scripts/HomologyAnalyzer.pl \
        -m 100 \
        -o analysis/homology/Sakai.100 \
        FASTA/Sakai.fasta \
        features/Sakai.features
rm analysis/homology/Sakai.100/*.hash

# Examine homology present within Sakai at 500 bp:
scripts/HomologyAnalyzer.pl \
        -m 500 \
        -o analysis/homology/Sakai.500 \
        FASTA/Sakai.fasta \
        features/Sakai.features
rm analysis/homology/Sakai.500/*.hash

# Examine homology present within MG1655 at 100 bp:
scripts/HomologyAnalyzer.pl \
        -m 100 \
        -o analysis/homology/MG1655.100 \
        FASTA/MG1655.fasta \
        features/MG1655.features
rm analysis/homology/MG1655.100/*.hash

# Examine homology present within MG1655 at 500 bp:
scripts/HomologyAnalyzer.pl \
        -m 500 \
        -o analysis/homology/MG1655.500 \
        FASTA/MG1655.fasta \
        features/MG1655.features
rm analysis/homology/MG1655.500/*.hash
