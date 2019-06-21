#!/bin/bash

# Create BLAST databases for Sakai, FRIK804, and MG1655:
makeblastdb -in FASTA/Sakai.fasta -dbtype nucl -title Sakai -out BLAST/Sakai
makeblastdb -in FASTA/FRIK804.fasta -dbtype nucl -title FRIK804 -out BLAST/FRIK804
makeblastdb -in FASTA/MG1655.fasta -dbtype nucl -title MG1655 -out BLAST/MG1655

printf "oriC:\n"

# BLAST locations of oriC:
printf "Sakai\n"
blastn -db BLAST/Sakai -query FASTA/oriC.fasta -outfmt "6 sstart send length"
printf "FRIK804\n"
blastn -db BLAST/FRIK804 -query FASTA/oriC.fasta -outfmt "6 sstart send length"
printf "MG1655\n"
blastn -db BLAST/MG1655 -query FASTA/oriC.fasta -outfmt "6 sstart send length"

printf "dif:\n"
# BLAST locations of dif:
printf "Sakai\n"
blastn -db BLAST/Sakai -query FASTA/dif.fasta -outfmt "6 sstart send length"
printf "FRIK804\n"
blastn -db BLAST/FRIK804 -query FASTA/dif.fasta -outfmt "6 sstart send length"
printf "MG1655\n"
blastn -db BLAST/MG1655 -query FASTA/dif.fasta -outfmt "6 sstart send length"

exit

printf "ISEc8:\n"

# BLAST locations of ISEc8:
printf "Sakai\n"
blastn -db BLAST/Sakai -query FASTA/ISEc8.fasta -outfmt "6 sstart send length"
printf "FRIK804\n"
blastn -db BLAST/FRIK804 -query FASTA/ISEc8.fasta -outfmt "6 sstart send length"
printf "MG1655\n"
blastn -db BLAST/MG1655 -query FASTA/ISEc8.fasta -outfmt "6 sstart send length"

printf "IS629:\n"
# BLAST locations of IS629:
printf "Sakai\n"
blastn -db BLAST/Sakai -query FASTA/IS629.fasta -outfmt "6 sstart send length"
printf "FRIK804\n"
blastn -db BLAST/FRIK804 -query FASTA/IS629.fasta -outfmt "6 sstart send length"
printf "MG1655\n"
blastn -db BLAST/MG1655 -query FASTA/IS629.fasta -outfmt "6 sstart send length"

exit

# BLAST locations of 5S rRNA:
blastn -db BLAST/Sakai -query FASTA/5S.fasta -outfmt "6 sstart send length"
blastn -db BLAST/FRIK804 -query FASTA/5S.fasta -outfmt "6 sstart send length"
blastn -db BLAST/MG1655 -query FASTA/5S.fasta -outfmt "6 sstart send length"

# BLAST locations of 16S rRNA:
blastn -db BLAST/Sakai -query FASTA/16S.fasta -outfmt "6 sstart send length"
blastn -db BLAST/FRIK804 -query FASTA/16S.fasta -outfmt "6 sstart send length"
blastn -db BLAST/MG1655 -query FASTA/16S.fasta -outfmt "6 sstart send length"

# BLAST locations of 23S rRNA:
blastn -db BLAST/Sakai -query FASTA/23S.fasta -outfmt "6 sstart send length"
blastn -db BLAST/FRIK804 -query FASTA/23S.fasta -outfmt "6 sstart send length"
blastn -db BLAST/MG1655 -query FASTA/23S.fasta -outfmt "6 sstart send length"
