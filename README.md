# EHEC genomic synteny repeat analysis and visualization

These scripts are associated with publications:
* Chronological set of E. coli O157:H7 bovine strains establishes a role for repeat sequences and mobile genetic elements in genome diversification (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06943-x)
* The impact of mobile genetic elements as drivers of genome diversification in bovine Escherichia coli O157:H7 (https://www.proquest.com/openview/75255941af6de228d499db701e4b5018/1?pq-origsite=gscholar&cbl=2026366&diss=y)

Bash scripts for generating information and figures in publications are located in the open (detailed at end). FASTA, FASTQ, GenBank, GFF, XML (omptical maps), and genome feature files needed to run these are all included in separate directories. Main scripts are located in /scripts directory.

Examples of the figures created with these scripts:  
<img src="https://github.com/eliotstanton/EHEC/blob/master/Fig1.large.png" width="480">  
Output of Synteny.pl (additional editing using Adobe Illustrator). Comparison of FRIK804 and Sakai chromosomes. Alignment of the FRIK804 (outer) and Sakai (inner) chromosome found disruption of synteny by large-scale structural alterations.  
  
  
<img src="https://github.com/eliotstanton/EHEC/blob/master/Fig6.FRIK804.png" width="480">
<img src="https://github.com/eliotstanton/EHEC/blob/master/Fig6.MG1655.png" width="480">  
Outputs of HomologyAnalyzer.pl (additional editing using Adobe Illustrator). Locations of direct and inverted contiguous repeats ≥75 bp in length in the chromosomes of FRIK804 and E. coli strain MG1655.  
  
  
<img src="https://github.com/eliotstanton/EHEC/blob/master/Fig3A.small.png" width="480">  
Output of Optical.pl (additional editing using Adobe Illustrator). NcoI restriction site maps of the chromosomes from FRIK804 (outer), FRIK1275 (middle), and FRIK1625 (inner).  
  
  
<img src="https://github.com/eliotstanton/EHEC/blob/master/Fig4.small.png" width="960">  
Output of AlignReads.pl, GrabORFs.pl, and WriteORFs.pl (additional editing using Adobe Illustrator). Detection of the boundaries of inter-prophage deletion (indel-4) in Φ804–9/Φ804–10 present in FRIK1275 and FRIK1625.  

--------------------------------------------------------------------------------

# Main scripts

## HomologyAnalyzer.pl
  Master script used for caluclating and visualising homology within a circular bacterial genome  

  HomologyAnalyzer.pl [FASTA] [Features]  
    -d Prohibit direct links from being drawn  
    -i Prohibit inverted links from being drawn  
    -m Minimum repeat length (default: 100)  
    -n nmer length for homology (default: 20 bp)  
    -o OUTPUT DIRECTORY (required)  
    -s Output files prefix (default: default)  
    
    Feature file format:  
    seqID	feature_type	feature_name	start	end  
    example:  
    0	prophage	prophage0	103894	163432  


## Synteny.pl
  Calculating and visualising related strains using Mauve and Circos  
  
  Synteny.pl [OPTIONS] [FASTA] [Features]  
    -m Minimum length for region alignment (default: 100)  
    -o Output directory (required)  
    -p Force progressiveMauve to run  
    -s Output files prefix (default: default)  
    
    Feature file format:  
    seqID	feature_type	feature_name	start	end  
    example:  
    0	prophage	prophage0	103894	163432  


## AlignReads.pl
  Used to align reads to a FASTA sequence or sequences using Bowtie.  

  AlignReads.pl [FASTA] [FASTQ1],[FASTQ2]  
    -c Scaling factor (default: 10)  
    -o Output directory (required)  
    -s Output files prefix (default: default)  
    -t Start coordinate (default: 1)  
    -u End coordinate (default: end of sequence)  

## GrabORFs.pl
  Used to pull ORF locations from a GFF file  

TODO: Document fully  


## HomologyLite.pl
  Used for determining homology shared between one or more short FASTA sequences  
  
  HomologyLite.pl [OPTIONS] [FASTA]  
    -m Minimum repeat length (default: 100)  
    -n nmer length for homology (default: 20 bp)  
    -o OUTPUT DIRECTORY (required)  
    -s Output files prefix (default: default)  


## ORFs.pl 
  Used for converting GenBank data over to SVG for use in figures  
  
  ORFs.pl [OPTIONS]  
    -a Mauve alignment file  
    -e Stem for file name  
    -g GenBank file (required)  
    -o Output directory (required)  
    -s Start location  
    -t Stop location  


## Optical.pl
  Used for visualising optical mapping  
    
  optical.pl [OPTIONS] [MAP1],[MAP2]  


## WriteORFs.pl
  Converts genomic features coordinates into SVF formatting for use in figures  
  
TODO: Document fully  


--------------------------------------------------------------------------------

# BASH scripts:
## run.sh
TODO: Document fully  


## run_ORFs.sh
TODO: Document fully  


## run_align.sh
TODO: Document fully  


## run_blast.sh
TODO: Document fully  


## run_homology.sh
TODO: Document fully  


## run_mapping.sh
TODO: Document fully  


## run_synteny.sh
TODO: Document fully  

## lineage2.sh
TODO: Document fully  
