#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;
#use Homology;

# ---------------------------------------------------------------------------- #

# File name:		AlignReads.pl
# Date created:		21 November, 2018
# Last modified:	04 December, 2018
# Created by:		Eliot Stanton

# Description:		This is a script for aligning reads to a FASTA sequence or
#					sequences using bowtie.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;
my @array_SAM;
my @array_align;
my %hash_FASTA;
my %hash_coverage;

# Import variables from the command-line provided by the user:
getopts('c:o:s:t:u:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_FASTQ		= $ARGV[1];
my $var_scale		= $hash_variables{c} || 10;
my $var_min			= $hash_variables{m} || 500;
my $var_nmer		= $hash_variables{n} || 20;
my $var_output		= $hash_variables{o};
my $var_prefix		= $hash_variables{s} || "default";
my $var_start		= $hash_variables{t} || 1;
my $var_end			= $hash_variables{u};
my $file_SAM		= "$var_output/$var_prefix.sam";
my $file_aligned	= "$var_output/$var_prefix.aligned.fastq";
my $var_reads		= 0;
my $var_length_all	= 0;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nAlignReads.pl [FASTA] [FASTQ1],[FASTQ2]

	-c Scaling factor (default: 10)
	-o Output directory (required)
	-s Output files prefix (default: default)
	-t Start coordinate (default: 1)
	-u End coordinate (default: end of sequence)\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output or $var_FASTA is missing, print $var_help and stop script:
unless ( $var_FASTA && $hash_variables{o} ) {

	print "$var_help\n";

	print "\tOutput directory required!\n" unless $hash_variables{o};

	print "\tFASTA file required!\n" unless $var_FASTA;

	print "\tFASTQ file(s) required!\n" unless $var_FASTQ;

	exit;

}

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;
my @array_FASTQ		= split /\,/, $var_FASTQ;

# Exit if coordinates are specified and there is more than one FASTA file:
if ( scalar @array_FASTA > 1 && $var_start && $var_end ) {

	print "$var_help\n";

	print "Start and end coordinates can only be used with one FASTA file\n";

	exit;

}

# Check if files in @array_FASTA and @array_FASTQ are present:
@array_FASTA		= @{ General::CheckFiles ( \@array_FASTA, \%hash_variables ) };
@array_FASTQ		= @{ General::CheckFiles ( \@array_FASTQ, \%hash_variables ) };

# Merge FASTA files back into a comma-separated list:
$var_FASTA	= join ",", @array_FASTA;

# Check if output directory is present, if it isn't, create it:
General::CheckDirectory ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import each FASTA sequence to @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# Store each FASTA sequence in %hash_FASTA:
for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

	my $var_reference 	= $array_FASTA[$i][1];
	my $var_length		= length $array_FASTA[$i][2];

	# Add $var_length to $var_length_all:
	$var_length_all		+= $var_length;

	$var_reference		= ( split /[>,\s]/, $var_reference)[1];

	$hash_FASTA { $var_reference }	= $var_length;

}

# Store default length for $var_end if required:
$var_end	= $var_length_all unless $var_end;

# ---------------------------------------------------------------------------- #

# Create Bowtie index unless it already exists:
unless ( -e "$var_output/$var_prefix.1.ebwt" ) {

	print "bowtie-build \\\n\t$var_FASTA \\\n\t$var_output/$var_prefix \\\n";
	print "\t--threads 2 \\\n\t--quiet";
	system ( "bowtie-build $var_FASTA $var_output/$var_prefix --threads 2 --quiet" );

}

# ---------------------------------------------------------------------------- #

# Run bowtie:
print "bowtie \\\n\t$var_output/$var_prefix \\\n\t-1 $array_FASTQ[0] \\\n\t";
print "-2 $array_FASTQ[1] \\\n\t-n 0 \\\n\t-I 500 \\\n\t-X 1000 \\\n\t";
print "-k 1 \\\n\t-S $file_SAM \\\n\t --al $file_aligned \\\n\t-p 2\n";
system ( "bowtie $var_output/$var_prefix -1 $array_FASTQ[0] -2 $array_FASTQ[1] -n 0 -I 500 -X 1000 -k 1 -S $file_SAM --al $file_aligned -p 2" );

# ---------------------------------------------------------------------------- #

# Open the $file_SAM:
open ( my $file_read, "<", "$file_SAM" ) or die "Unable to open $file_SAM!\n";

my $var_counter	= 3;

# Iterate through each line:
while (<$file_read>) {

	if ( $var_counter ) {

		$var_counter--;

		next;

	}

	my $var_line	= $_;

	chomp $var_line;

	# Split $var_line into an array:
	my @array_SAM	= split " ", $var_line;

	# Define variables holding information needed:
	my $var_ref		= $array_SAM[2];
	my $var_loc0	= $array_SAM[3];
	my $var_loc1	= $array_SAM[7];
	my $var_length	= $array_SAM[8];

	# Move on if there is no alignment or it's the second read:
	next unless $var_loc0;
	next unless $var_length > 0;

	# Move on if region is outside of $var_start and $var_end:
	next unless $var_loc0 >= $var_start;
	next unless $var_loc1 <= $var_end;

	# Iterate $var_reads by one:
	$var_reads++;

	my @array_temp	= ( $var_loc0, $var_loc1, $var_length );

	push( @{ $hash_coverage { $var_ref } }, \@array_temp); 

}

# Close $file_SAM:
close $file_read or die "Unable to close $file_SAM!\n";

# ---------------------------------------------------------------------------- #

# Create coverage information:
foreach my $var_ref ( keys %hash_coverage ) {

	# Open array:
	my @array_alignments	= @{$hash_coverage{$var_ref}};

	# Open FASTA sequence:
	my $var_length			= $hash_coverage{$var_ref};

	# Sort alignments:
	@array_alignments		= sort { $a -> [0] <=> $b -> [0] } @array_alignments;

	# Create array and populate with zeros:
	my @array_coverage		= (0)x$var_length;

	# Iterate through alignments and create coverage map:
	for ( my $i = 0; $i < scalar @array_alignments; $i++ ) {

		my $var_loc0	= $array_alignments[$i][0];
		my $var_loc1	= $array_alignments[$i][1];

		# Add coverage for each location:
		for ( my $j = $var_loc0; $j <= $var_loc1; $j++ ) { $array_coverage[$j]++; }

	}

	# Store @array_coverage in %hash_coverage:
	$hash_coverage { $var_ref }	= \@array_coverage;

}

# ---------------------------------------------------------------------------- #

# Iterate through %hash_coverage trimming and scaling:
foreach my $var_ref ( keys %hash_coverage ) {

	# Define array holding original coverage data:
	my @array_temp	= @{$hash_coverage{$var_ref}}[$var_start..$var_end];

	# Define location to export data to:
	my $file_out	= "$var_output/$var_prefix.$var_ref.out";

	# Define a second temporary array to hold scaled down sequence data:
	my @array_temp2;

	# Add file name to @array_temp2:
	push @array_temp2, "# $file_out";

	# Iterate through first temporary array by the number of elements define by
	# $var_scale:
	for ( my $i = 0; $i < scalar @array_temp; $i += $var_scale ) {

		# Define the value of the element at hand:
		my $var_value	= $array_temp[$i];

		# Push $var_value to second temporary array:
		push @array_temp2, $var_value;

	}

	# Print output file name to user:
	print "$file_out\n";

	General::ArrayToFile ( \@array_temp2, $file_out );

}

# ---------------------------------------------------------------------------- #

# Evaluate read coverage and report:
my $var_coverage	= $var_reads/$var_length_all;
$var_coverage		= sprintf ( "%.2f", $var_coverage );

# Print number of reads, total length of sequence, and total read coverage:
print "Total aligned reads:   $var_reads reads\n";
print "Total sequence length: $var_length_all bp\n";
print "Total coverage:        $var_coverage reads/bp\n\n";

# ---------------------------------------------------------------------------- #
