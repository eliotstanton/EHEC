#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;
use Homology;

# ---------------------------------------------------------------------------- #

# File name:		HomologyLite.pl
# Date created:		21 November, 2018
# Last modified:	30 November, 2018
# Created by:		Eliot Stanton

# Description:		This is a script for determining homology shared between
#					one or more short FASTA sequences.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('m:n:o:s:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_min			= $hash_variables{m} || 100;
my $var_nmer		= $hash_variables{n} || 20;
my $var_output		= $hash_variables{o};
my $var_prefix		= $hash_variables{s} || "default";
my $var_scale		= 10;
my $var_height		= 100;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nHomologyLite.pl [OPTIONS] [FASTA]
	-m Minimum repeat length (default: 100)
	-n nmer length for homology (default: 20 bp)
	-o OUTPUT DIRECTORY (required)
	-s Output files prefix (default: default)\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output or $var_FASTA is missing, print $var_help and stop script:
unless ( $var_FASTA && $hash_variables{o} ) {

	print "$var_help\n";

	print "\tOutput directory required!\n" unless $hash_variables{o};

	print "\tFASTA file(s) required!\n" unless $var_FASTA;

	exit;

}

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Check if files in @array_FASTA are present:
@array_FASTA		= @{ General::CheckFiles ( \@array_FASTA, \%hash_variables ) };

# Check if output directory is present, if it isn't, create it:
General::CheckDirectory ( \%hash_variables );

# Import data from @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# ---------------------------------------------------------------------------- #

# Generate a hash of each nmer present in all the FASTA sequences:
my $hash_nmer		= Homology::Nmer ( \%hash_variables, \@array_FASTA, $var_nmer );

# --------------------------------------

# Create formatted direct link data:
Homology::DirectLinks ( \%hash_variables, $hash_nmer, $var_nmer );

# Create formatted inverted link data:
Homology::InvertedLinks ( \%hash_variables, $hash_nmer, $var_nmer );

# Merge adjacent direct links:
my @array_direct	= @{Homology::MergeLinks ( \%hash_variables, "direct", $var_min )};

# Merge adjacent inverted links:
my @array_inverted	= @{Homology::MergeLinks ( \%hash_variables, "inverted", $var_min )};

# ---------------------------------------------------------------------------- #

# If entries are present in @array_direct, print them to user and to file:
if ( @array_direct ) {

	print "Direct Links:\n";

	for ( my $i = 0; $i < scalar @array_direct; $i++ ) {

		my @array_temp	= split " ", $array_direct[$i];

		my $var_seq0	= $array_temp[0];
		my $var_loc0	= $array_temp[1];
		my $var_loc1	= $array_temp[2];
		my $var_seq1	= $array_temp[3];
		my $var_loc2	= $array_temp[4];
		my $var_loc3	= $array_temp[5];

		my $var_length	= $var_loc1 - $var_loc0 + 1;

		print "$i: $var_seq0: $var_loc0 $var_loc1\t->\t$var_seq1: $var_loc2 $var_loc3\t$var_length bp\n";

		# Define length of region:
		$var_length		= $var_loc1 - $var_loc0;

		# Scale variables:
		$var_loc0		/= $var_scale;
		$var_loc2		/= $var_scale;
		$var_length		/= $var_scale;

		# Create SVG rectangles:
		my $var_string0	= General::RectangleSVG ( $var_loc0, 0, $var_length, 100, "(0,0,255)", 0 );
		my $var_string1	= General::RectangleSVG ( $var_loc2, 0, $var_length, 100, "(0,0,255)", 0 );

		push @{$array_FASTA[$var_seq0][3]}, $var_string0;
		push @{$array_FASTA[$var_seq1][3]}, $var_string1;

	}

}

# If entries are present in @array_inverted, print them to user and to file:
if ( @array_inverted ) {

	print "Inverted Links:\n";

	for ( my $i = 0; $i < scalar @array_inverted; $i++ ) {

		my @array_temp	= split " ", $array_inverted[$i];

		my $var_seq0	= $array_temp[0];
		my $var_loc0	= $array_temp[1];
		my $var_loc1	= $array_temp[2];
		my $var_seq1	= $array_temp[3];
		my $var_loc2	= $array_temp[4];
		my $var_loc3	= $array_temp[5];

		my $var_length	= $var_loc1 - $var_loc0 + 1;

		print "$i: $var_seq0: $var_loc0 $var_loc1\t->\t$var_seq1: $var_loc2 $var_loc3\t$var_length bp\n";

		# Define length of region:
		$var_length		= $var_loc1 - $var_loc0;

		# Scale variables:
		$var_loc0		/= $var_scale;
		$var_loc2		/= $var_scale;
		$var_length		/= $var_scale;

		# Create SVG rectangles:
		my $var_string0	= General::RectangleSVG ( $var_loc0, 0, $var_length, 100, "(255,0,0)", 0 );
		my $var_string1	= General::RectangleSVG ( $var_loc2, 0, $var_length, 100, "(255,0,0)", 0 );

		push @{$array_FASTA[$var_seq0][4]}, $var_string0;
		push @{$array_FASTA[$var_seq1][4]}, $var_string1;

	}

}

# ---------------------------------------------------------------------------- #

# Iterate through @array_FASTA and print SVG data to file:
for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

	# Define file path of FASTA sequence:
	my $var_file	= $array_FASTA[$i][0];

	# Define length of FASTA sequence:
	my $var_length	= length $array_FASTA[$i][2];

	# Define path for file to be printed:
	my $file_out	= "$var_output/$var_prefix.$i.svg";

	print "$i: $var_file -> $file_out\n";

	# Define an array to hold information to be printed:
	my @array_out;

	# Define a temporary array holding direct and inverted repeats:
	my @array_direct 	= @{$array_FASTA[$i][3]} if $array_FASTA[$i][3];
	my @array_inverted	= @{$array_FASTA[$i][4]} if $array_FASTA[$i][4];

	# Define the width of image using length of FASTA sequence:
	my $var_width		= $var_length / $var_scale;

	# Add header information:
	push @array_out, "<svg width=\"$var_width\" height=\"$var_height\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

	# Add background:
	push @array_out, "<g id=\"bg\">";
	my $var_string	= "<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
	$var_string		.= "height=\"$var_height\" style=\"fill:rgb(255,255,255)\"/>";
	push @array_out, $var_string;
	push @array_out, "</g>";

	# ----------------------------------

	for my $var_string ( @array_direct ) {

		push @array_out, $var_string;

	}

	for my $var_string ( @array_inverted ) {

		push @array_out, $var_string;

	}

	push @array_out, "</svg>";

	# ----------------------------------

	General::ArrayToFile ( \@array_out, $file_out );

}

# ---------------------------------------------------------------------------- #
