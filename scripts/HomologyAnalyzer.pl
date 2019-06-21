#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;
use Homology;
use Circos;

# ---------------------------------------------------------------------------- #

# File name:		HomologyAnalyzer.pl
# Date created:		27 October, 2018
# Last modified:	02 December, 2018
# Created by:		Eliot Stanton

# Description:		This is a master script for calculating and visualising
#					homology within a circular bacterial genome.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('dim:n:o:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_features	= $ARGV[1];
my $var_min			= $hash_variables{m} || 100;
my $var_nmer		= $hash_variables{n} || 20;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nHomologyAnalyzer.pl [FASTA] [Features]
	-d Prohibit direct links from being drawn
	-i Prohibit inverted links from being drawn
	-m Minimum repeat length (default: 100)
	-n nmer length for homology (default: 20 bp)
	-o OUTPUT DIRECTORY (required)
	-s Output files prefix (default: default)

	Feature file format:
	seqID	feature_type	feature_name	start	end

	example:
	0	prophage	prophage0	103894	163432\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output is missing, print $var_help and stop script:
unless ( $hash_variables{o} ) {

	print "$var_help\n";

	print "\tOutput directory required!\n" and exit;

}

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# If $var_FASTA or $var_features are missing, print $var_help and stop script:
unless ( scalar @ARGV == 2 ) {

	print "$var_help\n";
	
	print "\tFASTA and feature files are required!\n" and exit;

}

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Split $var_features into @array_features:
my @array_features	= split /\,/, $var_features;

if ( scalar @array_FASTA > 1 || scalar @array_features > 1 ) {

	print "$var_help\n";
	
	print "\tOnly one FASTA and Feature file allowed!\n" and exit;

}

# ---------------------------------------------------------------------------- #

# Check if Circos will run:
General::CheckCircos ( \%hash_variables );

# Check if files in @array_FASTA are present:
@array_FASTA		= @{ General::CheckFiles ( \@array_FASTA, \%hash_variables ) };

# Check if files in @array_features are present:
@array_features		= @{ General::CheckFiles ( \@array_features, \%hash_variables ) };

# Check if output directory is present, if it isn't, create it:
General::CheckDirectory ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import data from @array_FASTA:
@array_FASTA		= @{ General::ImportFASTA ( \@array_FASTA, $var_help ) };

# Import data from @array_features:
@array_features		= @{ General::FileToArray ( $array_features[0] ) };

# Catch and correct genome features that have coordinates in the wrong order:
for ( my $i = 0; $i < scalar @array_features; $i++ ) {

	my $var_loc0	= $array_features[$i][3];
	my $var_loc1	= $array_features[$i][4];

	if ( $var_loc1 < $var_loc0 ) {

		$array_features[$i][3]	= $var_loc1;
		$array_features[$i][4]	= $var_loc0;

	}

}

# ---------------------------------------------------------------------------- #

# Generate a hash of each nmer present in all the FASTA sequences:
my $hash_nmer		= Homology::Nmer ( \%hash_variables, \@array_FASTA, $var_nmer );

# --------------------------------------

# Create formatted direct link data:
Homology::DirectLinks ( \%hash_variables, $hash_nmer, $var_nmer );

# Create formatted inverted link data:
Homology::InvertedLinks ( \%hash_variables, $hash_nmer, $var_nmer );

# Merge adjacent direct links:
Homology::MergeLinks (  \%hash_variables, "direct", $var_min );

# Merge adjacent inverted links:
Homology::MergeLinks ( \%hash_variables, "inverted", $var_min );

# --------------------------------------

# Calculate overlap between links and features:
my $array_direct	= Homology::Overlap ( \%hash_variables, \@array_features, "direct" );
my $array_inverted	= Homology::Overlap ( \%hash_variables, \@array_features, "inverted" );

Temp ( \%hash_variables, \@array_features, "Direct" );
Temp ( \%hash_variables, \@array_features, "Inverted" );

exit;

# ---------------------------------------------------------------------------- #

# Create Circos config file:
Homology::Config ( \%hash_variables );

# Create Circos links files:
Homology::Links ( $array_direct, $array_inverted, \@array_features, \%hash_variables );

# Create Circos Karyotype file:
Homology::Karyotype ( \@array_FASTA, \%hash_variables );

# Create Circos Ideogram file:
Circos::Ideogram ( \%hash_variables );

# Create Circos highlights features:
Homology::Features ( \@array_features, \@array_FASTA, \%hash_variables );

# Create ticks file:
Homology::Ticks ( \@array_FASTA, \%hash_variables );

# Create Circos Image file:
Circos::Image ( \%hash_variables );

# Run Circos:
Circos::Run ( \%hash_variables );

# ---------------------------------------------------------------------------- #

sub Temp {

	# Arguments:
	my ( $hash_in, $array_features, $var_type )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_features	= @$array_features;
	my %hash_in;
	my @array_out;
	my @array_flat;
	my %hash_flat;

	# Variables:
	my $file_in;
	my $var_total		= 0;
	my $var_flat		= 0;

	# ----------------------------------

	# Set $var_variable equal to one (direct) or minus one (inverted):
	if ( $var_type eq "Direct" ) {

		$file_in		= $hash_out { file_direct_merged };

	}

	if ( $var_type eq "Inverted" ) {

		$file_in		= $hash_out { file_inverted_merged };

	}

	# Import data from $file_in, containing links, to @array_in:
	my @array_in		= @{General::FileToArray ( $file_in )};

	# ----------------------------------

	# Sort features from smallest to largest:
	@array_features	= sort { $a -> [5] <=> $b -> [5] } @array_features;

	# ----------------------------------

	# TODO: Iterate through locations of repeats and identify overlapping genomic
	# features:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0		= $array_in[$i][1];
		my $var_loc1		= $array_in[$i][2];

		my $var_loc2		= $array_in[$i][4];
		my $var_loc3		= $array_in[$i][5];

		my $var_length		= $var_loc1 - $var_loc0 + 1;

		$var_total			+= $var_length;

		# Add regions to %hash_in for flattening:
		$hash_in { "$var_loc0 $var_loc1" }	= 0;
		$hash_in { "$var_loc2 $var_loc3" }	= 0;

	}

	# ----------------------------------

	# Iterate through %hash_in and store locations (without repeats) in
	# @array_out:
	foreach my $var_locations ( keys %hash_in ) {

		my @array_temp	= split " ", $var_locations;

		push @array_flat, \@array_temp;

	}

	# Sort repeats in @array_flat by location:
	@array_flat	= sort { $a -> [0] <=> $b -> [0] || $b -> [1] <=> $a -> [1] } @array_flat;

	# ----------------------------------

	# Iterate down through @array_flat merging overlapping regions:
	for ( my $i = 0; $i < scalar @array_flat; $i++ ) {

		my $var_loc0		= $array_flat[$i][0];
		my $var_loc1		= $array_flat[$i][1];

		# Iterate down through subsequent entries:
		for ( my $j = $i+1; $j < scalar @array_flat; $j++ ) {

			my $var_loc2		= $array_flat[$j][0];
			my $var_loc3		= $array_flat[$j][1];

			last if $var_loc2 > $var_loc1;

			# If region is entirely inside, remove it and start over:
			if ( $var_loc0 <= $var_loc2 && $var_loc1 >= $var_loc3 ) {

				splice @array_flat, $j, 1;

				$i = -1;

				last;

			}

			# If region overlaps, merge them and start over:
			elsif ( $var_loc2 <= $var_loc1 && $var_loc3 >= $var_loc1 ) {

				$array_flat[$i][1]	= $var_loc3;

				splice @array_flat, $j, 1;

				$i = -1;

				last;

			}

		}

	}

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_flat; $i++ ) {

		my $var_loc0	= $array_flat[$i][0];
		my $var_loc1	= $array_flat[$i][1];

		# Iterate through @array_features:
		for ( my $j = 0; $j < scalar @array_features; $j++ ) {

			# Define feature type and its locations:
			my $var_type	= $array_features[$j][1];
			my $var_loc2	= $array_features[$j][3];
			my $var_loc3	= $array_features[$j][4];

			# Handle repeats that encompass right edge of element:
			if ( $var_loc0 <= $var_loc2 && $var_loc1 > $var_loc2 && $var_loc1 < $var_loc3 ) {

				my @array_temp	= ( $var_loc0, $var_loc2 - 1 );

				push @array_flat, \@array_temp;

				$array_flat[$i][0]	= $var_loc2;
				$array_flat[$i][1]	= $var_loc1;
				$array_flat[$i][2]	= $var_type;

				last;

			}

			# Handle repeats that encompass left edge of element:
			if ( $var_loc0 > $var_loc2 && $var_loc0 < $var_loc3 && $var_loc1 >= $var_loc3 ) {

				my @array_temp	= ( $var_loc3 + 1, $var_loc1 );

				push @array_flat, \@array_temp;

				$array_flat[$i][0]	= $var_loc0;
				$array_flat[$i][1]	= $var_loc3;
				$array_flat[$i][2]	= $var_type;

				last;

			}

			# Handle repeats that encompass the element:
			if ( $var_loc0 <= $var_loc2 && $var_loc1 >= $var_loc3 ) {

				my @array_temp0	= ( $var_loc0, $var_loc2 - 1 );
				my @array_temp1	= ( $var_loc3 + 1, $var_loc1 );

				push @array_flat, \@array_temp0;
				push @array_flat, \@array_temp1;

				$array_flat[$i][0]	= $var_loc2;
				$array_flat[$i][1]	= $var_loc3;
				$array_flat[$i][2]	= $var_type;

				last;

			}

			# Handle repeats nested within the element:
			if ( $var_loc0 > $var_loc2 && $var_loc1 < $var_loc3 ) {

				$array_flat[$i][0]	= $var_loc0;
				$array_flat[$i][1]	= $var_loc1;
				$array_flat[$i][2]	= $var_type;

				last;

			}

		}

	}

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_flat; $i++ ) {

		my $var_loc0	= $array_flat[$i][0];
		my $var_loc1	= $array_flat[$i][1];
		my $var_type	= $array_flat[$i][2] || "";

		my $var_length	= $var_loc1 - $var_loc0 + 1;

		if ( $var_type ) {

			$hash_flat { $var_type }	+= $var_length;

		}

		else {

			$hash_flat { "backbone" }	+= $var_length;

		}

		$var_flat		+= $var_length;

	}

	# ----------------------------------

	print "$var_type repeats:\n";

	print "\tTotal: $var_total\n";

	print "\tFlat: $var_flat\n";

	# Print length data for each type of genomic feature:
	foreach my $var_type ( sort { $a cmp $b } keys %hash_flat ) {

		my $var_total		= $hash_flat { $var_type } || 0;

		print "\t - $var_type $var_total bp\n";

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

sub Temp2 {

	# Arguments:
	my ( $var_loc0, $var_loc1, $array_features )	= @_;

	# Data structures:
	my @array_features	= @$array_features;

	# Variables:

	# ----------------------------------

	# ----------------------------------

	# End subroutine:
	return;

}

