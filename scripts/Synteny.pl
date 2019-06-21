#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;
use Synteny;
use Synteny;
use Circos;

# ---------------------------------------------------------------------------- #

# File name:		Synteny.pl
# Date created:		31 October, 2018
# Last modified:	19 February, 2019
# Created by:		Eliot Stanton

# Description:		This is a master script for calculating and visualising
#					synteny between closely related strains using Mauve and
#					Circos.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('m:o:ps:', \%hash_variables);

# Define variables used by this script:
my $var_FASTA		= $ARGV[0];
my $var_features	= $ARGV[1];
my $var_pM			= $hash_variables{p};

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nSynteny.pl [OPTIONS] [FASTA] [Features]
	-m Minimum length for region alignment (default: 100)
	-o Output directory (required)
	-p Force progressiveMauve to run
	-s Output files prefix (default: default)

	Feature file format:
	seqID	feature_type	feature_name	start	end

	example:
	0	prophage	prophage0	103894	163432

\n";

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

	print "\tFASTA and Feature file(s) are required!\n" and exit;

}

# Split $var_FASTA into @array_FASTA:
my @array_FASTA		= split /\,/, $var_FASTA;

# Split $var_features into @array_features:
my @array_features	= split /\,/, $var_features;

# ---------------------------------------------------------------------------- #

# Check if ProgressiveMauve will run:
General::CheckProgressiveMauve ( \%hash_variables);

# Check if files in @array_FASTA are present:
@array_FASTA	= @{ General::CheckFiles ( \@array_FASTA, \%hash_variables ) };

# Check if files in @array_features are present:
@array_features	= @{ General::CheckFiles ( \@array_features, \%hash_variables ) };

# Check if output directory is present, if it isn't, create it:
General::CheckDirectory ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import data from @array_features:
@array_features	= @{ General::ImportFiles ( \@array_features ) };

# ---------------------------------------------------------------------------- #

# Run progressiveMauve if forced to do so:
Synteny::progressiveMauve ( \%hash_variables, \@array_FASTA ) if $var_pM;

# Process $file_backbone created by progressiveMauve and split backbone and
# insert alignments into separate arrays:
%hash_variables	= %{Synteny::FormatAlignment ( \%hash_variables, \@array_FASTA )};

# Organise backbone regions and split transposed backbone regions into a
# separate array:
%hash_variables	= %{Synteny::OrganiseBackbone ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Add bracketing reference numbers to inserts and transposed regions, this
# subroutine all returns @array_backbone rearranged by sequence:
%hash_variables	= %{Synteny::BracketRefs ( \%hash_variables )};

# Split insert regions sitting within junctions of backbone regions into a
# separate array:
%hash_variables	= %{Synteny::JunctionInserts ( \%hash_variables)};

# Check for transposed regions in @array_inserts:
#%hash_variables	= %{TransposedInserts ( \%hash_variables )};
%hash_variables	= %{Synteny::TransposedInserts ( \%hash_variables )};

# Organise insert regions:
#%hash_variables	= %{OrganiseInserts ( \%hash_variables )};
%hash_variables	= %{Synteny::OrganiseInserts ( \%hash_variables )};

# Merge regions all together:
%hash_variables	= %{Synteny::MergeRegions ( \%hash_variables )};

# Check that all regions are present in array_merged:
Synteny::CheckRegions ( \%hash_variables );

# Assign adjusted global locations to all regions:
%hash_variables	= %{ Synteny::GlobalLocations ( \%hash_variables ) };

# ---------------------------------------------------------------------------- #

# Calculate overlap between different genomic features and 
#conserved/non-conserved regions:
Synteny::Overlap ( \@array_features, \@array_FASTA, \%hash_variables );

# Format genomic features:
%hash_variables	= %{ Synteny::FormatFeatures ( \@array_features, \%hash_variables ) };

# Import FASTA files to @array_FASTA:
@array_FASTA	= @{General::ImportFASTA ( \@array_FASTA )};

# Write conserved and non-conserved sequences to file in FASTA format:
Synteny::WriteSeq ( \%hash_variables, \@array_FASTA );

# Store final alignment data as files:
Synteny::StoreFinal ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Create SVG data for syntenic regions:
%hash_variables	= %{RegionsSVG ( \%hash_variables )};

# Create SVG data for features:
%hash_variables	= %{FeaturesSVG ( \%hash_variables )};

# Create SVG data for ticks:
%hash_variables	= %{TicksSVG ( \%hash_variables )};

# Merge SVG data and print to file:
MergeSVG ( \%hash_variables );

# ---------------------------------------------------------------------------- #

sub RegionsSVG {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{ array_merged }};
	my @array_out;

	# Variables:
	my $var_height		= $hash_out{var_span};
	my $var_factor		= $hash_out{var_factor};
	my $var_y			= 100;
	my $var_increment	= 200;
	my $var_scalar		= scalar @array_merged;
	my $var_stroke		= "0,0,0";

	# ----------------------------------

	# Iterate through @array_merged:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_copy		= $array_merged[$i][$j][4];
			my $var_length		= $array_merged[$i][$j][6];
			my $var_inverted	= $array_merged[$i][$j][8];
			my $var_loc0		= $array_merged[$i][$j][11];
			my $var_loc1		= $array_merged[$i][$j][12];

			$var_y += $var_height if $var_inverted;

			next if $var_loc0 	== 0;

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			my $var_colour	= Circos::ColoursNew ( "chromosome", $var_copy, $var_scalar );

			my $var_string	= General::RectangleSVG ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

			push @{$array_out[$i]}, $var_string;

			$var_y -= $var_height if $var_inverted;

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	# Store reference to @array_out in %hash_out:
	$hash_out { array_regions }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #


sub FeaturesSVG { 

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_features	= @{$hash_out{ array_features }};
	my @array_out;

	# Variables:
	my $var_height		= $hash_out{var_span};
	my $var_height2		= $hash_out{var_span}/2;
	my $var_factor		= $hash_out{var_factor};
	my $var_y			= 100;
	my $var_increment	= 200;
	my $var_scalar		= scalar @array_features;
	my $var_stroke		= "0,0,0";

	# ----------------------------------

	# Iterate through @array_merged:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		next unless $array_features[$i];

		# Iterate down through features for each sequence:
		for ( my $j = 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			# Define start/end locations of genomic feature:
			my $var_loc0		= $array_features[$i][$j][0];
			my $var_loc1		= $array_features[$i][$j][1];

			# Define the type of genomic feature, inverted status and copy:
			my $var_type		= $array_features[$i][$j][2];
			my $var_inverted	= $array_features[$i][$j][3];
			my $var_copy		= $array_features[$i][$j][4];
			my $var_orient		= $array_features[$i][$j][5] || "0";

			# Calculate length of genomic feature:
			my $var_length		= $var_loc1 - $var_loc0 + 1;

			# Change $var_height if feature is inverted:
			$var_y += $var_height if $var_inverted;

			# Change height if feature is reverse oriented:
			$var_y += $var_height2 if $var_orient eq "-";

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			# Calculate genomic feature colour:
			my $var_colour	= Circos::ColoursNew ( $var_type, $var_copy, $var_scalar );

			if ( $var_type eq "gene" ) {

				# Create string holding SVG formatted region:
				my $var_string	= General::RectangleSVG ( $var_loc0, $var_y, $var_length, $var_height2, $var_colour, $var_stroke );

				# Store $var_string in @array_out:
				push @{$array_out[$i]}, $var_string;

			}

			else {

				# Create string holding SVG formatted region:
				my $var_string	= General::RectangleSVG ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

				# Store $var_string in @array_out:
				push @{$array_out[$i]}, $var_string;

			}

			# Decrease $var_height if handling an inverted region:
			$var_y -= $var_height if $var_inverted;

			# Change height if feature is reverse oriented:
			$var_y -= $var_height2 if $var_orient eq "-";

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	# Store reference to @array_out in %hash_out:
	$hash_out { array_SVGfeatures }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #


sub TicksSVG { 

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{ array_merged }};
	my @array_out;

	# Variables:
	my $var_thick		= $hash_out{stroke_thickness}/2;
	my $var_height		= $hash_out{var_span};
	my $var_factor		= $hash_out{var_factor};
	my $var_increment	= 200;

	my $var_gap			= $hash_out{var_gap};
	my $var_y			= 100 - $var_gap;

	my $var_scalar		= scalar @array_merged;
	my $var_stroke		= "0,0,0";
	my $var_colour		= "0,0,0";

	$var_scalar			= scalar @{$array_merged[0]};
	my $var_length		= $array_merged[0][$var_scalar-1][2];

	my $var_1Mb_height		= $hash_out{var_1Mb_height};
	my $var_100kb_height	= $hash_out{var_100kb_height};
	my $var_10kb_height		= $hash_out{var_10kb_height};
	my $var_1kb_height		= $hash_out{var_1kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_1Mb				= 10**6;
	my $var_100kb			= 10**5;
	my $var_10kb			= 10**4;
	my $var_1kb				= 10**3;

	# ----------------------------------

	# Create horizontal lines:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_copy		= $array_merged[$i][$j][4];
			my $var_length		= $array_merged[$i][$j][6];
			my $var_inverted	= $array_merged[$i][$j][8];
			my $var_loc0		= $array_merged[$i][$j][11];
			my $var_loc1		= $array_merged[$i][$j][12];

			next if $var_loc0 == 0;

			$var_y += $var_height if $var_inverted;

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			my $var_string	= General::RectangleSVG ( $var_loc0, $var_y, $var_length, $var_thick, $var_colour, $var_stroke );

			push @{$array_out[$i]}, $var_string;

			$var_y -= $var_height if $var_inverted;

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	$var_y			= 100 - $var_gap;

	# Create vertical lines:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		my $var_length	= 0;

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_loc1	= $array_merged[$i][$j][2];

			$var_length = $var_loc1 if $var_length < $var_loc1;

		}

		# Iterate through by smallest increment for each individual tick:
		for ( my $j = 0; $j <= $var_length; $j+=$var_1kb ) {

			# Define variables for use outside of loop immediately below:
			my $var_diff		= 0;
			my $var_inverted	= 0;

			# Read through each region and find where to place tick:
			for ( my $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				# Define locations:
				my $var_loc0	= $array_merged[$i][$k][1];
				my $var_loc1	= $array_merged[$i][$k][2];
			 	$var_inverted	= $array_merged[$i][$k][8];
				my $var_loc2	= $array_merged[$i][$k][11];

				# If start of region contains next hit define difference between
				# local and global location values:
				if ( $var_loc0 < $j && $var_loc1 > $j ) {

					$var_diff	= $var_loc2 - $var_loc0;

					last;

				}

			}

			my $var_x		= ( $j + $var_diff )/$var_factor;

			my $var_tick_height	= $var_1kb_height;

			$var_tick_height	= $var_10kb_height if $j % ( $var_10kb ) == 0;
			$var_tick_height	= $var_100kb_height if $j % ( $var_100kb ) == 0;
			$var_tick_height	= $var_1Mb_height if $j % ( $var_1Mb ) == 0;

			$var_y += $var_height if $var_inverted;

			$var_y			-= $var_tick_height;

			my $var_string	= General::RectangleSVG ( $var_x, $var_y, $var_thick,
									 $var_tick_height, $var_colour, $var_stroke );

			$var_y			+= $var_tick_height;

			$var_y -= $var_height if $var_inverted;

			if ( $var_factor == 1000 ) {

				next if $var_tick_height == $var_1kb_height;

				push @{$array_out[$i]}, $var_string;

			}		

		}

		$var_y += $var_increment;

	}

	# ----------------------------------

	# Store reference to @array_out in %hash_out:
	$hash_out { array_SVGticks }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

sub MergeSVG {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_regions	= @{$hash_out{array_regions}};
	my @array_features	= @{$hash_out{array_SVGfeatures}};
	my @array_ticks		= @{$hash_out{array_SVGticks}};
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_out;

	# Variables:
	my $var_scalar			= scalar @array_regions;
	my $file_out			= $hash_out{file_svg};
	my $var_image_radius	= $hash_out{ var_image_radius };
	my $var_increment		= $hash_out { var_increment };
	my $var_factor			= $hash_out { var_factor };
	my $var_height			= 400;

	$var_scalar				= scalar @{$array_merged[0]};
	my $var_width			= ( $array_merged[0][$var_scalar-1][2] ) / $var_factor;

	# ----------------------------------

	# Add header information:
	push @array_out, "<svg width=\"$var_width\" height=\"$var_height\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

	# Add background:
	push @array_out, "<g id=\"bg\">";
	my $var_string	= "<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
	$var_string		.= "height=\"$var_height\" style=\"fill:rgb(255,255,255)\"/>";
	push @array_out, $var_string;
	push @array_out, "</g>";

	# ----------------------------------

	# Add regions from @array_regions:
	for ( my $i = 0; $i < scalar @array_regions; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_regions[$i]}; $j++ ) {

			push @array_out, $array_regions[$i][$j];

		}

	}

	# ----------------------------------

	# Add features from @array_features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		# TODO: Add tags:

		next unless $array_features[$i];

		for ( my $j = 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			push @array_out, $array_features[$i][$j];

		}

	}

	# ----------------------------------

	# Add ticks from @array_ticks:
	for ( my $i = 0; $i < scalar @array_ticks; $i++ ) {

		# TODO: Add tags:

		for ( my $j = 0; $j < scalar @{$array_ticks[$i]}; $j++ ) {

			push @array_out, $array_ticks[$i][$j];

		}

	}

	# ----------------------------------

	push @array_out, "</svg>";

	# ----------------------------------

	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Create Circos config file:
Synteny::Config ( \%hash_variables );

# Create Circos Karyotype file:
Synteny::Karyotype ( \%hash_variables );

# Create highlights files:
Synteny::Highlights ( \%hash_variables );

# Create Circos Ideogram file:
Circos::Ideogram ( \%hash_variables );

# Create Circos Image file:
Circos::Image ( \%hash_variables );

# Create ticks files:
Synteny::Ticks ( \%hash_variables );

# Run Circos:
Circos::Run ( \%hash_variables );

# ---------------------------------------------------------------------------- #

