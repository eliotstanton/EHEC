#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;
use Circos;

# ---------------------------------------------------------------------------- #

# File name:		ORFs.pl
# Date created:		04 January, 2019
# Last modified:	07 January, 2019
# Created by:		Eliot Stanton

# Description:		This is a script for converting ORFs in GenBank data into 
#					SVG image data.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;
my @array_genbank;
my @array_align;
my @array_out;

# Import variables from the command-line provided by the user:
getopts('a:e:g:n:o:s:t:', \%hash_variables);

# Define variables used by this script:
my $file_alignment	= $hash_variables{a};
my $var_stem		= $hash_variables{e} || 0;
my $file_genbank	= $hash_variables{g};
my $var_output		= $hash_variables{o};
my $var_start		= $hash_variables{s} || 0;
my $var_stop		= $hash_variables{t} || 0;

my $file_out		= "$var_output\/default.$var_stem.svg";
my $var_y			= 100;
my $var_scale		= 100;
my $var_width		= 1000;
my $var_height		= 400;

# ---------------------------------------------------------------------------- #

# Define help file:
my $var_help	= "\nORFs.pl [OPTIONS]

	-a Mauve alignment file
	-e Stem for file name
	-g GenBank file (required)
	-o Output directory (required)
	-s Start location
	-t Stop location\n";

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output is missing print $var_help and stop script:
unless ( $hash_variables{o} && $hash_variables{g} ) {

	print "$var_help\n";

	print "\tGenBank file required!\n" unless $hash_variables{g};

	print "\tOutput directory required!\n" unless $hash_variables{o};

	exit;

}

# Add more variables to %hash_variables:
#%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# TODO: Check if $file_genbank is present and contains data:

# Check if output directory is present, if it isn't, create it:
General::CheckDirectory ( \%hash_variables );

# ---------------------------------------------------------------------------- #

# Import GenBank data into @array_genbank:
@array_genbank		= @{ General::FileToArray ( $file_genbank ) };

# Import alignment data into @array_alignment:
if ( $file_alignment ) {

	@array_align		= @{ General::FileToArray ( $file_alignment ) };

}

# Reformat GenBank data to hold only coordinates and orientation:
@array_genbank		= @{ FormatGB ( \@array_genbank, $var_start, $var_stop ) };

# Align data with alignment file if required:
@array_out	= @{Align ( \@array_genbank, \@array_align )};

# Remove and shift locations to beginning if required:
if ( $var_start && $var_stop ) {

	@array_out	= @{ Shift ( \@array_out, $var_start, $var_stop ) };

	$var_width	= ( $var_stop - $var_start ) / $var_scale;

	print "$var_width\n";

}

# Translate location data into SVG formatting:
@array_out			= @{ Translate ( \@array_out, $var_y, $var_height, $var_scale ) };

# ---------------------------------------------------------------------------- #

unshift @array_out, "</g>";

# Add background:

my $var_string	= "<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
$var_string		.= "height=\"$var_height\" style=\"fill:rgb(255,255,255)\"/>";

unshift @array_out, $var_string;

unshift @array_out, "<g id=\"bg\">";

# Add header information:
unshift @array_out, "<svg width=\"$var_width\" height=\"$var_height\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

push @array_out, "</svg>";

# Print data in @array_out to file:
General::ArrayToFile ( \@array_out, $file_out );

print "$file_out\n";

# ---------------------------------------------------------------------------- #

# Subroutine name:	Translate
# Description:

# ----------------------------------------------------------

sub Translate {

	# Arguments:
	my ( $array_in, $var_y, $var_height, $var_scale )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# Variables:
	my $var_stroke	= "0,0,0";
	my $var_colour	= Circos::ColoursNew ( "chromosome", 1, 1 );
	$var_height		/= 4;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0		= $array_in[$i][0];
		my $var_loc1		= $array_in[$i][1];
		my $var_orientation	= $array_in[$i][2];

		my $var_length		= $var_loc1 - $var_loc0 - 1;

		my $var_x			= $var_loc0/$var_scale;
		$var_length			= $var_length/$var_scale;

		$var_y				+= $var_height if $var_orientation eq "-";

		my $var_colour		= "(0,0,0)";

		my $var_string		= General::RectangleSVG ( $var_x, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

		$var_y				-= $var_height if $var_orientation eq "-";

		push @array_out, $var_string;

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}
# ---------------------------------------------------------------------------- #

# Subroutine name:	Shift
# Description:

# ----------------------------------------------------------

sub Shift {

	# Arguments:
	my ( $array_in, $var_start, $var_stop )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# Variables:
	my $var_length	= $var_stop - $var_start;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

		$var_loc0 -= $var_start;
		$var_loc1 -= $var_start;

		$array_out[$i][0]	= $var_loc0;
		$array_out[$i][1]	= $var_loc1;
		$array_out[$i][2]	= $array_in[$i][2];

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Align
# Description:		

# ----------------------------------------------------------

sub Align {

	# Arguments:
	my ( $array_in, $array_align, $var_number )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_align	= @$array_align;
	my @array_out;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0	= $array_in[$i][0];
		my $var_loc1	= $array_in[$i][1];

		for ( my $j = 0; $j < scalar @array_align; $j++ ) {

			my $var_loc2	= $array_align[$j][0];
			my $var_loc3	= $array_align[$j][1];
			my $var_loc4	= $array_align[$j][2];
			my $var_loc5	= $array_align[$j][3];

			my $var_diff	= $var_loc4 - $var_loc2;

			if ( ( $var_loc2 <= $var_loc0 && $var_loc3 > $var_loc0 )
				&& ( $var_loc2 <= $var_loc1 && $var_loc3 > $var_loc1 ) ) {

				$var_loc1 += $var_diff;

				$var_loc0 += $var_diff;

			}

			if ( $var_loc2 > $var_loc1 ) {

				last;

			}

		}

		$array_out[$i][0]	= $var_loc0;
		$array_out[$i][1]	= $var_loc1;
		$array_out[$i][2]	= $array_in[$i][2];

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Translate
# Description:

# ----------------------------------------------------------

sub FormatGB {

	# Arguments:
	my ( $array_in, $var_start, $var_stop )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		if ( $array_in[$i][0] && $array_in[$i][0] =~ /^CDS/ && length $array_in[$i][0] == 3 ) {

			my @array_temp	= split /[\.,),(]/, $array_in[$i][1];

			my $var_orientation	= "+";

			if ( $array_temp[0] eq "complement" ) {

				$var_orientation	= "-";

			}

			my $var_loc0	= $array_temp[$#array_temp-2];
			my $var_loc1	= $array_temp[$#array_temp];

			@array_temp		= ( $var_loc0, $var_loc1, $var_orientation );

			if ( $var_start && $var_stop ) {

				next unless $var_loc0 < $var_stop;
				next unless $var_loc1 > $var_start;

				push @array_out, \@array_temp;

			}

			else {

				push @array_out, \@array_temp;

			}

		}

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

