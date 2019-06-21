#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;
use Homology;

# ---------------------------------------------------------------------------- #

# File name:		WriteORFs.pl
# Date created:		20 February, 2018
# Last modified:	20 February, 2018
# Created by:		Eliot Stanton

# Description:		This is script takes genomic features and writes them to
#					SVG format.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my @array_in	= @{General::FileToArray( $ARGV[0] )};
my @array_out;

# Define variable used by this script:
my $file_out	= $ARGV[1];
my $var_factor	= 10;
my $var_y		= 0;
my $var_y2		= 100;
my $var_height	= 200;
my $var_height2	= 100;
my $var_stroke	= 1;
my $var_colour	= "0,0,0";
my $var_width	= 10000;

# ---------------------------------------------------------------------------- #

for ( my $i = 0; $i < scalar @array_in; $i++ ) {

	# Define the type of genomic feature, locations, and orientation:
	my $var_loc0		= $array_in[$i][3];
	my $var_loc1		= $array_in[$i][4];
	my $var_type		= $array_in[$i][1];
	my $var_orient		= $array_in[$i][5] || "0";

	my $var_length		= $var_loc1 - $var_loc0 + 1;

	$var_length			/= $var_factor;
	$var_loc0			/= $var_factor;

	my $var_string;

	if ( $var_orient eq "+" ) {

		# Create string holding SVG formatted region:
		$var_string	= General::RectangleSVG ( $var_loc0, $var_y, $var_length, $var_height2, $var_colour, $var_stroke );

	}

	else {

		# Create string holding SVG formatted region:
		$var_string	= General::RectangleSVG ( $var_loc0, $var_y2, $var_length, $var_height2, $var_colour, $var_stroke );

	}

#	print "@{$array_in[$i]}\n";

#	print "$var_string\n";

	# Store $var_string in @array_out:
	push @array_out, $var_string;

}

# ---------------------------------------------------------------------------- #

# Print data to file:
for ( my $i = 0; $i < scalar @array_out; $i++ ) {

	my @array_temp;

	# Add header information:
	push @array_temp, "<svg width=\"$var_width\" height=\"$var_height\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">";

	# Add background:
	push @array_temp, "<g id=\"bg\">";
	my $var_string	= "<rect x=\"0\" y=\"0\" width=\"$var_width\" ";
	$var_string		.= "height=\"$var_height\" style=\"fill:rgb(255,255,255)\"/>";

	push @array_temp, "$var_string\n";

	push @array_temp, "</g>";

	for ( my $j = 0; $j < scalar @array_out; $j++ ) {

#		print "$array_out[$j]\n";

		push @array_temp, $array_out[$j];

	}

	push @array_temp, "</svg>";

	General::ArrayToFile ( \@array_temp, $file_out );

}


# ---------------------------------------------------------------------------- #

