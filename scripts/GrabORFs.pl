#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib './scripts';
use General;

# ---------------------------------------------------------------------------- #

# File name:		GrabORFs.pl
# Date created:		14 February, 2019
# Last modified:	15 February, 2019
# Created by:		Eliot Stanton

# Description:		This is a script for pulling ORF locations from a GFF file.

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my @array_in;
my @array_out;

# Import variables from the command-line provided by the user:
my $file_in		= $ARGV[0];
my $var_tag		= $ARGV[1];
my $file_out	= $ARGV[2];

# ---------------------------------------------------------------------------- #

# Open $file_in and store in @array_in:
@array_in	= @{General::FileToArray ( $file_in )};

# Iterate through @array_in:
for ( my $i = 0; $i < scalar @array_in; $i++ ) {

	my $var_loc0	= $array_in[$i][3];
	my $var_loc1	= $array_in[$i][4];

	my $var_orient	= $array_in[$i][6];

	my $var_string	= "0 gene $var_tag $var_loc0 $var_loc1 $var_orient";

	push @array_out, $var_string;

}

# ---------------------------------------------------------------------------- #

General::ArrayToFile ( \@array_out, $file_out );

# ---------------------------------------------------------------------------- #
