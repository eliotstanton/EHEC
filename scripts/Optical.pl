#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use lib "./scripts";
use General;
use Circos;
#use Svg;
use Optical;

# ---------------------------------------------------------------------------- #

# File name:        optical.pl
# Created by:       Eliot Stanton (estanton@wisc.edu)
# Created on:       18 December, 2018
# Last modified:	18 December, 2018

# Description:      This is a script for visualising optical mapping data.

# ---------------------------------------------------------------------------- #

# TODO: Add SVG option.z

# ---------------------------------------------------------------------------- #

# Define data-structures used by this script:
my %hash_variables;

# Import variables from the command-line provided by the user:
getopts('o:', \%hash_variables);

# Define variables for command-line options:
my $var_output		= $hash_variables {o};
my $var_input		= $ARGV[0];

# Define help file:
my $var_help	= "\noptical.pl [OPTIONS] [MAP1],[MAP2]";

# ---------------------------------------------------------------------------- #

# Store $var_help in %hash_variables:
$hash_variables	{ var_help }	= $var_help;

# If $var_output or $var_input is missing, print $var_help and stop script:
unless ( $var_input && $var_output ) {

	print "$var_help\n";

	print "\tInput map file required!\n" unless $var_input;

	print "\tOutput directory required!\n" unless $var_output;

	exit;

}

# Add more variables to %hash_variables:
%hash_variables		= %{General::Variables ( \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Check if output directory is present, if it isn't, create it:
General::CheckDirectory ( \%hash_variables );

# Split $var_input into @array_files:
my @array_files		= split /\,/, $var_input;

# Check input files:
@array_files	= @{General::CheckFiles ( \@array_files, \%hash_variables )};

# ---------------------------------------------------------------------------- #

# Import data from map file to @array_optical
my @array_optical	= @{General::ImportFiles ( \@array_files )};

# Convert data for each map into a location:
@array_optical		= @{Optical::Locations(\@array_optical)};

# ---------------------------------------------------------------------------- #

#  Create Circos configuration file:
Optical::Config ( \@array_optical, \%hash_variables );

# Create Circos highlights files:
Optical::Highlights ( \@array_optical, \%hash_variables );

# Create Circos karyotype files:
Optical::Karyotype ( \@array_optical, \%hash_variables );

# Create Circos ideogram file:
Circos::Ideogram ( \%hash_variables );

# Create Circos image file:
Circos::Image ( \%hash_variables );

# Create Circos ticks files:
Optical::Ticks ( \@array_optical, \%hash_variables );

# Run Circos:
Circos::Run ( \%hash_variables );

# ---------------------------------------------------------------------------- #
