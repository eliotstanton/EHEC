package Synteny;

use strict;
use warnings;
use Algorithm::NeedlemanWunsch;

# ---------------------------------------------------------------------------- #

# File name:		Synteny.pm
# Date created:		31 October, 2018
# Last modified:	30 November, 2018
# Created by:		Eliot Stanton

# Description:		This package contains subroutines specific to adjusting
#					synteny data.

# ---------------------------------------------------------------------------- #

# Subroutine name:	Alignment
# Description:		This subroutine handles alignment using the Needleman-Wunsch
#					algorithm.

# ----------------------------------------------------------

sub Alignment {

	# Arguments:
	my ($array_in0, $array_in1)	= @_;

	# Data structures:
	my @array_in0	= @$array_in0;
	my @array_in1	= @$array_in1;
	my @array_matrix;
	my @array_matrix2;

	# Variables:
	my $var_counter0	= -1;
	my $var_counter1	= -1;
	my $var_scalar0		= scalar @array_in0;
	my $var_scalar1		= scalar @array_in1;
	my $var_match		= 1;
	my $var_penalty		= 1;
	my $var_gap			= 1;

	# ----------------------------------

	# Add starting location values:
	$array_matrix[0][0]		= 0;
	$array_matrix2[0][0]	= "none";

	# Initialise the matrix:
	for ( my $i	= 1; $i < $var_scalar1 + 1; $i++ ) {

		$array_matrix[0][$i]	= $var_counter0;
		$array_matrix2[0][$i]	= "horizontal";

		$var_counter0--;

	}

	for ( my $i	= 1; $i <= $var_scalar0; $i++ ) {

		$array_matrix[$i][0]	= $var_counter1;
		$array_matrix2[$i][0]	= "vertical";

		$var_counter1--;

	}

	# ----------------------------------

	# Iterate down through matrix:
	for ( my $i = 0; $i < $var_scalar0 + 1; $i++ ) {

		for ( my $j	= 0; $j < $var_scalar1 + 1; $j++ ) {

			next if $i == 0 && $j == 0;

			my $var_a	= $array_in0[$i-1];
			my $var_b	= $array_in1[$j-1];

			unless ( $array_matrix[$i][$j] ) {

				my $diagonal	= $array_matrix[$i-1][$j-1];
				my $vertical	= $array_matrix[$i-1][$j] - $var_gap;
				my $horizontal	= $array_matrix[$i][$j-1] - $var_gap;

				if ( $var_a eq $var_b ) { 

					$diagonal	+= $var_match;

				}

				else { 

					$diagonal	-= $var_penalty;

				}

				if ( $diagonal >= $vertical ) {

					if ( $diagonal >= $horizontal ) {

						$array_matrix[$i][$j]	= $diagonal;
						$array_matrix2[$i][$j]	= "diagonal";

					}

					else {

						$array_matrix[$i][$j]	= $horizontal;
						$array_matrix2[$i][$j]	= "horizontal";

					}

				}

				else {

					if ( $vertical >= $horizontal ) {

						$array_matrix[$i][$j]	= $vertical;
						$array_matrix2[$i][$j]	= "vertical";

					}

					else {

						$array_matrix[$i][$j]	= $horizontal;
						$array_matrix2[$i][$j]	= "horizontal";

					}

				}

			}

		}

	}

	# ----------------------------------

	my @array_out0;
	my @array_out1;

	my $i	= $var_scalar0;
	my $j	= $var_scalar1;

	my $var_counter	= 20;

	# Trace-back through the matrix:
	while ($var_counter > 0) {

		last if $array_matrix2[$i][$j] eq "none";

		my $var_direction	= $array_matrix2[$i][$j];

		if ( $var_direction eq "horizontal" ) {

			unshift @array_out0, "$array_in1[$j-1]";
			unshift @array_out1, "-";

			$j--;

		}

		elsif ( $var_direction eq "vertical" ) {

			unshift @array_out0, "-";
			unshift @array_out1, "$array_in0[$i-1]";

			$i--;

		}

		elsif ( $var_direction eq "diagonal" ) {

			unshift @array_out0, "$array_in1[$j-1]";
			unshift @array_out1, "$array_in0[$i-1]";

			$i--;
			$j--;

		}

		$var_counter--;

	}

	# ----------------------------------

	return ( \@array_out0, \@array_out1 );

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketHash
# Description:		This subroutine stores regions in a hash of arrays using
#					pairs of bracketing reference IDs as keys.

# ----------------------------------------------------------

sub BracketHash {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my %hash_out;

	# ----------------------------------

	# Iterate through @array_in and store regions in %hash_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Iterate down through each region:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define bracketing reference IDs:
			my $var_bracket0	= $array_in[$i][$j][9];
			my $var_bracket1	= $array_in[$i][$j][10];

			# Move on if region is a blank:
			next if $var_bracket0 == 0 && $var_bracket1 == 0;

			# Add region from @array_in to %hash_merged:
			push ( @{$hash_out{"$var_bracket0 $var_bracket1"}}, $array_in[$i][$j]);

		}

	}

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketIDs
# Description:		This subroutine adds bracketIDs to each region:

# ----------------------------------------------------------

sub BracketIDs {

	# Arguments:
	my ( $array_in	)	= @_;

	# Data structures:
	my @array_out;

	# ----------------------------------

	# Sort regions by sequence and location:
	@array_out	= @{ SortBySequence ( $array_in ) };

	# ----------------------------------

	# Add bracket IDs to each region:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_out[$i]}; $j++ ) {

			# Define variables holding reference IDs of the current and former
			# region:
			my $var_refID0	= $array_out[$i][$j][3];
			my $var_refID1	= $array_out[$i][$j-1][3];

			# Store bracketing reference IDs to current and former region:
			$array_out[$i][$j][9]		= $var_refID1;
			$array_out[$i][$j-1][10]	= $var_refID0;

		}

	}

	# ----------------------------------

	# Sort regions by reference ID:
	@array_out	= @{ SortByReference ( \@array_out ) };

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketInserts
# Description:		This subroutine assigns bracketing reference IDs to each
#					insert region.

# ----------------------------------------------------------

sub BracketInserts {

	# Arguments:
	my ( $array_inserts, $array_backbone )	= @_;

	# Data structures:
	my @array_backbone		= @$array_backbone;
	my @array_inserts		= @$array_inserts;

	# ----------------------------------

	# Iterate through @array_inserts:
	for ( my $i = 0; $i < scalar @array_inserts; $i++ ) {

		# Iterate through regions:
		for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

			# Define sequence ID and initial location:
			my $var_seq		= $array_inserts[$i][$j][0];
			my $var_loc0	= $array_inserts[$i][$j][1];

			# Move on if region is blank:
			next if $var_loc0 == 0;

			# Iterate down through backbone regions:
			for ( my $k = 0; $k < scalar @{$array_backbone[$var_seq]}; $k++ ) {

				# Define initial backbone region location:
				my $var_loc1	= $array_backbone[$var_seq][$k][1];

				# If backbone region is larger define bracketing reference IDs
				# and end loop:
				if ( $var_loc1 > $var_loc0 ) {

					# Define reference IDs of the current and former backbone
					# regions:
					my $var_refID0		= $array_backbone[$var_seq][$k-1][3];
					my $var_refID1		= $array_backbone[$var_seq][$k][3];

					# Define inverted status of both regions:
					my $var_inverted0	= $array_backbone[$var_seq][$k-1][8];
					my $var_inverted1	= $array_backbone[$var_seq][$k][8];

					if ( $var_refID0 < $var_refID1 ) {

						# Store bracketing data in region:
						$array_inserts[$i][$j][9]	= $var_refID0;
						$array_inserts[$i][$j][10]	= $var_refID1;

					}

					else {

						# Store bracketing data in region:
						$array_inserts[$i][$j][9]	= $var_refID1;
						$array_inserts[$i][$j][10]	= $var_refID0;

					}

					# If $var_inverted0 and $var_inverted1 both are positive,
					# redefine region as inverted as well:
					if ( $var_inverted0 == 1 && $var_inverted1 == 1 ) {

						$array_inserts[$i][$j][8]	= 1;

					}

					last;

				}

			}

		}

	}

	# ----------------------------------

	# Return reference to \@array_inserts and end subroutine:
	return \@array_inserts;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	CheckRegions
# Description:		This subroutine checks that all regions are accounted for
#					and in the correct order.

# ----------------------------------------------------------

sub CheckRegions {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_all		= %{$hash_out{hash_all}};
	my @array_merged	= @{$hash_out{array_merged}};
	my %hash_merged;
	my %hash_inserts	= %{$hash_out{hash_inserts}};

	# ----------------------------------

	# Move regions in @array_merged to %hash_merged:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_seqID	= $array_merged[$i][$j][0];
			my $var_refID	= $array_merged[$i][$j][3];

			$hash_merged{"$var_seqID $var_refID"}	= $array_merged[$i][$j];

		}

	}

	# ----------------------------------

	# Check that all the regions in %hash_all are in %hash_merged:
	foreach my $var_entry ( keys %hash_all ) {

		unless ( $hash_merged { $var_entry } ) {

			my @array_temp	= @{$hash_all { $var_entry }};

			print "$var_entry: @array_temp\n";

		}

	}

	# ----------------------------------

	# Check that all the regions in @array_merged are in order:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		my $var_loc0	= $array_merged[$i][0][1];

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_loc1	= $array_merged[$i][$j][1];

			next if $var_loc1 == 0;

			if ( $var_loc0 > $var_loc1 ) {

				print "$i:\n";

				print "\t$j: @{$array_merged[$i][$j-1]} - $var_loc0\n";
				print "\t$j: @{$array_merged[$i][$j]} - $var_loc1\n";

			}

			$var_loc0 = $var_loc1;

		}

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Config
# Description:		This subroutine creates a config file for Circos

# ----------------------------------------------------------

sub Config {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_merged	= @{$hash_in{array_merged}};
	my @array_out;

	# Variables:
	my $file_out			= $hash_in { "file_config" };
	my $file_highlights		= $hash_in { "file_highlights" };
	my $file_ideogram		= $hash_in { "file_ideogram" };
	my $file_image			= $hash_in { "file_image" };
	my $file_karyotype		= $hash_in { "file_karyotype" };
	my $file_ticks			= $hash_in { "file_ticks" };
	my $var_scalar			= scalar @array_merged;
	my $var_increment		= $hash_in { var_increment };
	my $var_radius0			= $hash_in { var_radius0 };
	my $var_radius1			= $hash_in { var_radius1 };
	my $r					= "r";

	# ----------------------------------

	push @array_out, "#$file_out";
	push @array_out, "karyotype = $file_karyotype";
	push @array_out, "<<include $file_ideogram>>";
	push @array_out, "<ideogram>\nshow = no\n</ideogram>";

	push @array_out, "<highlights>";

	for ( my $i = 0; $i < $var_scalar; $i++ ) {

		# Add highlights:
		push @array_out, "\t<highlight>";
		push @array_out, "\t\tfile	= $file_highlights.$i\n\t\tideogram	= no";
		push @array_out, "\t\tr0	= $var_radius0$r\n\t\tr1	= $var_radius1$r";
	 	push @array_out, "\t</highlight>";

		# Add ticks:
		push @array_out, "\t<highlight>";
		push @array_out, "\t\tfile	= $file_ticks.$i\n\t\tideogram	= no";
	 	push @array_out, "\t</highlight>";
	 
	 	$var_radius0	-= $var_increment;
	 	$var_radius1	-= $var_increment;

	}
 
 	push @array_out, "</highlights>";
 
	push @array_out, "<<include $file_image>>";
	push @array_out, "<<include etc/colors_fonts_patterns.conf>>";
	push @array_out, "<<include etc/housekeeping.conf>>";

	# Write @array_out to file:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	BracketRefs
# Description:		This subroutine is a wrapper for adding bracketing reference
#					IDs to insert and transposed regions.

# ----------------------------------------------------------

sub BracketRefs {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my @array_inserts		= @{$hash_out{array_inserts}};
	my @array_transposed	= @{$hash_out{array_transposed}};

	# ----------------------------------

	# Sort regions in @array_backbone by sequence:
	@array_backbone		= @{Synteny::SortBySequence(\@array_backbone)};

	# Add bracketing reference IDs to @array_inserts:
	@array_inserts		= @{Synteny::BracketInserts(\@array_inserts,\@array_backbone)};

	# Add bracketing reference IDs to @array_transposed:
	@array_transposed	= @{Synteny::BracketInserts(\@array_transposed,\@array_backbone)};

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_backbone;

	# Store @array_inserts in %hash_out:
	$hash_out { array_inserts }		= \@array_inserts;

	# Store @array_inserts in %hash_out:
	$hash_out { array_transposed }	= \@array_transposed;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FeaturesSVG
# Description:		This subroutine formats genomic feature data for SVG.

# ----------------------------------------------------------

sub FeaturesSVG { 

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_features	= @{$hash_out{ array_features }};
	my @array_out;

	# Variables:
	my $var_height		= $hash_out{var_span};
	my $var_factor		= $hash_out{var_factor};
	my $var_y			= 100;
	my $var_increment	= $hash_out{var_increment} * $var_factor;
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

			# Calculate length of genomic feature:
			my $var_length		= $var_loc1 - $var_loc0 + 1;

			# Change $var_height if feature is inverted:
			$var_y += $var_height if $var_inverted;

			$var_length		/= $var_factor;
			$var_loc0		/= $var_factor;

			# Calculate genomic feature colour:
			my $var_colour	= Circos::ColoursNew ( $var_type, $var_copy, $var_scalar );

			# Create string holding SVG formatted region:
			my $var_string	= General::RectangleSVG ( $var_loc0, $var_y, $var_length, $var_height, $var_colour, $var_stroke );

			# Store $var_string in @array_out:
			push @{$array_out[$i]}, $var_string;

			# Decrease $var_height if handling an inverted region:
			$var_y -= $var_height if $var_inverted;

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

# Subroutine name:	FormatAlignment
# Description:		This subroutine takes output from Circos and formats data
#					for use by other subroutines.
# [0: SeqID][1: Loc0][2: Loc1][3: RefID][4: Copy][5: Length0][6: Length1]
# [7: Inverted0][8: Inverted1][9: Bracket0][10: Bracket1]
# [11: Global0][12: Global1]

# ----------------------------------------------------------

sub FormatAlignment {

	# Arguments:
	my ( $hash_in, $array_FASTA )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_in;
	my @array_inserts;
	my @array_inverted;
	my @array_backbone;
	my %hash_all;

	# Variables:
	my $file_backbone	= $hash_out { file_backbone };
	my $var_minimum		= $hash_out { var_minimum };

	# ----------------------------------

	# If $file_backbone does not exist, run progressiveMauve:
	progressiveMauve ( $hash_in, $array_FASTA ) unless -e $file_backbone;

	# Load $file_backbone into @array_in:
	@array_in	= @{ General::FileToArray ( $file_backbone ) };

	# ----------------------------------

	# Iterate through each alignment in @array_in:
	for ( my $i = 1; $i < scalar @array_in; $i++ ) {

		# Define variable to hold copy number:
		my $var_copy	= 0;

		# Define variable to hold global length:
		my $var_length0	= 0;

		# Define array to hold alignment data:
		my @array_alignment;

		# Define variable to hold inverted status:
		my $var_inverted0	= 0;

		# Define variable to hold sequence ID:
		my $var_seq		= 0;

		# Iterate through each pair of locations:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j+=2 ) {

			my $var_loc0		= $array_in[$i][$j];
			my $var_loc1		= $array_in[$i][$j+1];
			my $var_inverted1	= 0;

			if ( $var_loc0 < 0 ) {

				$var_loc0 *= -1;
				$var_loc1 *= -1;

				$var_inverted1	= 1;
				$var_inverted0	= 1;

			}

			$var_copy++ if $var_loc0 > 0;

			# Define length of region:
			my $var_length1	= $var_loc1 - $var_loc0 + 1;

			# If $var_length1 is greater than $var_length0, reassign 
			# $var_length0 to equal $var_length0:
			$var_length0 = $var_length1 if $var_length1 > $var_length0;

			# Create a temporary array to hold region data:
			my @array_temp	= (0)x13;

			$array_temp[0]	= $var_seq;
			$array_temp[1]	= $var_loc0;
			$array_temp[2]	= $var_loc1;
			$array_temp[3]	= $i;
			$array_temp[6]	= $var_length1;
			$array_temp[8]	= $var_inverted1;

			push @array_alignment, \@array_temp;

			$var_seq++;

		}
	 
		# Move on if $var_length0 is less than $var_minimum:
		next if $var_length0 < $var_minimum;

		# Store $var_length0, $var_copy and $var)inverted0 in each alignment:
		for ( my $j = 0; $j < scalar @array_alignment; $j++ ) {

			my $var_seqID			= $array_alignment[$j][0];
			my $var_loc0			= $array_alignment[$j][1];
			my $var_refID			= $array_alignment[$j][3];

			$array_alignment[$j][4] = $var_copy;
			$array_alignment[$j][5]	= $var_length0;
			$array_alignment[$j][7]	= $var_inverted0;

			next if $var_loc0 == 0;

			$hash_all{"$var_seqID $var_refID" }	= $array_alignment[$j];

		}

		# Add regions to appropiate array:
		if ( $var_copy == $var_seq ) {

			push @array_backbone, \@array_alignment;

		}

		else {

			push @array_inserts, \@array_alignment;

		}

	}

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_backbone;

	# Store @array_inserts in %hash_out:
	$hash_out { array_inserts }		= \@array_inserts;

	# Store %hash_all in %hash_out:
	$hash_out { hash_all }			= \%hash_all;

	# ----------------------------------

	# Return referene to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	FormatFeatures
# Description:		This subroutine create feature blocks with adjusted
#					locations and with feature type data.

# ----------------------------------------------------------

sub FormatFeatures {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_out;

	# ----------------------------------

	# Iterate down through genomic features for each sequence:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Go through each feature and find overlapping regions of the
		# chromosome:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define start, end, type, and orientation of feature:
			my $var_loc0	= $array_in[$i][$j][3];
			my $var_loc1	= $array_in[$i][$j][4];
			my $var_type	= $array_in[$i][$j][1];
			my $var_orient	= $array_in[$i][$j][5];

			# Iterate down through regions in @array_merged finding those region
			# blocks that overlap::
			for ( my $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				# Define start, end, and inverted status of each regions:
				my $var_loc2		= $array_merged[$i][$k][1];
				my $var_loc3		= $array_merged[$i][$k][2];
				my $var_copy		= $array_merged[$i][$k][4];
				my $var_inverted	= $array_merged[$i][$k][8];

				# Skip blanks:
				next if $var_loc2 == 0;

				# Define adjusted locations of region:
				my $var_loc4	= $array_merged[$i][$k][11];
				my $var_loc5	= $array_merged[$i][$k][12];

				# Define variables to hold adjusted locations of features:
				my $var_diff	= $var_loc4-$var_loc2;
				my $var_new0	= 0;
				my $var_new1	= 0;

				# If region straddles beginning of genomic feature:
				if ( $var_loc2 < $var_loc0 && $var_loc3 > $var_loc0 && $var_loc3 < $var_loc1 ) {

					$var_new0	= $var_loc0 + $var_diff;
					$var_new1	= $var_loc3 + $var_diff;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				# If region straddles end of genomic feature:
				elsif ( $var_loc2 > $var_loc0 && $var_loc2 < $var_loc1 && $var_loc3 > $var_loc1 ) {
 
					$var_new0	= $var_loc2 + $var_diff;
					$var_new1	= $var_loc1 + $var_diff;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				# If region straddles the entire genomic feature:
				if ( $var_loc2 < $var_loc0 && $var_loc3 > $var_loc1 ) {

					$var_new0	= $var_loc0 + $var_diff;
					$var_new1	= $var_loc1 + $var_diff;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				# If region straddles an interior portion of the feature:
				if ( $var_loc2 >= $var_loc0 && $var_loc3 <= $var_loc1 ) {

					$var_new0	= $var_loc2 + $var_diff;
					$var_new1	= $var_loc3 + $var_diff;

					my @array_temp	= ( $var_new0, $var_new1, $var_type, $var_inverted, $var_copy, $var_orient );

					push @{$array_out[$i]}, \@array_temp;

				}

				last if $var_loc2 > $var_loc1;

			}

		}

	}

	# ----------------------------------

	# Store @array_out in %hash_out:
	$hash_out { array_features }	= \@array_out;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	GappedAlignment
# Description:		This subroutine uses Needleman-Wunsch algorithm to identify
#					transposed regions.

# ----------------------------------------------------------

sub GappedAlignment {

	# Arguments:
	my ( $hash_in, $array_in0, $array_in1 )	= @_;

	# Data-structures:
	my %hash_out	= %$hash_in;
	my @array_in0	= @{$array_in0};
	my @array_in1	= @{$array_in1};
	my @array_out0;
	my @array_out1;

	# ----------------------------------

	# Scoring function:
	sub score_sub {

		# Gap penalty:
	  	if (!@_) { return -2; }

		# Mismatch scores -1, match +1
		return ($_[0] eq $_[1]) ? 1 : -1;

	} 

	# ----------------------------------

	my $var_matcher	= Algorithm::NeedlemanWunsch->new(\&score_sub);
	my $var_score	= $var_matcher->align( \@array_in0, \@array_in1, {
	   
		align => sub { 

			unshift @array_out0, $array_in0[$_[0]];
			unshift @array_out1, $array_in1[$_[1]];

		},

		shift_a => sub { 

			unshift @array_out0, $array_in0[$_[0]];
			unshift @array_out1, " ";
			$hash_out{$array_in0[$_[0]]}	= "1";

		},

		shift_b => sub { 

			unshift @array_out0, " ";
			unshift @array_out1, $array_in1[$_[0]];
			$hash_out{$array_in1[$_[0]]}	= "1";

		}
	} );

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	GlobalLocations
# Description		This subroutine assigns final locations for Circos.

# ----------------------------------------------------------

sub GlobalLocations {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_in	= @{ $hash_out { array_merged }};

	# Variables:
	my $var_start	= 1;
	my $var_end		= 0;

	# ----------------------------------

	# Iterate through @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variables to hold start and end locations:
		my $var_start	= 1;
		my $var_end		= 0;

		# Iterate down through regions assigning locations:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define variables holding global and location region lengths:
			my $var_length0	= $array_in[$i][$j][5];
			my $var_length1	= $array_in[$i][$j][6];

			$var_end	+= $var_length0;

			unless ( $var_length1 == 0 ) {

				my $var_diff	= int( ( $var_length0 - $var_length1 )/2 );

				my $var_loc0	= $var_start + $var_diff;
				my $var_loc1	= $var_end - $var_diff;

				$array_in[$i][$j][11]	= $var_loc0;
				$array_in[$i][$j][12]	= $var_loc1;

			}

			$var_start	= $var_end + 1;

		}

	}

	# ----------------------------------

	# Store @array_in back in %hash_out:
	$hash_out { array_merged }	= \@array_in;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Highlights
# Description:		This subroutine handles formatting highlights for Synteny.

# ----------------------------------------------------------

sub Highlights {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @{ $hash_in { array_merged } };
	my @array_features	= @{ $hash_in { array_features } };
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { file_highlights };
	my $var_thickness	= $hash_in { stroke_thickness };
	my $var_span		= $hash_in { var_span };
	my $var_scalar		= scalar @array_in;

	# ----------------------------------

	# Add highlights for regions:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Iterate down through @array_in formatting data for Circos:
		for ( my $j	= 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			my $var_copy		= $array_in[$i][$j][4];
			my $var_inverted	= $array_in[$i][$j][8];
			my $var_loc0		= $array_in[$i][$j][11];
			my $var_loc1		= $array_in[$i][$j][12];

			# Ignore blank regions:		 
			next if $var_loc0	== 0;

			my $var_border		= $hash_colours { chromosome_border };

			my $var_colour		= Circos::ColoursNew ( "chromosome", $var_copy, $var_scalar );

			my $var_string		= "chr0 $var_loc0 $var_loc1 ";

			$var_string			.= "offset=$var_span," if $var_inverted;
			$var_string			.= "stroke_color=$var_border,";
			$var_string			.= "stroke_thickness=$var_thickness,";
			$var_string			.= "fill_color=$var_colour";

			push @{$array_out[$i]}, $var_string;

		}

	}

	# ----------------------------------

	# Repeat the process for genomic features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		next unless $array_features[$i];

		# Go down through each genomic feature stored in @array_features:
		for ( my $j	= 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			# Define start and end locations of feature block:
			my $var_loc0		= $array_features[$i][$j][0];
			my $var_loc1		= $array_features[$i][$j][1];

			# Define feature type, inverted status and copy number:
			my $var_type		= $array_features[$i][$j][2];
			my $var_inverted	= $array_features[$i][$j][3];
			my $var_copy		= $array_features[$i][$j][4];

			# Define feature colours:
			my $var_border		= $hash_colours { "$var_type\_border" };

			my $var_colour		= Circos::ColoursNew ( $var_type, $var_copy, $var_scalar );

			my $var_string		= "chr0 $var_loc0 $var_loc1 ";

			$var_string			.= "offset=$var_span," if $var_inverted;
			$var_string			.= "stroke_color=$var_border," if $var_border;
			$var_string			.= "stroke_thickness=$var_thickness,";
			$var_string			.= "fill_color=$var_colour";

			push @{$array_out[$i]}, $var_string;

		}

	}

	# ----------------------------------

	# Write data from @array_out to file:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my @array_temp	= @{$array_out[$i]};

		unshift @array_temp, "# $file_out.$i";

		# Print @array_out to $file_out:
		General::ArrayToFile ( \@array_temp, "$file_out.$i" );

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	JunctionInserts
# Description:		This subroutine sequesters insert regions that sit in within
#					a junction of backbone regions into a separate array.

# ----------------------------------------------------------

sub JunctionInserts {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my @array_inserts		= @{$hash_out{array_inserts}};
	my %hash_junctions;
	my @array_junctions;

	# ----------------------------------

	# Identify $var_refIDs associated with junctions:
	for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

		# Define variable to hold previous reference ID:
		my $var_refID0		= $array_backbone[$i][-1][3];

		# Define variable to hold previous alignment inverted status:
		my $var_inverted0	= $array_backbone[$i][-1][7];

		# Iterate down through regions in @array_backbone:
		for ( my $j = 0; $j < scalar @{$array_backbone[$i]}; $j++ ) {

			# Define current alignment inverted status:
			my $var_inverted1	= $array_backbone[$i][$j][7];

			# Define current reference ID:
			my $var_refID1		= $array_backbone[$i][$j][3];

			# If current alignment and previous alignment don't agree on
			# inverted status store reference ID in %hash_junctionIDs:
			if ( $var_inverted0 != $var_inverted1 ) {

				# Store current and previous reference IDs in %hash_junctions:
				$hash_junctions { "$var_refID0 $var_refID1" } = 1;

			}

			# Store current alignments inverted status for next iteration:
			$var_inverted0		= $var_inverted1;

			# Store current reference ID for next iteration:
			$var_refID0			= $var_refID1;

		}

	}

	# ----------------------------------

	# Sequester inserts in @array_inserts associated with junctions and store
	# them in @array_junctions:
	for ( my $i = 0; $i < scalar @array_inserts; $i++ ) {

		# Iterate through regions in alignment:
		for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

			# Define variables to hold bracketing reference IDs:
			my $var_bracket0	= $array_inserts[$i][$j][9];
			my $var_bracket1	= $array_inserts[$i][$j][10];

			next if $var_bracket0	== 0;

			my $var_tag			= LowerHigher ( $var_bracket0, $var_bracket1 );

			if ( $hash_junctions {"$var_bracket0 $var_bracket1"} ) {

				push @array_junctions, $array_inserts[$i];

				splice @array_inserts, $i, 1;

				$i--;

				last;

			}

		}

	}

	# ----------------------------------

	# Store new version of @array_inserts in %hash_out:
	$hash_out { array_inserts }		= \@array_inserts;

	# Store @array_junctions in %hash_out:
	$hash_out { array_junctions }	= \@array_junctions;

	# Store %hash_junctionIDs in %hash_out:
	$hash_out { hash_junctions }	= \%hash_junctions;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Karyotype
# Description:		This subroutine produces a formatted Circos karyotype file.

# ----------------------------------------------------------

sub Karyotype {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @{$hash_in { array_merged }};
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { "file_karyotype" };
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 1 );

	# ----------------------------------

	my $var_scalar		= scalar @{$array_in[0]} - 1;
	my $var_length		= $array_in[0][$var_scalar][12];

	my $var_string	= "chr - chr0 0 0 $var_length $var_colour\n";

	push @array_out, $var_string;

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name: LowerHigher
# Description:		This subroutine takes two values and returns a string with
#					this lowere value first.

# ----------------------------------------------------------

sub LowerHigher {

	# Arguments:
	my ( $var_in1, $var_in2 )	= @_;

	# Variables:
	my $var_return;

	# ----------------------------------

	$var_return	= "$var_in1 $var_in2" if $var_in2 > $var_in1;

	$var_return	= "$var_in2 $var_in1" if $var_in1 > $var_in2;

	# ----------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	MergeSVG
# Description:		This subroutine merges various SVG data together and prints
#					it all to file.

# ----------------------------------------------------------

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
	my $var_height			= 100 + ( $var_increment * $var_factor * $var_scalar );

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

# Subroutine name:	MergeRegions
# Description:		This subroutine coordinates merging backbone, inserts, and
#					transposed regions together.

# ----------------------------------------------------------

sub MergeRegions {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my %hash_inserts		= %{$hash_out{hash_inserts}};
	my @array_merged;

	# ----------------------------------

	# Iterate through @array_backbone adding insert regions from %hash_inserts
	# and %hash_junctions:
	for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

		# Iterate down through regions:
		for ( my $j = 0; $j < scalar @{$array_backbone[$i]}; $j++ ) {

			# Add region from @array_backbone to @array_merged:
			push @{$array_merged[$i]}, $array_backbone[$i][$j];

			# Define bracketing reference IDs:
			my $var_bracket0	= $array_backbone[$i][$j][3];
			my $var_bracket1	= $array_backbone[$i][$j][10];

			my $var_tag	= LowerHigher ( $var_bracket0, $var_bracket1 );

			# If $var_tag is present in %hash_inserts, add the regions:
			if ( $hash_inserts { $var_tag } ) {

				# Make a temporary array holding regions:
				my @array_temp		= @{$hash_inserts{$var_tag}};

				next unless @array_temp;

				# Iterate through regions for specific sequence at hand:
				for ( my $k	= 0; $k < scalar @{$array_temp[$i]}; $k++ ) {

					# Add region from @array_backbone to @array_merged:
					push @{$array_merged[$i]}, $array_temp[$i][$k];

				}

			}

		}

	}

	# ----------------------------------

	# Store @array_merged in %hash_out:
	$hash_out { array_merged }		= \@array_merged;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	OrganiseInserts
# Description:		This subroutine organises different inserts together for
#					merging with backbone regions.

# ----------------------------------------------------------

sub OrganiseInserts {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_backbone		= @{$hash_out{array_backbone}};
	my @array_trans_inserts	= @{$hash_out{array_trans_inserts}};
	my %hash_junctionIDs	= %{$hash_out{hash_junctions}};

	# Variables:
	my $var_scalar			= scalar @array_backbone;

	# ----------------------------------

	# Organise insert regions into a hash based upon bracketing reference IDs:
	my %hash_inserts		= %{Synteny::BracketHash($hash_out{array_inserts})};

	# Organise transposed regions into a hash based upon bracketing reference IDs:
	my %hash_transposed		= %{Synteny::BracketHash($hash_out{array_transposed})};

	# Organise insert regions present in junction regions into a hash based upon
	# bracketing reference IDs:
	my %hash_junctions		= %{Synteny::BracketHash($hash_out{array_junctions})};

	# Organise transposed insert regions into a hash:
	my %hash_trans_inserts	= %{Synteny::BracketHash($hash_out{array_trans_inserts})};

	# ----------------------------------

	# Iterate through %hash_inserts organising inserts:
	foreach my $var_bracketIDs ( keys %hash_inserts ) {

		# Define array holding insert regions:
		my @array_regions	= @{$hash_inserts{$var_bracketIDs}};

		# Check if there are any regions in %hash_transposed:
		print "\t$var_bracketIDs\n" if $hash_transposed { $var_bracketIDs };

		# Check if there are any regions in %hash_junctions:
		print "\t$var_bracketIDs\n" if $hash_junctions { $var_bracketIDs };

		# Check if there are any regions in %hash_trans_inserts and add them
		# to @array_regions:
		if ( $hash_trans_inserts { $var_bracketIDs } ) {

			my @array_regions2	= @{$hash_trans_inserts{$var_bracketIDs}};

			delete $hash_trans_inserts { $var_bracketIDs };

			push (@array_regions, @array_regions2);

		}

		# Organise insert regions:
		@array_regions		= @{OrganiseRegions ( \@array_regions, $var_scalar )};

		# Store organised inserts back in %hash_inserts:
		$hash_inserts { $var_bracketIDs }	= \@array_regions;

	}

	# ----------------------------------

	# Iterate through %hash_transposed organising inserts:
	foreach my $var_bracketIDs ( keys %hash_transposed ) {

		# Define array holding insert regions:
		my @array_regions	= @{$hash_transposed{$var_bracketIDs}};

		# Check if there are any regions in %hash_transposed:
		print "\t$var_bracketIDs\n" if $hash_inserts{ $var_bracketIDs };

		# Check if there are any regions in %hash_junctions:
		print "\t$var_bracketIDs\n" if $hash_junctions { $var_bracketIDs };

		# Check if there are any regions in %hash_trans_inserts:
		print "\t$var_bracketIDs\n" if $hash_trans_inserts { $var_bracketIDs };

		# Organise insert regions:
		@array_regions		= @{OrganiseRegions ( \@array_regions, $var_scalar )};

		# Store organised inserts back in %hash_inserts:
		$hash_inserts { $var_bracketIDs }	= \@array_regions;

	}

	# ----------------------------------

	# Iterate through @array_backbone to identify groups of bracketing reference
	# IDs that represent each junction:
	for ( my $i = 0; $i < scalar @{$array_backbone[0]}; $i++ ) {

		# Define bracketing reference IDs:
		my $var_bracket0	= $array_backbone[0][$i][3];
		my $var_bracket1	= $array_backbone[0][$i][10];

		# Ensure lower reference ID comes first:
		my $var_tag			= Synteny::LowerHigher ( $var_bracket0, $var_bracket1 );

		# If bracketing reference IDs are present in %hash_junctionIDs, iterate
		# across sequences and gather all relevant regions:
		if ( $hash_junctionIDs { $var_tag } ) {

			# Define array to hold regions:
			my @array_regions;

			# Iterate down through @array_backbone:
			for ( my $j = 0; $j < scalar @array_backbone; $j++ ) {

				# Define bracketing reference IDs in the region at hand:
				my $var_bracket0	= $array_backbone[$j][$i][3];
				my $var_bracket1	= $array_backbone[$j][$i][10];

				# Ensure lower reference ID comes first:
				my $var_tag			= Synteny::LowerHigher ( $var_bracket0, $var_bracket1 );

				# If $var_tag is present in %hash_junctions add the data to
				# @array_regions:
				if ( $hash_junctions { $var_tag } ) {

					my @array_temp	= @{$hash_junctions{ $var_tag }};

					# Remove the regions from %hash_junctions:
					delete $hash_junctions { $var_tag };

					# Merge @array_regions and the temporary array:
					push ( @array_regions, @array_temp );

				}

			}

			# Organise regions:
			@array_regions	= @{OrganiseRegions ( \@array_regions, $var_scalar )};

			# Store regions for merging:
			for ( my $j = 0; $j < scalar @array_backbone; $j++ ) {

				my $var_bracket0	= $array_backbone[$j][$i][3];
				my $var_bracket1	= $array_backbone[$j][$i][10];

				my $var_tag			= Synteny::LowerHigher ( $var_bracket0, $var_bracket1 );

				$hash_inserts { $var_tag } = \@array_regions;

			}

		}

		$var_bracket0	= $var_bracket1;

	}

	# ----------------------------------

	# Iterate through any remaining regions in %hash_trans_inserts:
	foreach my $var_bracketIDs ( keys %hash_trans_inserts ) {

		# Define array holding insert regions:
		my @array_regions	= @{$hash_trans_inserts{$var_bracketIDs}};

		# Organise insert regions:
		@array_regions		= @{OrganiseRegions ( \@array_regions, $var_scalar )};

		# Store organised inserts back in %hash_inserts:
		$hash_inserts { $var_bracketIDs }	= \@array_regions;

	}

	# ----------------------------------

	# Store %hash_inserts in %hash_out:
	$hash_out { hash_inserts }	= \%hash_inserts;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	OrganiseRegions
# Description:		This subroutine takes an array containing common regions
#					and organises them for merging.

# ----------------------------------------------------------

sub OrganiseRegions {

	# Arguments:
	my ( $array_in, $var_scalar )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;
	my %hash_out;
	my %hash_regions;
	my @array_merged;
	my @array_final;

	# Variables:
	my $var_lowest	= $var_scalar;

	# ----------------------------------

	# Sort regions in @array_in by location:
	@array_in	= sort { $a -> [1] <=> $b -> [1] } @array_in;

	# Iterate through @array_in adding each regions reference IDs to @array_out 
	# organised by sequence and adding the region to %hash_regions organised
	# by sequence and reference IDs:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define sequene and reference IDs of region:
		my $var_seqID	= $array_in[$i][0];
		my $var_refID	= $array_in[$i][3];

		# Define the lowest sequence number used:
		$var_lowest	= $var_seqID if $var_seqID < $var_lowest;

		# Add regions to @array_out:
		push @{$array_out[$var_seqID]}, $var_refID;

		$hash_regions { "$var_seqID $var_refID" }	= $array_in[$i];

	}

	# Redefine @array_in to @array_out and undefine @array_out for further use:
	@array_in	= @array_out;
	undef @array_out;

	# ----------------------------------

	# Iterate through @array_in and perform sequence alignment using the lowest
	# sequence numbers and others:
	for ( my $i = $var_lowest+1; $i < scalar @array_in; $i++ ) {

		# Move on if sequence is empty:
		next unless $array_in[$i];

		my ( $array_align0, $array_align1 )	= Alignment ( $array_in[$var_lowest], $array_in[$i] );

		# Add guide alignment to @array_out:
		push @array_out, $array_align1;

		# Add alignment to be merged to @array_merged:
		push @array_merged, $array_align0;

	}

	# ----------------------------------

	# Move the first alignments from @array_out and @array_merged to
	# @array_final:
	$array_final[$var_lowest]	= $array_out[0];
	$array_final[$var_lowest+1]	= $array_merged[0];

	# Add first group of dashes to @array_out based on their presence in
	# the first guide already in @array_final:
	if ( $array_final[$var_lowest] ) {

		for ( my $i = 0; $i < scalar @{$array_final[$var_lowest]}; $i++ ) {

			my $var_chr0	= $array_final[$var_lowest][$i];

			if ( $var_chr0 eq "-" ) {

				for ( my $j = 1; $j < scalar @array_merged; $j++ ) {

					splice @{$array_merged[$j]}, $i, 0, "-";

				}

			}

		}

	}

	# Iterate through @array_out adding dashes as needed to complete alignment:
	for ( my $i = 1; $i < scalar @array_out; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_out[$i]}; $j++ ) {

			# Define the IDs of regions in @array_out and @array_merged:
			my $var_chr0	= $array_out[$i][$j];
			my $var_chr1	= $array_merged[$i][$j];

			# Add dashes:
			if ( $var_chr0 eq "-" && $var_chr1 ne "-" ) {

				for ( my $k = 0; $k < scalar @array_final; $k++ ) {

					splice @{$array_final[$k]}, $j, 0, "-";

				}

				for ( my $k = $i+1; $k < scalar @array_out; $k++ ) {

					splice @{$array_out[$k]}, $j, 0, "-";

				}

			}

		}

		# Add the finalised alignment to @array_final:
		push @array_final, $array_merged[$i];

	}

	# Undef @array_out:
	undef @array_out;

	# ----------------------------------

	# Add gaps to @array_final:
	for ( my $i = 0; $i < $var_scalar; $i++ ) {

		unless ( $array_in[$i] ) {

			next unless $array_final[$i];

			splice @array_final, $i, 0, "";

		}

	}

	# ----------------------------------

	# Iterate through @array_final and add regions to @array_out:
	for ( my $i = 0; $i < scalar @array_final; $i++ ) {

		next unless $array_final[$i];

		my $var_refID	= $array_final[$i];

		for ( my $j = 0; $j < scalar @{$array_final[$i]}; $j++ ) {

			my $var_refID		= $array_final[$i][$j];

			unless ( $var_refID eq "-" ) {

				my @array_region	= @{ $hash_regions { "$i $var_refID" } };

				$array_out[$i][$j]	= \@array_region;

				delete $hash_regions { "$i $var_refID" };

			}

		}

	}

	# ----------------------------------

	# If there are regions in %hash_regions, the number of regions to be aligned
	# is too simple and can just be added to @array_out:
	if ( scalar keys %hash_regions > 0 ) {

		# Add regions remaining in %hash_regions to @array_out:
		foreach my $var_regions ( keys %hash_regions ) {

				my @array_temp	= @{$hash_regions {$var_regions }};

				my $var_seqID	= $array_temp[0];

				push @{$array_out[$var_seqID]}, \@array_temp;

		}

		# Sort regions in each sequence by location.
		for ( my $i = 0; $i < scalar @array_out; $i++ ) {

			next unless $array_out[$i];

			@{$array_out[$i]}	= sort { $a -> [1] <=> $b -> [1] } @{$array_out[$i]};

		}

	}

	# ----------------------------------

	# Add blank regions to @array_out:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		next unless $array_out[$i];

		for ( my $j = 0; $j < scalar @{$array_out[$i]}; $j++ ) {

			next unless $array_out[$i][$j];

			my $var_length0	= 0;

			for ( my $k = 0; $k < scalar @array_out; $k++ ) {

				next unless $array_out[$k][$j];

				my $var_length1	= $array_out[$k][$j][5];

				$var_length0 = $var_length1 if $var_length1 > $var_length0;

			}

			for ( my $k = 0; $k < $var_scalar; $k++ ) {

				unless ( $array_out[$k][$j] ) {

					my @array_temp	= (0)x13;

					$array_temp[0]	= $k;

					$array_temp[5]	= $var_length0;

					$array_out[$k][$j]	= \@array_temp;

				}

				else {

					$array_out[$k][$j][5]	= $var_length0;

				}

			}

		}

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	OrganiseBackbone
# Description:		This subroutine 

# ----------------------------------------------------------

sub OrganiseBackbone {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_backbone	= @{$hash_out{array_backbone}};
	my @array_out;
	my @array_IDs;
	my @array_transposed;
	my %hash_IDs;

	# ----------------------------------

	# Sort regions by sequence and by location:
	@array_out	= @{SortBySequence ( \@array_backbone )};

	# Store reference IDs from each region as elements in @array_IDs:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_out[$i]}; $j++ ) {

			my $var_refID	= $array_out[$i][$j][3];

			push @{$array_IDs[$i]}, $var_refID;

		}

	}

	# ----------------------------------

	# Perform gapped alignment:
	for ( my $i = 1; $i < scalar @array_IDs; $i++ ) {

		%hash_IDs	= %{GappedAlignment( \%hash_IDs, $array_IDs[0], $array_IDs[$i] )};

	}

	# ----------------------------------

	# Remove transposed regions from @array_backbone:
	for ( my $i = 0; $i < scalar @array_backbone; $i++ ) {

		my $var_refID	= $array_backbone[$i][0][3];

		if ( $hash_IDs { $var_refID } ) {

			push @array_transposed, $array_backbone[$i];

			splice @array_backbone, $i, 1;

			$i--;

		}

	}

	# ----------------------------------

	# Add bracket IDs to @array_backbone:
	@array_backbone	= @{BracketIDs( \@array_backbone )};

	# ----------------------------------

	# Store @array_backbone in %hash_out:
	$hash_out { array_backbone }	= \@array_backbone;

	# Store @array_transposed in %hash_out:
	$hash_out { array_transposed }	= \@array_transposed;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Overlap
# Description:		This subroutine briefly calculates overlap between genomic
# 					features in each sequence and conserved/non-conserved
#					regions.

# ----------------------------------------------------------

sub Overlap {

	# Arguments:
	my ( $array_features, $array_FASTA, $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_features	= @$array_features;
	my @array_FASTA		= @$array_FASTA;
	my @array_out;

	# Variables:
	my $var_scalar	= scalar @array_merged;

	# ----------------------------------

	# Iterate through each collection of genomic features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		# Define temporary hashes to store data:
		my %hash_temp0;
		my %hash_temp1;

		# Iterate down through features:
		for ( my $j = 0; $j < scalar @{$array_features[$i]}; $j++ ) {

			my $var_type	= $array_features[$i][$j][1];
			my $var_loc0	= $array_features[$i][$j][3];
			my $var_loc1	= $array_features[$i][$j][4];

			# Iterate down through @array_merged identifying and tallying
			# overlapping regions:
			for ( my $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				my $var_loc2	= $array_merged[$i][$k][1];
				my $var_loc3	= $array_merged[$i][$k][2];
				my $var_copy	= $array_merged[$i][$k][4];

				next if $var_loc2 == 0;

				my $var_length	= 0;

				if ( $var_loc2 < $var_loc0 && $var_loc3 > $var_loc0 && $var_loc3 < $var_loc1 ) {

					$var_length	= $var_loc3 - $var_loc0 + 1;

				}

				if ( $var_loc2 > $var_loc0 && $var_loc2 < $var_loc1 && $var_loc3 > $var_loc1 ) {

					$var_length	= $var_loc1 - $var_loc2 + 1;

				}

				if ( $var_loc2 < $var_loc0 && $var_loc3 > $var_loc1 ) {

					$var_length	= $var_loc1 - $var_loc0 + 1;

				}

				if ( $var_loc2 > $var_loc0 && $var_loc3 < $var_loc1 ) {

					$var_length	= $var_loc3 - $var_loc2 + 1;

				}

				# Add conserved regions to %hash_temp0:
				if ( $var_copy == $var_scalar ) {

					$hash_temp0 { $var_type }	+= $var_length;

				}

				# Add non-conserved regions to %hash_temp1:
				else {

					$hash_temp1 { $var_type }	+= $var_length;

				}

				last if $var_loc2 > $var_loc0;

			}

		}

		# Store temporary hashes in @array_out:
		$array_out[$i][0]	= \%hash_temp0;
		$array_out[$i][1]	= \%hash_temp1;

	}

	# ----------------------------------

	# Iterate through @array_out reporting data to user:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		# Define temporary array holding regions:
		my @array_regions	= @{$array_merged[$i]};

		my $var_seq_length	= $array_regions[$#array_regions][2];

		print "$i: $array_FASTA[$i] - $var_seq_length bp\n";

		# Define temporary arrays holding length data:
		my %hash_temp0	= %{$array_out[$i][0]};
		my %hash_temp1	= %{$array_out[$i][1]};

		print "    Conserved:\n";

		# Print conserved data:
		for my $var_temp ( keys %hash_temp0 ) {

			my $var_length	= $hash_temp0{$var_temp};

			print "\t* $var_temp -\t$var_length bp\n";

		}

		print "    Non-conserved:\n";

		# Print non-conserved data:
		for my $var_temp ( keys %hash_temp1 ) {

			my $var_length	= $hash_temp0{$var_temp};

			print "\t* $var_temp -\t$var_length bp\n";

		}

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	progressiveMauve
# Description:		This subroutine handles submitting and running
#					progressiveMauve.

# ----------------------------------------------------------

sub progressiveMauve {

	# Arguments:
	my ( $hash_in, $array_FASTA )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_FASTA		= @$array_FASTA;

	# Variables:
	my $file_align		= $hash_in { file_align };
	my $file_backbone	= $hash_in { file_backbone };
	my $file_guide		= $hash_in { file_guide };
	my $var_help		= $hash_in { var_help };
	my $var_return;

	# ----------------------------------

	# Print details of running ProgressiveMauve to user:
	print "progressiveMauve \\\\
	--output=$file_align \\\\
	--output-guide-tree=$file_guide \\\\
	--backbone-output=$file_backbone \\\\
	@array_FASTA\n\n";

	# ----------------------------------

	$var_return	= system ( "progressiveMauve \\
	--output=$file_align \\
	--output-guide-tree=$file_guide \\
	--backbone-output=$file_backbone \\
	@array_FASTA" );

	# ----------------------------------

	# Remove .sslist files:
	foreach my $var_file ( @array_FASTA ) {

		$var_file	.= ".sslist";

		unlink $var_file;

	}

	# ----------------------------------

	# Print help and die if $var_return is non-zero:
	if ( $var_return ) {

		print "$var_help progressiveMauve failed to finish successfully!\n" and die;

		die;

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	RegionsSVG
# Description:		This subroutine formats region data for SVG.

# ----------------------------------------------------------

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
	my $var_increment	= $hash_out{var_increment} * $var_factor;
	my $var_scalar		= scalar @array_merged;
	my $var_stroke		= "0,0,0";

	# ----------------------------------

	# Iterate through @array_merged:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_copy		= $array_merged[$i][$j][4];
			my $var_length		= $array_merged[$i][$j][5];
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

# Subroutine name:	SortBySequence
# Description:		This subroutine sorts regions by sequence.

# ----------------------------------------------------------

sub SortBySequence {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	# Sort regions by sequence:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			push @{$array_out[$j]}, $array_in[$i][$j];

		}

	}

	# Sort regions in each sequence:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		@{$array_out[$i]}	= sort { $a -> [1] <=> $b -> [1] } @{$array_out[$i]};

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	SortByReference
# Description:		This subroutine sorts regions by their reference ID.

# ----------------------------------------------------------

sub SortByReference {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my %hash_out;
	my @array_out;

	# ----------------------------------

	# Iterate through each sequence in turn:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Iterate through regions in sequence:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define reference ID of region:
			my $var_refID	= $array_in[$i][$j][3];

			# Push region to %hash_out based upon reference ID:
			push @{ $hash_out{ $var_refID } }, $array_in[$i][$j];

		}

	}

	# ----------------------------------

	# Iterate through %hash_out moving alignments to @array_out:
	foreach my $var_refID ( keys %hash_out ) {

		# Create temporary array to hold alignment:
		my @array_temp	= @{$hash_out{$var_refID}};

		# Add alignment to @array_out:
		push @array_out, \@array_temp;

	}

	# Sort @array_out by locations in first region of each alignment:
	@array_out	= sort { $a -> [1] <=> $b -> [1] } @array_out;

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	StoreFinal
# Description:		This subroutine stores the original and adjusted locations
#					for each sequence to file:

# ----------------------------------------------------------

sub StoreFinal {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_merged	= @{$hash_out{ array_merged }};

	# Variables:
	my $file_out		= $hash_out { file_final };

	# ----------------------------------

	# Iterate through each sequence in @array_merged:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		# Define array to hold original and adjusted locations:
		my @array_out;

		# Iterate down through locations:
		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			# Define original and adjusted locaitons:
			my $var_loc0	= $array_merged[$i][$j][1];
			my $var_loc1	= $array_merged[$i][$j][2];
			my $var_loc2	= $array_merged[$i][$j][11];
			my $var_loc3	= $array_merged[$i][$j][12];

			# Move on if regions is a blank:
			next unless $var_loc0;

			# Define temporary array holding locations:
			my $var_string	= "$var_loc0 $var_loc1 $var_loc2 $var_loc3";

			# Store temporary array in @array_out:
			push @array_out, $var_string;

		}

		# Write @array_out to $file_out:
		General::ArrayToFile ( \@array_out, "$file_out.$i" );

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	TicksSVG
# Description:		This file creates ticks formatted for SVG.

# ----------------------------------------------------------

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
	my $var_increment	= $hash_out{var_increment} * $var_factor;

	my $var_gap			= $hash_out{var_gap};
	my $var_y			= 100 - $var_gap;

	my $var_scalar		= scalar @array_merged;
	my $var_stroke		= "0,0,0";
	my $var_colour		= "0,0,0";

	$var_scalar			= scalar @{$array_merged[0]};
	my $var_length		= $array_merged[0][$var_scalar-1][2];

	my $var_1Mb_height		= $hash_out{var_1Mb_height}/2;
	my $var_100kb_height	= $hash_out{var_100kb_height}/2;
	my $var_10kb_height		= $hash_out{var_10kb_height}/2;
	my $r					= "r";
	my $p					= "p";
	my $var_1Mb				= 10**6;
	my $var_100kb			= 10**5;
	my $var_10kb			= 10**4;

	# ----------------------------------

	# Create horizontal lines:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_copy		= $array_merged[$i][$j][4];
			my $var_length		= $array_merged[$i][$j][5];
			my $var_inverted	= $array_merged[$i][$j][8];
			my $var_loc0		= $array_merged[$i][$j][11];
			my $var_loc1		= $array_merged[$i][$j][12];

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

		for ( my $j = 0; $j <= $var_length; $j+=$var_10kb ) {

			my $k				= 0;
			my $var_diff		= 0;
			my $var_inverted	= 0;

			for ( $k = 0; $k < scalar @{$array_merged[$i]}; $k++ ) {

				my $var_loc0	= $array_merged[$i][$k][1];
				my $var_loc1	= $array_merged[$i][$k][2];
			 	$var_inverted	= $array_merged[$i][$k][8];
				my $var_loc2	= $array_merged[$i][$k][11];

				next if $var_loc0 == 0;

				if ( $var_loc0 < $j && $var_loc1 > $j ) {

					$var_diff	= $var_loc2 - $var_loc0;

					last;

				}

			}

			$var_y += $var_height if $var_inverted;

			if ( $j % $var_1Mb == 0 ) {

				my $var_x	= ( $j + $var_diff )/$var_factor;

				$var_y	-= $var_1Mb_height;

				my $var_string	= General::RectangleSVG ( $var_x, $var_y, $var_thick, $var_1Mb_height,
												 $var_colour, $var_stroke );

				push @{$array_out[$i]}, $var_string;

				$var_y	+= $var_1Mb_height;

			}

			elsif ( $j % $var_100kb == 0 ) {

				my $var_x	= ( $j + $var_diff )/$var_factor;

				$var_y	-= $var_100kb_height;

				my $var_string	= General::RectangleSVG ( $var_x, $var_y, $var_thick, $var_100kb_height,
												 $var_colour, $var_stroke );

				push @{$array_out[$i]}, $var_string;

				$var_y	+= $var_100kb_height;

			}

			else {

				my $var_x	= ( $j + $var_diff )/$var_factor;

				$var_y	-= $var_10kb_height;

				my $var_string	= General::RectangleSVG ( $var_x, $var_y, $var_thick, $var_10kb_height,
												 $var_colour, $var_stroke );

				push @{$array_out[$i]}, $var_string;

				$var_y	+= $var_10kb_height;

			}

			$var_y -= $var_height if $var_inverted;

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

# Subroutine name:	TransposedInserts
# Description:		This subroutine identifies inserts bracketed by different
#					backbone regions and inserts that belong with transposed
#					backbone regions.

# ----------------------------------------------------------

sub TransposedInserts {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_transposed	= @{$hash_out{array_transposed}};
	my @array_inserts		= @{$hash_out{array_inserts}};
	my %hash_transposed;
	my @array_trans_inserts;

	# ----------------------------------

	# Create a hash indexing bracketing reference numbers within 
	# @array_transposed:
	for ( my $i = 0; $i < scalar @array_transposed; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_transposed[$i]}; $j++ ) {

			my $var_bracket0	= $array_transposed[$i][$j][9];
			my $var_bracket1	= $array_transposed[$i][$j][10];

			next if $var_bracket0 == 0;

			my $var_tag			= Synteny::LowerHigher ( $var_bracket0, $var_bracket1 );

			$hash_transposed { $var_tag }	= 1;

		}

	}

	# ----------------------------------

	# Iterate through @array_inserts finding alignments associated with
	# transposed regions:
	for ( my $i = 0; $i < scalar @array_inserts; $i++ ) {

		for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

			my $var_bracket0	= $array_inserts[$i][$j][9];
			my $var_bracket1	= $array_inserts[$i][$j][10];

			next if $var_bracket0	== 0;

			my $var_tag			= Synteny::LowerHigher ( $var_bracket0, $var_bracket1 );

			# If the bracket IDs are the same as those stored in %hash_transposed,
			# store the alignment in @array_transposed:
			if ( $hash_transposed { $var_tag } ) {

				#  Create new alignments for each region:
				for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

					my $var_seqID	= $array_inserts[$i][$j][0];
					my $var_loc0	= $array_inserts[$i][$j][1];
					my $var_length	= $array_inserts[$i][$j][5];

					next if $var_loc0 == 0;

					my @array_temp2;

					for ( my $k = 0; $k < scalar @{$array_inserts[$i]}; $k++ ) {

						my @array_temp	= (0)x13;
						$array_temp[0]	= $k;
						$array_temp[5]	= $var_length;

						push @array_temp2, \@array_temp;

					}

					$array_temp2[$var_seqID]	= $array_inserts[$i][$j];

					push @array_transposed, \@array_temp2;

				}

				# Remove the alignment from @array_inserts:
				splice @array_inserts, $i, 1;

				$i--;

			}

		}

	}

	# ----------------------------------

	# Iterate through @array_inserts finding alignments with more than two
	# bracketing reference IDs:
	for ( my $i = 0; $i < scalar @array_inserts; $i++ ) {

		my %hash_bracketIDs;

		for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

			my $var_bracket0	= $array_inserts[$i][$j][9];
			my $var_bracket1	= $array_inserts[$i][$j][10];

			next if $var_bracket0	== 0;

			# Store bracketing reference IDs:
			$hash_bracketIDs { $var_bracket0 }	= 1;
			$hash_bracketIDs { $var_bracket1 }	= 1;

		}

		# If there is more than 2 keys in %hash_bracketIDs add the regions
		# to @array_trans_inserts and remove it from @array_inserts:
		if ( scalar keys %hash_bracketIDs > 2 ) {

			for ( my $j = 0; $j < scalar @{$array_inserts[$i]}; $j++ ) {

				my $var_seqID	= $array_inserts[$i][$j][0];
				my $var_loc0	= $array_inserts[$i][$j][1];
				my $var_length	= $array_inserts[$i][$j][5];

				next if $var_loc0 == 0;

				my @array_temp2;

				for ( my $k = 0; $k < scalar @{$array_inserts[$i]}; $k++ ) {

					my @array_temp	= (0)x13;
					$array_temp[0]	= $k;
					$array_temp[5]	= $var_length;

					push @array_temp2, \@array_temp;

				}

				$array_temp2[$var_seqID]	= $array_inserts[$i][$j];

				push @array_trans_inserts, \@array_temp2;

			}

			splice @array_inserts, $i, 1;

			$i--;

		}

	}

	# ----------------------------------

	# Store @array_transposed in %hash_out:
	$hash_out { array_transposed }	= \@array_transposed;

	# Store @array_trans_inserts in %hash_out:
	$hash_out { array_trans_inserts }	= \@array_trans_inserts;

	# Store updates @array_inserts in %hash_out:
	$hash_out { array_inserts }	= \@array_inserts;

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Ticks
# Description:		This subroutine takes region data and creates tick marks.

# ----------------------------------------------------------

sub Ticks {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_in		= @{$hash_out{ array_merged }};
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out		= $hash_out { file_ticks };
	my $var_thickness	= $hash_out { stroke_thickness };
	my $var_span0		= $hash_out { var_span };
	my $var_colour		= $hash_colours { ticks };

	my $var_gap			= $hash_out { var_gap };

	my $var_radius0		= $hash_out { var_radius0 };

	my $var_1Mb_height		= $hash_out{var_1Mb_height};
	my $var_100kb_height	= $hash_out{var_100kb_height};
	my $var_10kb_height		= $hash_out{var_10kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_increment		= $hash_out { var_increment };
	my $var_span			= $hash_out { var_span } * 2;

	# ----------------------------------

	# Iterate through @array_in adding horizontal and vertical marks for ticks:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variables for each size marker:
		my $var_1Mb		= 10**6;
		my $var_100kb	= 10**5;
		my $var_10kb	= 10**4;

		# Add first tick:
		my $var_string		= "chr0 1 1 ";
		$var_string			.= "stroke_color=$var_colour,";
		$var_string			.= "stroke_thickness=$var_thickness,";

		unless ( scalar @array_in == 2 && $i == 1 ) {

			$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
			$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

		}

		else {

			$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
			$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_1Mb_height$p-$var_gap$p";

		}

		push @{$array_out[$i]}, $var_string;

		# Iterate through regions in @array_in making up vertical ticks:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define original locations:
			my $var_loc0		= $array_in[$i][$j][1];
			my $var_loc1		= $array_in[$i][$j][2];

			# Move on if region is a blank:
			next if $var_loc0 == 0;

			# Define adjusted locations:
			my $var_loc2		= $array_in[$i][$j][11];
			my $var_loc3		= $array_in[$i][$j][12];

			# Define inverted status of region:
			my $var_inverted	= $array_in[$i][$j][8];

			# Define difference between original and adjusted locations:
			my $var_diff	= $var_loc2 - $var_loc0;

			# Handle ticks every 10 kb:
			if ( $var_10kb >= $var_loc0 && $var_10kb < $var_loc1 ) {

				my $var_temp	= $var_10kb + $var_diff;

				$var_10kb += 10**4;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "offset=$var_span0," if $var_inverted;
				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";

				unless ( scalar @array_in == 2 && $i == 1 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_10kb_height+$var_gap$p";

				}

				else {

					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_10kb_height$p-$var_gap$p";

				}

				push @{$array_out[$i]}, $var_string;

				$j--;

				next;

			}		 

			# Handle ticks every 100kb:
			if ( $var_100kb >= $var_loc0 && $var_100kb < $var_loc1 ) {

				my $var_temp	= $var_100kb + $var_diff;

				$var_100kb += 10**5;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "offset=$var_span0," if $var_inverted;
				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";

				unless ( scalar @array_in == 2 && $i == 1 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_100kb_height+$var_gap$p";

				}

				else {

					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_100kb_height$p-$var_gap$p";

				}

				push @{$array_out[$i]}, $var_string;

				$j--;

				next;

			}		 

			# Handle ticks every 1 Mb:
			if ( $var_1Mb >= $var_loc0 && $var_1Mb < $var_loc1 ) {

				my $var_temp	= $var_1Mb + $var_diff;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "offset=$var_span0," if $var_inverted;
				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";

				unless ( scalar @array_in == 2 && $i == 1 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

				}

				else {

					$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_1Mb_height$p-$var_gap$p";

				}

				push @{$array_out[$i]}, $var_string;

				$var_1Mb += 10**6;

				$j--;

			}

		}

		# Handle horizontal lines:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define original locations:
			my $var_loc0		= $array_in[$i][$j][1];
			my $var_loc1		= $array_in[$i][$j][2];

			# Move on if region is a blank:

			# Define adjusted locations:
			my $var_loc2		= $array_in[$i][$j][11];
			my $var_loc3		= $array_in[$i][$j][12];

			# Define inverted status of the region:
			my $var_inverted	= $array_in[$i][$j][8];

			my $var_string		= "chr0 $var_loc2 $var_loc3 ";

			$var_string			.= "offset=$var_span0," if $var_inverted;
			$var_string			.= "stroke_color=$var_colour,";
			$var_string			.= "stroke_thickness=$var_thickness,";

			unless ( scalar @array_in == 2 && $i == 1 ) {

				$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
				$var_string			.= "r1=$var_radius0$r+$var_gap$p";

			}

			else {

				$var_string			.= "r0=$var_radius0$r-$var_span$p-$var_gap$p,";
				$var_string			.= "r1=$var_radius0$r-$var_span$p-$var_gap$p";

			}

#			push @{$array_out[$i]}, $var_string;

		}

	 	$var_radius0	-= $var_increment;

		if ( scalar @array_in == 2 ) {

#			print "$var_span\n";

#			$var_radius0	-= $var_increment;

#			print "$var_radius0\n";

		}

	}

	# ----------------------------------

	# Iterate through @array_out and write data to @array_out:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my @array_temp	= @{$array_out[$i]};

		unshift @array_temp, "# $file_out.$i";

		# Print @array_out to $file_out:
		General::ArrayToFile ( \@array_temp, "$file_out.$i" );

	}

	# ----------------------------------

	# Return reference to %hash_out and end subroutine:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	WriteSeq
# Description:		This subroutine writes conserved and non-conserved
#					nucleotide sequences in FASTA format for each seqeunce.

# ----------------------------------------------------------

sub WriteSeq {

	# Arguments:
	my ( $hash_in, $array_FASTA )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_FASTA		= @$array_FASTA;
	my @array_merged	= @{$hash_out{array_merged}};
	my @array_conserved;
	my @array_nonconserved;

	# Variables:
	my $var_scalar			= scalar @array_merged;
	my $file_conserved		= $hash_out { file_conserved };
	my $file_nonconserved	= $hash_out { file_nonconserved };

	# ----------------------------------

	# Iterate down through merged synteny data:
	for ( my $i = 0; $i < scalar @array_merged; $i++ ) {

		# Define the FASTA sequence:
		my $var_seq	= $array_FASTA[$i][2];

		# Define the FASTA sequence header:
		my $var_header	= $array_FASTA[$i][1];

		# Add header data to conserved and non-conserved arrays:
		push @{$array_conserved[$i]}, "#$var_header";
		push @{$array_nonconserved[$i]}, "#$var_header";

		# Iterate through data for each sequence:
		for ( my $j = 0; $j < scalar @{$array_merged[$i]}; $j++ ) {

			my $var_loc0	= $array_merged[$i][$j][1];
			my $var_loc1	= $array_merged[$i][$j][2];
			my $var_refID	= $array_merged[$i][$j][3];
			my $var_copy	= $array_merged[$i][$j][4];
			my $var_length	= $array_merged[$i][$j][6];

			next if $var_loc0	== 0;

			my $var_subseq	= substr $var_seq, $var_loc0, $var_length;

			my $var_header	= ">$var_refID - $var_loc0, $var_length  bp, copies: $var_copy";

			# Add line breaks to $var_subseq:
			$var_subseq =~ s/(.{1,70})/$1\n/gs;

			if ( $var_copy == $var_scalar ) {

				push @{$array_conserved[$i]}, $var_header;
				push @{$array_conserved[$i]}, $var_subseq;

			}

			else {

				push @{$array_nonconserved[$i]}, $var_header;
				push @{$array_nonconserved[$i]}, $var_subseq;

			}

		}

	}

	# ----------------------------------

	# Print @array_conserved and @array_nonconserved to file:
	for ( my $i = 0; $i < $var_scalar; $i++ ) {

		General::ArrayToFile ( \@{$array_conserved[$i]}, "$file_conserved.$i.fasta" );
		General::ArrayToFile ( \@{$array_nonconserved[$i]}, "$file_nonconserved.$i.fasta" );

	}

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
