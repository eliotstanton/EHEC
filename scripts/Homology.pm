package Homology;

use strict;
use warnings;
use JSON::XS;

# ---------------------------------------------------------------------------- #

# File name:		Homology.pm
# Date created:		27 October, 2018
# Last modified:	18 November, 2018
# Created by:		Eliot Stanton

# Description:		This package contains subroutines specific to calculating
#					homology.

# ---------------------------------------------------------------------------- #

# Subroutine name:	Config
# Description:		This subroutine creates a config file for Circos

# ----------------------------------------------------------

sub Config {

	# Arguments:
	my ( $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_out;

	# Variables:
	my $file_out			= $hash_in { "file_config" };
	my $file_direct_links	= $hash_in { "file_direct_links" };
	my $file_features		= $hash_in { "file_features" };
	my $file_ideogram		= $hash_in { "file_ideogram" };
	my $file_inverted_links	= $hash_in { "file_inverted_links" };
	my $file_image			= $hash_in { "file_image" };
	my $file_karyotype		= $hash_in { "file_karyotype" };
	my $file_ticks			= $hash_in { "file_ticks" };
	my $var_radius0			= $hash_in { var_radius0 };
	my $var_radius1			= $hash_in { var_radius1 };
	my $var_gap				= $hash_in { var_gap };
	my $r					= "r";
	my $p					= "p";


	# ----------------------------------

	push @array_out, "#$file_out";

	push @array_out, "karyotype = $file_karyotype";
	push @array_out, "<<include $file_ideogram>>";
	push @array_out, "<ideogram>\nshow = no\n</ideogram>";

	push @array_out, "<highlights>";
	push @array_out, "\t<highlight>";
	push @array_out, "\t\tfile	= $file_features\nideogram	= no";
	push @array_out, "\t\tr0	= $var_radius0$r\n\t\tr1	= $var_radius1$r";
 	push @array_out, "\t</highlight>";
	push @array_out, "\t<highlight>";
	push @array_out, "\t\tfile	= $file_ticks\n\t\tideogram	= no";
	push @array_out, "\t</highlight>";
	push @array_out, "</highlights>";

	push @array_out, "<links>";
	push @array_out, "\tradius	= $var_radius1$r-$var_gap$p\nribbon	= yes";
 	push @array_out, "<link>\nfile	= $file_direct_links\n</link>" unless $hash_in{d};
 	push @array_out, "<link>\nfile	= $file_inverted_links\n</link>" unless $hash_in{i};
 	push @array_out, "</links>";

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

# Subroutine:		DirectLinks
# Description:		This subroutine creates an array of formatted direct links.

# ----------------------------------------------------------

sub DirectLinks {

	# Arguments:
	my ( $hash_in, $hash_nmer, $var_nmer )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_nmer		= %$hash_nmer if $hash_nmer;
	my @array_out;

	# Variables:
	my $file_out	= $hash_out { file_direct };

	# ----------------------------------

	# End subroutine early if the output file is already present:
	return if -e $file_out;

	# Iterate through each nmer and process those with more than one location:
	foreach my $var_seq ( keys %hash_nmer ) {

		my @array_tags	= @{$hash_nmer { $var_seq }};

		# Proceed if more than one location:
		if ( scalar @array_tags > 1 ) {

			# Iterate through locations and generate links:
			for ( my $i = 0; $i < scalar @array_tags - 1; $i++ ) {

				my $var_seq0	= $array_tags[$i][0];
				my $var_loc0	= $array_tags[$i][1];
				my $var_loc1	= $array_tags[$i][1] + $var_nmer - 1;

				for ( my $j = $i+1; $j < scalar @array_tags; $j++ ) {

					my $var_seq1	= $array_tags[$j][0];
					my $var_loc2	= $array_tags[$j][1];
					my $var_loc3	= $array_tags[$j][1] + $var_nmer - 1;

					# Define temporary array:
					my @array_temp;

					# Add lower sequence first:
					if ( $var_seq0 != $var_seq1 ) {

						if ( $var_seq0 < $var_seq1 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					else {

						# Add lower location first:
						if ( $var_loc0 < $var_loc2 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					# Store link in @array_out:
					push @array_out, \@array_temp;

				}

			}

		} 

	}

	# ----------------------------------

	# Sort elements in @array_out by location:
	@array_out	= sort { $a -> [0] <=> $b -> [0] 
						|| $a -> [1] <=> $b -> [1] 
						|| $a -> [3] <=> $b -> [3]
						|| $a -> [4] <=> $b -> [4] } @array_out;

	# ----------------------------------

	# Convert arrays in @array_temp into strings:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my $var_string	= join " ", @{$array_out[$i]};

		$array_out[$i]	= $var_string;

	}

	# ----------------------------------

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	my $var_scalar	= scalar @array_out;

	print "\t- $var_scalar direct links\n";

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Features
# Description:		This subroutine handles writing highlight files representing
#					genomic features.

# ----------------------------------------------------------

sub Features {

	# Arguments:
	my ( $array_in, $array_FASTA, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my @array_FASTA		= @$array_FASTA;
	my %hash_in			= %$hash_in;
	my %hash_colours	= %{Circos::Colours( )};
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { file_features };
	my $var_thickness	= $hash_in { stroke_thickness };
	my $var_length		= length $array_FASTA[0][2];
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 2 );
	my $var_stroke		= $hash_colours { "chromosome_border" };

	# ----------------------------------

	# Store file name at beginning of file:
	push @array_out, "#$file_out";

	# Create standing ideogram block(s):
	for ( my $i = 0; $i < scalar @array_FASTA; $i++ ) {

		my $var_length	= length $array_FASTA[$i][2];

		my $var_string	= "chr$i 1 $var_length fill_color=$var_colour";
		$var_string		.= ",stroke_color=$var_stroke";
		$var_string		.= ",stroke_thickness=$var_thickness";

		push @array_out, $var_string;

	}

	# Iterate through features in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_seq		= $array_in[$i][0];
		my $var_type	= $array_in[$i][1];
		my $var_loc0	= $array_in[$i][3];
		my $var_loc1	= $array_in[$i][4];

		my $var_colour	= Circos::ColoursNew ( $var_type, 1, 2 );
		my $var_stroke	= $hash_colours { "$var_type\_border" };

		my $var_string	= "chr$var_seq $var_loc0 $var_loc1 ";
		$var_string		.= "fill_color=$var_colour";
		$var_string		.= ",stroke_color=$var_stroke";
		$var_string		.= ",stroke_thickness=$var_thickness";

		push @array_out, $var_string;

	}

	# ----------------------------------

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:		InvertedLinks
# Description:		This subroutine creates an array of formatted inverted
#					links.

# ----------------------------------------------------------

sub InvertedLinks {

	# Arguments:
	my ( $hash_in, $hash_nmer, $var_nmer )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my %hash_nmer		= %$hash_nmer if $hash_nmer;
	my @array_out;

	# Variables:
	my $file_out	= $hash_out { file_inverted };

	# ----------------------------------

	# End the subroutine early if the output is already present:
	return if -e $file_out;

	# Iterate through each nmer and process those with more than one location:
	foreach my $var_seq ( keys %hash_nmer ) {

		# Generate reverse complement of $var_seq:
		my $var_reverse	= General::ReverseComplement ( $var_seq );

		# If $var_reverse is also present in %hash_in make links:
		if ( $hash_nmer { $var_reverse } ) {

			# Define both arrays of locations:
			my @array_tags1	= @{$hash_nmer{$var_seq}};
			my @array_tags2	= @{$hash_nmer{$var_reverse}};

			# Iterate through arrays making links:
			for ( my $i = 0; $i < scalar @array_tags1; $i++ ) {

				my $var_seq0	= $array_tags1[$i][0];
				my $var_loc0	= $array_tags1[$i][1];
				my $var_loc1	= $array_tags1[$i][1] + $var_nmer - 1;

				for ( my $j = 0; $j < scalar @array_tags2; $j++ ) {

					my $var_seq1	= $array_tags2[$j][0];
					my $var_loc2	= $array_tags2[$j][1];
					my $var_loc3	= $array_tags2[$j][1] + $var_nmer - 1;

					# Define temporary array:
					my @array_temp;

					# Add lower sequence first:
					if ( $var_seq0 != $var_seq1 ) {

						if ( $var_seq0 < $var_seq1 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					else {

						# Add lower location first:
						if ( $var_loc0 < $var_loc2 ) {

							@array_temp = ( $var_seq0, $var_loc0, $var_loc1, 
											$var_seq1, $var_loc2, $var_loc3 );

						}

						else {

							@array_temp = ( $var_seq1, $var_loc2, $var_loc3, 
											$var_seq0, $var_loc0, $var_loc1,);

						}

					}

					# Store link in @array_out:
					push @array_out, \@array_temp;

				}

			}

		}

		# Remove $var_seq from %hash_out:
		delete $hash_nmer { $var_seq };

	}

	# ----------------------------------

	# Sort elements in @array_out by location:
	@array_out	= sort { $a -> [0] <=> $b -> [0] 
						|| $a -> [1] <=> $b -> [1] 
						|| $a -> [3] <=> $b -> [3]
						|| $a -> [4] <=> $b -> [4] } @array_out;

	# ----------------------------------

	# Convert arrays in @array_temp into strings:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		my $var_string	= join " ", @{$array_out[$i]};

		$array_out[$i]	= $var_string;

	}

	# ----------------------------------

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	my $var_scalar	= scalar @array_out;

	print "\t- $var_scalar inverted links\n";

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Karyotype
# Description:		This subroutine produces a formatted Circos karyotype file.

# ----------------------------------------------------------

sub Karyotype {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @$array_in;
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { "file_karyotype" };
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 2 );

	# ----------------------------------

	push @array_out, "#$file_out";

	# Iterate through @array_FASTA writing data to @array_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_header		= $array_in[$i][0];

		my $var_length		= length $array_in[$i][2];

		my $var_string	= "chr - chr$i $i 0 $var_length $var_colour";

		push @array_out, $var_string;

	}

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Links
# Description:		This subroutine produces a file with Circos links data.

# ----------------------------------------------------------

sub Links {

	# Arguments:
	my ( $array_direct, $array_inverted, $array_features, $hash_in )	= @_;

	# Data structures:
	my @array_direct	= @$array_direct;
	my @array_inverted	= @$array_inverted;
	my @array_features	= @$array_features;
	my %hash_in			= %$hash_in;
	my %hash_colours	= %{Circos::Colours( )};
	my @array_direct_links;
	my @array_inverted_links;

	# Variables:
	my $file_direct		= $hash_in { "file_direct_links" };
	my $file_inverted	= $hash_in { "file_inverted_links" };
	my $var_thickness	= $hash_in { stroke_thickness_links };

	# ----------------------------------

	# Store file name at beginning of direct links file:
	push @array_direct_links, "#$file_direct";

	# Iterate through @array_direct writing links to @array_out:
	for ( my $i = 0; $i < scalar @array_direct; $i++ ) {

		my $var_seq0	= $array_direct[$i][0];
		my $var_loc0	= $array_direct[$i][1];
		my $var_loc1	= $array_direct[$i][2];

		my $var_seq1	= $array_direct[$i][3];
		my $var_loc2	= $array_direct[$i][4];
		my $var_loc3	= $array_direct[$i][5];

		my $var_type0	= $array_direct[$i][7];
		my $var_length0	= $array_direct[$i][8];
		my $var_type1	= $array_direct[$i][9];
		my $var_length1	= $array_direct[$i][10];

		my $var_colour;
		my $var_stroke;

		# If both ends of link have the same feature:
		if ( $var_type0 eq $var_type1 ) {

			$var_colour	= Circos::ColoursNew ( $var_type0, 1, 2 );
			$var_stroke	= $hash_colours { "$var_type0\_border" };

		}

		# Otherwise, use the smaller feature:
		else {

			if ( $var_length0 > $var_length1 ) {

				$var_colour	= $hash_colours { $var_type1 };

				$var_colour	= Circos::ColoursNew ( $var_type1, 1, 2 );
				$var_stroke	= $hash_colours { "$var_type1\_border" };

			}

			else {

				$var_colour	= $hash_colours { $var_type0 };

				$var_colour	= Circos::ColoursNew ( $var_type0, 1, 2 );
				$var_stroke	= $hash_colours { "$var_type0\_border" };

			}

		}

		# Create string holding data for Circos:
		my $var_string	= "chr$var_seq0 $var_loc0 $var_loc1 ";
		$var_string		.= "chr$var_seq1 $var_loc2 $var_loc3 ";
		$var_string		.= "color=$var_colour";

		# Push $var_string to @array_direct for writing to file:
		push @array_direct_links, $var_string;

	}

	# ----------------------------------

	# Store file name at beginning of inverted links file:
	push @array_inverted_links, "#$file_inverted";

	# Iterate through @array_inverted writing links to @array_out:
	for ( my $i = 0; $i < scalar @array_inverted; $i++ ) {

		my $var_seq0	= $array_inverted[$i][0];
		my $var_loc0	= $array_inverted[$i][1];
		my $var_loc1	= $array_inverted[$i][2];

		my $var_seq1	= $array_inverted[$i][3];
		my $var_loc2	= $array_inverted[$i][4];
		my $var_loc3	= $array_inverted[$i][5];

		my $var_type0	= $array_inverted[$i][7];
		my $var_length0	= $array_inverted[$i][8];
		my $var_type1	= $array_inverted[$i][9];
		my $var_length1	= $array_inverted[$i][10];

		my $var_colour;
		my $var_stroke;

		# If both ends of link have the same feature:
		if ( $var_type0 eq $var_type1 ) {

			$var_colour	= Circos::ColoursNew ( $var_type0, 1, 2 );
			$var_stroke	= $hash_colours { "$var_type0\_border" };

		}

		# Otherwise, use the smaller feature:
		else {

			if ( $var_length0 > $var_length1 ) {

				$var_colour	= Circos::ColoursNew ( $var_type1, 1, 2 );
				$var_stroke	= $hash_colours { "$var_type1\_border" };

			}

			else {

				$var_colour	= Circos::ColoursNew ( $var_type0, 1, 2 );
				$var_stroke	= $hash_colours { "$var_type0\_border" };

			}

		}

		# Create string holding data for Circos:
		my $var_string	= "chr$var_seq0 $var_loc0 $var_loc1 ";
		$var_string		.= "chr$var_seq1 $var_loc2 $var_loc3 ";
		$var_string		.= "color=$var_colour";

		# Add $var_string to @array_inverted_links for writing to file:
		push @array_inverted_links, $var_string;

	}

	# ----------------------------------

	# Write @array_direct_links to file:
	General::ArrayToFile ( \@array_direct_links, $file_direct );

	# Write @array_inverted_links to file:
	General::ArrayToFile ( \@array_inverted_links, $file_inverted );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:		MergeLinks
# Description:		This subroutine merges adjacent pairs of links.

# ----------------------------------------------------------

sub MergeLinks {

	# Arguments:
	my ( $hash_in, $var_variable, $var_min )	= @_;

	# Data structures:
	my %hash_out	= %$hash_in;
	my @array_in;
	my @array_out;

	# Variables:
	my $file_in;
	my $file_out;
	my $var_orientation	= $var_variable;

	# ----------------------------------

	# Set $var_variable equal to one (direct) or minus one (inverted):
	if ( $var_variable eq "direct" ) {

		$var_variable	= "1";

		$file_in		= $hash_out { file_direct };
		$file_out		= $hash_out { file_direct_merged };

	}

	if ( $var_variable eq "inverted" ) {

		$var_variable = "-1" ;

		$file_in		= $hash_out { file_inverted };
		$file_out		= $hash_out { file_inverted_merged };

	}

	@array_in		= @{General::FileToArray ( $file_in )};

	# ----------------------------------

	my $var_scalar2	= scalar @array_in;

	print "$var_scalar2 unmerged $var_orientation links\n";

	# Sort elements in @array_in by location:
	@array_in	= sort { 

		$a -> [1] <=> $b -> [1] || $a -> [3] <=> $b -> [3] 

	} @array_in;

	# ----------------------------------

	# Iterate down through each link in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_seq0	= $array_in[$i][0];
		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];

		my $var_seq1	= $array_in[$i][3];
		my $var_loc2	= $array_in[$i][4];
		my $var_loc3	= $array_in[$i][5];

		# Iterate downward from link at hand finding adjacent links and merging:
		for ( my $j = $i + 1; $j < scalar @array_in; $j++ ) {

			my $var_seq2	= $array_in[$j][0];
			my $var_loc4	= $array_in[$j][1];
			my $var_loc5	= $array_in[$j][2];

			my $var_seq3	= $array_in[$j][3];
			my $var_loc6	= $array_in[$j][4];
			my $var_loc7	= $array_in[$j][5];

			# End loop if locations are larger than the original:
			last if $var_loc4 > $var_loc0 + 1;

			# End loop is sequence ID varies:
			last if $var_seq0 != $var_seq2;

			# If links are adjacent merge the data:
			if ( $var_loc0 + 1 == $var_loc4 && $var_loc1 + 1 == $var_loc5 ) {

				if ( $var_loc2 + $var_variable == $var_loc6 && $var_loc3 + $var_variable == $var_loc7 ) {

					# Redefine locations in original element:
					$array_in[$i][2]	= $var_loc5;
					$array_in[$i][5]	= $var_loc7 if $var_variable == 1;
					$array_in[$i][4]	= $var_loc6 if $var_variable == -1;

					# Remove adjacent link from @array_in:
					splice @array_in, $j, 1;

					$j--;

					# Update location variables for next iteration:
					$var_loc0	= $var_loc4;
					$var_loc1	= $var_loc5;
					$var_loc2	= $var_loc6;
					$var_loc3	= $var_loc7;

				}

			}

		}

	}

	# ----------------------------------

	# Sort elements in @array_out by location:
	@array_in	= sort { $a -> [0] <=> $b -> [0] 
						|| $a -> [1] <=> $b -> [1] 
						|| $a -> [3] <=> $b -> [3]
						|| $a -> [4] <=> $b -> [4] } @array_in;

	# ----------------------------------

	# Remove links below minimum size and store them in @array_out:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0	= $array_in[$i][1];
		my $var_loc1	= $array_in[$i][2];

		my $var_length	= $var_loc1 - $var_loc0 + 1;

		if ( $var_length >= $var_min ) {

			my $var_string	= join " ", @{$array_in[$i]};

			push @array_out, $var_string;

		}

	}

	# ----------------------------------

	# Write @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	my $var_scalar	= scalar @array_out;

	print "\t- $var_scalar merged $var_orientation links\n";

	# ----------------------------------

	# Return @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine:		Nmer
# Description:		This subroutine stores each possible nmer length nucleotide
#					in a hash.

# ----------------------------------------------------------

sub Nmer {

	# Arguments:
	my ( $hash_in, $array_in, $var_nmer )	= @_;

	# Data structures:
	my %hash_in			= %$hash_in;
	my @array_in		= @{$array_in};
	my %hash_out;

	# Variables:
	my $file_hash		= $hash_in { file_hash };
	my $file_direct		= $hash_in { file_direct };
	my $file_inverted	= $hash_in { file_inverted };
	my $var_null		= ("N")x$var_nmer;

	# ----------------------------------

	return if -e $file_direct && -e $file_inverted;

	# If $file_hash exists load it into memory:
	if ( -e $file_hash ) {

		print "LOADING\n";

		# Open file containing stored hash.
		open ( my $file_read, '<', $file_hash ) or die " ERROR: UNABLE TO OPEN $file_hash!\n";

		my $hash_json = <$file_read>;

		%hash_out = %{decode_json( $hash_json )};

		close $file_read or die " ERROR: UNABLE TO CLOSE $file_hash!\n";

	}

	# ----------------------------------

	# Otherwise, create the hash from scratch:
	else {

		# Iterate through each FASTA sequence:
		for ( my $i = 0; $i < scalar @array_in; $i++ ) {

			# Define variable holding FASTA header:
			my $var_header		= $array_in[$i][1];

			# Define variable to hold sequence:
			my $var_sequence	= $array_in[$i][2];

			# Define variable to hold length of sequence:
			my $var_length		= length $var_sequence;

			print "$i: $var_header - $var_length bp\n";

			# Iterate through sequence character by character creating a hash of
			# sequences of $var_nmer length:
			for ( my $j = 0; $j < $var_length - $var_nmer + 1; $j++ ) {

				my $var_seq	= substr $var_sequence, $j, $var_nmer;

				my @array_temp	= ( $i, $j );

				next if $var_seq eq $var_null;

				push ( @{ $hash_out { $var_seq } }, \@array_temp );

			}

		}

	}

	# ----------------------------------

	# Store data to file if not already present:
	unless ( -e $file_hash ) {

		print "SAVING\n";

		# Encode %hash_out in JSON format:
		my $hash_json = encode_json( \%hash_out );

		# Open $file_hash to write data:
		open ( my $file_write, '>', $file_hash ) or die " ERROR: UNABLE TO OPEN $file_hash!\n";

		# Print the hash to file.
		print $file_write $hash_json;

		# Close $file_out:
		close $file_write or die " ERROR: UNABLE TO CLOSE $file_hash!\n";

	}

	# ----------------------------------

	# End subroutine and return %hash_out:
	return \%hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Overlap
# Description:		This subroutine calculates the overlap between links and
#					genomic features:

# ----------------------------------------------------------

sub Overlap {

	# Arguments:
	my ( $hash_in, $array_features, $var_orientation )	= @_;

	# Data structures:
	my %hash_out		= %$hash_in;
	my @array_features	= @$array_features;
	my %hash_features;
	my %hash_total;
	my @array_out;
	my @array_flat;

	# Variables:
	my $var_total		= 0;
	my $var_flat		= 0;
	my $file_in;

	# ----------------------------------

	# Set $var_variable equal to one (direct) or minus one (inverted):
	if ( $var_orientation eq "direct" ) {

		$file_in		= $hash_out { file_direct_merged };

	}

	if ( $var_orientation eq "inverted" ) {

		$file_in		= $hash_out { file_inverted_merged };

	}

	# Import data from $file_in, containing links, to @array_in:
	my @array_in		= @{General::FileToArray ( $file_in )};

	# ----------------------------------

	# Calculate and store the length of each feature in @array_features:
	for ( my $i = 0; $i < scalar @array_features; $i++ ) {

		# Define start and end locations:
		my $var_loc0	= $array_features[$i][3];
		my $var_loc1	= $array_features[$i][4];

		# Store length:
		$array_features[$i][5]	= $var_loc1 - $var_loc0 + 1;

	}

	# Sort features from smallest to largest:
	@array_features	= sort { $a -> [5] <=> $b -> [5] } @array_features;

	# ----------------------------------

	# Sort repeats in @array_in by length smallest to largest:
#	@array_flat	= sort { $a -> [1] <=> $b -> [1] } @array_in;

	# Create a flattened array of links:
#	for ( my $i = 0; $i < scalar @array_flat; $i++ ) {

#		my $var_loc0		= $array_flat[$i][1];
#		my $var_loc1		= $array_flat[$i][2];

#		print "$i: $var_loc0 $var_loc1\n";

#		for ( my $j = $i+1; $j < scalar @array_flat; $j++ ) {

#			my $var_loc2		= $array_flat[$j][1];
#			my $var_loc3		= $array_flat[$j][2];

#			print "\t$j: $var_loc2 $var_loc3\n";

			# If region is entirely inside, remove it and start over:
#			if ( $var_loc0 <= $var_loc2 && $var_loc1 >= $var_loc3 ) {

#				print "\t$j: $var_loc2 $var_loc3\n";

#				splice @array_flat, $i, 1;

#				$i = -1;

#				last;

#			}

			# If region overlaps, merge them and start over:
#			elsif ( $var_loc2 <= $var_loc1 && $var_loc3 >= $var_loc1 ) {

#				print "\t$j: $var_loc2 $var_loc3\n";

#				$array_flat[$i][2]	= $var_loc3;

#				splice @array_in, $i, 1;

#				$i = -1;

#				last;

#			}

#		}

#	}

	# Add length information to each link in @array_flat:
#	for ( my $i = 0; $i < scalar @array_flat; $i++ ) {

#		my $var_loc0		= $array_flat[$i][1];
#		my $var_loc1		= $array_flat[$i][2];
#		my $var_length		= $var_loc1 - $var_loc0 + 1;

		# Store information:
#		$array_flat[$i][6]	= $var_length;

		# Add to $var_total:
#		$var_flat			+= $var_length;

#	}

	# ----------------------------------

	# Add length information to each link in @array_in:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_loc0		= $array_in[$i][1];
		my $var_loc1		= $array_in[$i][2];
		my $var_length		= $var_loc1 - $var_loc0 + 1;

		# Store information:
		$array_in[$i][6]	= $var_length;

		# Add to $var_total:
		$var_total			+= $var_length;

	}

	# Sort repeats in @array_in by length smallest to largest:
	@array_in	= sort { $b -> [6] <=> $a -> [6] } @array_in;

	# ----------------------------------

	# Iterate through @array_in and identify overlapping genomic features:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define sequence ID and start/end locations of first portion of link:
		my $var_seq0			= $array_in[$i][0];
		my $var_loc0			= $array_in[$i][1];
		my $var_loc1			= $array_in[$i][2];

		# Define sequence ID and start/end locations of second portion of link:
		my $var_seq1			= $array_in[$i][3];
		my $var_loc2			= $array_in[$i][4];
		my $var_loc3			= $array_in[$i][5];

		# Iterate down through @array_features and identify all features that
		# overlap with first and second portions of the link:
		for ( my $j = 0; $j < scalar @array_features; $j++ ) {

			my $var_type	= $array_features[$j][1];
			my $var_loc4	= $array_features[$j][3];
			my $var_loc5	= $array_features[$j][4];
			my $var_length0	= $array_features[$j][5];

			# First portion overlapping left-hand region of the feature:
			if ( $var_loc0 < $var_loc4 && $var_loc1 > $var_loc4 && $var_loc1 < $var_loc5 ) {

				my $var_length	= $var_loc1 - $var_loc4 + 1;

				unless ( $array_in[$i][7] && $array_in[$i][8] ) {

					$array_in[$i][7]	= $var_type;
					$array_in[$i][8]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# Second portion overlapping left-hand region of feature:
			if ( $var_loc2 < $var_loc4 && $var_loc3 > $var_loc4 && $var_loc3 < $var_loc5 ) {

				my $var_length	= $var_loc3 - $var_loc4 + 1;

				unless ( $array_in[$i][9] && $array_in[$i][10] ) {

					$array_in[$i][9]	= $var_type;
					$array_in[$i][10]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# First portion overlapping right-hand region of the feature:
			if ( $var_loc0 > $var_loc4 && $var_loc0 < $var_loc5 && $var_loc1 > $var_loc5 ) {

				my $var_length	= $var_loc5 - $var_loc0 + 1;

				unless ( $array_in[$i][9] && $array_in[$i][10] ) {

					$array_in[$i][9]	= $var_type;
					$array_in[$i][10]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# Second portion overlapping right-hand region of the feature:
			if ( $var_loc2 > $var_loc4 && $var_loc2 < $var_loc5 && $var_loc3 > $var_loc5  ) {

				my $var_length	= $var_loc5 - $var_loc2 + 1;

				unless ( $array_in[$i][7] && $array_in[$i][8] ) {

					$array_in[$i][7]	= $var_type;
					$array_in[$i][8]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# First portion of link is completely within the feature:
			if ( $var_loc0 > $var_loc4 && $var_loc1 < $var_loc5 ) {

				my $var_length	= $var_loc1 - $var_loc0 + 1;

				unless ( $array_in[$i][7] && $array_in[$i][8] ) {

					$array_in[$i][7]	= $var_type;
					$array_in[$i][8]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# Second portion of link is completely within the feature:
			if ( $var_loc2 > $var_loc4 && $var_loc3 < $var_loc5 ) {

				my $var_length	= $var_loc3 - $var_loc2 + 1;

				unless ( $array_in[$i][9] && $array_in[$i][10] ) {

					$array_in[$i][9]	= $var_type;
					$array_in[$i][10]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# First portion of link completely covers feature:
			if ( $var_loc0 < $var_loc4 && $var_loc1 > $var_loc5 ) {

				my $var_length	= $var_loc5 - $var_loc4 + 1;

				unless ( $array_in[$i][7] && $array_in[$i][8] ) {

					$array_in[$i][7]	= $var_type;
					$array_in[$i][8]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

			# Second portion of the link entirely overlapped by feature:
			if ( $var_loc0 < $var_loc4 && $var_loc1 > $var_loc5 ) {

				my $var_length	= $var_loc5 - $var_loc4 + 1;

				unless ( $array_in[$i][9] && $array_in[$i][10] ) {

					$array_in[$i][9]	= $var_type;
					$array_in[$i][10]	= $var_length0;

				}

				$hash_total { $var_type }	+= $var_length;

			}

		}

		unless ( $array_in[$i][7] && $array_in[$i][8] ) {

			$array_in[$i][7]	= "chromosome";
			$array_in[$i][8]	= 0;

		}

		unless ( $array_in[$i][9] && $array_in[$i][10] ) {

			$array_in[$i][9]	= "chromosome";
			$array_in[$i][10]	= 0;

		}

	}

	# ----------------------------------

	# Print link orientation:
#	print "$var_orientation links:\n";

	# Print length data for each type of genomic feature:
#	foreach my $var_type ( sort { $a cmp $b } keys %hash_total ) {

#		my $var_total		= $hash_total { $var_type } || 0;

#		print "\t - $var_type $var_total bp\n";

#	}

	# Print total length of links:
#	print "\t - Total: $var_total bp\n";

	# Print total flat length of links:
#	print "\t - Flat: $var_flat bp\n";

	# ----------------------------------

	@array_out	= @array_in;

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Ticks
# Description:		This subroutine takes region data and creates tick marks.

# ----------------------------------------------------------

sub Ticks {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my %hash_out			= %$hash_in;
	my @array_in			= @$array_in;
	my %hash_colours		= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $var_length0			= 1;
	my $var_length1			= 0;
	my $file_out			= $hash_out { file_ticks };
	my $var_thickness		= $hash_out { stroke_thickness };
	my $var_colour			= $hash_colours { ticks };
	my $var_gap				= $hash_out { var_gap };
	my $var_radius0			= $hash_out { var_radius0 };
	my $var_1Mb_height		= $hash_out{var_1Mb_height};
	my $var_100kb_height	= $hash_out{var_100kb_height};
	my $var_10kb_height		= $hash_out{var_10kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_1Mb				= 10**6;
	my $var_100kb			= 10**5;
	my $var_10kb			= 10**4;

	# ----------------------------------

	# Add file name to beginning of file:
	push @array_out, "#$file_out";

	# Add ring highlight(s):
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_length	= length $array_in[$i][2];

		my $var_string	= "chr$i 1 $var_length ";
		$var_string		.= "stroke_color=$var_colour,";
		$var_string		.= "stroke_thickness=$var_thickness,";
		$var_string		.= "r0=$var_radius0$r+$var_gap$p,";
		$var_string		.= "r1=$var_radius0$r+$var_gap$p";

		push @array_out, $var_string;

	}

	# ----------------------------------

	# Add first tick:
	my $var_string		= "chr0 1 1 ";
	$var_string			.= "stroke_color=$var_colour,";
	$var_string			.= "stroke_thickness=$var_thickness,";
	$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
	$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

	push @array_out, $var_string;

	# Iterate through @array_in adding vertical marks for ticks:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_length2	= length $array_in[$i][2];

		$var_length1 	+= $var_length2;

		my $var_loc		= 1;

		for ( my $j = $var_length0; $j <= $var_length1; $j++ ) {

			if ( $j % $var_10kb == 0 ) {

				my $var_string	= "chr$i $var_loc $var_loc ";
				$var_string		.= "stroke_color=$var_colour,";
				$var_string		.= "stroke_thickness=$var_thickness,";

				if ( $j % $var_1Mb == 0 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

				} 

				elsif ( $j % $var_100kb == 0 ) {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_100kb_height+$var_gap$p";

				} 

				else {

					$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
					$var_string			.= "r1=$var_radius0$r+$var_10kb_height+$var_gap$p";

				} 

				push @array_out, $var_string;

			}

			$var_loc++;

		}

		$var_length0	= $var_length1 + 1;

	}
	# ----------------------------------

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, "$file_out" );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
