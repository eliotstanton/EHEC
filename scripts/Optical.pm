package Optical;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# File name:		Optical.pm
# Date created:		18 Devember, 2018
# Last modified:	18 December, 2018
# Created by:		Eliot Stanton

# Description:		This package contains subroutines specific to generating
#					images of aligned optical maps.

# ---------------------------------------------------------------------------- #

# Subroutine name:	Config
# Description:		This subroutine generates a a configuration file for use by
# 					Circos.

# ----------------------------------------------------------

sub Config {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my %hash_in			= %$hash_in;
	my @array_out;

	# Variables:
	my $file_out			= $hash_in { "file_config" };
	my $file_highlights		= $hash_in { "file_highlights" };
	my $file_ideogram		= $hash_in { "file_ideogram" };
	my $file_image			= $hash_in { "file_image" };
	my $file_karyotype		= $hash_in { "file_karyotype" };
	my $file_ticks			= $hash_in { "file_ticks" };
	my $var_scalar			= scalar @array_in;
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

# Subroutine name:	Highlights
# Description:		This subroutine generates highlights files representing each
#					restriction fragment and for use in Circos

# ----------------------------------------------------------

sub Highlights {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my %hash_in			= %$hash_in;
	my @array_out;

	# Variables:
	my $file_out		= $hash_in { file_highlights };
	my $var_thickness	= $hash_in { stroke_thickness };
	my $var_span		= $hash_in { var_span };
	my %hash_colours	= %{Circos::Colours()};
	my $var_scalar		= scalar @array_in;

	# ----------------------------------

	# Add highlights for regions:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Iterate down through @array_in formatting data for Circos:
		for ( my $j	= 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			my $var_loc0	= $array_in[$i][$j][0];
			my $var_loc1	= $array_in[$i][$j][1];

			my $var_border		= $hash_colours { ticks };
			my $var_colour		= Circos::ColoursNew ( "chromosome", $var_scalar, $var_scalar );
			my $var_string		= "chr0 $var_loc0 $var_loc1 ";

			$var_string			.= "stroke_color=$var_border,";
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

# Subroutine name:	Karyotype
# Description:		Generate a karyotype file for use by Circos.

# ----------------------------------------------------------

sub Karyotype {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my @array_in		= @$array_in;
	my %hash_in			= %$hash_in;
	my %hash_colours	= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $var_length		= 0;
	my $file_out		= $hash_in { "file_karyotype" };
	my $var_colour		= Circos::ColoursNew ( "chromosome", 1, 1 );

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		my $var_max		= scalar @{$array_in[$i]} - 1;

		my $var_length0	= $array_in[$i][$var_max][1];

		$var_length		= $var_length0 if $var_length0 > $var_length;

	}

	my $var_string	= "chr - chr0 0 0 $var_length $var_colour\n";

	push @array_out, $var_string;

	# Print @array_out to $file_out:
	General::ArrayToFile ( \@array_out, $file_out );

	# ----------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Locations
# Description:		This subroutine takes fragment sizes and formats them into
#					adjusted locations.

# ----------------------------------------------------------

sub Locations {

	# Arguments:
	my ( $array_in )	= @_;

	# Data structures:
	my @array_in	= @$array_in;
	my @array_out;

	# ----------------------------------

	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define variable to hold start and end locations of each fragment:
		my $var_start	= 1;
		my $var_end		= 0;

		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			my $var_length	= $array_in[$i][$j][0];

			if ( $var_length < 0 ) {

				$var_length	= abs($var_length);

				$var_end += $var_length if $var_length > 0;

				$var_start += $var_length if $var_length > 0;

			}

			else {

				$var_end += $var_length if $var_length > 0;

#				print "\t$j: $var_length - $var_start $var_end\n";

				my @array_temp	= ( $var_start, $var_end, $var_length );

				push @{$array_out[$i]}, \@array_temp;

				$var_start += $var_length if $var_length > 0;

			}

		}

	}

	# ----------------------------------

	# Return reference to @array_out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine name:	Ticks
# Description:		Generate ticks that match locations of fragments.

# ----------------------------------------------------------

sub Ticks {

	# Arguments:
	my ( $array_in, $hash_in )	= @_;

	# Data structures:
	my @array_in			= @$array_in;
	my %hash_in				= %$hash_in;
	my %hash_colours		= %{Circos::Colours()};
	my @array_out;

	# Variables:
	my $file_out			= $hash_in { file_ticks };
	my $var_thickness		= $hash_in { stroke_thickness };
	my $var_span			= $hash_in { var_span };
	my $var_colour			= $hash_colours { ticks };

	my $var_gap				= $hash_in { var_gap };

	my $var_radius0			= $hash_in { var_radius0 };

	my $var_1Mb_height		= $hash_in {var_1Mb_height};
	my $var_100kb_height	= $hash_in {var_100kb_height};
	my $var_10kb_height		= $hash_in {var_10kb_height};
	my $r					= "r";
	my $p					= "p";
	my $var_increment		= $hash_in { var_increment };

	# ----------------------------------

	# Iterate through @array_in adding horizontal and vertical marks for ticks:
	for ( my $i = 0; $i < scalar @array_in; $i++ ) {

		# Define marker lengths:
		my $var_1Mb		= 10**6;
		my $var_100kb	= 10**5;
		my $var_10kb	= 10**4;

		# Add first tick:
		my $var_string		= "chr0 1 1 ";
		$var_string			.= "stroke_color=$var_colour,";
		$var_string			.= "stroke_thickness=$var_thickness,";
		$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
		$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

		push @{$array_out[$i]}, $var_string;

		my $var_start	= 1;
		my $var_end		= 0;

		# Iterate through regions in @array_in creating original locations:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define locations and length:
			my $var_loc0			= $array_in[$i][$j][0];
			my $var_loc1			= $array_in[$i][$j][1];
			my $var_length			= $array_in[$i][$j][2];

			$var_end				+= $var_length;

			$array_in[$i][$j][2]	= $var_start;
			$array_in[$i][$j][3]	= $var_end;

			$var_start				+= $var_length;

		}

		# Iterate through regions in @array_in making up vertical ticks:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define original locations:
			my $var_loc0		= $array_in[$i][$j][2];
			my $var_loc1		= $array_in[$i][$j][3];

			# Define adjusted locations:
			my $var_loc2		= $array_in[$i][$j][0];
			my $var_loc3		= $array_in[$i][$j][1];

			my $var_diff		= $var_loc2 - $var_loc0;

			# Handle ticks every 10 kb:
			if ( $var_10kb >= $var_loc0 && $var_10kb < $var_loc1 ) {

				my $var_temp	= $var_10kb + $var_diff;

				$var_10kb += 10**4;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";
				$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
				$var_string			.= "r1=$var_radius0$r+$var_10kb_height+$var_gap$p";

				push @{$array_out[$i]}, $var_string;

				$j--;

				next;

			}

			# Handle ticks every 100kb:
			if ( $var_100kb >= $var_loc0 && $var_100kb < $var_loc1 ) {

				my $var_temp	= $var_100kb + $var_diff;


				$var_100kb += 10**5;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";
				$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
				$var_string			.= "r1=$var_radius0$r+$var_100kb_height+$var_gap$p";

				push @{$array_out[$i]}, $var_string;

				$j--;

				next;

			}

			# Handle ticks every 1 Mb:
			if ( $var_1Mb >= $var_loc0 && $var_1Mb < $var_loc1 ) {

				my $var_temp	= $var_1Mb + $var_diff;

				my $var_string		= "chr0 $var_temp $var_temp ";

				$var_string			.= "stroke_color=$var_colour,";
				$var_string			.= "stroke_thickness=$var_thickness,";
				$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
				$var_string			.= "r1=$var_radius0$r+$var_1Mb_height+$var_gap$p";

				push @{$array_out[$i]}, $var_string;

				$var_1Mb += 10**6;

				$j--;

			}

		}

		# Handle horizontal lines:
		for ( my $j = 0; $j < scalar @{$array_in[$i]}; $j++ ) {

			# Define locations:
			my $var_loc0		= $array_in[$i][$j][0];
			my $var_loc1		= $array_in[$i][$j][1];

			my $var_string		= "chr0 $var_loc0 $var_loc1 ";

			$var_string			.= "stroke_color=$var_colour,";
			$var_string			.= "stroke_thickness=$var_thickness,";
			$var_string			.= "r0=$var_radius0$r+$var_gap$p,";
			$var_string			.= "r1=$var_radius0$r+$var_gap$p";

#			push @{$array_out[$i]}, $var_string;

		}

	 	$var_radius0	-= $var_increment;

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

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #

1;
