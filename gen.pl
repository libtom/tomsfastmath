#!/usr/bin/perl -w
#
# Generates a "single file" you can use to quickly
# add the whole source without any makefile troubles
#
# SPDX-License-Identifier: Unlicense
use strict;

open( OUT, ">mpi.c" ) or die "Couldn't open mpi.c for writing: $!";
print OUT "#define TFM_PRE_GEN_MPI_C\n";
foreach my $filename (glob "src/*/*fp_*.c") {
   next if ($filename eq "fp_sqr_comba_generic.c");
   open( SRC, "<$filename" ) or die "Couldn't open $filename for reading: $!";
   print OUT "/* Start: $filename */\n";
   print OUT while <SRC>;
   print OUT "\n/* End: $filename */\n\n";
   close SRC or die "Error closing $filename after reading: $!";
}
print OUT "\n/* EOF */\n";
close OUT or die "Error closing mpi.c after writing: $!";
