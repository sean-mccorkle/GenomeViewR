#!/usr/bin/env perl
# Program:      build_index.pl
# Programmer:   Sean R. McCorkle
#               Biology Departmnet, Brookhaven National Laboratory
# Language:     perl
#
# Description:  Creates an index of blocks of SAMfile for seek() operations,
#               to allow rapid retreival of 1 megabyte regions surrounding
#               a given chromosome/position.
#
# Usage:        build_index.pl <SAMfile> >output_index
#
#               where <SAMfile> is a sorted SAM-format filey
#
#               output (on stdout) is an ascii, where each record contains
#
#               <chrom> <beginning_position_of_block> <line number> <file position>
#
#
use strict;

my $samfile = shift;

my $blocksize = 1000000;
my $last_chrom = "";
my $last_chrom_pos = 0;

open( SAM, $samfile ) || die "can't open $samfile: $!\n";

my $line_num = 0;           # for now, first line is numbered 0
my $file_pos = tell( SAM );

while ( $_ = <SAM> )
   {
    # Reading the input SAM file, its assumed that chromosome 
    # name is in column 2 (counting from 0), position is column 3

    my ( $chrom, $chrom_pos ) = (split( /\t/ ))[2,3];

    if ( ($chrom ne $last_chrom) )
          {
           $last_chrom = $chrom;
           $last_chrom_pos = 0;
          }
    while ( $last_chrom_pos <= $chrom_pos )
       {
        print "$chrom $last_chrom_pos $line_num $file_pos\n";
        $last_chrom_pos += $blocksize;
       }
    $file_pos = tell( SAM );
    $line_num++;
   }
close( SAM );
