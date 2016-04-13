#!/usr/bin/env perl
#
# create index of samfile for seek() operations
# 
#  chromosome name is in column 2 (counting from 0), position is column 3
#
use strict;

my $samfile = shift;

my $blocksize = 1000000;
my $last_chrom = "";
my $last_chrom_pos = 0;

open( SAM, $samfile ) || die "can't open $samfile: $!\n";

my $file_pos = tell( SAM );

while ( $_ = <SAM> )
   {
    my ( $chrom, $chrom_pos ) = (split( /\t/ ))[2,3];

    if ( ($chrom ne $last_chrom) )
          {
           $last_chrom = $chrom;
           $last_chrom_pos = 0;
          }
    while ( $last_chrom_pos <= $chrom_pos )
       {
        print "$chrom $last_chrom_pos $file_pos\n";
        $last_chrom_pos += $blocksize;
       }
    $file_pos = tell( SAM );
   }
close( SAM );
