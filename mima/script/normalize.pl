#!/usr/bin/perl -w
use strict;
BEGIN {
    unshift @INC, "/srv/scratch/mrcbio/scripts/";
}
#
use Matrix;
die "perl $0 <matrix.table>  <normalize to number> \n" unless(@ARGV == 2);
#
normalized($ARGV[0], $ARGV[1]);
