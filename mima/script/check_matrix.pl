#!/usr/bin/perl -w
use strict; 

die "$!\n" unless open(I, "$ARGV[0]"); 
my $col = 0;
my $row = 0;

my $head = <I>; chomp($head); 
my @h = split("\t", $head);
$col = $#h +1; 
while(<I>){
	chomp;
	my @m = split 	/\t/;
	if($col != $#m +1){ die "inconsistent $col\t",$#m+1,"\n"," $head\n$_\n"; }  
	$row ++; 
}
close I; 
print "row $row\tcol $col\n";
