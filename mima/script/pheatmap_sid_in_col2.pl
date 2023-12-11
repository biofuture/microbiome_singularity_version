#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $usage="perl $0 -i <table.txt> -m <metadata.txt> -c <group1,group2..> -o <output_dir>
-i <absolute path of table.txt>    Required. table.txt without annotation, such as otu_table.txt without 'Taxonomy', taxa.txt or any other tables.
-o <absolute path of output_dir>        Requried.
-m <absolute path of mapping file>	Requried. SampleID is in the second column.
-c <names of groups> 	Requried.
-n [numbers of groups that samples divide into]	default 1.
-h [help]
";

use vars qw($opt_i $opt_o $opt_m $opt_c $opt_n $opt_h);
getopts('i:o:m:c:n:h');

if ($opt_h){       
	print $usage,"\n";
        exit(1);
}

unless ($opt_i and $opt_o and $opt_c and $opt_m){
        print "$usage\n";
        exit(1);
}

if (-e $opt_o){
        print "output_dir already exists\n";
        exit(1);
}
`mkdir $opt_o`;

if(not $opt_n){
	$opt_n = 1;
}

my @meta = split(/,/,$opt_c);
my %meta;
foreach (@meta){
	$meta{$_} = 1;
}
open MAP,"$opt_m" or die "can't open $opt_m\n";
open OUT,">$opt_o/annotation.map.txt" or die "can't create $opt_o/annotation.map.txt.\n";
chomp(my $map_first = <MAP>);
my @map_first = split(/\t/,$map_first);
print OUT $map_first[1];
foreach my $i(0..$#map_first){
	if ($meta{$map_first[$i]}){
		print OUT "\t$map_first[$i]";
	}
}
print OUT "\n";
while(<MAP>){
	chomp;
	my @map_line = split(/\t/);
	print OUT $map_line[1];
	foreach my $i(0..$#map_first){
        	if ($meta{$map_first[$i]}){
                	print OUT "\t$map_line[$i]";
        	}
	}
	print OUT "\n";
}
close OUT;
close MAP;


#### R
open R,">$opt_o/pheatmap.r" or die "can't create pheatmap.r\n";
print R "

#library(psych)
library(pheatmap)

# Ignore \"#\" while reading
read.table.x <- function(filename, ...){
  lines <- readLines(filename)
  n <- grep(\"^#\", lines)
  if(length(n) > 0){start <- n[length(n)]}else{start <- 1}
  end <- length(lines)
  x <- read.table(text=lines[start:end],header=T,sep=\"\\t\",comment.char = '',check.names=F,...)
}

taxa <- read.table.x(\"$opt_i\",row.names = 1)
anno <- read.table.x(\"$opt_o/annotation.map.txt\",row.names = 1)

## heatmap
pheatmap(taxa,clustering_method=\"ward.D\",color = colorRampPalette(c(\"steelblue\",\"snow\",\"tomato3\"))(1000),cutree_cols = $opt_n,
         cellheight=6,cellwidth=6,fontsize=8,
         annotation_col = anno,border_color=NA,filename=\"$opt_o/pheatmap.pdf\")
## heatmap

";
`module load R/3.5.1`;
`R --no-save < $opt_o/pheatmap.r`;
#### R

