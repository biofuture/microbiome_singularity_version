#!/usr/bin/perl -w
use strict;

die "perl $0 <meta_data_online.txt> <x_Group> <output_dir>\n" unless @ARGV==3;

mkdir $ARGV[2];#R plot, R plot script and the input file for plot will be saved here.

# Create R script for plot -------------------------------------------------
open R,">$ARGV[2]/box_jitter.plot.r" or die "can't create $ARGV[2]/box_jitter.plot.r\n";
print R "
# Load the packages needed
library(ggplot2);library(grid)

# Read the data
read.table.meta = function(filename, ...){
	lines <- readLines(filename)
	n <- grep(\"^\#\", lines)
	if(length(n) > 0){start <- n[length(n)]}else{start <- 1}
	end <- length(lines)
	x <- read.table(text=lines[start:end],header=T,sep=\"\\t\",comment.char = '',...)
}
map <- read.table.meta(\"$ARGV[0]\")
map\$copynumb <- map\$X.of16Sreads/map\$CellNumber

if(length(levels(map\$$ARGV[1]))>2){
  r <- kruskal.test(map\$copynumb~map\$$ARGV[1])
}else{
  r <- wilcox.test(map\$copynumb~map\$$ARGV[1])
}

# Colors for choose
#colors = c(\"#80B1D3\",\"#B3DE69\",\"#FFFFB3\",\"#8DD3C7\",\"#4daf4a\",\"#377eb8\",\"#BEBADA\",\"#FB8072\",\"#FDB462\",\"#FCCDE5\",\"#BC80BD\",\"#CCEBC5\",\"#FFED6F\",\"#CD4F39\",\"#BC41A4\",\"#4F94CD\",\"#E41A1C\",\"#00CD66\",\"#CD3278\",\"#CD8A96\",\"#00C5CD\",\"#CDCD00\",\"#CD85CD\",\"#CD853F\",\"#8B5A2B\",\"#5CACEE\",\"#EE5C42\",\"#00EE76\",\"#EE4A8C\",\"#EED8AE\",\"#00E5EE\",\"#EEEE00\",\"#EED2EE\",\"#EE9A49\",\"#E41A1C\",\"#377EB8\",\"#FF6A6A\",\"#87CEFA\",\"#6E8B3D\",\"#FFEBCD\",\"#B2DFEE\")



# Plot use ggplot2
  ggplot(map,aes(x=$ARGV[1],y=copynumb))+stat_boxplot(geom = \"errorbar\",width=0.3) + geom_boxplot(notch = T, outlier.shape = NA) + geom_jitter(alpha = 0.5,height = 0) +theme_bw()+ylab(\"Average 16S copy numbers\")+xlab(paste(\"P = \",r\$p.value,sep =\"\"))
  
  ggsave(filename = paste(\"$ARGV[2]\",\"/\",\"16S_copynum_$ARGV[1].pdf\",sep = \"\"),height = 4,width = 4)


";
# Create R script for plot -----------------------------------------finished

`R --no-save < $ARGV[2]/box_jitter.plot.r`;

