#!/usr/bin/perl -w
use strict;
use Getopt::Std;
#use List::MoreUtilsi qw/first_index/;

my $usage="perl $0 -d <dm.txt> -m <map.txt> -t <colname:groupname> -o <output_file>
This program calculates the mean distance to group A for all samples in metadata and selects and output the results of  samples that you are interested in.
-d <path of distance matrices file>    Required. Such as weighted_unifrac_dm.txt
-o <output_dir>        Requried.
-m <path of mapping file>      Requried. Must contains the group A.
-t [group A colname:groupname]	e.g. if calculate mean distance to healthy control --> group:control. Or calculate the mean distance to patients with T2DM --> group:T2DM 
-l [path of list of target samples]	Or provide a list of samples instead of colname:groupname 
-g [colname:groupname of output samples]	Indicate the samples that you are interested in.
-s [path of list of output samples]	Or provied a list of output samples instead of colname:groupname
-h [help]
";

use vars qw($opt_d $opt_o $opt_m $opt_t $opt_l $opt_g $opt_s $opt_h);
getopts('d:o:m:t:l:g:s:h');

if ($opt_h){
        print $usage,"\n";
        exit(1);
}

unless ($opt_o and $opt_m and $opt_d){
        print "$usage\n";
        exit(1);
}

unless ($opt_t or $opt_l){
	print "$usage\n";
	print "Error: target samples should be provided using either -t or -l\n";
	exit(1);
}
mkdir $opt_o;
## Read in the target samples
my %target;
my $d_name;
if($opt_t){
	my @target_group = split(/\:/,$opt_t);
	open MAP,"$opt_m" or die "can't open $opt_m \n";
	chomp(my $header = <MAP>);
	my @header = split(/\t/,$header);
	my $index;
	if (grep /^$target_group[0]$/, @header){
		$index = first_index($target_group[0],@header);
	}else{
		die "Please check if $opt_m contains $target_group[0]\n";
	}
	while(<MAP>){
		chomp;
		my @line = split(/\t/);
		$target{$line[0]} = 0 if $line[$index] eq $target_group[1];
	}
	close MAP;
	$d_name = "mean_distance2$target_group[0]\_$target_group[1]";
}elsif ($opt_l){
	open LIST,"$opt_l" or die "can't open $opt_l \n";
	while(<LIST>){
		chomp;
		$target{$_} = 0;
	}
	close LIST;
	$d_name = "mean_distance";
}

## Read in the target samples, finished

## sum distance in %sum_dm
open DM,"$opt_d" or die "can't open $opt_d\n";
chomp(my $dmhead=<DM>);
my @dmhead=split(/\t/,$dmhead);

my @target_index;
foreach my $t (keys %target){
	if (grep /^$t$/, @dmhead){
		push @target_index, first_index($t,@dmhead);
	}else{
		print "distance metrics does not contain $t\n ";
	}
}

my %sum_dm;
my %mean_dm;
while(<DM>)
{
        chomp;
        my @line=split(/\t/);
        foreach my $i(@target_index){
                $sum_dm{$line[0]}+=$line[$i];
        }
	$mean_dm{$line[0]} = $sum_dm{$line[0]}/($#target_index+1);
}
close DM;
## sum distance in %sum_dm, finished

## Add mean distance to metadata
open OUT,">$opt_o/metadata.txt" or die "can't create $opt_o/metadata.txt\n";
open MAP,"$opt_m" or die "can't open $opt_m\n";
chomp(my $maphead=<MAP>);
print OUT "$maphead\t$d_name\n";
while(<MAP>)
{
        chomp;
        my @line=split(/\t/);
	if (exists $mean_dm{$line[0]}){
	        print OUT "$_\t$mean_dm{$line[0]}\n";
	}else{
		print OUT "NA\n";
	}
}
close MAP;
close OUT;
## Add mean distance to metadata, finished


## filter
if($opt_g){
	my @output_list;
        my @output_group = split(/\:/,$opt_g);
        open MAP,"$opt_m" or die "can't open $opt_m \n";
        chomp(my $header = <MAP>);
        my @header = split(/\t/,$header);
        my $index;
        if (grep /^$output_group[0]$/, @header){
                $index = first_index($output_group[0],@header);
        }else{
                die "Please check if $opt_m contains $output_group[0]\n";
        }
        while(<MAP>){
                chomp;
                my @line = split(/\t/);
                push @output_list,$line[0] if $line[$index] eq $output_group[1];
        }
        close MAP;
	filter_samples("$opt_o/metadata.txt","$opt_o/filtered_metadata.txt",@output_list);
}elsif ($opt_s){
	my @output_list;
        open LIST,"$opt_s" or die "can't open $opt_s \n";
        while(<LIST>){
                chomp;
                push @output_list,$_;
        }
        close LIST;
        filter_samples("$opt_o/metadata.txt","$opt_o/filtered_metadata.txt",@output_list);
}


########## sub declarations come here
sub filter_samples{
	my ($input_meta,$output_meta,@list) = @_;
	

	open OUT,">$output_meta" or die "can't create $output_meta\n";
	open MAP,"$input_meta" or die "can't open $input_meta\n";
	my $maphead=<MAP>;
	print OUT $maphead;
	while(<MAP>)
	{
	        my @line=split(/\t/);
	        if (grep /^$line[0]$/, @list){
	                print OUT $_;
	     	}
	}
	close MAP;
	close OUT;
}

sub first_index{
	my ($element,@list) = @_;
	my $index = -1;
	for my $i (0..$#list){
		$index = $i if $list[$i] eq $element;
	}
	return $index;
}
