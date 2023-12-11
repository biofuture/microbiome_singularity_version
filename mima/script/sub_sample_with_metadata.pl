#!/usr/bin/perl
use strict;

#This script subset a tsv feature table with definition in the meta data

die "$!\nperl $0 <matrix.tsv> <meta_data.txt> <output_dir>  <oprefix> <Columnname:Group1,Group2>\n

perl  $0 matrix.tsv meta_data.txt <outputdir> <oprefix> Cater1:Group1,Group2,GroupN Cater2:Group1,Group2,Group3\n" unless (@ARGV == 5 || @ARGV == 6);

unless(-d "$ARGV[2]"){
	`mkdir $ARGV[2]`;
} 
my $oprefix = "$ARGV[3]";
my $filsid= "$ARGV[2]/$oprefix.filtered_sid.txt";
my $filmetadata = "$ARGV[2]/$oprefix.filtered_metadata.txt";
my %filterids;
if(@ARGV == 5){	
	my ($cname, $groups) = split(/\:/, $ARGV[4]);
	my $headh = "";
	if($groups eq ""){
		##keep all ids
		die "$!\n" unless open(META, "$ARGV[1]");
		$headh = <META>; chomp($headh);
		while(<META>){
			chomp;
			my @tem = split /\t/;
		       	$filterids{$tem[0]} = $_;	
		}
		close META;		
	}else{
		my @detailgroups = split(/\,/, $groups);
		my %gcater = ();
		for my $k(@detailgroups){ $gcater{$k} = 1; }
		die "$!\n" unless open(META, "$ARGV[1]");
		my $index = 0; 
		my $head = <META>; chomp($head); $headh = $head; 
		my @headcol = split("\t", $head);
		for(my $i=1; $i <= $#headcol; $i++){
			if($headcol[$i] eq $cname){
				$index = $i;
				last;
			}
		}

		while(<META>){
			chomp;
			my @tem = split /\t/;
			if(exists $gcater{$tem[$index]}){
				#keep this sample
				$filterids{$tem[0]} = $_;
			}
		}
		close META;

	}	
	die "$!\n" unless open(METAO, ">$filmetadata");
	die "$!\n" unless open(METASID, ">$filsid");
	print METAO "$headh\n";
	for my$id (keys %filterids){
		print METASID "$id\n";
		print METAO "$filterids{$id}\n";
	}
	close METAO;
	
}elsif(@ARGV == 6){
	my ($cname1, $groups1) = split(/\:/, $ARGV[4]);
	my ($cname2, $groups2) = split(/\:/, $ARGV[5]);
	my $headh = "";
	if($groups1 eq "" && $groups2 eq ""){
		##keep all ids
		die "$!\n" unless open(META, "$ARGV[1]");
		$headh = <META>; chomp($headh);
		while(<META>){
			chomp;
			my @tem = split /\t/;
		       	$filterids{$tem[0]} = 1;	
		}
		close META;		
	}else{
		my @detailgroups1 = split(/\,/, $groups1);
		my @detailgroups2 = split(/\,/, $groups2);
		my %gcater1 = ();
		my %gcater2 = ();
		for my $k(@detailgroups1){ $gcater1{$k} = 1; }
		for my $k(@detailgroups2){ $gcater2{$k} = 1; }
		die "$!\n" unless open(META, "$ARGV[1]");
		my $index1 = 0;
		my $index2 = 0;
		my $head = <META>; chomp($head); $headh = $head;
		my @headcol = split("\t", $head);
		for(my $i=1; $i <= $#headcol; $i++){
			if($headcol[$i] eq $cname1){
				$index1 = $i;
			}
			if($headcol[$i] eq $cname2){
				$index2 = $i;
			}
		}

		while(<META>){
			chomp;
			my @tem = split /\t/;
			if(exists $gcater1{$tem[$index1]} && $gcater2{$tem[$index2]}){
				#keep this sample
				$filterids{$tem[0]} = $_;
			}
		}
		close META;

	}	
	die "$!\n" unless open(METAO, ">$filmetadata");
	die "$!\n" unless open(METASID, ">$filsid");
	print METAO  "$headh\n";
	for my$id (keys %filterids){
		print METAO "$filterids{$id}\n";
		print METASID "$id\n";
	}
	close METAO;
}

##fetch the table with selected IDs 
#
my $osh = "$ARGV[2]/$oprefix.shell_transform.sh";
die "$!\n" unless open(SH, ">$osh");
my $tablebiom = "$ARGV[2]/$oprefix.raw.biom";
my $tableselbiom = "$ARGV[2]/$oprefix.sel.biom";
my $otable = "$ARGV[2]/$oprefix.output_table.txt";
print SH "source activate /srv/scratch/z3524677/bashbin/miniconda3/envs/qiime2-2020.2\n";
print SH "biom convert -i $ARGV[0] -o $tablebiom --table-type \"Taxon table\" --to-hdf5\n";
print SH "biom subset-table -i $tablebiom -a sample -s $filsid -o $tableselbiom\n";
print SH "biom convert -i $tableselbiom -o $otable --to-tsv\n";
`sh $osh`;




