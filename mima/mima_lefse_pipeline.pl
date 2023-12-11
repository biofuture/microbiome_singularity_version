#!/usr/bin/perl -w
use strict; 

#This is for automatic LEfSe for many groups of comparisions

#dependence LEfSe, biom 

##This script compare all groups in a given enviroment 
#Input1: A big talbe with all samples as column and row as features
#Input2: A meta data file used to do LEfSe analysis
#Input3: A vector indicating the groups chosed to do the analysis: this verctor can be a small fraction of the samples
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);

##Generalize dir for this program
our (@dirset,$scriptpath);
BEGIN {
    @dirset = split(/\//,$Bin);
    $scriptpath = join("/", @dirset);
    unshift @INC, "$scriptpath";
}

our ($opt_h, $opt_i, $opt_o, $opt_g, $opt_m, $opt_f, $opt_d, $opt_x, $opt_y) = "";
my  $usage = <<USE;
        Author: JIANG Xiaotao
	Modify: 2020-07-28
	-i table.tsv, the table used to do LEfSe
	-m meta datafile, a table with rows as sample and column as differnt meta data
	-g group vector should be column_name:group1,group2,group3 #the detail groups in one column to compare 
	-o output dir/ string 
	-f prefix to the output file,  string format 
	-d LDA schore forcut off default 2 
	-x width of figure  default 8 inch
	-y height of figure default 8 inch 
	-h print this help information 
USE

getopts('i:o:g:f:m:d:x:y:h');
if($opt_h  || !$opt_i ){
	die "$usage\n";
}

#die "perl $0 <table.tsv> <meta.data.tsv> <Group_Vectors_to_Compare> <output_dir> <outputprefix>\n" unless (@ARGV == 5);
#
##read the vector to decide which samples are chosen 
my $VEC = $opt_g; ##The vector should metadatacolumn:Group1:Group2:Group3 
my $META = $opt_m;
my $TABLE = $opt_i;
my $oprefix = $opt_f;
$opt_d ||= 2; 
$opt_x ||= 8;
$opt_y ||= 10; 
my $lda = $opt_d;
my $w = $opt_x;
my $h = $opt_y;
unless(-d "$opt_o"){
	`mkdir $opt_o`;
}

my $ODIR = $opt_o;
my $fout1 = "$ODIR/$oprefix.select.table.txt";
my %sidg;
my $olog = "$ODIR/log.txt";
die "$!\n" unless open(LOG, ">$olog");

fetech_samples_group($META,$VEC,$fout1,\%sidg);
my $foutforlefse = "$ODIR/$oprefix.select.lefse.txt";
prepare_LEfSe_tables($TABLE, $fout1, \%sidg, $foutforlefse, $ODIR);
lefse($foutforlefse, $ODIR, $oprefix);
my $signaturetable = "$ODIR/$oprefix.signature.table.forbetadiv.txt";
my $filtersig = "$ODIR/lefse/$oprefix.filter.res";
signature_table($filtersig, $foutforlefse, $signaturetable);

sub signature_table{
	#
	my ($fs, $ft, $fo) = @_;
	#Store all the signatures name and group information 
	my %sig = ();
	die "$!\n" unless open(LEFSE, "$fs");
	die "$!\n" unless open(FO, ">$fo");
	while(<LEFSE>){
		chomp;
		my @tem = split /\s+/;
		#	my @formatname = split("_", $tem[0]);
		#	my $repname = join("-", $formatname[0], $formatname[1]);
		$sig{$tem[0]} = 1;
	}
	close LEFSE;

	die  "$!\n" unless open(FTAB, "$ft");
	my $head = <FTAB>;
	print FO $head;
	while(<FTAB>){
		chomp;
		my @m = split /\t/;
		$m[0] =~ s/\|/\./g;
		#my @n = split(/\:/, $m[0]);
		if(exists $sig{$m[0]}){
			print FO "$_\n";
			delete($sig{$m[0]});
		}
	}
	close FTAB;
	
	my @check = keys %sig;
	if($#check != 0 ){
		for my $k (keys %sig){
			print "$k\t$sig{$k}\n";
		}
	
	} 

}

##run LEfSe analysis for this pair 
sub lefse{
	#
	my ($input, $odir) = @_;
	unless(-d "$odir/lefse"){ `mkdir $odir/lefse`;} 
	my $fmt = "$odir/lefse/$oprefix.fmt";
	my $res = "$odir/lefse/$oprefix.res";
	my $filres = "$odir/lefse/$oprefix.filter.res";
	my $resplot = "$odir/lefse/$oprefix.res.pdf";
	my $cladoplot = "$odir/lefse/$oprefix.clado.pdf";
	my $zipfeature = "$odir/lefse/$oprefix.feature.zip";
	
	print LOG "lefse-format_input.py $input $fmt -c 2 -u 1\n";
	print LOG "run_lefse.py $fmt $res\n";
	`lefse-format_input.py $input $fmt -c 2 -u 1 -o 1000000`;
	#`lefse-format_input.py $input $fmt -c 2 -u 1`;
	`run_lefse.py $fmt $res -l $lda`;
	die "$!\n" unless open(TEM, "$res");
	die "$!\n" unless open(TOUT, ">$filres");
	while(<TEM>){
		chomp;
		my @m = split /\t/;
		my $flag = 1;
		for  my $element (@m){
			if($element eq "" || $element eq "-"){
				$flag =0;
			}
		}

		if($flag == 1){
			$m[0] =~ s/\./dot/g; #replace all . with dot  firstly 
			#print $m[0],"\tdot\n";
			$m[0] =~ s/__/\|/g;
			$m[0] =~ s/-//g; ##remove all dash 
                        $m[0] =~ s/\_p/\.p/g;
                        $m[0] =~ s/\_c/\.c/g;
                        $m[0] =~ s/\_o/\.o/g;
                        $m[0] =~ s/\_f/\.f/g;
                        $m[0] =~ s/\_g/\.g/g;
                        $m[0] =~ s/\_s/\.s/g;
			#$m[0] =~ s/\_/underline/g;
                        $m[0] =~ s/\|/__/g;
			
			#formating $m[0] if it is a 16S output 
			#output the name with nearest name for taxa like s__ g__

			my @tem = split(/\./, $m[0]);
			my @temp = ();
			for(my $i=$#tem; $i >= 0; $i--){
				if($tem[$i] =~ m/\S+__$/){
					unshift  @temp, $tem[$i];
				}else{
					unshift @temp, $tem[$i];
					my $mlast = join("_", @temp);
					my @otemps = (); 
					for(my $j =0; $j <= $i -1; $j ++){
						push @otemps, $tem[$j];
					}
					push @otemps, $mlast;
				       	$m[0] = join(".", @otemps);	
					last;
				}
			}
				
			my $ona = join("\t", @m);
			print TOUT "$ona\n";
		}
	}
	close TEM;close TOUT;
	print LOG "lefse-plot_res.py $filres $resplot --format svg --left_space 0.35\n";
	print LOG "lefse-plot_cladogram.py $filres $cladoplot --format pdf --max_lev 7\n";
	print LOG "lefse-plot_features.py $fmt $filres $zipfeature --archive zip -f diff --width 12 --height 10 --dpi 300 --format png\n";
	`lefse-plot_res.py $filres $resplot --format pdf --width $w --height $h`;
	`lefse-plot_cladogram.py $filres $cladoplot --format pdf --max_lev 11`;
	#`lefse-plot_features.py $fmt $filres $zipfeature --archive zip -f diff --width 15 --height 10 --dpi 300 --format png`;
	##generating the signature table by the names 
}

##convert table into biom and use biom to subset the table and then adding the categry information to the target table for LEfSe analysis 
sub prepare_LEfSe_tables{

	##using biom 
	#conda activate ampliconanalysis 
	
	my ($tsvtable, $selsampleid, $sidg, $tablelefse, $odir) = @_;
	#convert table into biom 
	#die "$odir\n";
	unless(-d "$odir/temp"){`mkdir $odir/temp`;}
	my $tablebiom = "$odir/temp/table.biom";
	my $tableselbiom = "$odir/temp/tablesel.biom";
	my $tablesel = "$odir/temp/tablesel.txt";
	#biom convert -i merge_taxonomy.txt -o test.biom --table-type "Taxon table"
	`biom convert -i $tsvtable -o $tablebiom --table-type \"Taxon table\" --to-hdf5`;
	`biom subset-table -i $tablebiom -a sample -s $selsampleid -o $tableselbiom`;
	`biom convert -i $tableselbiom -o $tablesel --to-tsv`;
	##rm the temp files 
	##adding the group information to LEfSe
	die "$!\n" unless open(ADDMETA, "$tablesel");
	die "$!\n" unless open(LEFSE, ">$tablelefse");
		
	<ADDMETA>; ##shift the header #
	my $head = <ADDMETA>;chomp($head);
	my @head = split(/\t/, $head);
	my @cater = ();
	for(my $i = 1; $i <= $#head; $i++){
		if(exists $$sidg{$head[$i]}){
			push @cater, $$sidg{$head[$i]};
		}
	}
	my $catrow = join("\t", "Cater", @cater);
	
	print LEFSE "$head\n$catrow\n";
	print LEFSE <ADDMETA>;

	close ADDMETA;
	close LEFSE;
	#`rm  $tablebiom; rm $tableselbiom; $rm $tablesel;`;
}

sub fetech_samples_group{
	#This function processing meta data file to fectch a small subset of samples ids by category information	
	#column_lable:Group1:Group2
	#Meta data file
	#The OUTPUT of this function is a file with the selected Sample IDs and a hash with sampleid->group information 

	my ($meta, $vector, $sidfile,$sidgroupmap) = @_;
	
	my ($columnid, $group) = split(/\:/, $vector);
	my @group = split(",", $group); 
	#print "$columnid\n";
	my %grouphash = ();
	for my $k (@group){
		#put all groups info into hash
		#print "$k\n";
		$grouphash{$k} =1;
	}
	die "$!\n" unless open(META, "$meta");
	my $head = <META>;
	chomp($head);
	my @h = split("\t", $head);
	my $index = 0; 
	for(my $i = 0; $i<= $#h; $i++){
		if($h[$i] eq $columnid){
			$index = $i;
			last;
		}
	}
	#print "$index\n";
	die "$!\n" unless open(SAMPLE, ">$sidfile");
	while(<META>){
		chomp;
		my @tem = split /\t/;
		if(exists $grouphash{$tem[$index]}){
			##select this sample
			print SAMPLE "$tem[0]\n";
			${$sidgroupmap}{$tem[0]} = $tem[$index];
		}
	}
	close META;
	close SAMPLE;
}
