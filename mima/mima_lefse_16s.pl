#!/usr/bin/perl -w
#use strict;

#this pipeline could process a subset of samples in a study by the selected meta data file and investigate on specific factor one time 
#developed by Dr. Xiao-Tao JIANG
die "perl $0 <table.qza> <taxonomy.qza> <Investigate_category> <select_meta.tsv> <output_dir>
        <table.qza>      qiime2 feature output table
        <taxonomy.qza>   qiime2 feature taxonomy annotation
        <Investigate_category>  the column name in the meta data file chosed to do comparison
	<Selected_meta data> <meta data file>
        <output_dir>     output directory include all the annotation files

\n" unless (@ARGV==5);

##in order to use this pipeline
#the qiime should be installed 
#ampliconanalysis env should be there, which including lefse and maAslin and others software

##filter samples for differential abundance analysis, processing LEfSe, output all LEfSe combination 
#qiime feature-table filter-samples --i-table ../table.qza --m-metadata-file sample.metadata.tsv --output-dir filterstool
#for both taxonomy and feature level 

#using qiime to export the data 
my $odir = $ARGV[4];
unless(-d "$odir"){ `mkdir $odir`;}
my $reptaxa = $ARGV[1];
my $smeta = $ARGV[2];
my $featuretable = $ARGV[0];
my $selectmetadata   = $ARGV[3];
my $filtertable = "$odir/filter.table.qza";
my $taxaqzv = "$odir/filter.taxa.qzv";

my $qiime = "$odir/qiime.sh";
die "$!\n" unless open(QIM, ">$qiime");
#print QIM "source activate qiime2-2018.11\nexport PYTHONPATH=/srv/scratch/z3524677/bin/miniconda3/envs/qiime2-2018.8/lib/python3.5/site-packages\n";
#print QIM "source activate /srv/scratch/mrcbio/bin/miniconda3/envs/qiime2-2018.11\n";
print QIM "qiime feature-table filter-samples --i-table $featuretable --m-metadata-file $selectmetadata --o-filtered-table $filtertable\n";
print QIM "qiime taxa barplot   --i-table $filtertable --i-taxonomy $reptaxa --m-metadata-file $selectmetadata --o-visualization $taxaqzv\n";
print QIM "qiime tools export --input-path $taxaqzv --output-path $odir/taxarank\n";
#print QIM "source deactivate\n";
close QIM;
`sh  $qiime`;

##testing all ranks from kingdom to species 
for(my $i=1; $i<=7; $i++){
	lefse_table("$odir/taxarank/level-$i.csv", $smeta, "$odir/taxa-rank-$i.txt");
	my $outdir = "$odir/lefse-rank$i";
	unless(-d $outdir){`mkdir $outdir`;}
	##pay special attention to the directory structure 
}
#lefse("$odir/taxa-rank-$i.txt", $outdir);
##testing at feature level for LEfSe
`mkdir $odir/lefse_allrank`;
my $ofile = "$odir/lefse_allrank/merge_taxonomy.csv";

die "$!" unless open(MG, ">$ofile");
die "$!" unless open(TEMP, "$odir/taxa-rank-1.txt");
my $h1 = <TEMP>; my $h2 = <TEMP>;
print MG $h1;
print MG $h2;
close TEMP;

for(my $i=2; $i<=7; $i++){
	die "$!" unless open(TMP, "$odir/taxa-rank-$i.txt");
	<TMP>; <TMP>;
	while(<TMP>){
		chomp;
		my @tem = split /\t/;
		$tem[0] =~ s/\;/\./g;
		my $o = join("\t", @tem);
		print MG "$o\n";
	}
	close TMP;
}

#lefse("$odir/lefse_allrank/merge_taxonomy.csv", "$odir/lefse_allrank");


=head 
##preform MaAsLin analysis for the meta data and tables
##generating the pcl file and the read.config file 
my $mergetaxa = "$odir/merge_level2-7.tsv";
`head -2 $odir/taxa-rank-1.txt | tail -1 > $mergetaxa`;
for(my $i=2; $i<=7; $i++){
	my $onormal = "$odir/taxa-normalize.txt";
	normalize("$odir/taxa-rank-$i.txt", $onormal);
	`tail -n+3 $onormal >> $mergetaxa`;
}

my $o_pcl = "$odir/maAslin.pcl";
my $o_config = "$odir/maAslin.read.config";
maAsLin_table($mergetaxa, $selectmetadata, $o_pcl, $o_config);
my $maaslinsh = "$odir/maaslin.sh";
my $maaslinr = "$odir/maaslin.R";
my $maaslinodir = "$odir/maaslin_output";
die "$!\n" unless open(MAS, ">$maaslinsh");
die "$!\n" unless open(RS, ">$maaslinr"); 
#print MAS "source activate ampliconanalysis\n";
#\nexport PYTHONPATH=/home/z3524677/.local/lib/python2.7/site-packages\n";
print RS  "library(Maaslin)\n";
print RS  "Maaslin(\'$o_pcl\',\'$maaslinodir\',strInputConfig=\'$o_config\')\n";
close RS;
print MAS "R CMD BATCH  $maaslinr\n";
close MAS;
`sh $maaslinsh`;
=cut 

##process class level and testing significant difference for different meta data 
#generate_lefsetable("level-3.csv", "Detail_TYPE", "ofile.txt");
#lefse_command("ofile.txt", $odir);

sub lefse_command {
	##this function run all the program
	#in order to run this function, the ampliconanalysis  should be installed which including the lefse init, in the ampliconanalysis package, the lefse, maAslin, vegan, ggplot2, phyloseq packages were installed.  
	#user should  enter into the ampliconanalysis env 
	##Reference for this too;: 
	#Genome Biol. 2011; 12(6): R60.
	
	my ($infile, $dir) = @_;
	my @m = split("/", $infile);
	my $prefix = join("/", $dir, $m[-1]);
	$prefix =~s/\.txt$//;	

	my $scrip = "$dir/runlefse.sh";
	die "$!\n" unless open(TEM, ">$scrip");
	#print TEM "source activate ampliconanalysis\nexport PYTHONPATH=/home/z3524677/.local/lib/python2.7/site-packages\n";
	#print TEM ". activate lefse-conda\n";
	print TEM "lefse-format_input.py $infile $prefix.format -c 1 -u 2 -o 1000000\n";
	print TEM "run_lefse.py $prefix.format $prefix.res\n";
	print TEM "lefse-plot_res.py $prefix.res $prefix.pdf --format pdf\n";
	print TEM "lefse-plot_cladogram.py $prefix.res $prefix.clade.pdf --format pdf\n"; 
	print TEM "lefse-plot_features.py $prefix.format $prefix.res $prefix.feature.zip --format pdf --archive zip\n";			
=head 	
	#print TEM "source deactivate\n";
	`. activate lefse-conda && lefse-format_input.py $infile $prefix.format -c 1 -u 2 -o 1000000`;
	`. activate lefse-conda && run_lefse.py $prefix.format $prefix.res`;
	`. activate lefse-conda && lefse-plot_res.py $prefix.res $prefix.pdf --format pdf`;
	`. activate lefse-conda && lefse-plot_cladogram.py $prefix.res $prefix.clade.pdf --format pdf`; 
	`. activate lefse-conda && lefse-plot_features.py $prefix.format $prefix.res $prefix.feature.zip --format pdf --archive zip`;			
=cut
	close TEM;	
	#`. activate lefse-conda && sh  $scrip`;
}


sub lefse_table {

	my ($table, $meta, $ofile) = @_;
	##read the qiime2 output taxonomy table and generate lefse format input table with class and id 
	die "$!\n" unless open(TEM, "$table");

	my $head = <TEM>; chomp($head);
	my @name = split(",", $head);
	my $index = 1;
	for(my $i=2; $i <= $#name; $i++){ if($name[$i] =~ m/^k__/){ $index++;}else{ last;}  }
	
	my @metaselect = ();
	my $sid;
	for(my $i=1; $i <= $#name; $i++){ if($name[$i] eq $meta){ $sid = $i; push @metaselect, $name[$i]; }  }
	
	my @rows = ();
	my @transposed = ();
	#print "$sid\t$index\n";
	
	for(my $i =0; $i <= $index; $i++){
		$rows[0][$i] = $name[$i];
	}
	my $inx = 1;
	while(<TEM>){
		chomp;
		my @tem = split(",", $_);
		for(my $i =0; $i <= $index; $i++){
			$rows[$inx][$i] = $tem[$i];
		}
		$inx ++;	
		die "$_\n" unless(exists $tem[$sid]);
		push @metaselect, $tem[$sid];
	}
	close TEM;
	
	for my $row (@rows) {
 		for my $column (0 .. $#{$row}) {
    			push(@{$transposed[$column]}, $row->[$column]);
  		}
	}
		
	die "$!\n"  unless open(OF, ">$ofile");
	#die "file exists\n" if(-e $ofile);
	my $ohead = join("\t", @metaselect);
	print OF  "$ohead\n";
	for my $line (@transposed){
		#for my $element (0..$#{$line}){
			#print OF $line->[$element],"\t";
		#}
		$$line[0] =~ s/\; /\./g;
		my $oline = join("\t", @$line);
		print OF "$oline\n";
	}
	close OF;
}##This function merge 


sub maAsLin_table{
	##This function is to generate input tables for maAsLin 
	#Refer the reference paper for this file Morgan XC, Tickle TL, Sokol H, Gevers D, Devaney KL, Ward DV, Reyes JA, Shah SA, LeLeiko N, Snapper SB, Bousvaros A, Korzenik J, Sands BE, Xavier RJ, Huttenhower C. Dysfunction of the intestinal microbiome in inflammatory bowel disease and treatment. Genome Biol. 2012 Apr 16;13(9):R79.
	my ($ftable, $meta, $opcl, $oconfig) = @_;	
		
	die "$!\n" unless open(PCL, ">$opcl");
	die "$!\n" unless open(CONFIG, ">$oconfig");
	die "$!\n" unless open(TABLE, "$ftable");
	die "$!\n" unless open(META, "$meta");
	
	##As the table is already selected, hence use all the samples IDs in the read.config
	##for PCL
	my $thead = <TABLE>; chomp($thead);
	my @sid = split(/\t/, $thead);
	my $mhead = <META>; chomp($mhead);
	$mhead =~ s/\r//g;
	my @mid = split(/\t/, $mhead);
	my %meta2sid;
	while(<META>){
		chomp;
		s/\r//g;
		my @tem = split /\t/;
		for(my $i =1; $i <= $#tem; $i++){
			$meta2sid{$mid[$i]}{$tem[0]} = $tem[$i];
		}
	}		
	close META;
	
	##output meta data to pcl
	print PCL "$thead\n";
	
	my @mids = sort keys %meta2sid;
	for my $mta (@mids){
		print PCL "$mta";
		for(my $i=1; $i<=$#sid; $i++){
			die "$mta\n" unless(exists $meta2sid{$mta}{$sid[$i]});
			print PCL "\t",$meta2sid{$mta}{$sid[$i]};
		}
		print PCL "\n";
	}

	while(<TABLE>){
		print PCL "$_";
	}
	close TABLE;
	close PCL;

	print CONFIG "Matrix: Metadata\n";
	shift @sid;
	my $idslist = join(",", @sid);
	print CONFIG "Read_PCL_Columns: $idslist\n"; 
	print CONFIG "Read_PCL_Rows: -$mid[1]\n";	
	
	print CONFIG "\n\n";
	print CONFIG "Matrix: Abundance\n";
	print CONFIG "Read_PCL_Columns: $idslist\n";
	print CONFIG "Read_PCL_Rows: $mid[-1]-\n";	
	close CONFIG;
}##maAslin pcl config table 

sub normalize {
	my ($infile, $ofile) = @_;
	
	die "$!\n" unless open(INF, "$infile");
	die "$!\n" unless open(OUT, ">$ofile");
	
	my @sum = ();$sum[0] = 0;
	<INF>; <INF>;
	while(<INF>){
		chomp;	
		my @m = split /\t/;
		for(my $i =1; $i <=$#m ;$i++){
			$sum[$i] += $m[$i];
		}
	
	}
	close INF;

	die "$!\n" unless open(INF, "$infile");
	my $lin1 = <INF>; 
	my $lin2 = <INF>;
	print OUT "$lin1$lin2";
	while(<INF>){
		chomp;	
		my @m = split /\t/;
		print OUT "$m[0]";
		for(my $i =1; $i <=$#m ;$i++){
			print OUT "\t", $m[$i]/$sum[$i];
		}
		print OUT "\n";
	}
	close INF;

}##normalize 

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
			print TOUT "$_\n";
		}
	}
	close TEM;close TOUT;
	print LOG "lefse-plot_res.py $filres $resplot --format svg --left_space 0.35\n";
	print LOG "lefse-plot_cladogram.py $filres $cladoplot --format pdf --max_lev 7\n";
	print LOG "lefse-plot_features.py $fmt $filres $zipfeature --archive zip -f diff --width 12 --height 10 --dpi 300 --format png\n";
	`lefse-plot_res.py $filres $resplot --format pdf --width $w --height $h`;
	`lefse-plot_cladogram.py $filres $cladoplot --format pdf --max_lev 7`;
	`lefse-plot_features.py $fmt $filres $zipfeature --archive zip -f diff --width 15 --height 10 --dpi 300 --format png`;

	##generating the signature table by the names
}
