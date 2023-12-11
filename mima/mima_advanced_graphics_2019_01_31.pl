#!/usr/bin/perl -w
use strict;

##generate the selected groups files 
die "perl $0 <taxa.tsv> <meta.data.tsv> <Group_Vectors_to_Compare> <output_dir> <outputprefix> <top_N>
	taxa.tsv tab	seperated taxa abundance table to draw stack barplot
	meta.dta.tsv	meta data file for samples, row as samples; columm as meta data
	Group_Vectors_to_Compare	The Groups to show taxa; format syntaxa as Column_name:group1,group2,groop3
	output_dir	the output direcotry of the PDF format output
	top_N	Draw the top N taxa; summup the left taxa as Others\n" unless (@ARGV == 6);

##read the vector to decide which samples are chosen 
my $VEC = $ARGV[2]; ##The vector should metadatacolumn:Group1:Group2:Group3 
my $META = $ARGV[1];
my $TABLE = $ARGV[0];
my $oprefix = $ARGV[4];
unless(-d "$ARGV[3]"){
        `mkdir $ARGV[3]`;
}

my $ODIR = $ARGV[3];
my $taxatable = "$ODIR/$oprefix.select.table.txt";
my $taxaaveragetable = "$ODIR/$oprefix.select.average.table.txt";
my $taxatopN = "$ODIR/$oprefix.select.table_top_$ARGV[5].txt";
my $taxaaveragetopN = "$ODIR/$oprefix.select.average.top_$ARGV[5].table.txt";
my $taxatopNseperate = "$ODIR/$oprefix.select.table_top_$ARGV[5]_seperategroup.txt";
my $taxaaveragetopNseperate = "$ODIR/$oprefix.select.average.top_$ARGV[5]_seperategroup.table.txt";
my $metatable = "$ODIR/$oprefix.select.meta.txt";

format_taxa_table($META,$VEC,$metatable,$TABLE, $taxatable, $taxaaveragetable);
top_N($taxatable, $ARGV[5], "$taxatopN");
top_N($taxaaveragetable, $ARGV[5], "$taxaaveragetopN");
top_N_seperate($taxaaveragetable, $taxatable, $ARGV[5],$taxaaveragetopNseperate, $taxatopNseperate );

my $abar = "$ODIR/$oprefix.select.average.top_$ARGV[5]_seperategroup.barchart.pdf";
my $sidbar = "$ODIR/$oprefix.select.top_$ARGV[5]_seperategroup.barchart.pdf";
my $merabar = "$ODIR/$oprefix.select.average.top_$ARGV[5]_.barchart.pdf";
my $mertaxbar = "$ODIR/$oprefix.select.top_$ARGV[5]_.barchart.pdf";
draw_stackbar($taxaaveragetopNseperate, $ARGV[5], $abar);
draw_stackbar($taxatopNseperate, $ARGV[5], $sidbar);
draw_stackbar($taxaaveragetopN, $ARGV[5], $merabar);
draw_stackbar($taxatopN, $ARGV[5], $mertaxbar);

sub draw_stackbar{

	my ($infile, $n, $output) = @_;

my $rscript = "$infile.$n.R";
die "$!\n" unless open(RS, ">$rscript");
my $text =<<RSP;
library(ggplot2)
library(reshape2)
otable <- read.table(file="$infile", header = TRUE, sep = "\\t")
otable <- setNames(melt(otable), c('rows', 'vars', 'values'))
order <- aggregate(otable\$values,by=list(otable\$rows),sum)
order <- order[order(order\$x),]
otable\$rows <- factor(otable\$rows,levels = order\$Group.1)
otable <- otable[order(otable\$rows),]
colors<- c("#B2DFEE","#FFEBCD","#6E8B3D","#87CEFA","#FF6A6A","#377EB8","#E41A1C","#EE9A49","#EED2EE","#EEEE00","#00E5EE","#EED8AE","#EE4A8C","#00EE76","#EE5C42","#5CACEE","#8B5A2B","#CD853F","#CD85CD","#CDCD00","#00C5CD","#CD8A96","#CD3278","#00CD66","#E41A1C","#4F94CD","#BC41A4","#CD4F39","#FFED6F","#CCEBC5","#BC80BD","#FCCDE5","#FDB462","#FB8072","#BEBADA","#377eb8","#4daf4a","#8DD3C7","#FFFFB3","#B3DE69","#80B1D3","#000000","#fbb040","#be1e2d","#39b54a","#ffe600","#ee2a7b", "#00aeef", "#662d91")
col <- rev(rev(colors)[1:length(order\$Group.1)])
pdf("$output", width=12,  height=12)
ggplot(otable, aes(x=vars, y = values, fill =rows)) + geom_bar(position = "stack", stat = "identity") +  xlab("Samples") + ylab("Percentage") + theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(size=0,color="white"),panel.grid.minor=element_line(size=0,color="white"),legend.key = element_rect(fill = "white"), axis.text.x = element_text(angle = 90)) + scale_colour_manual(name="Microbial Taxa", values=col,aesthetics = c("fill"))  
dev.off()
RSP
print RS $text; 

#  + scale_fill_brewer(palette="Set3") + 
`R CMD BATCH $rscript`;
}

sub format_taxa_table{

        #This function processing meta data file to fectch a small subset of samples ids by category information        
        #column_lable:Group1:Group2
        #Meta data file
        #The OUTPUT of this function is a file with the selected Sample IDs and a hash with sampleid->group information 

        my ($meta, $vector, $smetas, $table, $allselect, $averageselect) = @_;

        my ($columnid, @group) = split(/\:/, $vector);
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
	my %sid2group;
        die "$!\n" unless open(SMETA, ">$smetas");
	print SMETA "$head\n";
        while(<META>){
                chomp;
                my @tem = split /\t/;
                if(exists $grouphash{$tem[$index]}){
                        ##select this sample
			print SMETA "$_\n";
			$sid2group{$tem[0]} = $tem[$index];
                }
        }
        close META;
        close SMETA;


	die "$!\n" unless open(TABLE, "$table");
	my $shead = <TABLE>; chomp($shead);
	my @shead = split("\t", $shead);

	my @selectsamples = (); ##this is to label samples to be selected 
	my @selectidname = (); 
	my %index2group; #record the index id to group  information 

	my %average;  # taxa -> group -> SID  this hash store the whole selected matrix into hash style to generate the group average abundance table 
	for(my $i=1; $i <= $#shead; $i++){		
		if(exists $sid2group{$shead[$i]}){
			push @selectsamples, $i;
			push @selectidname, $shead[$i];
				
			die "line 85 \n" unless(exists  $sid2group{$shead[$i]});
			$index2group{$i} = $sid2group{$shead[$i]}; ##
		}
				
	}
	my $selectidtablehead =  join("\t", "", @selectidname);
	my @gaverhead = sort keys %grouphash;
	my $averhead = join("\t", "", @gaverhead);

	my $otable1 = $allselect;
	my $otable2 = $averageselect;
	die "$!\n" unless open(OT1,">$otable1" );
	die "$!\n" unless open(OT2,">$otable2" );
	print OT1 "$selectidtablehead\n";
	print OT2 "$averhead\n";

	while(<TABLE>){
		chomp;
		my @tem = split("\t", $_);
		print OT1 "$tem[0]";
		for(my $i=0; $i <= $#selectsamples; $i++){
			my $indexid = $selectsamples[$i];
			my $sampleid  = $shead[$indexid];
			print OT1 "\t$tem[$indexid]";
			my $groupinfo = $index2group{$indexid};
			#print "$groupinfo\n";
			$average{$tem[0]}{$groupinfo}{$sampleid} = $tem[$indexid]; ##taxa -> group -> sid = abundance 
		}
		print OT1 "\n";
	}
	close TABLE;
	close OT1;
	
	#output average for each group, select top N taxa 	
	for my $tax (sort keys %average){
		print OT2 $tax;
		for my $group (sort keys %{$average{$tax}}){
				
			my $sum = 0;
			my $num = 0;
			for my $sid (keys %{$average{$tax}{$group}}){
				$sum += $average{$tax}{$group}{$sid};
				$num ++;
			}
			
			die "$!\n" if($num == 0);
			my $ave = $sum / $num; 
			print OT2  "\t$ave";
		}
		print OT2 "\n";
	}	
	close OT2;
}

sub top_N{  ##overal average abundance top N 
	#the input is a data matrix
	#the output is the top N abundant feature 
	
	my ($infile, $n, $ofile) = @_; 	
	die "$!\n" unless open(IN, "$infile"); 
	die "$!\n" unless open(OT, ">$ofile"); 
	
	my %averageabu;
	my $head = <IN>;
	print OT $head;
	while(<IN>){
		chomp;
		my @tem = split /\t/;
		my $sum = 0;
		for(my $i=1;$i<=$#tem; $i++){ $sum += $tem[$i];}
		$averageabu{$tem[0]} = $sum;
	}
	close IN;
	
	my %topN;
	my $index =1;
	for my $k (sort {$averageabu{$b} <=> $averageabu{$a}} keys %averageabu ){
		if($index <= $n){
			$topN{$k} = 1; #select for the top N 	
		}
		$index ++;
	}

	##generate the others taxa 

	die "$!\n" unless open(INPUT, "$infile"); 
	my @sumarray = ();
	$sumarray[0] = "Others";
	<INPUT>;
	while(<INPUT>){
		chomp;
		my @tem = split /\t/;
		if(exists $topN{$tem[0]}){
			print OT  "$_\n";
			for(my $i=1; $i <= $#tem; $i++){
				$sumarray[$i] += $tem[$i];
			}
		}
	}
	close INPUT;

	##calculate the Others as 1 - sum of all the other taxa relative abundance 
	for(my $i=1; $i <= $#sumarray; $i++){
		$sumarray[$i] = 1 - $sumarray[$i];
	}
	my $other = join("\t", @sumarray);
	print OT "$other\n";
	close OT;
}

sub top_N_seperate{
	##This is to generate a top N for each group, the input is the average table for each group and the overal table 
	#the output is the top_N_seperate to average and single samples 
	
	my ($avetable, $sidtable, $n, $avertop, $sidtop) = @_;

	my %selecttaxa;
	my %storematrix; #group->taxa = value
	die "$!\n" unless open(IN, "$avetable");
	my $head = <IN>; chomp($head);
	my @headarray  = split("\t", $head);

	while(<IN>){
		chomp;
		my @tem = split /\t/;
		for(my $i=1; $i <= $#tem; $i++){
			$storematrix{$headarray[$i]}{$tem[0]} = $tem[$i];
		}
	
	}
	close IN;

	for my $groupname (sort keys %storematrix){
		my $index =1;
		for my $tax (sort {$storematrix{$groupname}{$b} <=> $storematrix{$groupname}{$a} }  keys %{$storematrix{$groupname}}){  ##sort by value and keep the top N taxa 
			if($index <= $n){
				$selecttaxa{$tax} = 1; ##select this taxa as it is the top N in this group 
			}
			$index ++;
		}
	}

	##selecttaxa store all the taxa need to keep 
	
	die "$!\n" unless open(IN1, "$avetable");
	die "$!\n" unless open(IN2, "$sidtable");
	die "$!\n" unless open(OUT1, ">$avertop");
	die "$!\n" unless open(OUT2, ">$sidtop");


	my @sumarray1 = ();
	$sumarray1[0] = "Others";
	my @sumarray2 = ();
	$sumarray2[0] = "Others";

	my $hIN1 = <IN1>;
	print OUT1 $hIN1;
	while(<IN1>){
		chomp;
		my @tem = split /\t/;
		if(exists $selecttaxa{$tem[0]}){
			print OUT1 "$_\n";		
			for(my $i=1; $i <= $#tem; $i++){
				$sumarray1[$i] += $tem[$i];
			}
		}
	}
	close IN1;

	my $hIN2 = <IN2>;
	print OUT2 $hIN2;
	while(<IN2>){
		chomp;
		my @tem = split /\t/;
		if(exists $selecttaxa{$tem[0]}){
			print OUT2 "$_\n";		
			for(my $i=1; $i <= $#tem; $i++){
				$sumarray2[$i] += $tem[$i];
			}
		}
	}
	close IN2;
	
	##calculate the Others as 1 - sum of all the other taxa relative abundance 
	for(my $i=1; $i <= $#sumarray1; $i++){
		$sumarray1[$i] = 1 - $sumarray1[$i];
	}
	my $other1 = join("\t", @sumarray1);
	print OUT1 "$other1\n";
	close OUT1;
	
	for(my $i=1; $i <= $#sumarray2; $i++){
		$sumarray2[$i] = 1 - $sumarray2[$i];
	}
	my $other2 = join("\t", @sumarray2);
	print OUT2 "$other2\n";
	close OUT2;
}

