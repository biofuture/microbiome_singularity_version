#!/usr/bin/perl -w
#Author: Huimin Zheng
#Date: 2019.12.16
#Contact: zhenghuimin91@126.com
#Copyright: Hong-Wei Zhou
use strict;
use Getopt::Std;

my $usage="perl $0 -m <metadata.txt> -i <table.txt> -g <group> -o <output_dir> -a <add-pseudocount>

e.g. perl $0 -m metadata.txt -i table.txt -g Group -a 1 -o output
e.g. perl $0 -m metadata.txt -i table.txt -g Group -a 1 -c Control -o output
e.g. perl $0 -m metadata.txt -i table.txt -g Group -a 1 -s Group:CaseA,Control -c Control -o output
e.g. perl $0 -m metadata.txt -i table.txt -g Group -a 1 -s Age:12w,20w -c Control -o output

        -m <path of metadata.txt>    Required. Each row represent a sample. Sample IDs should be in first column. 
        -i <path of table.txt>        Required. Each column represent a sample. 
        -o <path of output_dir> Required.
        -g <group or treatment>	Required.	Do statistic analysis to compare among groups. For two groups, wilcoxon rank sum test will be applied and output volcano plot. For more than two groups, kruskal wallis test and dunn test will be applied. 
        -a <add-pseudocount>	Required. To deal with 0 when calculate log2 fold of change. Please do not use '1' if input table is relative abundance instead of sequence counts.
        -s [select samples]	Select samples according to metadata. e.g. If only want to compare case and control in certain age. Use '-s Age:12w,20w' to keep samples where the value for 'Age' in the mapping file is '12w','20w'.
        -c [group of control]	log2 fold of change will be calculated based on defined control group. default by ASCII order of group name.
        -p [suppress output of boxplots]	Only 'All', 'Nonsig', 'None'. Use 'All' to suppress all boxplots or use 'Nonsig' to suppress boxplots for non-significant features according to original p values. default Nonsig.
        -q [FDR threshold]	range (0,1]. default 0.05.
        -e [Threshold of effect size for barplot]	Only features with Hodges-Lehmann > threshold will be keeped in barplot. default 0.
        -l [min prevalence]	Range [0,1]. Filter features based on min prevalence. Require at least total samples * min prevalence values. default 0.
        -n [features annotations]	Add features annotations to statistic summary tables. e.g. taxonomy for OTUs. File should contains 2 columns, first column is Features ID, second column is annotations. 
        -h [help]

Author: Huimin Zheng, E-mail: zhenghuimin91\@126.com or 328093402\@qq.com

Version: v1.0

Update: 2019/12/16

Notes:
Required softwares: R
Required R packages: ggplot2, ggpubr, EnhancedVolcano, dunn.test
";

use vars qw($opt_m $opt_i $opt_o $opt_g $opt_a $opt_s $opt_c $opt_p $opt_q $opt_e $opt_l $opt_n $opt_h);
getopts('m:i:o:g:a:s:c:p:q:e:l:n:h');

if ($opt_h){       
	print $usage,"\n";
        exit(1);
}

unless ($opt_m and $opt_i and $opt_g and $opt_o and $opt_a){
        print "$usage\n";
        exit(1);
}

if (not $opt_c){
        $opt_c = "Control";
}

if (not $opt_l){
        $opt_l = 0;
}

`mkdir $opt_o`;

my @select;
my @selected_groups;

if (not $opt_s){
	$opt_s = "select all";
}else{
	@select = split(/:/,$opt_s);
	@selected_groups = split(/,/,$select[1]);
	$opt_s = "filter";
}

if (not $opt_p){
        $opt_p = "Nonsig";
}elsif($opt_p ne "None" and $opt_p ne "Nonsig" and $opt_p ne "All"){
	print "$usage\nPlease make sure value of -p is None, All or Nonsig.\n";
        exit(1);
}

if (not $opt_q){
        $opt_q = 0.05;
}

if (not $opt_e){
        $opt_e = 0;
}

#Get the location of this script, since it has used other scripts, we need to get the location of all of them
use File::Spec;

my $path_script=File::Spec->rel2abs(__FILE__);
my @path_script=split(/\//,$path_script);
pop @path_script;
my $dir_script=join("/",@path_script);
#Get the location of this script finished


## main #######################################
my %taxa;# $taxa{sampleID}{taxa} = value;
my %taxa_name;# $taxa{taxa_name} = simplized taxa name;
my %taxa_name_r;# $taxa{simplized taxa name} = taxa name; for deduplicate simplized taxa name
my %map;# $map{$sampleID} = class/group;

### read metadata and select samples needed
open M,"$opt_m" or die "can't open $opt_m\n";
#### check groups to compare or select are in metadata
my $index = -1;
my $select_index = -1;
chomp(my $map_header = <M>);
my @map_header = split(/\t/,$map_header);
for my $i (1..$#map_header){
        if ($map_header[$i] eq $opt_g){
                $index = $i;
        }
	if($opt_s eq "filter" and $map_header[$i] eq $select[0]){
		$select_index = $i;
	}
}

if ($index == -1){
    die "Please make sure $opt_g is in $opt_m\n";
}

if ($opt_s eq "filter" and $select_index == -1){
	die "Please check parameter of '-s' and make sure $select[0] is in $opt_m\n";
}
#### store samples selected in %map
while(<M>){
        chomp;
        my @map_line = split(/\t/);
	if ($opt_s eq "select all" or (grep {$_ eq $map_line[$select_index]} @selected_groups)){
	        $map{$map_line[0]} = $map_line[$index];
	}
}
close M;


### read features table
open IN, "$opt_i" or die "can't open $opt_i\n";

my $matched = 0; my $check_lines = 0;
my $taxa_header;my @taxa_header;
while($matched == 0 and $check_lines <5){
	chomp($taxa_header = <IN>);
	@taxa_header = split(/\t/,$taxa_header);
	## count matched samples to find the sample ID line. stop if it's not in 5 lines
	foreach my $id (1..$#taxa_header){
		$matched += 1  if exists $map{$taxa_header[$id]};
	}
	$check_lines++;
#	print "$check_lines\t$matched";
}

while(<IN>){
	chomp;
	my @taxa_line = split(/\t/);
	for my $i(1..$#taxa_header){
		$taxa{$taxa_header[$i]}{$taxa_line[0]} = $taxa_line[$i];
	}

	## simplize taxa name
	my @taxa_name = split(/;/,$taxa_line[0]);
	#my $taxa_phylum = $taxa_name[1];
	my $taxa_last = $taxa_name[-1];
	while($taxa_name[-1] =~ /__$/ or $taxa_name[-1] =~ /Other$/){
		pop @taxa_name;
	}
	my $taxa_name = pop @taxa_name;
	if ($taxa_name ne $taxa_last){
		$taxa_name = "$taxa_name\.$taxa_last";
	}
	## my $taxa_name = $taxa_line[0]; ## or original taxa name 
	if(exists $taxa_name_r{$taxa_name}){
		my $upper_level = pop @taxa_name;
		$taxa_name = "$upper_level\.$taxa_name";
	}
	$taxa_name{$taxa_line[0]} = $taxa_name;
	$taxa_name_r{$taxa_name} = $taxa_line[0];
	
}
close IN;

## write out table as input of R script
open OUT ,">$opt_o/$opt_g.input_table.txt" or die "can't create file $opt_o/$opt_g.input_table.txt\n";

print OUT "SID\tGroup";
foreach my $taxa_keys(sort keys %taxa_name){
	print OUT "\t$taxa_name{$taxa_keys}";
}
print OUT "\n";


foreach my $SID(keys %map){
	if (exists $taxa{$SID}){
		print OUT "$SID\t$map{$SID}";
		foreach my $taxa_keys(sort keys %taxa_name){
			print OUT "\t$taxa{$SID}{$taxa_keys}";
		}
		print OUT "\n";
	}
}
close OUT;


## write out taxa name mapping file. To trace original taxa name.
if (not $opt_n){
	open OUT2, ">$opt_o/features_annotation.txt" or die "can't create $opt_o/features_annotation.txt\n";
	print OUT2 "FeaturesID\tAnnotation\n";
	foreach my $taxa_keys(sort keys %taxa_name){
        	print OUT2 "$taxa_name{$taxa_keys}\t$taxa_keys\n";
	}
	close OUT2;
}else{
	`cp $opt_n $opt_o/features_annotation.txt`;
}

## Write R code
open R,">$opt_o/Rplot.r" or die "can't create";
print R "
## Readin 
dta <- read.table(\"$opt_o/$opt_g.input_table.txt\",header = T,row.names = 1,sep = \"\\t\",quote=\"\",comment.char=\"\")
anno <- read.table(\"$opt_o/features_annotation.txt\",header = T,row.names = 1,sep = \"\\t\",quote=\"\",comment.char=\"\")
anno <- t(data.frame(t(anno)))

## filter features based on min prevalence
min_samples <- nrow(dta)*$opt_l
dta <- dta[,c(TRUE,colSums(dta[,c(2:ncol(dta))] > 0) > min_samples)]
write.table(dta, \"$opt_o/$opt_g.filterd_table.txt\", sep=\"\\t\",quote = F,row.names=T,col.names=NA)#write table

## Kruskal wallis
KW<-lapply(dta[,c(2:ncol(dta))], function(x) kruskal.test(x~dta\$Group))
kw_pval <-sapply(KW, function(x) x[[\"p.value\"]][[1]][1])
FDR<- p.adjust(kw_pval, \"fdr\")
pval_fdr<-cbind(kw_pval,FDR)#merge pvalue and fdr
mn <- data.frame(t(data.frame(aggregate(.~Group,data=dta,FUN=mean),row.names = 1)))
p_mn <- merge(pval_fdr,mn,by=\"row.names\")
p_mn <- p_mn[order(p_mn\$kw_pval),]
wp_mn <- merge(p_mn,anno,by.x=\"Row.names\",by.y=\"row.names\",all.x=T)
write.table(wp_mn, \"$opt_o/$opt_g.kw.stat_summary.txt\", sep=\"\\t\",quote = F,row.names=F)#write table 

if (length(levels(dta\$Group))>2){
## Dunn test
library(dunn.test)
# list taxa in otu table  
tCol = colnames(dta[,c(2:ncol(dta))])
succeedDunnTest = list()
#loop to run dunn test
out <- matrix(NA, nrow=length(tCol), ncol=length(combn(levels(as.factor(dta\$Group)),2,list))) 
i <- 1
for (t in tCol){ ## overwrites errors and continues loop
  tryCatch(
    {
      output = dunn.test(dta[,t],as.factor(dta\$Group), method = \"BH\", table =TRUE)
      out[i,] <- output\$P.adjusted # prints matrix for pvalues for all treatment comparisons (i) 
      colnames(out) <- output\$comparisons # column names are comparisons
      succeedDunnTest[[i]] <- t
      i <- i + 1 
    },
    error=function(cond) {
      write(sprintf(\"Error processing %s\", t), file = stderr())
    },
    warning=function(cond) {
    },
    finally={
    }
  )    
}
out <- out[!apply(is.na(out),1,all),] # remove empty rows (NA) from failed runs
rownames(out) <- succeedDunnTest # names rows with taxa names
write.table(out, \"$opt_o/$opt_g.pairwise_dunn.test.txt\",quote = F,sep = \"\\t\")
}


wil_plot  <- function(dta,outdir,pseudocount,controlgroup){
  if (length(levels(dta\$Group))==2){
    wil_p_mn = data.frame()
    g = levels(dta\$Group)  
  
    wil=lapply(dta[,c(2:ncol(dta))], function(x) tryCatch({wilcox.test(x~dta\$Group, conf.int=TRUE) },error=function(cond) {r <- wilcox.test(x~dta\$Group)
      r\$estimate <- NA
    return(r)}))
	
    wil_pval = sapply(wil, function(x) x[[\"p.value\"]][[1]][1])
    HodgesLehmann = sapply(wil, function(x) x[[\"estimate\"]][[1]][1])
    wil_FDR = p.adjust(wil_pval, \"fdr\")
    p_fdr = cbind(wil_pval,wil_FDR,HodgesLehmann)#merge pvalue and fdr
    mn = data.frame(t(data.frame(aggregate(.~Group,data=dta,FUN=mean),row.names = 1)))
    wil_p_mn = merge(p_fdr,mn,by=\"row.names\")

    if (names(wil_p_mn)[length(names(wil_p_mn))-1] == controlgroup){
      wil_p_mn\$d = wil_p_mn[,length(names(wil_p_mn))] - wil_p_mn[,length(names(wil_p_mn))-1]
      wil_p_mn\$Log2FC = log2(wil_p_mn[,length(names(wil_p_mn))-1] / wil_p_mn[,length(names(wil_p_mn))-2])
      wil_p_mn\$adjustedLog2FC = log2((wil_p_mn[,length(names(wil_p_mn))-2]+pseudocount) / (wil_p_mn[,length(names(wil_p_mn))-3]+pseudocount))
    }else if (names(wil_p_mn)[length(names(wil_p_mn))] == controlgroup){
      wil_p_mn\$d = wil_p_mn[,length(names(wil_p_mn))-1] - wil_p_mn[,length(names(wil_p_mn))]
      wil_p_mn\$Log2FC = log2(wil_p_mn[,length(names(wil_p_mn))-2] / wil_p_mn[,length(names(wil_p_mn))-1])
      wil_p_mn\$adjustedLog2FC = log2((wil_p_mn[,length(names(wil_p_mn))-3]+pseudocount) / (wil_p_mn[,length(names(wil_p_mn))-2]+pseudocount))
    }else{
      warning(\"Warning message: log2 fold of change will be calculated by default direction\")
      wil_p_mn\$d = wil_p_mn[,length(names(wil_p_mn))] - wil_p_mn[,length(names(wil_p_mn))-1]
      wil_p_mn\$Log2FC = log2(wil_p_mn[,length(names(wil_p_mn))-1] / wil_p_mn[,length(names(wil_p_mn))-2])
      wil_p_mn\$adjustedLog2FC = log2((wil_p_mn[,length(names(wil_p_mn))-2]+pseudocount) / (wil_p_mn[,length(names(wil_p_mn))-3]+pseudocount))
    }

    ## order
    wil_p_mn = wil_p_mn[order(wil_p_mn\$wil_pval),]

    ## correct the direction of HodgesLehmann
    if (sum(sign(wil_p_mn\$HodgesLehmann*wil_p_mn\$d),na.rm = T)<0){
      wil_p_mn\$HodgesLehmann = wil_p_mn\$HodgesLehmann*-1
    }

    ## directions
    wil_p_mn\$HL_direction =  ifelse(sign(wil_p_mn\$HodgesLehmann)>0,\"Increase\",\"Decrease\")
    wil_p_mn\$d_direction =  ifelse(sign(wil_p_mn\$d)>0,\"Increase\",\"Decrease\")
    wwil_p_mn <- merge(wil_p_mn,anno,by.x=\"Row.names\",by.y=\"row.names\",all.x=T)
    write.table(wwil_p_mn, paste(outdir,\"/\",g[1],\"_\",g[2],\".wil.stat_summary.txt\",sep = \"\"), sep=\"\\t\",quote = F,row.names=F)#write table

    ## volcano plot 
    library(EnhancedVolcano)

    EnhancedVolcano(wil_p_mn,lab = wil_p_mn\$Row.names,x = 'd',y = 'wil_pval',pCutoff = 0.05,FCcutoff = 0,xlab = 'change of abundance',
                    transcriptPointSize = 1.5,transcriptLabSize = 3.0,col=c('black', 'black', '#BE1E2D', '#BE1E2D'))
    ggsave(paste(outdir,\"/\",g[1],\"_\",g[2],\".volcanoplot.d.p.pdf\",sep = \"\"))

    EnhancedVolcano(wil_p_mn,lab = wil_p_mn\$Row.names,x = 'd',y = 'wil_FDR',pCutoff = $opt_q,FCcutoff = 0,xlab = 'change of abundance',ylab = '-Log10 adjusted P',
                    transcriptPointSize = 1.5,transcriptLabSize = 3.0,col=c('black', 'black', '#BE1E2D', '#BE1E2D'))
    ggsave(paste(outdir,\"/\",g[1],\"_\",g[2],\".volcanoplot.d.FDR$opt_q.pdf\",sep = \"\"))

    EnhancedVolcano(wil_p_mn,lab = wil_p_mn\$Row.names,x = 'adjustedLog2FC',y = 'wil_FDR',pCutoff = $opt_q,FCcutoff = 1.0,
                    transcriptPointSize = 1.5,transcriptLabSize = 3.0,ylab = '-Log10 adjusted P')
    ggsave(paste(outdir,\"/\",g[1],\"_\",g[2],\".volcanoplot.Log2FC.FDR$opt_q.pdf\",sep = \"\"))

    EnhancedVolcano(wil_p_mn,lab = wil_p_mn\$Row.names,x = 'adjustedLog2FC',y = 'wil_pval',pCutoff = 0.05,FCcutoff = 1.0,
                    transcriptPointSize = 1.5,transcriptLabSize = 3.0)
    ggsave(paste(outdir,\"/\",g[1],\"_\",g[2],\".volcanoplot.Log2FC.p.pdf\",sep = \"\"))

    library(ggplot2)
    p = subset(wil_p_mn,HodgesLehmann!=\"NA\"&wil_FDR<$opt_q&(HodgesLehmann>$opt_e|HodgesLehmann < (-$opt_e)))
    p\$Row.names = factor(p\$Row.names,levels = p\$Row.names[order(p\$HodgesLehmann)])
    ggplot(p,aes(x=Row.names,y=HodgesLehmann,fill=HL_direction))+geom_bar(stat = \"identity\")+coord_flip()+xlab(\"\")+ylab(\"Difference in location (Hodges-Lehmann estimator)\")+theme_classic()+scale_fill_manual(values=c(\"Decrease\"=\"#39b54a\",\"Increase\"=\"#BE1E2D\"))
    ggsave(paste(outdir,\"/\",g[1],\"_\",g[2],\".barplot.FDR$opt_q.HL$opt_e.pdf\",sep = \"\"))

    p = subset(wil_p_mn,HodgesLehmann!=\"NA\"&wil_pval<0.05)
    p\$Row.names = factor(p\$Row.names,levels = p\$Row.names[order(p\$HodgesLehmann)])
    ggplot(p,aes(x=Row.names,y=HodgesLehmann,fill=HL_direction))+geom_bar(stat = \"identity\")+coord_flip()+xlab(\"\")+ylab(\"Difference in location (Hodges-Lehmann estimator)\")+theme_classic()+scale_fill_manual(values=c(\"Decrease\"=\"#39b54a\",\"Increase\"=\"#BE1E2D\"))
    ggsave(paste(outdir,\"/\",g[1],\"_\",g[2],\".barplot.p.pdf\",sep = \"\"))
  }
}

comp <- combn(levels(as.factor(dta\$Group)),2,list)
for (pair in comp){
  sub <- subset(dta,Group %in% pair)
  sub <- droplevels.data.frame(sub)
  if (\"$opt_c\" %in% pair){
    wil_plot(sub,\"$opt_o\",$opt_a,\"$opt_c\")
  }else{
    wil_plot(sub,\"$opt_o\",$opt_a,pair[1])
  }
}


";
if ($opt_p ne "All"){
`mkdir $opt_o/plots`;
`mkdir $opt_o/plots/sig_with_FDR$opt_q`;
`mkdir $opt_o/plots/sig_with_FDR$opt_q/barplot`;
`mkdir $opt_o/plots/sig_with_FDR$opt_q/dotplot`;
`mkdir $opt_o/plots/sig_with_FDR$opt_q/violin`;
`mkdir $opt_o/plots/sig_only_without_FDR_correction`;
`mkdir $opt_o/plots/sig_only_without_FDR_correction/barplot`;
`mkdir $opt_o/plots/sig_only_without_FDR_correction/dotplot`;
`mkdir $opt_o/plots/sig_only_without_FDR_correction/violin`;


if($opt_p ne "Nonsig"){
	`mkdir $opt_o/plots/nonsig`;
	`mkdir $opt_o/plots/nonsig/barplot`;
	`mkdir $opt_o/plots/nonsig/dotplot`;
	`mkdir $opt_o/plots/nonsig/violin`;
}

print R "
## boxplots
library(ggpubr)

sig <- subset(p_mn,p_mn\$FDR<$opt_q)
sig_no_correction <- subset(p_mn,p_mn\$kw_pval<0.05)
nonsig <- subset(p_mn,p_mn\$kw_pval>=0.05)

comp <- combn(levels(as.factor(dta\$Group)),2,list)
tCol = colnames(dta[,c(2:ncol(dta))])

col <- c(\"#B2DFEE\",\"#FFEBCD\",\"#6E8B3D\",\"#87CEFA\",\"#FF6A6A\",\"#377EB8\",\"#E41A1C\",\"#EE9A49\",\"#EED2EE\",\"#EEEE00\",\"#00E5EE\",\"#EED8AE\",\"#EE4A8C\",\"#00EE76\",\"#EE5C42\",\"#5CACEE\",\"#8B5A2B\",\"#CD853F\",\"#CD85CD\",\"#CDCD00\",\"#00C5CD\",\"#CD8A96\",\"#CD3278\",\"#00CD66\",\"#E41A1C\",\"#4F94CD\",\"#BC41A4\",\"#CD4F39\",\"#FFED6F\",\"#CCEBC5\",\"#BC80BD\",\"#FCCDE5\",\"#FDB462\",\"#FB8072\",\"#BEBADA\",\"#377eb8\",\"#4daf4a\",\"#8DD3C7\",\"#FFFFB3\",\"#B3DE69\",\"#80B1D3\",\"#000000\",\"#fbb040\",\"#be1e2d\",\"#39b54a\",\"#ffe600\",\"#ee2a7b\", \"#00aeef\", \"#662d91\")

col <- rev(col)

for (t in tCol){
  a <- strsplit(t,'.',fixed = T)
  
  if (t %in% sig\$Row.names){
    outdirprefix <- \"$opt_o/plots/sig_with_FDR$opt_q/\"
  }else if (t %in% sig_no_correction\$Row.names){
    outdirprefix <- \"$opt_o/plots/sig_only_without_FDR_correction/\"
  }else if (t %in% nonsig\$Row.names){
    outdirprefix <- \"$opt_o/plots/nonsig/\"
  }

  if ((!t %in% nonsig\$Row.names)|(t %in% nonsig\$Row.names & \"$opt_p\" != \"Nonsig\")){
  ggbarplot(dta, x = \"Group\", y = t,error.plot = \"upper_errorbar\",
            add = c(\"mean_se\", \"jitter\"),add.params = list(color=\"black\"),
            fill = \"Group\", palette = col)+
  rotate_x_text(angle = 45)+ # Add horizontal line at base mean
    stat_compare_means(comparisons = comp,label = \"p.signif\", method = \"wilcox.test\",  p.adjust.method = \"BH\")# +  # Pairwise comparison against all
  #stat_compare_means(method = \"kruskal.test\")        # Add global annova p-value
  ggsave(paste(outdirprefix,\"barplot/\",t,\".barplot.png\",sep=\"\"), width = 15, height = 15, units = c(\"cm\"), dpi = 320)

  ggviolin(dta, x = \"Group\", y = t, color = \"black\",fill = \"Group\", font.label = list(size = 30),palette = col,
            add = \"boxplot\", add.params = list(width = 0.2,fill=\"white\"))+
    rotate_x_text(angle = 45)+ # Add horizontal line at base mean
    stat_compare_means(comparisons = comp,label = \"p.signif\", method = \"wilcox.test\",  p.adjust.method = \"BH\")# +  # Pairwise comparison against all
  #stat_compare_means(method = \"kruskal.test\")        # Add global annova p-value
  ggsave(paste(outdirprefix,\"violin/\",t,\".violin.png\",sep=\"\"), width = 15, height = 15, units = c(\"cm\"), dpi = 320)

  ggboxplot(dta, x = \"Group\", y = t, color = \"Group\", font.label = list(size = 30),palette = col,
            add = \"dotplot\", legend = \"none\")  + xlab(\"\")+ylab(a[[1]][length(a[[1]])])+
    rotate_x_text(angle = 45)+ # Add horizontal line at base mean
    stat_compare_means(comparisons = comp,label = \"p.signif\", method = \"wilcox.test\",  p.adjust.method = \"BH\")# +  # Pairwise comparison against all
    #stat_compare_means(method = \"kruskal.test\")        # Add global annova p-value
  ggsave(paste(outdirprefix,\"dotplot/\",t,\".dotplot.png\",sep=\"\"), width = 15, height = 15, units = c(\"cm\"), dpi = 320)
  }

}
";
}

close R;

`Rscript $opt_o/Rplot.r`;
