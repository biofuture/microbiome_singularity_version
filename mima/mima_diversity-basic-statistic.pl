#!/usr/bin/perl -w 

use strict;
use Getopt::Std;
##This is a designed flow for Miseq PE amplicon data analysis, special designed for V3V4. Users Also can change the database 
##pre-requirement qiime2 installed. 
##source activate qiime2-2018.8
##prepare manefest file
##prepare meta data file
##put the decontaminated demux.qza table.qza rep.qza file into output dir

my $usage = <<USG;
 Author Xiao-tao JIANG
    Email: biofuture.jiang\@gmail.com
    Date 17/06/18
    Modified: 10/12/2018    
	
    perl $0 -m sample.metatdata.tsv -c clasifier.qza -o output_dir -n [threads number] -h [verbose] 
    
    -m  sample_metadata.tsv
    -c  The trianed classifer for your 16S region
    -n  [Threads number used when split file default 6]
    -d  depth of default [5000]	
    -o  output directory
    -h  print this help information

     NOTICE: before running this flow, the table.qza and rep_seqs.qza should be placed uner output_dir and the depth should be also determined
USG

##
our($opt_m, $opt_c, $opt_n, $opt_o, $opt_d, $opt_h)="";
getopts('m:c:n:o:d:h');


if($opt_h || ($opt_m eq ""))
{
    die $usage;
}

my $commandsh = "beta-diversity-commands.sh";
my $loginfo = "log.txt";
die "$!\n" unless open(SHELL, ">$commandsh"); ##store all the shell script to run 
#die "$!\n" unless open(LOG, ">$loginfo"); ##store all the shell script to run 

##define files 
my $samplemetadata = "$opt_m";
my $classifier = "$opt_c";
my $odir = "$opt_o";

##As a standard pipeline, all files in a project keep as the same file name and won't change 

my $demuxqza = "$odir/demux.qza";
my $demuxqzv = "$odir/demux.qzv";
my $repseq = "$odir/rep-seqs.qza";
my $statsdada2 = "$odir/stats-dada2.qza";
my $table = "$odir/table.qza";
my $tablesum = "$odir/table-summarize.qzv";
my $reptaxa = "$odir/taxonomy.qza";
my $taxabar = "$odir/taxa-bar-plots.qzv";
my $alignrep = "$odir/aligned-rep-seqs.qza";
my $maskedalignrep = "$odir/masked-aligned-rep-seqs.qza";
my $unrootedtree = "$odir/unrooted-tree.qza";
my $rootedtree = "$odir/rooted-tree.qza";
my $alphararfac = "$odir/alpha-rarefaction.qzv";
##
my $rarefytable = "$odir/rarefy.table.qza";
my $normtaxabar = "$odir/normalized.taxa-bar-plots.qzv";
my $coretable = "$odir/core.table.qza";

unless (-d $odir){ `mkdir $odir`;} 

=head 
##decide the parameters of trimming ##this shoudl be decided during the qc steps
my $tf = 290;
my $tr = 220;
my $pf = 17;
my $pr = 21;
=cut 

##decide the depth of normalization 
$opt_d ||= 5000; ##This parameters is determined by the QC process 
my $depthnorm = $opt_d;

##--------------------------------------------------------------------------------------------------------Overall Basic analysis for all samples in a study----------------------------------------------------------------------------##
#unless (-d "$odir"){ `mkdir $odir`;}
#Step 0 Enter int QIIME2
#source activate qiime2-2018.4
print SHELL "#Step 0 activate the QIIME2 conda environment\n";
#print SHELL "source activate qiime2-2018.11\n";
#print SHELL "export /srv/scratch/z3524677/bin/miniconda3/envs/qiime2-2018.8/lib/python3.5/site-packages\n"; #The default PYTHONPATH for qiime2 is here with python3 
#print SHELL "#decide the parameters of trimming\n";
##Step 1 import your data 
#print SHELL "\n#Step 1 import pair-end illumina fastq data (For Miseq data, the quality format is Sanger Phred33)\n";
#print SHELL "qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $manifest --output-path $demuxqza --source-format PairedEndFastqManifestPhred33\n";

##Step 2 denoizing with dada2
#print SHELL "\n#Step 2 denoising and merging with dada2\n";
#print SHELL "qiime dada2 denoise-paired --i-demultiplexed-seqs $demuxqza --p-trunc-len-f $tf --p-trunc-len-r $tr --p-trim-left-f $pf --p-trim-left-r $pr  --p-n-threads $opt_n  --o-representative-sequences $repseq  --o-denoising-stats $statsdada2 --o-table $table\n";

##sum up taxonomy after normalization of 
print SHELL "qiime feature-table summarize --i-table $table --m-sample-metadata-file $samplemetadata --o-visualization $tablesum\n";
##summarize feature table 
#print SHELL "qiime demux summarize  --i-data $demuxqza  --o-visualization $demuxqzv\n";

##Step 3 clasify the representative features 
print SHELL "\n#Step 3 taxonomy assignment of feature representative sequences\n";
print SHELL "qiime feature-classifier classify-sklearn  --i-classifier $classifier --i-reads $repseq --o-classification  $reptaxa\n";

##adding taxonomy information to meta data file 
##rarify table.qza firstly and then summarize  
print SHELL "qiime taxa barplot   --i-table $table --i-taxonomy $reptaxa --m-metadata-file $samplemetadata --o-visualization $taxabar\n";
print SHELL "qiime feature-table rarefy --i-table $table --p-sampling-depth $depthnorm --o-rarefied-table $rarefytable\n";
print SHELL "qiime taxa barplot   --i-table $rarefytable --i-taxonomy $reptaxa --m-metadata-file $samplemetadata --o-visualization $normtaxabar\n";

##Step 4 constructing phylogenetic tree of the features 
print SHELL "\n#Step 4 constructing phylogenetic tree\n";
print SHELL "qiime alignment mafft  --i-sequences $repseq --o-alignment $alignrep\n";

print SHELL "qiime alignment mask  --i-alignment $alignrep  --o-masked-alignment $maskedalignrep\n";

print SHELL "qiime phylogeny fasttree  --i-alignment $maskedalignrep  --o-tree $unrootedtree --p-n-threads $opt_n\n";

print SHELL "qiime phylogeny midpoint-root --i-tree  $unrootedtree --o-rooted-tree $rootedtree\n";

##Step 5 core diversity analysis
print SHELL "\n#Step 5  core diversity analysis\n";

##Select samples from the sample.meta.data.tsv
##update the table for analysis
print SHELL "qiime feature-table filter-samples --i-table $table --m-metadata-file $samplemetadata --o-filtered-table $coretable\n";
print SHELL "qiime diversity core-metrics-phylogenetic  --i-phylogeny  $rootedtree  --i-table $coretable --p-sampling-depth $depthnorm --m-metadata-file $samplemetadata   --output-dir $odir/core-metrics-results\n";

##Step 6 alpha diversity analysis, rarefaction curves analysis, alpha group significance analysis 
print SHELL "\n#Step 6 alpha diversity analysis, rarefaction curves analysis, alpha group significance analysis\n";
print SHELL "qiime diversity alpha-rarefaction   --i-table $table   --i-phylogeny $rootedtree   --p-max-depth  $depthnorm   --m-metadata-file $samplemetadata   --o-visualization $alphararfac\n";

unless(-d "$odir/alpha-group"){ `mkdir $odir/alpha-group`;} 

print SHELL "qiime diversity alpha-group-significance   --i-alpha-diversity $odir/core-metrics-results/faith_pd_vector.qza   --m-metadata-file $samplemetadata   --o-visualization $odir/alpha-group/faith-pd-group-significance.qzv\n";
print SHELL "qiime diversity alpha-group-significance   --i-alpha-diversity $odir/core-metrics-results/evenness_vector.qza   --m-metadata-file $samplemetadata   --o-visualization $odir/alpha-group/evenness-group-significance.qzv\n";
print SHELL "qiime diversity alpha-group-significance   --i-alpha-diversity $odir/core-metrics-results/observed_features_vector.qza   --m-metadata-file $samplemetadata   --o-visualization $odir/alpha-group/observed_features_vector.qzv\n";
print SHELL "qiime diversity alpha-group-significance   --i-alpha-diversity $odir/core-metrics-results/shannon_vector.qza --m-metadata-file $samplemetadata   --o-visualization $odir/alpha-group/shannon_vector.qzv\n";

##alpha correlation, correlating numerical meta data with alpha diversity 
#
##Step 7 beta diversity analysis, group significance analysis 

print SHELL "\n#Step 7 beta diversity analysis, group significance analysis\n";

##Generate batch commands for beta diversity group significance analysis 

my $nam = `head -1 $samplemetadata`;chomp($nam);
my @name = split(/\t/, $nam); shift @name;  
unless(-d "$odir/beta-group"){ `mkdir $odir/beta-group`;} 

for my $k (@name){
	my $o = $k;
	$o =~  s/[ \/\(\)\=\;\:]//g;  ## remove special characters not appropriate for file name
	print SHELL "qiime diversity beta-group-significance   --i-distance-matrix $odir/core-metrics-results/bray_curtis_distance_matrix.qza --m-metadata-file $samplemetadata  --m-metadata-column \"$k\"   --o-visualization $odir/beta-group/$o.bray_curtis.qzv  --p-pairwise\n";	
	print SHELL "qiime diversity beta-group-significance   --i-distance-matrix $odir/core-metrics-results/jaccard_distance_matrix.qza --m-metadata-file $samplemetadata  --m-metadata-column \"$k\"   --o-visualization $odir/beta-group/$o.jaccard.qzv  --p-pairwise\n";	
	print SHELL "qiime diversity beta-group-significance   --i-distance-matrix $odir/core-metrics-results/unweighted_unifrac_distance_matrix.qza  --m-metadata-file $samplemetadata  --m-metadata-column \"$k\"   --o-visualization $odir/beta-group/$o.unweighted_unifrac.qzv  --p-pairwise\n";	
	print SHELL "qiime diversity beta-group-significance   --i-distance-matrix $odir/core-metrics-results/weighted_unifrac_distance_matrix.qza --m-metadata-file $samplemetadata  --m-metadata-column \"$k\"   --o-visualization $odir/beta-group/$o.weighted_unifrac.qzv  --p-pairwise\n";	

}

##unzip qza and qzv format file into detail text information of the results
##-------------------------------------------------------------------------------Seperating Analysis-------------------------------------------------------------------------------------------------------------------##
##Seperate analysis: in the following part, it will automatically select the chosen part of data from the above generated results according to users'need to generate relevant  analysis 
##Column chosen, groups chosen,
#For example,  in the sample.meta.data.tsv, just choose one column of detail types and from this column, only choose two subgroups for analysis: Detail_TYPE:F_10_BXSB_Wild:F_20_BXSB_Wild

#Step 8, biomarker identification with lefse 
#print SHELL "\n#Step 8, biomarker identification with lefse or with genisis in QIIME2\b"; 
##input the otu table and meta data file identify all the signinficantly different signature with ANCOM algorithms 

##quit qiime2 conda environment  
#print SHELL "source deactivate\n";
#unless(-d "$odir/lefse"){`mkdir $odir/lefse`} ;
#print SHELL "source activate ampliconanalysis";
##As the default python version  is python2.7, for the lefse used python the correct numpy packages is under this directory 
#print SHELL "export PYTHONPATH=export PYTHONPATH=/home/z3524677/.local/lib/python2.7/site-packages"

##generate from feature to phylum automatically for all catergory


##perform LEfSe analysis to identify those differential abundant taxa/genes  also with ANCOM 
#citation:
##
#print SHELL "perl lefse_qiime.pl -i $odir/lefse/rarefiled.tab.tsv -m $samplemetadata -c $cat -o $catdirout -t pair-wise\n";
#`nohup bash $commandsh &`;
##unrap all the qzv and some of the qza 

#MaAsLin analysis to associate sOTUs with alpha diversity and sOTU 
##prepare the pcl file
#Set the default p-value for MaAsLin 0.05 and q-value as 0.2
#unless(-d "$odir/maaslin"){`mkdir $odir/masslin`} ;

##Step 9 
##Using R to draw vector graphics for alph diversity boxplot, taxa stackplot (Phylum, class, genus, species level), PcoA, RDA/CCA plot, variable partition, correlation network analysis 


__END__
