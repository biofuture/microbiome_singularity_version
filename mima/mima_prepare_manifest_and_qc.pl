#!/usr/bin/perl -w
use strict;

die "perl $0 <Absolute path to store fastq files> <output.manifest> <total_1.fq> <total_2.fq> <outputdir> <# of threads> <Singularity_sif>
	<inputdir_abs_path>  The absolute directory to store the pair-end fastq files
	<output.manifest>    The manifest file that qiime needed when running dada2
	<Merged_1.fq>        The overall fastq file pair 1
	<Merged_2.fq>        The overall fastq file pair 2
	<Output_dir>         Output directory store all the QC files and dada2 output
	<#ofthreads>	     Deafault 1
	<Singularity sif img> the singualrity image to use 	

	##This script is to generate manifest file for qiime2 
	##developed by Dr. Xiaotao Jiang 
	#Email: xiaotao.jiang\@unsw.edu.au
	#copyright: 2018-2021
\n" unless (@ARGV == 7);
my $singularitysif = $ARGV[6]; 

manifest($ARGV[0], $ARGV[1]);
unless(-d "$ARGV[4]"){ `mkdir $ARGV[4]`;} ##create the output dir
`cp $ARGV[1] $ARGV[4]/`;
overall_quality($ARGV[2], $ARGV[3], $ARGV[1], $ARGV[4]);

sub overall_quality{
	#using fastp to generate the quality report for the whole run data
	#this function will submit the job to katana with qsub
	
	my ($fq1, $fq2,$manifest, $odir) = @_;
	
	my $temp = "$odir/qc.pbs";
	my $ofq1 = "$odir/ALL_R1.fq";
	my $ofq2 = "$odir/ALL_R2.fq";
	my $outreport = "$odir/fastp.outreport.html";
	die "$!\n" unless open(OT, ">$temp");
	print OT "#!/bin/bash\n#PBS -l nodes=1:ppn=$ARGV[5]\n#PBS -l mem=80gb\n#PBS -l walltime=100:00:00\nexport LC_ALL=en_AU.utf8\nexport LANG=en_AU.utf8\n";
	#print OT "source activate /srv/scratch/mrcbio/bin/miniconda3/envs/qiime2-2019.1\n";
	#print OT "export PYTHONPATH=/srv/scratch/z3524677/bin/miniconda3/envs/qiime2-2018.8/lib/python3.5/site-packages/\n";		
	print OT "cd $odir\n";
	print OT "singularity exec $singularitysif  bash -c \' fastp -i $fq1 -I $fq2 -o $ofq1 -O $ofq2 -h $outreport\'\n";
	print OT "singularity exec $singularitysif  bash -c \'. activate qiime2-2020.8 && qiime tools import --type \'SampleData[PairedEndSequencesWithQuality]\' --input-path $manifest --output-path demux.qza --input-format PairedEndFastqManifestPhred33\'\n";
	print OT "singularity exec $singularitysif  bash -c \'. activate qiime2-2020.8 && qiime dada2 denoise-paired  --i-demultiplexed-seqs demux.qza  --p-trunc-len-f 295 --p-trunc-len-r 220 --p-trim-left-f 17  --p-trim-left-r 21  --p-n-threads $ARGV[5] --o-representative-sequences rep-seqs.qza --o-table table.qza  --o-denoising-stats stats-dada2.qza\'\n";
	close OT;
}

sub manifest{
	#input as the absolute directory
	#all pair-end files should be placed into the directory with clear R1 as the forward file and R2 as the reverse file
	
	##output of the function is the manifest file in the ofile 

	my ($path, $ofile) = @_; 
	open (MAF, ">$ofile") or die "$!";
	opendir (DIR, $path) or die $!;
	print MAF "sample-id,absolute-filepath,direction\n";
	while (my $file = readdir(DIR)) {
		#check whether it is the first file
		my @m = split("/", $file);
		$m[-1] =~ /^(\S+)\_R1\S+$/;
		my $id = $1;
        	if($file =~ /_R1/){
			#check whether there is the second file
			my $f2 = $file;
			$f2 =~ s/_R1/_R2/;
			die "$path/$f2\n$!" unless(-e "$path/$f2");
			print MAF "$id,$path/$file,forward\n$id,$path/$f2,reverse\n";
		}
	}
	closedir(DIR); 
	close MAF;
}
