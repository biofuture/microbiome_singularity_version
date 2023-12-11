#!/usr/bin/perl -w

die "perl $0 <rep_seqs.qza> <table.qza> <outputdir> <host_genome.fa_bowtie2_index> <threads>
	<rep_seqs.qza> qiime format represntative sequences output from dada2
	<table.qza>    qimme format features table 
	<outputdir>    output directory to store the decontaminated updated rep_seqs.qza and table.qza 
	<host_genome.fa_bowtie2_index> bowtie2 index of the host genome 
	<threads>      number of threads used in bowtie2 mapping 
	
	This script use bowtie2 to maping the rep_seqs.qza aganist host genome to identify the contaminated sequences and update the sequences and table with qiime2. The updated files can be used for downstream analysis 

\n" unless (@ARGV == 5);

unless(-d "$ARGV[2]"){
	`mkdir $ARGV[2]`;
}
$ARGV[4] ||= 4;
#my $qiime = $ARGV[4];
my $odir = $ARGV[2];
my $mousebowtieindex = "$ARGV[3]"; # mouse.udb mouse database build for usearch index 
my $pbs = "$ARGV[2]/decontaminationhost.sh";
die "$!\n" unless open(PBS, ">$pbs");
#print PBS "source activate $qiime\n";
print PBS "qiime tools export --input-path $ARGV[0] --output-path $odir/rep-seqsold\n";
#die "error qiime tools export --input-path\n" unless (-e "$odir/rep-seqs/dna-sequences.fasta");
#alignment with host database 
my $dnaseqs = "$odir/rep-seqsold/dna-sequences.fasta";
my $samout = "$odir/rep-seqsold/align-hostgenome.sam.txt"; 
my $matechedfasta = "$odir/rep-seqsold/matched_dnaseq.fasta";
#$mousebowtieindex ||= "/srv/scratch/z3524677/db/mousegenome/mouse.genome.bowtie2";
#print PBS "usearch -usearch_local  $odir/rep-seqs/dna-sequences.fasta  -db $hostudb -id 0.6 -evalue 1e-5 -blast6out $bblastm6  -strand both -maxaccepts 1  -matched matched_dnaseq.fasta -query_cov 0.6\n";
print PBS "bowtie2 -x $mousebowtieindex -U $dnaseqs -f --sensitive-local  -p $ARGV[4] --un unaligned_bowtie2.fasta --al $matechedfasta > $samout\n";
my $metadata = "$odir/metadata.table.txt"; 
print PBS "grep \">\" $matechedfasta  | perl -ne 'chomp; s/^>//; print \"\$_\\t1\\n\";' > $metadata ; echo  -e \"feature-id\\tSIDs\"  | cat - $metadata > temp.out | mv temp.out $metadata\n";
my $table = "$ARGV[1]";
my $otable = "$odir/table.qza";
print PBS "qiime feature-table filter-features --i-table $table --m-metadata-file $metadata  --p-where \"SIDS==1\" --p-exclude-ids --o-filtered-table $otable\n"; 
my $filteredseqs = "$odir/filtered_dnaseqs.fasta";
my $cleandnaseq = "$odir/rep-seqs.qza";
print PBS "perl -e 'open(I, \"$metadata\"); <I>;while(<I>){ chomp; \@m = split /\\t/; \$hash{\$m[0]} =1;  }   open(II, \"$dnaseqs\"); while(\$id = <II>){ \$seq = <II>; chomp(\$id); \$id =~ s/^>//; if(exists \$hash{\$id}){}else{print \">\$id\\n\$seq\";} }  ' > $filteredseqs\n";
print PBS "qiime tools import --type 'FeatureData[Sequence]' --input-path $filteredseqs --output-path $cleandnaseq\n";
print PBS "qiime tools export --input-path $otable --output-path $odir/otable\n";
print PBS "biom convert -i $odir/otable/feature-table.biom -o $odir/otable/feature.txt --to-tsv\n";
my $tablesummarize = "$odir/tablesummarize.qzv";
print PBS "qiime feature-table summarize --i-table $otable  --o-visualization  $tablesummarize\n";
print PBS "qiime tools export --input-path  $tablesummarize --output-path $odir/otablesummarize\n";
close PBS;
