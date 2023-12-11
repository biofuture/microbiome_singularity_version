package  Matrix;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(normalized add_groupinfo);

#This package is mainly to process matrix
#such as normalize matrix by column
sub normalized{
    #F P  is matrix file

    my @tm = @_;
    my $numb = $tm[1];
$numb ||= 1;
    die "This File does not exits $tm[0] $!\n" unless open(I, "$tm[0]") ;
    my $head = <I>; chomp($head);
    my @nam = split("\t",$head);
    my $num = $#nam;
    my %sum;
    while(<I>)
    {
        my @tem = split("\t",$_);
        for(my $i =1; $i <= $num; $i++)
        {
            $sum{$i}+=$tem[$i];
        }
    }
    close I;

    die "Normalized $!\n" unless open(I,"$tm[0]");
    my $out = "$tm[0].normalized_$numb.txt";
    die "OutPut file $out\n" unless open(T,">$out");
    print T "$head\n";
    <I>;
    while(<I>)
    {
        my @tem = split("\t",$_);
        print T  "$tem[0]";
        for(my $i =1; $i <= $num; $i++)
        {
            my $on = int(($tem[$i] / $sum{$i}) * $numb);
            print T "\t$on";
        }
        print T "\n";
    }
    close I;
    close T;
}#normanized


sub add_groupinfo{
    #1. the matrix with first row with sample ID
    #2. group files

    my @ar = @_;
    my $fline = `head -1 $ar[0]`; chomp($fline);
    die "$! \n" unless open(I1, "$ar[1]");
    my $ot = "$ar[0].addgroup";
    die "$! \n" unless open(T1, ">$ot");
    die "$! \n" unless open(I0, "$ar[0]");

    my %ghash;
    <I1>;
    while(<I1>)
    {
        chomp;
        my @ts = split(/\t/, $_);
        $ghash{$ts[0]} = $ts[1];
    }
    close I1;
    
    my @tm = split(/\t/,$fline);
    my @otm=();
    if($tm[0] eq ""){

        $tm[0] = "SID";    
    }
    push @otm, "Class";
    for(my $i= 1; $i <= $#tm; $i++)
    {
        if(exists $ghash{$tm[$i]})
        {
            push @otm, $ghash{$tm[$i]};
        }
        else
        {
            die "wrong\n";    
        }
    }
    
    my $o = join("\t", @otm);
    my $osi = join("\t", @tm);
    print T1  "$o\n$osi\n";
     <I0>;
     while(<I0>)
     {
         #print T1 $_;
         s/\;/\|/g;
         print T1 $_;
    }
    close I0; close T1;

}#addgroup

1;
__END__
