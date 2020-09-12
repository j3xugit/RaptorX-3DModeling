#!/usr/bin/perl


(@ARGV >=5 ) or print <<"VEND";
Usage: GenFeat_nodssp.pl [sequence file] [pssm file] [feature matrix] [output file] [home directory]

VEND

require "$ARGV[4]/bin/BioDSSP.pm";

my $resorder="ARNDCQEGHILKMFPSTWYV";
my %resfeat;
my $seq;
open fhSEQ,"<$ARGV[0]";
$seq=<fhSEQ>;
chomp($seq);

open OFEAT,">$ARGV[3]";
open MATFILE,"<$ARGV[2]";

# Read feature matrix
my $i=0;
while(<MATFILE>)
{
    chomp;
    $_=~s/^\s+//;
    @parts=split/\s+/;
    $resfeat{substr($resorder,$i,1)}=$_;
    $i++;
    $resfeat{"X"}="0 "x scalar(@parts); 

}
close MATFILE;
## End of reading


$nsample=0;
%aaorder=qw(A 0 R 1 N 2 D 3 C 4 Q 5 E 6 G 7 H 8 I 9 L 10  K 11 M 12  F 13  P 14 S 15  T 16  W 17  Y 18  V 19);

#Parse PSSM
my @tmprst=ParsePSSM("$ARGV[1]");
my @pssmSeq=@{$tmprst[0]};
my @pssmFeat=@{$tmprst[1]};

my $pssmSequence=join "", @pssmSeq;
die "PSSM contains gap or sequence inconsistent with the input sequence!\n$seq\n$pssmSequence\n" if ($seq ne $pssmSequence);
#add identity matrix
my @allfeats;
for($j=0;$j<scalar(@pssmSeq);$j++)
{
    $allfeats[$j]=$pssmFeat[$j]." ".$resfeat{$pssmSeq[$j]};
}

#Write all features to output file
print OFEAT "1\n";
print OFEAT scalar(@allfeats),"\n";
print OFEAT join("\n",@allfeats);
print OFEAT "\n";
print OFEAT "0\n"x scalar(@pssmSeq);# all label are set to zeros
close OFEAT;

close fhSEQ;

