#!/usr/bin/perl
$vvv=0;
$raptorx=$ENV{"HOME"}."/GitBucket/casp13_pipeline/TGT_Package";

#main directory
$CNFHOME=$raptorx."/util/SS8_Predict";
$PSIBLASTEXE=$raptorx."/util/BLAST/bin/blastpgp";
$NR=$raptorx."/database/NR_new/nr90";
$outputdir="$CNFHOME";

require "$CNFHOME/bin/BioDSSP.pm";
use File::Temp qw/ tempfile tempdir /;
use List::Util qw[min max];

#check arguments
(@ARGV >= 1) or die <<"VEND";
Usage: run_raptorx-ss8.pl [sequence file]
sequence file should be a one-line file containing a peptide sequence.
VEND

    @p=@ARGV[1..@ARGV-1];
%par=@p;
if(defined $par{-outdir}){
    $outputdir=$par{-outdir};
}
if(defined $par{-pssm}){
    $FNTMPPSSM=$par{-pssm};
}
#Deal with the filename of sequence.
@paths=split /\//,$ARGV[0];
$pathsnum=scalar(@paths)-1;
my @tmp=split /\./,$paths[$pathsnum]; #./
my $numtmp=scalar(@tmp)-2;
$numtmp = 0 if($numtmp<0);
$seqPrefix=join(".",@tmp[0..$numtmp]);

#Read sequence in FASTA format and check sequence
open FSEQ,"<$ARGV[0]" or die $!;
@seqlines=<FSEQ>;
$seq="";
for($i=0;$i<scalar(@seqlines);$i++)
{
    if($seqlines[$i]=~/^>/){}
    else{
	chomp($seqlines[$i]);
	$seqlines[$i]=~s/\s+//g;
	$seq=$seq.$seqlines[$i];
    }
}
$seq=uc $seq;
$seq=~m/[ARNDCQEGHILKMFPSTWYVX]*/g;
#$seq=~/[BDAC]*/g;
#print $seq," ","  $& ","\n";
die "The sequence contains unknown amino acids other than ARNDCQEGHILKMFPSTWYVX \n" if ($seq ne $&);
close FSEQ;

#make tmp sequence file and test working directory.
$jobtmpdir=tempdir(DIR => "./" );
chomp($jobtmpdir);
$jobid=int(rand(99999999));
$jobid="$seqPrefix.$jobid";
$FNTMPSEQ="$jobtmpdir/cnfsseight.$jobid.seq";
$FNTINPUT="$jobtmpdir/cnfsseight.$jobid.input";
open FTMPSEQ,">$FNTMPSEQ";
print FTMPSEQ "$seq\n";
close FTMPSEQ;
if(!defined $par{-pssm}){
    $FNTMPPSSM="$jobtmpdir/cnfsseight.$jobid.seq.pssm";
    #Run psiblast to generate PSSM
    #print STDERR "$PSIBLASTEXE -d $NR -j 1 -h 0.001 -i $FNTMPSEQ -Q $FNTMPPSSM\n";
    print STDERR "Running PSIBLAST...\n" if $VVV==1;
    $tmp=`$PSIBLASTEXE -d $NR -j 5 -h 0.001 -i $FNTMPSEQ -Q $FNTMPPSSM `;
}
#print STDERR $tmp;
#exit;
#Generate test input file
print STDERR "Generating features...\n" if $VVV==1;
$tmp=`$CNFHOME/bin/GenFeat_nodssp.pl $FNTMPSEQ $FNTMPPSSM $CNFHOME/data/p.p2.c.unit $FNTINPUT $CNFHOME`;
#print STDERR $tmp;
#Run CNF_predict
$FNTMPRST="$jobtmpdir/cnfsseight.$jobid.result";
$FNTMPNULL="$jobtmpdir/cnfsseight.$jobid.null";
open FTMPNULL,">$FNTMPNULL";
print FTMPNULL "0\n";
close FTMPNULL;
#exit;
print STDERR "Predicting 8-class secondary structure using CNF model...\n" if $VVV==1;
#$cmd="$CNFHOME/bin/bcnf_mpitp CONF $CNFHOME/data/CNF.ax.norm.conf ACT null PREDICT1 1 FEATMASK 0 ALLRST all.rst.csv BIEN 1 TRAIN $FNTMPNULL TEST $FNTINPUT RESULT $FNTMPRST RESUME $CNFHOME/data/model.p0-900.1561189 &> $jobtmpdir/cnfsseight.$jobid.stderr";
$cmd="$CNFHOME/bin/bcnf_mpitp CONF $CNFHOME/data/CNF.ax.norm.conf ACT null PREDICT1 1 FEATMASK 0 ALLRST all.rst.csv BIEN 1 TRAIN $FNTMPNULL TEST $FNTINPUT RESULT $FNTMPRST RESUME $CNFHOME/data/model.p0-900.1561189 ";
$cmdout=`$cmd 2>&1`;

#`./bcnf_mpitp CONF CNF.ax.norm.conf ACT null PREDICT1 1 FEATMASK 0 ALLRST all.rst.csv BIEN 1 TRAIN train.null TEST output RESULT result.cnfss8  `;

#Parse the output of cnf ss eight prediction
open fhCNFOutput,"<$FNTMPRST" or die $!;
$rstline=<fhCNFOutput>;
$rstline=~s/^\s+//;
chomp($rstline);
@p=split/\s+/, $rstline;
my @alllab=split //, $p[6];
my @prob=@p[7..scalar(@p)-1];
close fhCNFOutput;

@alllab=convert_label_to_letter(@alllab);
#Output the data in letter form, each line is 
# [number] [amino acid] [ss in eight letters]  [eight probability]
print STDERR "Formating results...\n" if $VVV==1;


open fhRESULT,">$outputdir/$seqPrefix.ss8";
@allseq=split //, $seq;
print fhRESULT "#RaptorX-SS8: eight-class secondary structure prediction results\n";
print fhRESULT "#probabilities are in the order of H G I E B T S L(loops), the 8 secondary structure types used in DSSP\n\n";
for($i=0;$i<@allseq;$i++)
{
    
    print fhRESULT sprintf("%4d %s %s   ",$i+1,$allseq[$i], $alllab[$i]);
#    print fhRESULT "$allseq[$i]\t $alllab[$i]\t";
    for($k=0;$k<8;$k++)
    {
	print fhRESULT sprintf("%.3f ",$prob[$i*8+$k]);
#	print fhRESULT $prob[$i*8+$k],"\t";
    } 
    print fhRESULT "\n";
} 
close fhRESULT;

open fhHO,">$outputdir/$seqPrefix.ss8.horiz";
print fhHO "#RaptorX-SS8: eight-class secondary structure prediction results\n\n";
for(my $i=0;$i<@allseq;$i=$i+60)
{
    my $end;
    if($i+60>@allseq)
    {
	$end=@allseq-1;
    }
    else{
	$end=$i+60-1;
    }
    print fhHO "Conf: ",join("",@hoConf[$i..$end]),"\n";
    print fhHO "Pred: ",join("",@alllab[$i..$end]),"\n";
    print fhHO "  AA: ",join("",@allseq[$i..$end]),"\n";
    print fhHO "      ";
    for(my $j=$i+10;$j<=$end+1;$j=$j+10){
	print fhHO sprintf("%10d",$j);
    }
    print fhHO "\n"x3;
}
close fhHO;
print STDERR "DONE ss8.\n" if $VVV==1;
`rm -rf output*`;
if(!defined $ARGV[1]){
`rm -f cnfsseight.$jobid.*`;
}

print STDERR "DONE.\n" if $VVV==1;
`rm -rf output?*`;
if(!defined $par{-keep}){
`rm -f cnfsseight.$jobid.*`;
`rm -rf $jobtmpdir`;
}

