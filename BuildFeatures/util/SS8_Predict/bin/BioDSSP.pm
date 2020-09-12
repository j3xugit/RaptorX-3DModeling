#!/usr/bin/perl -w


sub getdssp {
#return @seq, @label, @marks
    open FDSSP,"<$_[0]" or die $!;
    while(<FDSSP>)
    {
	last if(/RESIDUE AA STRUCTURE/);
    }
    my @allaa=qw();
    my @alllab=qw();
    my @allbe=qw();
    my $offset=-1;
    my @mark=qw();
    
#    push @mark,0;
    my $markat=0;
    $lastaa="!";
    while(<FDSSP>)
    {
	chomp;
#	next if(/GAP/);
	last if(/END/);
	my $l=substr($_,16,1);
	if($l eq ' '){$l='L';}
	my $aa=substr($_,13,1);
	my $pos=substr($_,5,5);
	if($aa ne '!')
	{
	    if($aa eq lc($aa))
	    {
		$aa='C'; #s-s bond
	    }
#	    push @label,$labelmap{$l};
	    push @alllab, $l;
	    push @allaa, $aa;
	    $markat++;
	}
	if($lastaa eq '!' && $aa ne '!' )
	{
	    push @mark,$markat-1;
	}
	if($lastaa ne '!' && $aa eq '!' )
	{
	    push @mark,$markat-1;
	}
	$lastaa=$aa;
	last if(eof(FDSSP));
    }
    if($lastaa ne '!' )
    {
	push @mark,$markat;
    }
    close FDSSP;
    return (\@allaa, \@alllab ,  \@mark);
}

sub getdssp_short {
#return @seq, @label, @marks from simplified dssp file like ~/db/dssp
    open FDSSP,"<$_[0]" or die $!;
    while(<FDSSP>)
    {
	last if(/RS NUM SS/);
    }
    my @allaa=qw();
    my @alllab=qw();
    my @allbe=qw();
    my $offset=-1;
    my @mark=qw();
    
#    push @mark,0;
    my $markat=0;
    $lastaa="!";
    while(<FDSSP>)
    {
	chomp;
	next if(/GAP/);
	last if(/END/);
	my $l=substr($_,14,1);
	if($l eq 'L'){$l='C';}
	my $aa=substr($_,6,1);
	my $pos=substr($_,7,5);
	if($aa ne '!')
	{
	    if($aa eq lc($aa))
	    {
		$aa='C'; #s-s bond
	    }
#	    push @label,$labelmap{$l};
	    push @alllab, $l;
	    push @allaa, $aa;
	    $markat++;
	}
	if($lastaa eq '!' && $aa ne '!' )
	{
	    push @mark,$markat-1;
	}
	if($lastaa ne '!' && $aa eq '!' )
	{
	    push @mark,$markat-1;
	}
	$lastaa=$aa;
	last if(eof(FDSSP));
    }
    if($lastaa ne '!' )
    {
	push @mark,$markat;
    }
    close FDSSP;
    return (\@allaa, \@alllab ,  \@mark);
}



sub get_psipred {
#return @seq, @label, @marks
    open FSS2,"<$_[0]" or die $!;
    <FSS2>;
    <FSS2>;
    my @allaa=qw();
    my @alllab=qw();
    my @allprob=qw();
    while(<FSS2>)
    {
	chomp;
	s/^\s+//;
	@p=split/\s+/;
	die 'Woops! $_' if(scalar(@p)!=6);
	push @allaa,$p[1];
	push @alllab,$p[2];
	push @allprob,@p[3..5];
    }
    close FSS2;
    return (\@allaa, \@alllab ,  \@allprob);
}

sub convert_label83 {
#letter will be converted to letters
#numbers will be converted to numbers
#all input should be upper case.
    my %map83=qw(H H G H B E E E C C S C T C I C L C);
    if(scalar(@_)==1)
    {
	my @rst;
	#deal with a string
	for($i=0;$i<length($_[0]);$i++)
	{
	    push @rst, $map83{substr($_[0],$i,1)};
	}
	return join("",@rst);
    }
    else
    {
	#deal with an array of chars.
	my @rst;
	#deal with a string
	for($i=0;$i<scalar(@_);$i++)
	{
	    push @rst, $map83{uc($_[$i])};
	}
	return @rst;
    }
}

sub convert_label_to_letter {
#letter will be converted to letters
#numbers will be converted to numbers
#all input should be upper case.
    my %map=qw(0 h 1 g 2 i 3 e 4 b 5 t 6 s 7 l);
    if(scalar(@_)==1)
    {
	my @rst;
	#deal with a string
	for($i=0;$i<length($_[0]);$i++)
	{
	    push @rst, uc $map{substr($_[0],$i,1)};
	}
	return join("",@rst);
    }
    else
    {
	#deal with an array of chars.
	my @rst;
	#deal with a string
	for($i=0;$i<scalar(@_);$i++)
	{
	    push @rst, uc $map{$_[$i]};
	}
	return @rst;
    }
}

sub convert_label_to_letter3 {
#letter will be converted to letters
#numbers will be converted to numbers
#all input should be upper case.
    my %map=qw(0 h 1 e 2 c);
    if(scalar(@_)==1)
    {
	my @rst;
	#deal with a string
	for($i=0;$i<length($_[0]);$i++)
	{
	    push @rst, uc $map{substr($_[0],$i,1)};
	}
	return join("",@rst);
    }
    else
    {
	#deal with an array of chars.
	my @rst;
	#deal with a string
	for($i=0;$i<scalar(@_);$i++)
	{
	    push @rst, uc $map{$_[$i]};
	}
	return @rst;
    }
}



sub HammingDistanceSeqs{
    $s1=$_[0];
    $s2=$_[1];
    @sp1=split(//,$s1);
    @sp2=split(//,$s2);
#    $maxlen=scallar(@sp1);
#    $maxlen=@sp2 if(@sp2>$maxlen);
    my $corr;
    for(my $ii=0;$ii<@sp1;$ii++)
    {
	$corr++ if($sp1[$ii] eq $sp2[$ii]);
    }
    if(@sp1 != 0)
    {
	return $corr/@sp1;
    }
    else
    {
	return -1;
    }
}

sub HammingDistanceArrays{
    $s1=$_[0];
    $s2=$_[1];
    @sp1=@{$s1};
    @sp2=@{$s2};
#    $maxlen=scallar(@sp1);
#    $maxlen=@sp2 if(@sp2>$maxlen);
    my $corr;
    for(my $ii=0;$ii<@sp1;$ii++)
    {
	$corr++ if($sp1[$ii] eq $sp2[$ii]);
    }
    if(@sp1 != 0)
    {
	return $corr/@sp1;
    }
    else
    {
	return -1;
    }
}

sub ParsePSSM{
#    print STDERR $_[0],"\naaa\n";
    open fhPSSM,"<$_[0]" or die $!;
    <fhPSSM>;<fhPSSM>;<fhPSSM>;
    my @allseq;
    my $allfeat="";
    my $len=0;
    my @allfeats=qw();
    while(<fhPSSM>)#for($k=0;$k<$len;$k++)
    {
	my $f1=$_;
	chomp($f1);
	$f1=~s/^\s+//;
	my @parts=split(/\s+/,$f1);
#    die "Parse PSSM file error!\n" if(@parts!=44);
	last if(scalar(@parts)!=44);
	push @allseq, $parts[1];
	$len++;
	my $neff=$parts[42]." "; # Neff value;
#normalize psp
        my $neffpos=0;
        for($i=22;$i<42;$i++)
        {
            if($parts[$i]!=0){
                $neffpos = $neffpos - $parts[$i]/100 * log($parts[$i]/100);
            }
        }
        $neffpos=exp($neffpos);
        $globalNeff=$globalNeff+$neffpos;

	for($ipssm=2;$ipssm<=21;$ipssm++){
	    $parts[$ipssm]=sprintf("%.4f",($parts[$ipssm]+7)/14.0);
	}
	my $f2=join(" ",@parts[2..21]);
	$f2=sprintf("%.3f",$neffpos/10)." ".$f2;
	push @allfeats,$f2;
    }
    close fhPSSM;
#    print STDERR join("",@allseq);
    return (\@allseq,\@allfeats);
}


#### the end  #####
1;
