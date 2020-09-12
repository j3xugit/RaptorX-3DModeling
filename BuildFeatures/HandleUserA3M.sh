#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 " fastaFile outDir "
	echo "	This script handles user-provided MSA in a3m format"
	echo "	fastaFile: a user-provided a3m file"
	echo "	outDir: the output directory. the resultant files will be saved to $outDir/$target_OUT/"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

#-------- helper ----#
meff_cdhit=$cmdDir/Meff_CDHIT/meff_cdhit

#-------- input arguments ------#
userMSA=$1
if [ ! -f $userMSA ]; then
	echo "ERROR: invalid user-provided MSA file: $userMSA"
	exit 1
fi

fulnam=`basename $userMSA`
target=${fulnam%.*}

if [ ! -f $target.a3m ]; then
	cp $userMSA $target.a3m
fi

out_folder=$2
outdir=$out_folder/${target}_OUT
mkdir -p $outdir

#-------- generate a valid fasta file----#
head -2 $userMSA > $outdir/$target.seq

input_fasta=$target.fasta
head -2 $userMSA > $target.fasta

#-> 1. for threading
mkdir -p $outdir/${target}_thread
if [ ! -f "$outdir/${target}_thread/$target.a3m" ] ||
   [ ! -f "$outdir/${target}_thread/$target.tgt" ]
then

	echo "Generating threading tgt file..."
	tmp_root=$target"_tmp"/
	mkdir -p $tmp_root
	$cmdDir/util/A3M_To_TGT -i $outdir/$target.seq -I $target.a3m -o $target.tgt -t $tmp_root
	cp $target.a3m $outdir/${target}_thread
	mv $target.tgt $outdir/${target}_thread
	rm -rf $tmp_root
fi


#-> 2. for prediction of inter-atom or inter-residue relationship
DESDIR=$outdir/${target}_contact
mkdir -p $DESDIR/${target}_user

if [ ! -f "$DESDIR/${target}_user/$target.a3m" ]; then
	echo "Generating user a3m file for contact prediction..."
	mv $target.a3m $DESDIR/${target}_user/
	cp $outdir/$target.seq $DESDIR/${target}_user/
	mv $target.fasta $DESDIR/${target}_user/${target}.fasta_raw

	#-> calculate meff
	$meff_cdhit -i $DESDIR/${target}_user/$target.a3m > $DESDIR/${target}_user/$target.meff
fi
