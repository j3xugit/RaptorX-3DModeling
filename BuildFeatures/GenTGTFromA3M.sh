#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 a3mFile [ResDir]
	exit 1
fi

a3mfile=$1

ResDir=`pwd`
if [ $# -ge 2 ]; then
	ResDir=$2
	if [ ! -d $ResDir ]; then
		mkdir -p $ResDir
	fi
fi

if [[ -z "${DistFeatureHome}" ]]; then
        echo "Please set environmental variable DistFeatureHome to the installation folder of this program, e.g., $HOME/3DModeling/BuildFeatures "
        exit 1
fi


## get target name
fulnam=`basename $a3mfile`
relnam=${fulnam%.*}

# ---- generate A2M file from A3M file ---------#
a2mfile=$ResDir/$relnam.a2m
$DistFeatureHome/util/A3M_To_PSI $a3mfile $a2mfile.tmp
grep -v "ss_pred\|ss_conf" $a2mfile.tmp | awk '{print substr($0,34,length($0)-32) }' > $a2mfile
rm -f $a2mfile.tmp

# ---- generate FASTA file from A2M file -------#
fastafile=$ResDir/$relnam.seq
echo ">$relnam" > $fastafile
head -n1 $a2mfile >> $fastafile

# ---- generate TGT file from A3M file ---------#
tgtfile=$ResDir/$relnam.tgt
tmp=$ResDir/${relnam}-contact
mkdir -p $tmp
numLines=`wc -l $a3mfile | cut -f1 -d' '`
if [ $numLines -gt 80000 ]; then
	echo "WARNING: converting A3M to TGT by sampling 40000 seqs from $a3mfile"
	a3mfile2=$ResDir/$relnam.a3m.sampled
	python $DistFeatureHome/Helpers/SampleA3MByNumber.py $a3mfile 40000 $a3mfile2
        $DistFeatureHome/util/A3M_To_TGT -i $fastafile -I $a3mfile2 -o $tgtfile -t $tmp 1> $relnam.ws1 2> $relnam.ws2
else
        $DistFeatureHome/util/A3M_To_TGT -i $fastafile -I $a3mfile -o $tgtfile -t $tmp 1> $relnam.ws1 2> $relnam.ws2
fi

OUT=$?
if [ $OUT -ne 0 ]; then
	echo "ERROR: failed to run util/A3M_To_TGT -i $fastafile -I $a3mfle -o $tgtfile -t $tmp"
        exit 1
fi
rm -rf $tmp
rm -f $relnam.ws1
rm -f $relnam.ws2
