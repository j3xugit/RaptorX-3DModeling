#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 SeqFile ModelFolder [cutoff-method]
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=$DistanceFoldingHome/Scripts4SPICKER
program=$cmdDir/GenInputInfo4SPICKER.py

seqFile=$1
target=`basename $seqFile .fasta`
target=`echo $target | cut -f1 -d'.' `

modelFolder=$2

cutoff=-1
if [ $# -ge 3 ]; then
	cutoff=$3
fi

python $program -s ${target}-Spicker -c $cutoff $seqFile $modelFolder
