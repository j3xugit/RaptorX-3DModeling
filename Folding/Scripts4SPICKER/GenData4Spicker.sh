#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 targetList SeqDir MetaFolder [cutoff-method]
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/GenInputInfo4SPICKER.py

targets=`cat $1`
SeqDir=$2
MetaFolder=$3

cutoff=-1
if [ $# -ge 4 ]; then
	cutoff=$4
fi

for target in $targets
do
	seqFile=$SeqDir/${target}.fasta

	for modelFolder in $MetaFolder/${target}-*
	do
		python $program -s ${target}-Spicker -c $cutoff $seqFile $modelFolder &
	done
done
