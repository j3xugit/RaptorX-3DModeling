#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 targetList SeqDir MetaFolder [cutoff]
	exit 1
fi

targets=`cat $1`
SeqDir=$2
MetaFolder=$3

cutoff=-1
if [ $# -ge 4 ]; then
	cutoff=$4
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for target in $targets
do
	for modelFolder in $MetaFolder/${target}-*
	do
		sbatch -p contrib-cpu -J Data4Spicker -o Spicker-${target}.out -e Spicker-${target}.out $cmdDir/GenOneTargetData4Spicker.sh $SeqDir/${target}.fasta $modelFolder $cutoff
	done
done
