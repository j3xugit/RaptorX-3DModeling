#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 proteinListFile tgtDir tplpklDir [ResDir]
	exit 1
fi

proteinList=$1
TGTDir=$2
TPLPKLDir=$3

ResDir=`pwd`
if [ $# -ge 4 ]; then
	ResDir=$4
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/TPLGenerator/TGTpkl_To_TPL2

if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

for i in `cat $proteinList`
do
	$program -i $TPLPKLDir/${i}.tpl.pkl -I $TGTDir/${i}.tgt -o $ResDir/${i}.tpl -H $cmdDir/TPLGenerator/
done
