#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 'proteinList inputDir [outDir]'
	echo '	This script simplifies a list of tpl.pkl files'
	exit 1
fi

listFile=$1
if [ ! -f $listFile ]; then
	echo 'ERROR: invalid protein list file $listFile'
	exit 1
fi

inDir=$2
if [ ! -d $inDir ]; then
	echo 'ERROR: invalid folder for input tpl.pkl files'
	exit 1
fi

outDir=`pwd`
if [ $# -ge 3 ]; then
	outDir=$3
fi

program=$ModelingHome/Alignment/SimplifyTPLPKL.py
if [ ! -f $program ]; then
	echo "ERROR: invalid program file $program"
	exit 1
fi

mkdir -p $outDir

for i in `cat $listFile `
do
	python $program $inDir/${i}.tpl.pkl $outDir
done
