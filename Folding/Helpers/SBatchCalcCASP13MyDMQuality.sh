#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 metaFolder [NativeDir [ DM2MyDM_file [flag] ] ]
	echo "metaFolder: a folder containing a set of subfolders. Each subfolder has name target-FoldResults that contains the decoys of one target"
	echo "NativeDir: the folder for experimental stucture files"
	echo "DM2MyDM_file: a file that mapps official domain names to my domain names. default CASP13FM2MyDM.txt"
	exit 1
fi

metaFolder=$1

nativeDir=CASP13DM-Native
if [ $# -ge 2 ]; then
	nativeDir=$2
fi

mappingFile=CASP13FM2MyDM.txt
if [ $# -ge 3 ]; then
	mappingFile=$3
fi

flag=''
if [ $# -ge 4 ]; then
	flag=$4
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/CalcMyDMQuality.sh

while read p; do
  	echo "$p"
	DMname=`echo $p | cut -f1 -d':'`
	MyDMname=`echo $p | cut -f3 -d':'`
	for decoyFolder in $metaFolder/${MyDMname}-*Results/
	do
		if [ -z "$flag" ]; then
			sbatch -p contrib-cpu -J Quality4$DMname -o $DMname.quality.out -e $DMname.quality.out $program $decoyFolder $nativeDir/$DMname.pdb
		else
			sbatch -p contrib-cpu -J Quality4$DMname -o $DMname.quality.out -e $DMname.quality.out $program $decoyFolder $nativeDir/$DMname.pdb $flag
		fi
	done

done <${mappingFile}
