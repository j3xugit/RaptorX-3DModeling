#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "proteinListFile inputfolder [savefolder]"
	echo "	This script derives distance information for protein threading from predicted dist matrix"
	echo "	inputfolder: the folder for predicted dist matrix in PKL format"
	echo "	savefolder: the folder for derived dist information, including files ending with .pkl and .epad_prob, default current work directory"
	exit 1
fi

proteinListFile=$1

if [ ! -f $proteinListFile ]; then
	echo "ERROR: invalid protein list file $proteinListFile"
	exit 1
fi

inputFolder=$2
if [ ! -d $inputFolder ]; then
	echo "ERROR: invalid input folder for predicted dist informaiton: $inputFolder"
	exit 1
fi

savefolder=`pwd`
if [ $# -ge 3 ]; then
	savefolder=$3
	if [ ! -d $savefolder ]; then
		mkdir -p $savefolder
	fi
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/DeriveDistInfo4Threading.sh

numAllowedJobs=15
keywords=`basename $program`
myself=`basename $0 `

for target in `cat $proteinListFile `
do
        while true
        do
                numRunningJobs=`ps -x | grep ${keywords} | grep -v $myself | wc -l`
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
			$program $inputFolder/${target}.predictedDistMatrix.pkl $savefolder & 
			sleep 1
                        break
                else
                        a=`expr $RANDOM % 3 `
                        sleep $a
                fi
        done
done

wait
