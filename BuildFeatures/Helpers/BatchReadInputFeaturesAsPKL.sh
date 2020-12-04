#!/bin/bash

if [[ -z "${DL4DistancePredHome}" ]]; then
        echo "ERROR: Please set environmental variable DL4DistancePredHome to the installation folder of DL4DistancePrediction4 "
        exit 1
fi

if [ $# -lt 2 ]; then
	echo $0 proteinListFile inputFeatureMetaFolder [savefolder]
	echo "	This script reads input features for a list of proteins and save the features as individual PKL files"
	echo "	inputFeatureMetaFolder: the umbrella folder containing a list of subfolders, each with name feat_proteinName_contact"
	exit 1
fi

proteinListFile=$1
if [ ! -f $proteinListFile ]; then
	echo "ERROR: invalid file for a list of proteins $proteinListFile"
	exit 1
fi

inputFolder=$2
if [ ! -d $inputFolder ]; then
	echo "ERROR: invalid input folder for protein features $inputFolder"
	exit 1
fi

savefolder=`pwd`
if [ $# -ge 3 ]; then
	savefolder=$3
	if [ ! -d $savefolder ]; then
		mkdir -p $savefolder
	fi
fi

program=$DL4DistancePredHome/ReadSingleInputFeature.py

numAllowedJobs=15
keywords=`basename $program`
myself=`basename $0 `

for target in `cat $proteinListFile `
do
        while true
        do
                ## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | grep -v $myself | wc -l`
                if [ $numRunningJobs -lt $numAllowedJobs  ]; then
			python $program $target $inputFolder/feat_${target}_contact $savefolder & 
                        break
                else
                        a=`expr $RANDOM % 3 `
                        sleep $a
                fi

        done

done

wait
