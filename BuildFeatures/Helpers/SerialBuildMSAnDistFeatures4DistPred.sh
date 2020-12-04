#!/bin/bash

if [ -z "$DistFeatureHome" ]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the folder of BuildFeatures/ "
        exit 1
fi

if [ $# -lt 1 ]; then
	echo "$0 proteinListFile [DBID]"
	echo "	This script builds MSA and distance prediction input features for a list of proteins"
	echo "	The MSA and features will be built one-by-one protein"
	echo "	for each protein, one MSA will be generated"
	echo "   DBID: 2015, 2016, 2017 and 2018(default)"
	exit 1
fi

list=$1
GPU=-1

DBID=2018
if [ $# -eq 2 ]; then
	DBID=$2
fi

numCPUs=2

MSAprogram=$DistFeatureHome/HHblitsWrapper/BuildMSA4DistPred.sh
FEATprogram=$DistFeatureHome/GenDistFeaturesFromMSA.sh 
SEQDIR=SEQ/

DB="uniclust30_2017_10/uniclust30"
if [ $DBID == "2015" ]; then
	DB="uniprot20_2015_06/uniprot20"
elif [ $DBID == "2016" ]; then
	DB="uniprot20_2016_05/uniprot20"
elif [ $DBID == "2018" ]; then
	DB="uniclust30_2018_08/uniclust30"
fi

DB=$DistFeatureHome/HHblitsWrapper/databases/${DB}

MSADir=MSA_${DBID}_E001
mkdir -p $MSADir

DestDir=Features4Train_${DBID}_E001
mkdir -p $DestDir

for i in `cat $list `
do
	$MSAprogram -d $DB -c $numCPUs -o $MSADir/ $SEQDIR/$i.fasta
	if [ ! -f $MSADir/${i}.a3m ]; then
		echo "WARNING: invalid MSA file $MSADir/${i}.a3m"
		continue
	fi
	$FEATprogram -o $DestDir/feat_${i}_contact -g $GPU $MSADir/${i}.a3m
done
