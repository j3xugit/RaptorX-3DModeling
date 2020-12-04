#!/bin/bash

if [ -z "$DistFeatureHome" ]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the folder of BuildFeatures/ "
        exit 1
fi

if [ $# -lt 2 ]; then
	echo "$0 proteinListFile SeqDir [DBID]"
	echo "	This script builds MSAs for a list of proteins, one MSA for each protein"
	echo "	GNU parallel will be used to build MSAs for multiple proteins simultaneously"
	echo "	SeqDir: the folder for input sequence files"
	echo "	DBID: identification of uniclust database version: 2015, 2016, 2017 (default) and 2018"
	echo "	the resultant files will be saved to a folder named after MSA_DBID_E001"
	exit 1
fi

list=$1

SEQDIR=$2
if [ ! -d $SEQDIR ]; then
	echo $SEQDIR for sequence files not exist
	exit 1
fi

DBID=2017
if [ $# -eq 3 ]; then
	DBID=$3
fi

## number of CPUs for each protein
numCPUs=2

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
if [ ! -d $MSADir ]; then
	mkdir -p $MSADir
fi

program=$DistFeatureHome/HHblitsWrapper/BuildMSA4DistPred.sh
if [ ! -x $program ]; then
	echo invalid or non-executable program $program
	exit 1
fi

cat $list | parallel --memfree 20G --load 98% --delay 1 $program -d $DB -c $numCPUs -o $MSADir/ $SEQDIR/{}.fasta
