#!/bin/bash

if [ -z "$DistFeatureHome" ]; then
	echo "ERROR: please set environmental variable DistFeatureHome to the folder of BuildFeatures/ "
	exit 1
fi

if [ $# -lt 2 ]; then
	echo "$0 proteinListFile SeqDir [DBID]"
	echo "  SeqDir: the folder for input sequence files"
        echo "  DBID: identification of uniclust database version: 2015, 2016, 2017 and 2018(default)"
        echo "  the resultant files will be saved to a folder named MSA4Threading_{DBID}_E001"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/BuildMSA4Threading.sh
if [ ! -x $program ]; then
        echo "ERROR: invalid or non-executable program $program"
        exit 1
fi

list=$1
if [ ! -f $list ]; then
        echo "ERROR: invalid protein list file $list"
        exit 1
fi

SEQDIR=$2
if [ ! -d $SEQDIR ]; then
        echo "ERROR: invalid folder for protein sequences $SEQDIR"
        exit 1
fi


DBID=2018
if [ $# -eq 3 ]; then
	DBID=$3
fi

DB="uniclust30_2017_10/uniclust30"
if [ $DBID == "2015" ]; then
	DB="uniprot20_2015_06/uniprot20"
elif [ $DBID == "2016" ]; then
	DB="uniprot20_2016_05/uniprot20"
elif [ $DBID == "2018" ]; then
	DB="uniclust30_2018_08/uniclust30"
fi
DB=$DistFeatureHome/HHblitsWrapper/databases/${DB}

MSADir=MSA4Threading_${DBID}_E001
if [ ! -d $MSADir ]; then
	mkdir -p $MSADir
fi

cat $list | parallel $program -s $MSADir -d $DB $SEQDIR/{}.fasta
