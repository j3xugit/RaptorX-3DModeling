#!/bin/sh

if [ -z "$DistFeatureHome" ]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the folder of BuildFeatures/ "
        exit 1
fi

numAllowedJobs=4
numCPUs=2

DBID=2017
DB="uniclust30_2017_10/uniclust30"

ResultDir=""

if [ $# -lt 2 ]; then
	echo "$0 [ -d DBID | -n numJobs | -c numCPUs | -o ResultDir ] proteinListFile SeqDir"
	echo "	This script builds MSAs for a list of proteins, one MSA for each protein"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	SeqDir: the folder for input sequence files"
	echo "	DBID: identification of uniclust database version: 2015, 2016, 2017 (default), 2018 and 2020"
	echo "	numJobs: the number of protein sequences to be run simultaneously, default $numAllowedJobs"
	echo "	numCPUs: the number of CPUs to be used for each protein sequence, default $numCPUs"
	echo "	ResultDir: the folder for result saving"
	echo "		by default, the resultant files will be saved to a folder MSA_DBID_E001"
	exit 1
fi

while getopts ":o:d:c:n:" opt; do
        case ${opt} in
                o )
                  ResultDir=$OPTARG
                  ;;
                c )
                  numCPUs=$OPTARG
                  ;;
                n )
                  numAllowedJobs=$OPTARG
                  ;;
		d )
		  DBID=$OPTARG
		  ;;
                \? )
                  echo "Invalid Option: -$OPTARG" 1>&2
                  exit 1
                  ;;
                : )
                  echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                  exit 1
                  ;;
        esac
done
shift $((OPTIND -1))

if [ $# -ne 2 ]; then
        Usage
        exit 1
fi


list=$1

SEQDIR=$2
if [ ! -d $SEQDIR ]; then
	echo $SEQDIR for sequence files not exist
	exit 1
fi

if [ $DBID == "2015" ]; then
	DB="uniprot20_2015_06/uniprot20"
elif [ $DBID == "2016" ]; then
	DB="uniprot20_2016_05/uniprot20"
elif [ $DBID == "2018" ]; then
	DB="uniclust30_2018_08/uniclust30"
elif [ $DBID == "2020" ]; then
	DB="uniclust30/uniclust30"
fi
DB=$DistFeatureHome/HHblitsWrapper/databases/${DB}

if [ -z $ResultDir ]; then
	MSADir=MSA_${DBID}_E001
else
	MSADir=$ResultDir
fi

if [ ! -d $MSADir ]; then
	mkdir -p $MSADir
fi

program=$DistFeatureHome/HHblitsWrapper/BuildMSA4DistPred.sh
if [ ! -x $program ]; then
	echo invalid or non-executable program $program
	exit 1
fi

keywords=`basename $program`
myself=`basename $0 `

for i in `cat $list`
do
        while true
        do
                ## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -le $numAllowedJobs ]; then
			$program -d $DB -c $numCPUs -o $MSADir/ $SEQDIR/$i.fasta &
                        break
                else
                        a=`expr $RANDOM % 4 `
                        sleep $a
                fi
        done
done

wait
