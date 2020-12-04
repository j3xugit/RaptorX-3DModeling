#!/bin/bash

savefolder=`pwd`
numAllowedJobs=15
host=`hostname`
if [ "$host" == "raptorx10.uchicago.edu" ]; then
	numAllowedJobs=50
fi

function Usage 
{
	echo $0 "[ -n numJobs | -d savefolder ] proteinListFile folder4predictedDistProbMatrix"
	echo "	This script derives distance potential from predicted dist information for protein threading"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	folder4predictedDistMatrix: a folder for predicted dist matrix files, each ending with .predictedDistMatrix414C.pkl"
	echo "	numJobs: the number of jobs to be simultaneously run, default $numAllowedJobs"
	echo "	savefolder: the folder for result saving, default current work directory"
}

while getopts ":n:d:" opt; do
        case ${opt} in
                n )
                  numAllowedJobs=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
                  ;;
                g )
                  GPU=$OPTARG
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

proteins=$1
if [ ! -f $proteins ]; then
	echo "ERROR: invalid protein list file $proteins"
	exit 1
fi

inputFolder=$2
if [ ! -d $inputFolder ]; then
	echo "ERROR: invalid folder for predicted dist prob distribution "
	exit 1
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

#echo $proteins
#echo $inputFolder
#echo $savefolder

program=GenPairwisePotentialFromPrediction.py

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

keywords=$program
myself=`basename $0 `

for i in `cat $proteins`
do
        while true
        do
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l  `
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
			distFile=$inputFolder/${i}.predictedDistMatrix414C.pkl
                        python $cmdDir/$program -a CbCb -s $savefolder $distFile &
			sleep 1
                        break
		else
                	a=`expr $RANDOM % 3 `
                	sleep $a
                fi
        done
done

wait

