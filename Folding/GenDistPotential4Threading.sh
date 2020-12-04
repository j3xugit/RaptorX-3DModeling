#!/bin/bash

savefolder=`pwd`

function Usage 
{
	echo $0 "[ -d savefolder ] predictedDistProbMatrix"
	echo "	This script derives distance potential from predicted dist information for protein threading"
	echo "	predictedDistMatrix: a predicted dist matrix file, each ending with .predictedDistMatrix414C.pkl"
	echo "	savefolder: the folder for result saving, default current work directory"
}

while getopts ":d:" opt; do
        case ${opt} in
                d )
                  savefolder=$OPTARG
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

if [ $# -ne 1 ]; then
	Usage
	exit 1
fi

distFile=$1
if [ ! -f $distFile ]; then
	echo "ERROR: invalid file for predicted dist prob distribution $distFile"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/GenPairwisePotentialFromPrediction.py
if [ ! -f $program ]; then
	echo "ERROR: invalid file for program $program"
	exit 1
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

target=`basename $distFile | cut -f1 -d'.' `
potFile=$savefolder/${target}.pairPotential.DFIRE16.pkl
python $program -a CbCb -s $potFile $distFile 

#echo $target
#rawFile=$savefolder/${target}.pairPotential.DFIRE.18.1.61.Wt4D.pkl

#mv $rawFile $potFile
