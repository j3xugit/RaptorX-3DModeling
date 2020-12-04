#!/bin/bash

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome to installtion folder of Folding"
        exit 1
fi

savefolder=`pwd`
filter=''
alpha=1.61

function Usage 
{
	echo $0 "[ -d savefolder | -f filter4modelNames | -a alpha ] initModelFolder predictedPairInfo [predictedPropertyInfo]"
	echo "	This script prints out a list of jobs for relaxation"
	echo "	initModelFolder: a folder for initial models"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance and orientation"
	echo "	predictedPropertyInfo: could be empty, a string 'cst' or a PKL file for predicted Phi/Psi angles"
	echo "	    when empty or 'cst', predictedPairInfo shall be a Rosetta constraint file instead of a PKL file"
	echo "	    Otherwise, both predictedPairInfo and predictedPropertyInfo shall be a PKL file"
	echo "	-f: a filter string used to select initial models for relaxation, default empty, i.e., relax all models"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-a: alpha used in DFIRE potential, default $alpha. If larger than 20, one value in [1.57, 1.63] will be randomly chosen."
}

while getopts ":d:f:a:" opt; do
        case ${opt} in
                f )
                  filter=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
                  ;;
		a )
		  alpha=$OPTARG
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

if [ $# -lt 2 -o $# -gt 3 ]; then
	Usage
	exit 1
fi

initModelFolder=$1
if [ ! -d $initModelFolder ]; then
	echo "ERROR: invalid folder for initial models $initModelFolder"
	exit 1
fi

pairMatrixFile=$2
if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: invalid file for predicted distance/orientation info: $pairMatrixFile"
	exit 1
fi

propertyFile="cst"
if [ $# -ge 3 ]; then
        propertyFile=$3
fi
if [ "$propertyFile" != "cst" ]; then
        if [ ! -f $propertyFile ]; then
                echo "ERROR: invalid file for predicted property info: $propertyFile"
                exit 1
        fi
fi


program=${DistanceFoldingHome}/Scripts4Rosetta/GenPotentialNFoldRelax.sh
if [ ! -x $program ]; then
	echo "ERROR: invalid or non-executable program $program"
	exit 1
fi

for model in $initModelFolder/*${filter}*.pdb
do
	bname=`basename $model .pdb`
	if [ `ls -1 $savefolder/$bname.relaxed.* 2>/dev/null | wc -l ` -gt 0 ]; then
		continue
	fi
        echo "$program -r 2 -d $savefolder -a $alpha $model $pairMatrixFile $propertyFile"
done
