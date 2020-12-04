#!/bin/bash

#if [[ -z "${DistanceFoldingHome}" ]]; then
#        echo "ERROR: Please set the environmental variable DistanceFoldingHome to installation folder of Folding/"
#        exit 1
#fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

savefolder=`pwd`
alpha=1.61
filter=''
## at most this number of jobs to be run simultaneously
numAllowedJobs=`nproc --all`

function Usage 
{
	echo $0 "[ -f filter4initiModels | -d savefolder | -c numCPUs | -a alpha ] initModelFolder predictedPairInfo [predictedPropertyInfo]"
	echo "	This script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then relax 3D models"
	echo "	initModelFolder: a folder for initial models; all model files shall end with .pdb"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance and orientation information"
        echo "	predictedPropertyInfo (optional): a string cst or a PKL file for predicted Phi/Psi angles"
        echo "	     when it is empty or equal to cst, predictedPairInfo shall be a Rosetta constraint file instead of a PKL file"
        echo "	-c: the number of CPUs to be used, default $numAllowedJobs"
        echo "	-d: the folder for result saving, default current work directory"
	echo "	-f: a filter string for the selection of intial models to be relaxed, default empty, i.e., all PDB files in initModelFolder will be relaxed"
	echo "	-a: alpha value for DFIRE potential, default 1.61. If > 20, a random value will be used"
}

while getopts ":c:d:f:a:" opt; do
        case ${opt} in
                c )
                  numAllowedJobs=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
                  ;;
		f )
		  filter=$OPTARG
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

#program=${DistanceFoldingHome}/Scripts4Rosetta/GenPotentialNFoldRelax.sh
program=${cmdDir}/GenPotentialNFoldRelax.sh
if [ ! -x $program ]; then
        echo "ERROR: invalid model relaxation program $program"
        exit 1
fi
keywords=`basename $program`

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

for model in $initModelFolder/*${filter}*pdb
do
	if [ ! -f $model ]; then
		continue
	fi

	bname=`basename $model .pdb`
	# if there is already one relaxed model, skip
	if [ `ls -1 $savefolder/$bname.relaxed.* 2>/dev/null | wc -l ` -gt 0 ]; then
                continue
        fi

        while true
        do
                numRunningJobs=`ps -x | grep ${keywords} | wc -l  `
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
                        $program -r 2 -d $savefolder -a $alpha $model $pairMatrixFile $propertyFile &
			sleep 2
                        break
		else
                	a=`expr $RANDOM % 5 `
                	sleep $a
                fi
        done
done

wait
