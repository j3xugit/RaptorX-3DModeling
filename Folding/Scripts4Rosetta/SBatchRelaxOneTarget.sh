#!/bin/bash

savefolder=`pwd`

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome to installation folder of Folding"
        exit 1
fi

filter=''
partition=contrib-cpu

function Usage 
{
	echo $0 "[ -f filter4initiModels | -d savefolder | -q partition ] initModelFolder predictedPairInfo predictedPropertyInfo"
	echo "	This script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then relax 3D models"
	echo "	initModelFolder: a folder for initial models; all model files shall end with .pdb"
	echo "  predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance and orientation information"
        echo "  predictedPropertyInfo: a string cst or a PKL file for predicted Phi/Psi angles"
        echo "          when it is cst, it indicates that predictedPairInfo is a Rosetta constraint file instead of a PKL file"
        echo "  -d: the folder for result saving, default current work directory"
	echo "	-f: a filter string for the selection of intial models to be relaxed, default empty, i.e., all PDB files in initModelFolder will be relaxed"
	echo "	-q: compute node queue name in a slurm cluster, default $partition"
}

while getopts ":d:f:q:" opt; do
        case ${opt} in
                d )
                  savefolder=$OPTARG
                  ;;
		f )
		  filter=$OPTARG
		  ;;
		q )
		  partition=$OPTARG
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

if [ $# -lt 3 ]; then
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

propertyFile=$3
if [ "$propertyFile" != "cst" ]; then
	if [ ! -f $propertyFile ]; then
		echo "ERROR: invalid file for predicted property info: $propertyFile"
		exit 1
	fi
fi

program=${DistanceFoldingHome}/Scripts4Rosetta/GenPotentialNFoldRelax.sh
if [ ! -x $program ]; then
        echo "ERROR: invalid program for relaxation $program"
        exit 1
fi
keywords=`basename $program`

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

for model in $initModelFolder/*filter*pdb
do
	bname=`basename $model .pdb`
        sbatch -p $partition -J ${bname}R -o $bname.out -e $bname.out -C avx $program -r 2 -d $savefolder $model $pairMatrixFile $propertyFile
done
