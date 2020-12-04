#!/bin/bash

savefolder=`pwd`
filter=""
alpha=1.61

function Usage 
{
	echo $0 "[ -d savefolder | -f filter4models | -a alpha] proteinListFile metaFolder4initialModels cstFolder [propertyFolder]"
	echo "	This script prints out relaxation jobs for a list of proteins to be submitted to a cluster or GNU parallel"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	metaFolder4initModels: the meta folder for initial 3D models, in which each target shall have a subfolder target-FoldResults/"
	echo "	cstFolder: the folder for predicted dist/orientation files (which have names like proteinName.predictedDistMatrix.pkl) or Rosetta constraint files (which have names like proteinName.pairPotential4Rosetta.SPLINE.txt)"
	echo "	propertyFolder: empty or a folder. When empty, cstFolder shall be the folder for Rosetta constraint files. Otherwise cstFolder and propertyFolder shall be folders for predicted distance/orientation and property files, respectively"
	echo "	-d: the folder for the final 3D models, default current work directory"
	echo "	-f: a filter string used to select models (through file names) to be relaxed, default empty, i.e., all models are selected"
	echo "	-a: alpha used in DFIRE potential, default $alpha; if larger than 20, one value between 1.57 and 1.63 will be randomly chosen."
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

if [ $# -lt 3 -o $# -gt 4 ]; then
	Usage
	exit 1
fi

proteinList=$1
if [ ! -f $proteinList ]; then
        echo "ERROR: invalid protein list file $proteinListFile"
        exit 1
fi

inDir=$2
if [ ! -d $inDir ]; then
        echo "ERROR: invalid folder for input files $inDir"
        exit 1
fi

cstFolder=$3
if [ ! -d $cstFolder ]; then
        echo "ERROR: invalid folder for predicted distance/orientation or Rosetta constraint files: $cstFolder"
        exit 1
fi

propertyFolder='cst'
if [ $# -eq 4 ]; then
	propertyFolder=$4
fi
if [ "$propertyFolder" != "cst" ]; then
	if [ ! -d $propertyFolder ]; then
		echo "ERROR: invalid folder for predicted property files: $propertyFolder"
		exit 1
	fi
fi

if [ ! -d $savefolder ]; then
        mkdir -p $savefolder
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for target in `cat $proteinList`
do
	if [ "$propertyFolder" == "cst" ]; then
                cstFile=$cstFolder/${target}.pairPotential4Rosetta.SPLINE.txt
                propertyFile=""
        else
                cstFile=$cstFolder/${target}.predictedDistMatrix.pkl
                propertyFile=${propertyFolder}/${target}.predictedProperties.pkl
                if [ ! -f $propertyFile ]; then
                        seqName=`echo $target | cut -f1 -d'-' `
                        propertyFile=${propertyFolder}/${seqName}.predictedProperties.pkl
                fi
        fi

	if [ -z "${filter}" ]; then
		$cmdDir/PrintJob4RelaxOneTarget.sh -a $alpha -d $savefolder/${target}-RelaxResults $inDir/${target}-FoldResults $cstFile $propertyFile
	else
		$cmdDir/PrintJob4RelaxOneTarget.sh -a $alpha -d $savefolder/${target}-RelaxResults -f $filter $inDir/${target}-FoldResults $cstFile $propertyFile
	fi
done
