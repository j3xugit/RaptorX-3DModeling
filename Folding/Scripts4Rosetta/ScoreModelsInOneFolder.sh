#!/bin/bash

alpha=1.61
savefolder=`pwd`

function Usage 
{
	echo $0 "[-a alpha | -d savefolder ] ModelFolder cstFolder [propertyFolder]"
	echo "	This script scores unrelaxed 3D models in one folder"
	echo "	-a: alpha used in DFIRE potential, default $alpha. Not used when cstFolder is the folder for Rosetta constraint files"
	echo "		used only when cstFolder is the folder for predicted dist/orientation"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	ModelFolder: a folder for 3D decoys. A model folder shall have name like XXX-FoldResults where XXX is a target name or an alignment name"
	echo "	When propertyFolder is empty, cstFolder shall be a folder for Rosetta constraints. Otherwise cstFolder and propertyFolder shall be folders for predicted distance and Phi/Psi angles"
	echo "	The scores of all models in ModelFolder will be saved to a file in savefolder with name XXX-FoldResults.score where XXX is the basename of cstFolder"
}

while getopts ":a:d:" opt; do
        case ${opt} in
                a )
                  alpha=$OPTARG
                  ;;
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

if [ $# -lt 2 -o $# -gt 3 ]; then
	Usage
	exit 1
fi

modelFolder=$1
if [ ! -d $modelFolder ]; then
	echo "ERROR: invalid folder for 3D models: $modelFolder"
	exit 1
fi

cstFolder=$2
if [ ! -d $cstFolder ]; then
	echo "ERROR: invalid folder for predicted distance information or Rosetta distance constraints: $cstFolder "
	exit 1
fi

propertyFolder='cst'
if [ $# -eq 3 ]; then
	propertyFolder=$3
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

ScoreProgram=$cmdDir/ScoreModels.py
if [ ! -f $ScoreProgram ]; then
	echo "ERROR: invalid model scoring program $ScoreProgram"
	exit 1
fi

bname=`basename $modelFolder -FoldResults`
savefile=$savefolder/${bname}-FoldResults.score

if [ "$propertyFolder" == "cst" ]; then
	cstfile=$cstFolder/${bname}.pairPotential4Rosetta.SPLINE.txt
	if [ ! -f $cstfile ]; then
		echo "ERROR: invalid cstfile for scoring models in $modelFolder:", $cstfile
		exit 1
	fi
else
	propertyFile=${propertyFolder}/${bname}.predictedProperties.pkl
	if [ ! -f $propertyFile ]; then
		target=`echo $bname | cut -f1 -d'-'`
		propertyFile=${propertyFolder}/${target}.predictedProperties.pkl
	fi
	pairMatrixFile=${cstFolder}/${bname}.predictedDistMatrix.pkl

	cstfile="$pairMatrixFile -I $propertyFile,alpha=$alpha"

fi

python $ScoreProgram $modelFolder $cstfile -s $savefile
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run python $ScoreProgram $modelFolder $cstfile -s $savefile"
	exit 1
fi
