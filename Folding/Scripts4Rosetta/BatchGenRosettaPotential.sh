#!/bin/bash

savefolder=`pwd`

atomType="CbCb+CaCa+NO+TwoROri"

PotentialType=DFIRE
Alpha=1.61
MaxDist=18
seqSep=1

## TopRatio * SeqLen is the number of top orientation constraints to be used
TopRatio=25

## weight for Phi/Psi potential
w4phipsi=1

## 
querySeqFolder=/dev/null

function Usage {
	echo "$0 [ -A atomType | -a alpha | -c distCutoff | -t TopRatio | -w w4phipsi | -s seqSep | -q querySeqFolder | -d savefolder ] targetListFile folder4predictedPairInfo folder4predictedProperty "
	echo "	This script generates rosetta constraints for a list of proteins from predicted distance/orientation and property information"
        echo "	targetListFile: a file for a list of protein names, each in one row"
        echo "	folder4predictedPairInfo: a folder for predicted distance and orientation files in PKL format, each ending with .predictedDistMatrix.pkl"
        echo "	folder4predictedProperty: a folder for predicted Phi Psi angles in PKL format, each ending with .predictedProperty.pkl"
	echo " "
	echo "	-A: the type of atom pairs and orientation to be used in generating potential, default $atomType"
        echo "	-a: alpha for DFIRE distance reference state, default $Alpha"
        echo "	-c: the maximum distance cutoff for DFIRE distance potential, default $MaxDist"
        echo "	-t: the number of orientation constraints to be used per residue, default $TopRatio"
        echo "	-w: the weight for backbone Phi/Psi angle information, default $w4phipsi"
        echo "	-s: sequence separation for distance potential (default $seqSep), i.e., using two atoms only if their residue indices i, j satisfy abs(i-j)>= this value"
	echo " "
	echo "	-q: the folder for the sequence files, default empty. When not empty, sequence files will be used to check sequence consistency with the predicted distance/angle files"
        echo "	-d: the folder for result saving, default current work directory"
	echo "		one file and one subfolder will be generated under this folder for each protein"
}

while getopts ":A:a:c:t:w:s:q:d:" opt; do
        case ${opt} in
		A )
		  atomType=$OPTARG
		  ;;
		a )
                  Alpha=$OPTARG
                  ;;
                c )
                  MaxDist=$OPTARG
                  ;;
		t )
		  TopRatio=$OPTARG
		  ;;
                w )
                  w4phipsi=$OPTARG
                  ;;
                s )
                  seqSep=$OPTARG
                  ;;
		q )
		  querySeqFolder=$OPTARG
		  if [ ! -d $querySeqFolder ]; then
			echo "ERROR: invalid folder for query sequences: $querySeqFolder"
			exit 1
		  fi
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

if [ $# -ne 3 ]; then
	Usage
	exit 1
fi

targetList=$1
if [ ! -f $targetList ]; then
        echo "ERROR: invalid target list file $targetList"
        exit 1
fi

pairFolder=$2
if [ ! -d $pairFolder ]; then
        echo "ERROR: invalid folder for predicted distance/orientation files $pairFolder"
        exit 1
fi

propertyFolder=$3
if [ ! -d $propertyFolder ]; then
        echo "ERROR: invalid folder for predicted phi/psi angles $propertyFolder"
        exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/GenRosettaPotential.sh
if [ ! -x $program ]; then
        echo "ERROR: non-existing or non-excutable program $program"
        exit 1
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

options="-A $atomType -a $Alpha -c $MaxDist -t $TopRatio -w $w4phipsi -s $seqSep -d $savefolder"

for i in `cat $targetList`
do
	propertyFile=$propertyFolder/${i}.predictedProperties.pkl
	if [ ! -f $propertyFile ]; then
		seqName=`echo $i | cut -f1 -d'-'`
		propertyFile=$propertyFolder/${seqName}.predictedProperties.pkl
	fi

	if [ -d $querySeqFolder ]; then
		querySeqFile=$querySeqFolder/${i}.fasta
		if [ ! -f $querySeqFile ]; then
			querySeqFile=$querySeqFolder/${i}.seq
		fi
		if [ -f querySeqFile ]; then
			tmpOptions=$options" -q $querySeqFile "
		else
			tmpOptions=$options
		fi
	fi
		
	$program $tmpOptions $pairFolder/${i}.predictedDistMatrix.pkl $propertyFile &

	sleep 1
done

wait
