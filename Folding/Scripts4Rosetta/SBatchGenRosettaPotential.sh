#!/bin/bash

avefolder=`pwd`
w4phipsi=1
seqSep=1
MaxDist=18
Alpha=1.61
TopRatio=25

host=`hostname`

if [ "$host" == "raptorx3.uchicago.edu" ]; then
	QUEUE=cpu
else
	QUEUE=contrib-cpu
fi

function Usage {
	echo "$0 [-a Alpha | -c distCutoff | -t TopRatio | -w w4phipsi | -s seqSep | -d savefolder | -q jobqueue ] targetListFile folder4predictedPairInfo folder4predictedProperty "
	echo "	This script runs on a slurm cluster to generate rosetta constraints for a list of proteins from predicted distance/orientation and property information"
	echo "	targetListFile: a file for a list of protein names, each in one row"
	echo "	folder4predictedPairInfo: a folder for predicted distance and orientation files in PKL format, each ending with .predictedDistMatrix.pkl"
        echo "	folder4predictedProperty: a folder for predicted Phi Psi angles in PKL format, each ending with .predictedProperty.pkl"
	echo "	-a: alpha for DFIRE distance reference state, default $Alpha"
        echo "	-c: the maximum distance cutoff for DFIRE distance potential, default $MaxDist"
        echo "	-t: the number of orientation constraints to be used per residue, default $TopRatio"
        echo "	-w: the weight for backbone Phi/Psi angle information, default $w4phipsi"
        echo "	-s: sequence separation for distance potential (default $seqSep), i.e., using two atoms only if their residue indices i, j satisfy abs(i-j)>= this value"
        echo "	-d: the folder for result saving, default current work directory"
	echo "	-q: the job queue in a slurm cluster, default $QUEUE. You may use contrib-cpu-long for large proteins."
}

while getopts ":a:c:t:w:s:d:q:" opt; do
        case ${opt} in
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
                d )
                  savefolder=$OPTARG
                  ;;
		q )
		  QUEUE=$OPTARG
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
if [ ! -f $program ]; then
	echo "ERROR: invalid path info for potential generation program $program"
	exit 1
fi


for i in `cat $targetList`
do
	propertyFile=$propertyFolder/${i}.predictedProperties.pkl
        if [ ! -f $propertyFile ]; then
                seqName=`echo $i | cut -f1 -d'-'`
                propertyFile=$propertyFolder/${seqName}.predictedProperties.pkl
        fi

	if [ "$host" == "raptorx3.uchicago.edu" ]; then
		sbatch -p $QUEUE -J $i-potential -o $i-potential.out -e $i-potential.out $program -a $Alpha -c $MaxDist -t $TopRatio -w $w4phipsi -s $seqSep -d $savefolder $pairFolder/${i}.predictedDistMatrix.pkl $propertyFile 
	else
		sbatch -p $QUEUE -J $i-potential -o $i-potential.out -e $i-potential.out -C avx -x cpu1,cpu23 $program -a $Alpha -c $MaxDist -t $TopRatio -w $w4phipsi -s $seqSep -d $savefolder $pairFolder/${i}.predictedDistMatrix.pkl $propertyFile 
	fi
done
