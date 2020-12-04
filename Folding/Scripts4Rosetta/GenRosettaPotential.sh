#!/bin/bash

atomType="CbCb+CaCa+NO+TwoROri"

PotentialType=DFIRE
MaxDist=18
Alpha=1.61
seqSep=1

## TopRatio * SeqLen is the number of top orientation constraints to be used
TopRatio=25

## weight for Phi/Psi potential
w4phipsi=1

querySeqFile=""
savefolder=`pwd`

function Usage {
	echo "$0 [-A atomType | -a alpha | -c distCutoff | -t topRatio | -w w4phipsi | -s seqSep | -q querySeqFile | -d saveFolder] predictedPairPKLfile predictedPropertyPKLfile"
	echo "	This script generates Rosetta constraints from predicted inter-atom distance/orientation and backbone phi/psi angles"
	echo "		The generated distance potential is built upon the DFIRE reference state"
	echo "	predictedPairFile: a PKL file for predicted inter-residue distance and orientation information"
	echo "	predictedPropertyFile: a PKL file for predicted backbone phi/psi angles"
	echo "	-q: the protein sequence file in FASTA format. When provided, check the sequence consistency between querySeqFile and predictedPairPKLfile/predictedPropertyPKLfile"
	echo "	-d: the folder for result saving, default current work directory"
	echo "		The resulant file has name XXX.pairPotential4Rosetta.SPLINE.txt where XXX is the base name of predictedPairPKLfile"
	echo " "
	echo "	-A: the atom pair type for distance and orientation information, default $atomType"
	echo "	-a: alpha for DFIRE distance reference state, default $Alpha"
	echo "	-c: the maximum distance cutoff for DFIRE distance potential, default $MaxDist"
	echo "	-t: the number of top orientation constraints to be used per residue, default $TopRatio"
	echo "	-w: the weight for backbone Phi/Psi angle constraint, default $w4phipsi"
	echo "	-s: sequence separation for inter-atom distance potential (default $seqSep), i.e., considering two atoms only if their residue index difference is at least this value"
}

if [[ -z "$DistanceFoldingHome" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome to the installation folder of Folding"
        exit 1
fi

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
		  querySeqFile=$OPTARG
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

if [ $# -ne 2 ]; then
	Usage
	exit 1
fi

pairMatrixFile=$1
if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: invalid file for predicted distance/orientation info: $pairMatrixFile"
	exit 1
fi

target=`basename $pairMatrixFile`
target=`echo $target | cut -f1 -d'.'`
#echo target=$target

propertyFile=$2
if [ ! -f $propertyFile ]; then
	echo "ERROR: invalid file for predicted property info: $propertyFile"
	exit 1
fi

date; echo "Generating Rosetta potential from $pairMatrixFile and $propertyFile ..."

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

PotentialDir=$savefolder


## generate pairwise potential
PairPotentialFile=$PotentialDir/${target}.pairPotential.$PotentialType.$MaxDist.$Alpha.pkl
python $DistanceFoldingHome/GenPairwisePotentialFromPrediction.py -s $PairPotentialFile -a $atomType -r $PotentialType+$MaxDist+$Alpha $pairMatrixFile 

## the resultant CST file has name ${target}.pairPotential4Rosetta.SPLINE.txt
cstfile=${PotentialDir}/${target}.pairPotential4Rosetta.SPLINE.txt

options="-t $TopRatio -s $seqSep -d $PotentialDir -a $atomType"
if [ ! -z "$querySeqFile" ]; then
	options=$options" -q $querySeqFile "
fi

python $DistanceFoldingHome/Scripts4Rosetta/GeneratePairPotential4Rosetta.py $options $PairPotentialFile
if [ $? -ne 0 -o ! -f $cstfile ]; then
        echo "ERROR: failed to run python $DistanceFoldingHome/Scripts4Rosetta/GeneratePairPotential4Rosetta.py $options $PairPotentialFile"
        exit 1
fi

## add Phi/Psi potential
#propertyFile2=$PotentialDir/${target}.predictedProperty.pkl
#ln -s `readlink -f $propertyFile` $propertyFile2

PhiPsiPotentialFile=$PotentialDir/${target}.PhiPsi4AMBERPERIODIC.w${w4phipsi}.txt
options="  -w ${w4phipsi} -s $PhiPsiPotentialFile "
if [ ! -z "$querySeqFile" ]; then
	options=$options" -q $querySeqFile "
fi

python $DistanceFoldingHome/GenPropertyPotential4Rosetta.py $options $propertyFile
if [ $? -ne 0 -o ! -f $PhiPsiPotentialFile ]; then
        echo "ERROR: failed to run python $DistanceFoldingHome/GenPropertyPotential4Rosetta.py $options $propertyFile"
        exit 1
fi
cat $PhiPsiPotentialFile >> $cstfile

rm -f $PairPotentialFile
rm -f $PhiPsiPotentialFile
#rm -f $propertyFile2

date; echo "The resultant Rosetta constraint files are saved at at $cstfile"
