#!/bin/bash


if [ "$#" -lt 4 ]; then
	echo $0 proteinName SeqFolder DistanceFolder AngleFolder [ atomPairType [ boundType [ cutoff ] ] ] 
	echo "  DistanceFolder: the folder for predicted distance file"
	echo "  AngleFolder: the folder for predicted angle and secondary structure file"
	echo "  atomPairType: the atom pair type, e.g., All, CaCa, CbCb, CaCa+CbCb+CgCg, default All "
	echo "  boundType specifies the type of restraints, negsative indicates contact restraints, default is 3"
	echo "  cutoff specifies the distance cutoff for restraints (<=15). When contacts are used, it is Zscore for contact selection."
	echo "	   the default value for cutoff is 15 (for distance restraints) and 3 (for contact selection)" 
	echo "  Example: " $0 " T0859-D1 CASP12FM_DistancePrediction CASP12FM_PropertyPrediction CbCb 3 15"
	exit 1
fi

if [ -z "${DistanceFoldingHome}" ]; then
	echo "Please set the environmental variable DistanceFoldingHome to the installation folder of the Folding module, e.g., $HOME/3DModleing/Folding"
	exit 1
fi

FoldHome=$DistanceFoldingHome
program=$FoldHome/scripts/FoldByPredictedRestraints.sh

target=$1
SeqDir=$2
DistMatrixDir=$3
AngleDir=$4
ResultDir=FoldingResult-`basename ${DistMatrixDir}`-`basename ${AngleDir}`/

dihedralEnergyConstant=5.0
atomPairTypes='All'
boundType=3
defaultDistCutoff=17
defaultZCutoff=3


distFile=${DistMatrixDir}/${target}.predictedDistMatrix.pkl
if [ ! -f ${distFile} ]; then
	echo "file does not exist: " ${distFile}
	exit 1
fi

angleFile=${AngleDir}/${target}.predictedProperties.pkl
if [ ! -f ${angleFile} ]; then
	echo "file does not exist: " ${angleFile}
	exit 1
fi


if [ ! -d ${ResultDir} ]; then
	mkdir -p  ${ResultDir}
fi

if [ "$#" -ge 5 ]; then
	atomPairTypes=$5
fi

if [ "$#" -ge 6 ]; then
	boundType=$6
fi

if [ ${boundType} -lt 0 ]; then
	cutoff=$defaultZCutoff
else
	cutoff=$defaultDistCutoff
fi

if [ "$#" -ge 7 ]; then
	cutoff=$7
fi

if (( $(echo "$cutoff < 0.1 " | bc -l) )); then
	echo 'Error: the cutoff value is too small'
	exit 1
fi


if (( $(echo "$cutoff > 3.5 " | bc -l) )); then
	if [ ${boundType} -lt 0 ]; then
		echo 'Error: the cutoff for contact Zscore may be too big.'
		exit 1
	fi
fi

if (( $(echo "$cutoff > 25 " | bc -l) )); then
	if [ ${boundType} -ge 0 ]; then
		echo 'Error: the cutoff for distance restraints may be too big.'
		exit 1
	fi
fi

OutDir=${ResultDir}/${target}_${atomPairTypes}_t${boundType}_c${cutoff}_A${dihedralEnergyConstant}/
##check if the 5 models have been generated or not
if [ -e ${OutDir}/${target}_model5.pdb ]; then
	echo 'There are already 5 models generated for ' $target
	exit 0
fi

#echo $OutDir

seqFile=${SeqDir}/${target}.fasta
if [ ! -f $seqFile ]; then
	seqFile=${SeqDir}/${target}.seq
fi

if [ ! -f $seqFile ]; then
	echo "ERROR: cannot locate the input sequence file!"
	exit 1
fi

cmdStr="${program} -i $seqFile -d ${distFile} -p ${angleFile} -e ${dihedralEnergyConstant} -t $boundType -a ${atomPairTypes} -c $cutoff -o ${OutDir} -m 20:5"
$cmdStr
