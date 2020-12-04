#!/bin/bash

numModels=100

if [ $# -lt 3 ]; then
	echo $0 "alignmentFile PDBDir ResDir [numModels]"
	echo "	This script runs RosettaCM on an alignment file which may have multiple templates"
	echo "	alignemntFile: a fasta file for query-template alignment. The file shall have name Query-XXX.fasta where Query is protein under prediction and XXX represents one or multiple templates"
	echo "		A template shall have a name like yyyyZ where yyyy is the PDB ID and Z is the chain ID"
	echo "	PDBDir: the folder for template structure files"
	echo "	numModels: the number of models to be generated, default $numModels"
	echo "	The result will be saved to ResDir/alnName-RosettaResults where alnName is the basename of alignmentFile"
	exit 1
fi

alnFile=$1
if [ ! -f $alnFile ]; then
	echo "ERROR: invalid alignment file $alnFile"
	exit 1
fi

PDBDir=$2
if [ ! -d $PDBDir ]; then 
	echo "ERROR: invalid folder for template structure file"
	exit 1
fi

ResDir=$3
if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

if [ $# -ge 4 ]; then
	numModels=$4
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

## parse the name of alignment file into query

bname=`basename $alnFile .fasta`
queryName=`echo $bname | cut -f1 -d'-'`

ResDir=$ResDir/${bname}-RosettaResults
mkdir -p $ResDir

alnFile2=$ResDir/$bname.fasta
cp $alnFile $alnFile2

RosettaHome=/mnt/data/RaptorXCommon/RaptorX-Threading/BuildModel_Package/
RosettaCmd=$RosettaHome/oneline_rosetta.sh

$RosettaCmd -i $alnFile2 -q $queryName -d $PDBDir -o $ResDir -m $numModels -H $RosettaHome

if [ $? -ne 0 ]; then
        echo "ERROR: failed to run RosettaCM for $alnFile"
        exit 1
fi

mv $ResDir/${queryName}.rosetta_pdb $ResDir/TMP_ROSETTA_${queryName}/
mv $ResDir/TMP_ROSETTA_${queryName}/${queryName}.models/*.pdb $ResDir/
