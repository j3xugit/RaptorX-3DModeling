#!/bin/bash

numModels=100

if [ $# -lt 3 ]; then
	echo $0 "alignmentFile PDBDir ResDir [numModels]"
	echo "	alignemntFile: a fasta file for a pairwise query-template alignment. Only a single template is supported"
	echo "		This file shall have name like QueryName-{*}-templateName.fasta where a template name could be XXXX_Y where XXXX is the PDB ID and Y is the chain ID"
	echo "	PDBDir: the folder for template structure file"
	echo "	numModels: the number of models to be generated, default $numModels"
	echo "	The result will be saved to ResDir/alnName-RelaxResults where alnName is the basename of alignmentFile"
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

## parse the name of alignment File into query and template

bname=`basename $alnFile .fasta`
queryName=`echo $bname | cut -f1 -d'-'`
templateName=`echo $bname | rev | cut -d'-' -f1 | rev`

#echo $queryName $templateName
ResDir=$ResDir/${bname}-RelaxResults
mkdir -p $ResDir

## revise the template name for Rosetta
newTPLName=`echo $templateName | sed -r 's/_//g' `

#echo $newTPLName

alnFile2=$ResDir/$bname.fasta
cp $alnFile $alnFile2

sed -i s/$templateName/$newTPLName/g $alnFile2

if [ -f $PDBDir/${templateName}.pdb ]; then
	cp $PDBDir/${templateName}.pdb $ResDir/${newTPLName}.pdb

elif [ -f $PDBDir/${templateName}.cif ]; then
	echo "WARNING: need to convert a .cif file to a .pdb file using CIF2PDB.py "
	$cmdDir/CIF2PDB.py $PDBDir/${templateName}.cif $ResDir/${newTPLName}.pdb
	#exit 1
else
	echo "ERROR: invalid PDB file in $PDBDir"
	exit 1
fi

RosettaHome=/mnt/data/RaptorXCommon/RaptorX-Threading/BuildModel_Package/
RosettaCmd=$RosettaHome/oneline_rosetta.sh

$RosettaCmd -i $alnFile2 -q $queryName -d $ResDir -o $ResDir -m $numModels -H $RosettaHome

if [ $? -ne 0 ]; then
        echo "ERROR: failed to run RosettaCM for $alnFile"
        exit 1
fi

mv $ResDir/${queryName}.rosetta_pdb $ResDir/TMP_ROSETTA_${queryName}/
mv $ResDir/TMP_ROSETTA_${queryName}/${queryName}.models/*.pdb $ResDir/
