#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 "groupFile ListDir SpecFile [NoTemplate]"
	echo "	This script generates metaData files for model training"
	echo "	groupFile: a file ending with .group.txt. When templates are not used, this file is just a list of proteins for training and validation, each in one row"
	echo "		When templates are not used, each line in this file starts with a training/validation protein and then followed by some seq-template pairs"
	echo "	ListDir: a folder containing some list files. Each contains a list of proteins (one in each row) for training or validation"
	echo "	SpecFile: this file contains path information for needed seq, feature and ground truth files"
	echo "	NoTemplate: if specified, not use templates; default use templates"
	exit 1
fi

groupFile=$1
if [ ! -f $groupFile ]; then
	echo "ERROR: invalid group file $groupFile"
	exit 1
fi

ListDir=$2
if [ ! -d $ListDir ]; then
	echo "ERROR: invalid folder for list files $ListDir"
	exit 1
fi

SpecFile=$3
if [ ! -f $SpecFile ]; then
	echo "ERROR: invalid specification file $SpecFile"
	exit 1
fi

NoTemplate=0
if [ $# -ge 4 ]; then
	NoTemplate=1
fi

groupName=`basename $groupFile .group.txt`

## generate group.txt files for each list in ListDir
for i in $ListDir/*.list
do
	b=`basename $i | cut -f2 -d'.' `
	join $i $groupFile > $groupName.$b.group.txt &
done

wait

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for gf in $groupName.train?.group.txt
do
	if [ $NoTemplate -eq 1 ]; then
		python $cmdDir/GenerateMetaData.py $gf $SpecFile &
	else
		python $cmdDir/GenerateMetaData.py $gf $SpecFile 2 &
	fi
done

for gf in $groupName.valid?.group.txt
do
	if [ $NoTemplate -eq 1 ]; then
		python $cmdDir/GenerateMetaData.py $gf $SpecFile 1 &
	else
		python $cmdDir/GenerateMetaData.py $gf $SpecFile 3 &
	fi
done

for gf in $groupName.small?.group.txt
do
	if [ $NoTemplate -eq 1 ]; then
		python $cmdDir/GenerateMetaData.py $gf $SpecFile 1 &
	else
		python $cmdDir/GenerateMetaData.py $gf $SpecFile 3 &
	fi
done

wait

rename group.txt.m m *.group.txt.metaData.json
