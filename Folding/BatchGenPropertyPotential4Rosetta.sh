#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 proteinListFile inputFolder [savefolder]
	echo	"proteinListFile: a list file for proteins"
	echo	"inputFolder: the folder containing predicted property files for the proteins in the list. Each file shall end with .predictedProperties.pkl "
	echo	"savefolder: the folder for result saving, default current work directory"
	exit 1
fi

proteinListFile=$1
inputFolder=$2

saveFolder=`pwd`
if [ $# -ge 3 ]; then
	saveFolder=$3
	if [ ! -d $saveFolder ]; then
		mkdir -p $saveFolder
	fi
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/GenPropertyPotential4Rosetta.py

for i in `cat $proteinListFile`
do
	python $program -s $saveFolder $inputFolder/${i}.predictedProperties.pkl
done
