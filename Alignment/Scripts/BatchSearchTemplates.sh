#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 proteinListFile queryDir [ResultDir]
	exit 1
fi

proteinListFile=$1
queryDir=$2

ResultDir=`pwd`
if [ $# -ge 3 ]; then
	ResultDir=$3
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

mkdir -p $ResultDir

for i in `cat $proteinListFile`
do
	$cmdDir/SearchTemplates.sh $queryDir/$i.hhm $ResultDir
done
