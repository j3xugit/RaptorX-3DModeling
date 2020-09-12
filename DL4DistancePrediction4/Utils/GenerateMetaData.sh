#!/bin/sh

if [ $# -lt 3 ]; then
	echo $0 groupFile ListDir SpecFile
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
	python $cmdDir/GenerateMetaData.py $gf $SpecFile 2 &
done

for gf in $groupName.valid?.group.txt
do
	python $cmdDir/GenerateMetaData.py $gf $SpecFile 3 &
done

for gf in $groupName.small?.group.txt
do
	python $cmdDir/GenerateMetaData.py $gf $SpecFile 3 &
done

wait

rename group.txt.m m *.group.txt.metaData.json
