#!/bin/bash

if [ $# -lt 2 ]; then
	echo "Usage: $0 DMName ModelFolder [NativeFolder]"
	echo "	DMName: the official domain name"
	echo "	MyDMName: the domain name used in folding the target "
	echo "	ModelFolder: the folder for all predicted 3D models"
	echo "	NativeFolder: the folder for ground truth, default CASP13DM-Native/"
	exit 1
fi

target=$1
MyDMName=$1
modelDir=$2

NativeDir=CASP13DM-Native/
if [ $# -ge 3 ]; then
	NativeDir=$3
fi

PID=$?

OUT=$(mktemp /tmp/${target}-${PID}.XXXXXXXXXX) || { echo "Failed to create temp file"; exit 1; }
cp /dev/null ${modelDir}/${target}-quality.txt

for i in ${modelDir}/${MyDMName}*.pdb 
do 
	DeepScore $i ${NativeDir}/${target}.pdb -n -2 > $OUT
	quality=`grep -A6 "1st input protein" $OUT | tail -1 | cut -f5- -d' '`
	b=`basename $i .pdb`
	echo $b $quality >> ${modelDir}/${target}-quality.txt
done

cut -f1,6- -d' ' ${modelDir}/${target}-quality.txt | sort -k4,4 -rn > ${modelDir}/${target}-quality.txt.sorted
rm -f  ${modelDir}/${target}-quality.txt

rm -f $OUT

