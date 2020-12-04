#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "target Folder (the folder for the folding job results)"
	exit 1
fi

target=$1
folder=$2

if [ ! -d $folder ]; then
	echo $folder does not exist
	exit 1
fi

currDir=`pwd`

cd ${folder}

rm -f *_${target}_sub_embed_*.pdb

stats=${target}.noe_statis
cp /dev/null $stats

for i in *_${target}_*.pdb
do
	bname=`basename $i`
	score=`grep "REMARK noe" $i | cut -f2 -d'=' `
	if [ -z $score ]; then
		echo $i does not contain NOE score. Skip it
		continue
	fi
	echo $bname $score >>  $stats
done

sort -k2,2 -g $stats > ${stats}_sort

top5models=`head -5 ${stats}_sort | cut -f1 -d' '`

i=1
for model in $top5models
do
	cp $model ${target}_model${i}.pdb
	i=`expr $i + 1`
done

cd $currDir
