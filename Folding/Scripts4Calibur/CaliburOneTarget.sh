#!/bin/sh

if [ $# -lt 1 ]; then
	echo $0 modelFolder [savefolder]
	exit 1
fi

modelFolder=$1
bname=`basename $modelFolder`

savefolder=`pwd`
if [ $# -ge 2 ]; then
	savefolder=$2
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

saveFile=$savefolder/$bname.calibur

listFile=$(mktemp $savefolder/$bname.models.XXXXXX)
ls $modelFolder/*.pdb > $listFile
calibur $listFile > $saveFile
rm -f $listFile

i=0
centers=`grep Centroid $saveFile | cut -f2 -d' '`
for center in $centers
do
	cfile=`readlink -f $center `
	ln -s $cfile $savefolder/center${i}.pdb
	i=` expr $i + 1 `
done
