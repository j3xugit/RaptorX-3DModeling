#!/bin/sh

if [ $# -lt 3 ]; then
	echo $0 targetListFile ModelMetaFolder [savefolder]
	exit 1
fi

targets=`cat $1`
ModelMetaDir=$2

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

savefolder=`pwd`
if [ $# -ge 3 ]; then
	savefolder=$3
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

for target in $targets
do
	echo "Clustering decoys for $target by calibur..."
	$cmdDir/CaliburOneTarget.sh $ModelMetaDir/${target}-RelaxResults $savefolder/${target}-CaliburResults &
done
