#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 MyDMname MetaFolder NativeFile
	exit 1
fi

target=$1
metaFolder=$2
native=$3

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for modelFolder in $metaFolder/${target}-*Results
do
	$cmdDir/CalcModelQuality.sh $modelFolder $native &
done
wait

