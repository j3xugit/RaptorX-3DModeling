#!/bin/sh

if [ $# -lt 3 ]; then
	echo "$0 [-s] targetListFile SeqDir ModelMetaFolder [savefolder]"
	echo "-s: when specified, select models by energy"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

selectModel=false
while getopts ":s" opt; do
       case ${opt} in
        s )
          selectModel=true
          ;;
        \? )
          echo "Invalid Option: -$OPTARG" 1>&2
                 exit 1
                  ;;
        : )
          echo "Invalid Option: -$OPTARG requires an argument" 1>&2
          exit 1
                  ;;
       esac
done
shift $((OPTIND -1))

if [ $# -lt 3 ]; then
       Usage
       exit 1
fi

targets=`cat $1`
SeqDir=$2
ModelMetaDir=$3

savefolder=`pwd`
if [ $# -ge 4 ]; then
	savefolder=$4
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

for target in $targets
do
	echo "Clustering decoys for $target by Spicker..."
	decoyFolder=$ModelMetaDir/${target}-RelaxResults
	if [ ! -d $decoyFolder ]; then
		decoyFolder=$ModelMetaDir/${target}-FoldResults
	fi

	if [ ! -d $decoyFolder ]; then
		echo "WARNING: cannot find decoy folder for $target "
		continue
	fi

	if [ "$selectModel" = true ]; then
		$cmdDir/SpickerOneTarget.sh -s $SeqDir/${target}.fasta $decoyFolder $savefolder/${target}-SpickerResults &
	else
		$cmdDir/SpickerOneTarget.sh $SeqDir/${target}.fasta $decoyFolder $savefolder/${target}-SpickerResults &
	fi
done
