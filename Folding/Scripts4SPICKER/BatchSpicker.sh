#!/bin/bash

if [ $# -lt 3 ]; then
	echo "$0 [-s] targetListFile SeqDir ModelMetaFolder [savefolder]"
	echo "	This script runs SPICKER on a list of proteins, each has a set of decoys"
	echo "	-s: when specified, select models by potential, otherwise not"
	echo "	SeqDir: a folder for all sequence files, each in FASTA format"
	echo "	ModelMetaFolder: a meta folder containing a list of subfolders, each subfolder has the decoys of one protein"
	echo "	savefolder: a folder for result saving. One subfolder in this folder will be created for each protein"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

selectModel=0
while getopts ":s" opt; do
       case ${opt} in
        s )
          selectModel=1
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

	options=" -d $savefolder/${target}-SpickerResults "
	if [ $selectModel -eq 1 ]; then
		options=$options" -s "
	fi

	$cmdDir/SpickerOneTarget.sh $options $SeqDir/${target}.fasta $decoyFolder &
done

wait
