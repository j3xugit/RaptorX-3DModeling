#!/bin/bash

alpha=1.61
savefolder=`pwd`

function Usage 
{
	echo $0 "[-a alpha | -d savefolder ] MyDMName MetaFolder cstFolder [propertyFolder]"
	echo "	This script scores unrelaxed models of one protein"
	echo "	The models may be dispersed in several folders, generated from different constraints"
	echo "	-a: alpha used in DFIRE potential, default $alpha. Not used if cstFolder is the folder for Rosetta constraint files"
	echo "		used only when cstFolder is the folder for predicted dist/orientation"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	When propertyFolder is empty, then cstFolder shall be folder for constraints. Otherwise cstFolder and propertyFolder shall be folders for predicted distance and Phi/Psi angles"
}

while getopts ":a:d:" opt; do
        case ${opt} in
                a )
                  alpha=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
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

if [ $# -lt 3 -o $# -gt 4 ]; then
	Usage
	exit 1
fi

target=$1
metaFolder=$2
if [ ! -d $metaFolder ]; then
	echo "ERROR: invalid folder for 3D models: $metaFolder"
	exit 1
fi

cstFolder=$3
if [ ! -d $cstFolder ]; then
	echo "ERROR: invalid folder for predicted distance information or Rosetta distance constraints: $cstFolder "
	exit 1
fi

propertyFolder='cst'
if [ $# -eq 4 ]; then
	propertyFolder=$4
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/ScoreModelsInOneFolder.sh

if [ ! -x $program ]; then
	echo "ERROR: non-existing or non-excetuable program $program"
	exit 1
fi

for modelFolder in $metaFolder/${target}-*FoldResults
do
	$program -d $savefolder $modelFolder $cstFolder $propertyFolder &
done
wait

