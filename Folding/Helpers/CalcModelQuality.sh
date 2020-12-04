#!/bin/bash

function Usage 
{
	echo "$0 [ -f filter4decoyName | -d savefolder ] ModelFolder NativeFile"
	echo "	This script calculates the quality of 3D models of a one protein. Please make sure DeepScore is searchable by your bin path"
	echo "	ModelFolder: the folder for predicted 3D models. All model files shall end with .pdb"
	echo "	NativeFile: the experimental structure file. It shall end with .pdb"
	echo "	-f: a filter string used to select a subset of decoys in ModelFolder, default empty, i.e., all files ending with .pdb will be selected"
	echo "	-d: the folder for the resultant model quality file, default empty. When empty, the result file will be saved to ModelFolder, i.e., ModelFolder shall be writable"
	echo "	    a proteinDomainName-quality.txt.sorted file will be generated containing decoy file name, native length, RMSD, TM, MaxSub, GDT, GHA"
}

savefolder=''
filter=''

while getopts ":f:d:" opt; do
        case ${opt} in
                f )
                  filter=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
		  if [ ! -d $savefolder ]; then
			mkdir -p $savefolder
		  fi
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

if [ $# -lt 2 ]; then
	Usage
	exit 1
fi

modelDir=$1
if [ ! -d $modelDir ]; then
	echo "ERROR: invalid decoy folder $modelDir"
	exit 1
fi

nativeFile=$2
if [ ! -f $nativeFile ]; then
	echo "ERROR: invalid experimental structure file $nativeFile"
	exit 1
fi
target=`basename $nativeFile .pdb`

if [ -z "${savefolder}" ]; then
	if [ ! -w $modelDir ]; then
		echo "ERROR: you are trying to write the result to $modelDir, but it is not writable"
		exit 1
	fi
	qualityFile=${modelDir}/${target}-quality.txt
else
	qualityFile=${savefolder}/${target}-quality.txt
fi

cp /dev/null $qualityFile

PID=$?
OUT=$(mktemp /tmp/${target}-${PID}.XXXXXXXXXX) || { echo "Failed to create temp file"; exit 1; }

for i in ${modelDir}/*${filter}*pdb 
do 
	if [ ! -f $i ]; then
		continue
	fi
	DeepScore $i ${nativeFile} -n -2 > $OUT
	quality=`grep -A6 "1st input protein" $OUT | tail -1 | cut -f5- -d' '`
	b=`basename $i .pdb`
	echo $b $quality >> $qualityFile
done

if [ -s "$qualityFile" ]; then
	cut -f1,6- -d' ' $qualityFile | sort -k4,4 -rn > ${qualityFile}.sorted
fi

rm -f $qualityFile
rm -f $OUT
