#!/bin/bash

nativeDir=CASP13DM-Native
mappingFile=CASP13FM2MyDM.txt

function Usage 
{
	echo $0 "[ -f filter4decoyFiles | -d savefolder ] proteinListFile/domainMappingFile decoyFolder [NativeFolder]"
	echo "	This script calculates the quality (TM, GDT) of decoys of a list of proteins. Please make sure DeepScore is searchable by your bin path"
	echo "	proteinListFile/domainMappingFile: a file for a list of proteins or a list of domain mapping, e.g., $mappingFile"
	echo "	     when it is a protein list file, it shall have only one column, i.e., each row has one protein name"
	echo "	     when it is a protein domain mapping file, it shall have three columns separated by :. The columns are: native domain name, domain boundary and target name used in model building"
	echo "	decoyFolder: a folder containing a set of subfolders. Each subfolder, with name like targetName-*Results, contains 3D models of one protein"
	echo "	NativeDir: the folder for experimentally-solved stucture files, default $nativeDir"
	echo "	     all decoy files and native files shall end with .pdb"
	echo "	-f: a filter string used to select a subset of decoys. Only decoy file name containing this filter will be evaluated; default empty, i.e., all decoys will be evaluated"
	echo "	-d: the folder for result saving, default empty. When empty, the results will be saved to decoyFolder, i.e., decoyFolder shall be writable "
	echo "	     for each protein or domain, a proteinDomainName-quality.txt.sorted file will be generated containing decoy file name, native length, RMSD, TM, MaxSub, GDT, GHA"
}

filter=''
savefolder=''

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

mappingFile=$1
if [ ! -f $mappingFile ]; then
	echo "ERROR: invalid protein list file or domain mapping file $mappingFile"
	exit 1
fi

metaFolder=$2
if [ ! -d $metaFolder ]; then
	echo "ERROR: invalid folder for decoys/models to be evaluated"
	exit 1
fi

if [ $# -ge 3 ]; then
	nativeDir=$3
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/CalcModelQuality.sh

if [ ! -x $program ]; then
	echo "ERROR: invalid executable program $program"
	exit 1
fi

## count the min number of columns in mappingFile
minNumCols=` awk 'BEGIN {FS = ":"} ; {print NF}' $mappingFile | sort -nu | head -1 `
maxNumCols=` awk 'BEGIN {FS = ":"} ; {print NF}' $mappingFile | sort -nu | tail -n 1 `
#echo $minNumCols $maxNumCols
if [ $minNumCols -ne $maxNumCols ]; then
	echo "ERROR: all rows in the protein list file or domain mapping file shall have the same number of columns"
	exit 1
fi

numCols=$maxNumCols
if [ $numCols -ne 1 -a $numCols -ne 3 ]; then
	echo "ERROR: the protein list file or domain mapping file shall have either 1 or 3 columns"
	exit 1
fi

while read p; do
  	echo "$p"
	if [ $numCols -eq 3 ]; then
		DMname=`echo $p | cut -f1 -d':'`
		MyDMname=`echo $p | cut -f3 -d':'`
	else
		DMname=$p
		MyDMname=$p
	fi

	for decoyFolder in $metaFolder/${MyDMname}-*Results/
	do
		if [ ! -d $decoyFolder ]; then
			continue
		fi

		if [ -z "$filter" ]; then
			if [ -z "$savefolder" ]; then
				$program $decoyFolder $nativeDir/$DMname.pdb &
			else
				$program -d $savefolder $decoyFolder $nativeDir/$DMname.pdb &
			fi
		else
			if [ -z "$savefolder" ]; then
				$program -f $filter $decoyFolder $nativeDir/$DMname.pdb &
			else
				$program -f $filter -d $savefolder $decoyFolder $nativeDir/$DMname.pdb &
			fi
		fi
	done

done <${mappingFile}

wait
