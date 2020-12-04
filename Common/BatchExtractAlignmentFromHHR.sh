#!/bin/bash

E="0.00001,0.001,0.1"
MaxSeqID=0.4
maxNumTemplates=15

HHRDir=HHR_2017MSA2PDB70Apr2020
HHMDir=MSA_2017_E001
PDBClusterFile=/mnt/data/RaptorXCommon/PDB-Clusters/bc-40.out
templates=TPL-PDB70_29Apr2020.list

savefolder=''
groupFile=''

function Usage
{
	echo $0 "[-E evalue | -I MaxSeqID | -d savefolder | -s savefile | -c PDBClusterFile | -t templateList | -n maxNumTemplates ] targetListFile HHMDir HHRDir"
	echo "	This script extracts pairwise alignments from .hhr files, mainly for training"
	echo "	targetListFile: a file for a list of query proteins, each in one row"
	echo "	HHRDir: the folder for .hhr files, e.g., $HHRDir. Eah file in this folder shall end with .pdb70.hhr"
	echo "	HHMDir: the folder for .hhm files of query proteins, e.g., $HHMDir. Each file in this folder shall end with .hhm"
	echo "	This script will generate one group.txt file and also one folder containing all pairwise alignments"
	echo " "
	echo "	-d: the folder for result saving. The default folder name depends on Evalue and MaxSeqID"
	echo "	-s: the .group.txt file for result saving. By default, this file name consists of Evalue, MaxSeqID and basename of targetListFile"
	echo "	-E: the list of Evalues for cutoff, default $E"
	echo "	-I: the maximum allowed seqID (scale 0-1) between query and template, default $MaxSeqID"
	echo "	-c: the file for removing redundant templates, e.g., bc-40.out downloaded from PDB website, default $PDBClusterFile"
	echo "	-t: the file for a list of allowed templates, default $templates"
	echo "	-n: the maximum number of selected templates for each query sequence, default $maxNumTemplates"
}

while getopts ":E:I:d:c:t:s:n:" opt; do
        case ${opt} in
                E )
                  E=$OPTARG
                  ;;
                I )
                  MaxSeqID=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
                  ;;
		c )
		  PDBClusterFile=$OPTARG
		  ;;
		t )
		  templates=$OPTARG
		  ;;
		s )
		  groupFile=$OPTARG
		  ;;
		n )
		  maxNumTemplates=$OPTARG
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

targets=$1
if [ ! -f $targets ]; then
	echo "ERROR: invalid protein list file $targets"
	exit 1
fi

if [ ! -d $HHRDir ]; then
	echo "ERROR: invalid folder for the .hhr files $HHRDir"
	exit 1
fi

if [ ! -d $HHMDir ]; then
	echo "ERROR: invalid folder for the .hhm files $HHMDir"
	exit 1
fi

if [ ! -f $PDBClusterFile ]; then
	echo "ERROR: invalid PDB cluster file $PDBClusterFile"
	exit 1
fi

if [ ! -f $templates ]; then
	echo "ERROR: invalid file for template list $templates"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/ExtractAlignmentFromHHR.py
if [ ! -f $program ]; then
	echo "ERROR: invalid program $program"
	exit 1
fi

if [ -z "$savefolder" ]; then
	savefolder=Align4S35ToNewPDB70-ID${MaxSeqID}-E${E}/
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

if [ -z "$groupFile" ]; then
	groupFile=group-S35-NewPDB70-E${E}ID${MaxSeqID}-`basename $targets`
fi
cp /dev/null $groupFile

for target in `cat $targets`
do
	#echo $target
	python $program $HHRDir/${target}.pdb70.hhr $HHMDir/${target}.hhm -E $E -I $MaxSeqID -c $PDBClusterFile -t $templates -d $savefolder -n $maxNumTemplates >> $groupFile
done
