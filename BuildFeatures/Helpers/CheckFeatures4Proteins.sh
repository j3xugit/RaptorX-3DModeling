#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 proteinList MetaFolder
	echo "	This scripts check the distance/orientation features for a list of proteins"
	echo "	MetaFolder: the folder contains a list of subfolders XXX_OUT where XXX represents one protein"
	echo "	    XXX_OUT shall contain XXX_contact, which in turn shall contain subfolders XXX_YYY and feat_XXX_YYY where YYY is a MSA generation method"
	echo "	This script will output a list of proteins for which some feature files are not correctly generated"
	exit 1
fi

listFile=$1
if [ ! -f $listFile ]; then
	echo "ERROR: invalid protein list file $listFile"
	exit 1
fi

metaFolder=$2
if [ ! -d $metaFolder ]; then
	echo "ERROR: invalid meta folder $metaFolder"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for protein in `cat $listFile`
do
	featureFolder=$metaFolder/${protein}_OUT/
	if [ ! -d $featureFolder ]; then
		echo $protein invalid $featureFolder
		continue
	fi
	$cmdDir/CheckFeatures4OneProtein.sh $protein $featureFolder
done
