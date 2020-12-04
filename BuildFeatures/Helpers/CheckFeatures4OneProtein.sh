#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 proteinName OutFolder
	echo "	This scripts check the distance/orientation features for one protein"
	echo "	OutFolder: it shall contain a subfolder XXX_contact, which in turn shall contain subfolders XXX_YYY and feat_XXX_YYY where YYY is a MSA generation method"
	echo "	This script will output *proteinName invalid* as the message if some features are not correctly generated"
	exit 1
fi

protein=$1

OutFolder=$2
if [ ! -d $OutFolder ]; then
	echo "ERROR: invalid out folder $OutFolder"
	exit 1
fi

featureFolder=$OutFolder/${protein}_contact
if [ ! -d $featureFolder ]; then
	echo $protein invalid $featureFolder
	exit 1
fi

for a3mfolder in $featureFolder/${protein}_*
do
	bname=`basename $a3mfolder`
	ffolder=$featureFolder/feat_${bname}
	if [ ! -d $ffolder ]; then
		echo "$protein invalid $ffolder"
		continue
	fi

	inputFeatureFile=$ffolder/${protein}.inputFeatures.pkl
	extraCCMFile=$ffolder/${protein}.extraCCM.pkl
	a2mFile=$ffolder/${protein}.a2m

	if [ ! -f $inputFeatureFile -o ! -f $extraCCMFile -o ! -f $a2mFile ]; then
		echo "$protein invalid at least one of $inputFeatureFile $extraCCMFile $a2mFile"
		continue
	fi
done
