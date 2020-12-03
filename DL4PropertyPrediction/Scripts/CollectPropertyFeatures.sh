#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "proteinName rootDir [resultDir [MSA_mode] ]"
	echo "	rootDir: a folder containing information for property prediction, e.g., 1pazA_OUT/. It shall contain subfolders proteinName_contact/ and proteinName_thread/, which in turn shall contain .a3m or .hhm files"
	echo "	resultDir: a folder for the resultant feature file named as proteinName.propertyFeatures.pkl, default RootDir/PropertyPred/"
	echo "	MSA_mode: which MSAs to be used, all (default), uce3, ure5 or user?"
	exit 1
fi

target=$1

RootDir=$2
if [ ! -d $RootDir ]; then
	echo "ERROR: the root folder for feature files does not exist: $RootDir"
	exit 1
fi

resDir=${RootDir}/PropertyPred/
if [ $# -ge 3 ]; then
	resDir=$3
fi
if [ ! -d $resDir ]; then
	mkdir -p $resDir
fi

MSA="all"
if [ $# -ge 4 ]; then
	MSA=$4
fi

if [ -z "$DL4PropertyPredHome" ]; then
        echo "ERROR: please set environmental variable DL4PropertyPredHome to the instllation directory of DL4PropertyPrediction"
        exit 1
fi

if [ -z "$DistFeatureHome" ]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the instllation directory of BuildFeatures"
        exit 1
fi

a3mFile=${RootDir}/${target}_thread/${target}.a3m
if [ ! -f $a3mFile ]; then
	echo "ERROR: the a3m file does not exist: $a3mFile"
	exit 1
fi

hhmFile=${RootDir}/${target}_thread/${target}.hhm
if [ ! -f $hhmFile ]; then
	$DistFeatureHome/util/hhmake -i $a3mFile -o $hhmFile
	if [ $? -ne 0 ]; then
		echo "ERROR: failed to run $DistFeatureHome/util/hhmake -i $a3mFile -o $hhmFile"
		exit 1
	fi
fi

if [ "$MSA" == "all" ]; then
	methods="uce3 ure5 user"
else
	methods=$MSA
fi

hhmFiles=$hhmFile
for method in $methods
do
	featureFolder=${RootDir}/${target}_contact/feat_${target}_${method}
	hhmf=$featureFolder/${target}.hhm
	if [ -f $hhmf ]; then
		hhmFiles=$hhmFiles" "$hhmf
	else
		a3mf=$featureFolder/${target}.a3m
		if [ -f $a3mf ]; then
			$DistFeatureHome/util/hhmake -i $a3mf -o $hhmf
			hhmFiles=$hhmFiles" "$hhmf
		fi
	fi
done

python $DL4PropertyPredHome/GenPropertyFeaturesFromMultiHHMs.py ${target} $hhmFiles $resDir
if [ $? -ne 0 ]; then
	echo "ERROR: failed to collect property feature files for $target !"
	exit 1
fi

inputFeature=${resDir}/${target}.propertyFeatures.pkl
if [ ! -f $inputFeature ]; then
        echo "ERROR: failed to generate the basic property feature file: " $inputFeature
        exit 1
fi

echo "The input features for property prediction are saved in $inputFeature"
