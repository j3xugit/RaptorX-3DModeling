#!/bin/bash


if [ $# -lt 2 ]; then
	echo $0 "proteinName rootDir [resultDir]"
	echo "   rootDir: the folder containing files for property prediction, e.g., 1pazA_OUT/. This folder shall contain subfolders proteinName_contact/ and proteinName_thread/, which in turn shall contain .tgt files for input features"
	echo "   resultDir: the result folder for generated feature file named as proteinName.propertyFeatures.pkl (default RootDir/PropertyPred/)"
	##echo "   MSA: specify file which MSAs are used to derive features (default all).
	exit 1
fi

target=$1
RootDir=$2
resDir=${RootDir}/PropertyPred/

if [ $# -ge 3 ]; then
	resDir=$3
fi

MSA="all"
if [ $# -ge 4 ]; then
	MSA=$4
fi

if [ -z "${DL4PropertyPredHome}" ]; then
        echo "Please set environmental variable DL4PropertyPredHome, e.g., $HOME/3DModeling/DL4PropertyPrediction/ "
        exit 1
fi

tgtFile=${RootDir}/${target}_thread/${target}.tgt
if [ ! -f $tgtFile1 ]; then
	echo "ERROR: the threading tgt file $tgtFile does not exist"
	exit 1
fi

tgtFiles=$tgtFile

if [ "$MSA" == "all" ]; then

	for method in uce3 ure5 user
	do
		featureFolder=${RootDir}/${target}_contact/feat_${target}_${method}
		tgtf=$featureFolder/${target}.tgt
		if [ -f $tgtf ]; then
			tgtFiles=$tgtFiles" "$tgtf
		fi
	done

	#echo "tgtFiles= $tgtFiles"
	python $DL4PropertyPredHome/GenPropertyFeaturesFromMultiTGTs.py ${target} $tgtFiles

	if [ $? -ne 0 ]; then
		echo "ERROR: Failed to collect input feature files for property prediction for $target !"
		exit 1
	fi
fi

inputFeature=${target}.propertyFeatures.pkl
if [ ! -f $inputFeature ]; then
        echo "ERROR: Failed to generate the basic property feature file: " $inputFeature
        exit 1
fi

if [ ! ${resDir} -ef ./ ]; then
	mkdir -p ${resDir}
	mv $inputFeature ${resDir}/
	if [ $? -ne 0 ]; then
		echo "ERROR: Failed to save the input feature file for property prediction to ${resDir}"
		exit 1
	fi

fi


echo "The input features for property prediction are saved in ${resDir}/$inputFeature"
