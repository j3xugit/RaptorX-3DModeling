#!/bin/bash


if [ $# -lt 2 ]; then
	echo $0 "targetName infolder [ ResDir [MSA-options] ]"
	echo "   infolder (e.g., 1pazA_OUT/) shall contains feature folders such as target_contact/feat_targetName_uce0, target_contact/feat_targetName_uer3"
	echo "   possible MSA options, currently only the default value is supported"
	echo "   result files are saved to $ResDir or infolder/DistancePred/"
	exit 1
fi

if [[ -z "${DL4DistancePredHome}" ]]; then
	echo "Please set environmental variable DL4DistancePredHome, e.g., $HOME/3DModeling/DL4DistancePrediction3/ "
	exit 1
fi

## collect features
proteinName=$1
rootDir=$2

resDir=${rootDir}/DistancePred/

if [ $# -ge 3 ]; then
	resDir=$3
fi

MSA="all"
if [ $# -ge 4 ]; then
	MSA=$4
fi

#echo "Use ${MSA} MSAs ..."


folder1=${rootDir}/${proteinName}_contact/feat_${proteinName}_uce0/
folder2=${rootDir}/${proteinName}_contact/feat_${proteinName}_uce3/
folder3=${rootDir}/${proteinName}_contact/feat_${proteinName}_ure3/
folder4=${rootDir}/${proteinName}_contact/feat_${proteinName}_ure5/

folder5=${rootDir}/${proteinName}_contact/feat_${proteinName}_user/

if [ "$MSA" == "all" ]; then
	folders=""
	for i in $folder1 $folder2 $folder3 $folder4 $folder5
	do
		if [ -d $i ]; then
			folders=$folders" "$i
		fi
	done
	#echo "folders= $folders"

	python $DL4DistancePredHome/ReadOneProteinFeatures.py $proteinName $folders

	#if [ ! -d $folder5 ]; then
	#	python $DL4DistancePredHome/ReadOneProteinFeatures.py $proteinName $folder1 $folder2 $folder3 $folder4
	#else
	#	python $DL4DistancePredHome/ReadOneProteinFeatures.py $proteinName $folder5
	#fi

	if [ $? -ne 0 ]; then
                echo "Failed to collect input feature files for distance prediction!"
                exit 1
        fi

else
	echo "unsupported MSA option" $MSA
	exit 1
fi


inputFeature=${proteinName}.distanceFeatures.pkl

if [ ! -f $inputFeature ]; then
	echo "Failed to generate an input feature file for " proteinName
	exit 1
fi

if [ ! ${resDir} -ef ./ ]; then
	mkdir -p ${resDir}
	mv $inputFeature ${resDir}/
	if [ $? -ne 0 ]; then
        	echo "Failed to save the input feature file for distance prediction to ${resDir}"
        	exit 1
	fi

fi

echo "The input feature file is saved to" ${resDir}
