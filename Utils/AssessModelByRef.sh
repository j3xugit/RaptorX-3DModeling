#!/bin/bash

if [ $# -lt 2 ]; then
        echo $0 "targetModel modelFolder1 modelFolder2 ..."
	echo "This script assesses the quality of a model by comparing it to some models in provided folders"
        echo "	targetModel: the model for which you would like to predict local and global quality"
        echo "	modelFolder: a folder for reference models. Multiple folders can be used"
	echo "	A new PDB file will be generated and saved in the same folder as the targetModel, so please make sure that the folder of the targetModel is writable"
        exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

targetModel=$1
if [ ! -f $targetModel ]; then
        echo "ERROR: invalid target model file" $targetModel
        exit 1
fi

## needs revision later
ResDir=`dirname $targetModel`

modelFile=`basename $targetModel .pdb`
refListFile=$(mktemp -t RXrefList4${modelFile}-XXXXXXXXXX)

modelFolders=${@:2}
for modelFolder in $modelFolders
do
	ls $modelFolder/*.pdb | shuf | head -40 >> $refListFile
done

if [ ! -f $refListFile ]; then
        echo "ERROR: failed to generate reference model file for $targetModel"
        exit 1
fi

qualityFile=$ResDir/Quality-${modelFile}.txt
CAqualityFile=$ResDir/CA.Quality-${modelFile}.txt
python ${cmdDir}/EstimateModelError.py -o $qualityFile $targetModel $refListFile > /dev/null

rm -f $refListFile

grep -v REMARK $qualityFile | grep -w CA | awk '{print $NF}' > $CAqualityFile

## revise the model file
$cmdDir/AddErrorEstimate2PDB $targetModel $CAqualityFile $ResDir/$modelFile.quality.pdb

