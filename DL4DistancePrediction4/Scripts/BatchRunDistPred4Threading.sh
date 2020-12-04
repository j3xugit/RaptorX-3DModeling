#!/bin/bash

#modelName=EC47C31C16CL99PDB30V2020CbCbTwoRModels
modelName=EC47C31C16CL99PDB30CbCbTwoRV7C

if [ $# -lt 2 ]; then
	echo $0 proteinListFile inputFeatureFolder [gpu]
	exit 1
fi

proteinListFile=$1
inputFolder=$2

if [[ "$proteinListFile" != *.list ]]; then
	echo "ERROR: the protein list file shall end with .list"
	exit 1
fi

GPU=cuda1
if [ $# -ge 3 ]; then
	GPU=$3
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

$cmdDir/PredictPairRelation4MultiInputs.sh -m $modelName -g $GPU $proteinListFile $inputFolder

inputFolderBaseName=`basename $inputFolder`
ListBaseName=`basename $proteinListFile`
ListBaseName=${ListBaseName%.*}
resultDir=LocalResults/Dist_${ListBaseName}-${inputFolderBaseName}_${modelName}

reduceProgram=$DL4DistancePredHome/ReducePredictedDistMatrix.py

numAllowedJobs=15
keywords=`basename $reduceProgram`
myself=`basename $0 `

for target in `cat $proteinListFile `
do
        while true
        do
                ## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | wc -l`
                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
                        #$program $SeqDir/$target.fasta $ResDir $MSAmode &
			python $reduceProgram -s $resultDir $resultDir/${target}.predictedDistMatrix.pkl & 
			python $reduceProgram -r CbCb_Discrete12C -s $resultDir $resultDir/${target}.predictedDistMatrix.pkl &
                        break
                else
                        a=`expr $RANDOM % 3 `
                        sleep $a
                fi

        done

done

wait

