#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "targetName infolder [numJobsPerModel (5, 10 or 20, default 20) ]"
	echo "  infolder (e.g., 1pazA_OUT/) contains subfolders for predicted distance and property information"
	echo "  results are also saved to infolder"
	exit 1
fi

target=$1
rootDir=$2

numJobsPerModel=20
if [ $# -ge 3 ]; then
	numJobsPerModel=$3
fi

if [ $numJobsPerModel -lt 5 -o $numJobsPerModel -gt 20 ] ; then
        echo "incorrect value for numJobsPerModel: " $numJobsPerModel
        exit 1
fi

if [ -z "${DistanceFoldingHome}" ]; then
	echo "Please set the environmental variable DistanceFoldingHome to the installation folder of the Folding module, e.g., $HOME/3DModeling/Folding/ "
	exit 1
fi


numJobs=`expr 40 / $numJobsPerModel `
Wrapper=$DistanceFoldingHome/scripts/FoldingWrapper${numJobsPerModel}.sh

#listFile=`mktemp $target.list.XXXXXX`
#echo $target > $listFile

currDir=`pwd`
cd $rootDir

#SeqDir=$rootDir
#PredDistDir=$rootDir/DistancePred/Dist_Server_EC52CSet10820N11410Atom18Models/
#PredPhiPsiSSDir=$rootDir/PropertyPred/PhiPsiSS8_Server_SeqSet10820Models/

SeqDir=./
PredDistDir=DistancePred/Dist_Server_EC25CSet10820N11410Atom15Models/
PredPhiPsiSSDir=PropertyPred/PhiPsiSS8_Server_SeqSet10820Models/

#$DistanceFoldingHome/scripts/DistributeJobs.sh $numJobs $listFile $Wrapper $SeqDir $PredDistDir $PredPhiPsiSSDir All
$Wrapper $target $SeqDir $PredDistDir $PredPhiPsiSSDir All 3 15

if [ $? -ne 0 ]; then
	echo "Errors happen in the folding process for $target "
	exit 1
fi
