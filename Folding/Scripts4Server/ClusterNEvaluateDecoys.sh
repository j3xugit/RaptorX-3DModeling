#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 "targetName (e.g. T0859-D1) [ DataName (e.g., CASP12DM) ]"
	exit -1
fi

target=$1
DataName=CASP12DM

if [ $# -ge 2 ]; then
	DataName=$2
fi

DataName2=$DataName
if [ $DataName == "PDB10820Test1718" ]; then
	DataName2=PDB10820TestData
fi



DecoyDir1=ModelingResult-Dist_${DataName}_EC52CSet10820Atom8Models-PhiPsiSS8_${DataName2}_SeqSet10820/
DecoyDir2=ModelingResult-Dist_${DataName}_EC52CSet11410Atom10Models-PhiPsiSS8_${DataName2}_SeqSet10820/

DecoyDir3=ModelingResult-Dist_${DataName}_EC25CSet10820Atom8Models-PhiPsiSS8_${DataName2}_SeqSet10820/
DecoyDir4=ModelingResult-Dist_${DataName}_EC25CSet11410Atom7Models-PhiPsiSS8_${DataName2}_SeqSet10820/

DecoyDir5=ModelingResult-Dist_${DataName}_EC52CSet11410MEAtom9Models-PhiPsiSS8_${DataName2}_SeqSet10820/
DecoyDir6=ModelingResult-Dist_${DataName}_EC52CSet10820MEAtom8Models-PhiPsiSS8_${DataName2}_SeqSet10820/

DecoyDir7=ModelingResult-Dist_${DataName}_EC25CSet10820N11410Atom15Models-PhiPsiSS8_${DataName2}_SeqSet10820/
DecoyDir8=ModelingResult-Dist_${DataName}_EC52CSet10820N11410Atom18Models-PhiPsiSS8_${DataName2}_SeqSet10820/ 


Method1=CbCb_t3_c15_A5.0
Method2=All_t3_c15_A5.0
Method3=CaCa+CbCb+CgCg_t3_c15_A5.0

## collect models
cp /dev/null $target.models.list
#ls ${DecoyDir1}/${target}_All_t3_c15_A5.0/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir1}/${target}_CaCa+CbCb+CgCg_t2_c15_A5.0/${target}_model?.pdb >> $target.models.list
ls ${DecoyDir1}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list

#ls ${DecoyDir2}/${target}_All_t3_c15_A5.0/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir2}/${target}_All_t-1_c3_A5.0/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir2}/${target}_CaCa+CbCb+CgCg_t3_c15_A5.0/${target}_model?.pdb >> $target.models.list
ls ${DecoyDir2}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list

#ls ${DecoyDir3}/${target}_All_t3_c15_A5.0/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir3}/${target}_CaCa+CbCb+CgCg_t3_c15_A5.0/${target}_model?.pdb >> $target.models.list
ls ${DecoyDir3}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir3}/${target}_All_t3_c11_A5.0/${target}_model?.pdb >> $target.models.list

#ls ${DecoyDir4}/${target}_All_t3_c15_A5.0/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir4}/${target}_CaCa+CbCb+CgCg_t2_c15_A5.0/${target}_model?.pdb >> $target.models.list
ls ${DecoyDir4}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir4}/${target}_All_t3_c13_A5.0/${target}_model?.pdb >> $target.models.list

#ls ${DecoyDir5}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir6}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list
ls ${DecoyDir7}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list
#ls ${DecoyDir8}/${target}_${Method2}/${target}_model?.pdb >> $target.models.list

native=${DataName2}-PDB/${target}.pdb
cp /dev/null $target.deepscore

## run DeepScore to calculate model quality
for m in `cat $target.models.list`
do
	scores=`DeepScore $m $native | python scripts/CollectModelQuality.py `
	echo $scores $m >> $target.deepscore
done

## run maxcluster to do consensus analysis
maxcluster64bit -l $target.models.list  > $target.maxcluster

## evaluate the consensus analysis result
python scripts/AnalyzeMaxCluster.py $target.maxcluster $target.deepscore | tee $target.consensus

## clean up
rm -f $target.deepscore $target.maxcluster $target.models.list

