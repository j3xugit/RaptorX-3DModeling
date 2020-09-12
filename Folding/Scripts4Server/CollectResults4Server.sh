#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "target target_OUT"
	exit 1
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
	echo "please set the environmental variable DistanceFoldingHome to the installation folder of the Folding module, e.g., $HOME/3DModeling/Folding"
	exit 1
fi

target=$1
SrcDir=$2
ResDir=${target}_Results
mkdir -p $ResDir

cp $DistanceFoldingHome/0README_overall $ResDir/

## copy sequence files
cp $SrcDir/${target}.seq $ResDir/
cp $SrcDir/${target}.seq $ResDir/${target}.fasta

##copy MSA and ccmpred
srcdir=${SrcDir}/${target}_contact/feat_${target}_user

if [ ! -d $srcdir ]; then
	srcdir=$SrcDir/${target}_contact/feat_${target}_uce3
fi

cp $srcdir/${target}.ss3 $ResDir/

## for display purpose. Since drawing a large MSA is slow, here we display up to 5000 sequences
head -2  $srcdir/${target}.a3m > ${target}.a3m.head
seqLen=`tail -1 ${target}.a3m.head | wc -c | cut -f1 -d' ' `
tail -n +3 $srcdir/${target}.a3m | paste -sd '*\n' | shuf -n 5000 | sed 's/*/\n/g' >> ${target}.a3m.head
cut -c1-$seqLen ${target}.a3m.head > $ResDir/${target}.a3m
rm -f ${target}.a3m.head

for i in user uce0 uce3 ure3 ure5
do
	srcdir=${SrcDir}/${target}_contact/feat_${target}_$i
	if [ -d $srcdir ]; then
		gzip < $srcdir/${target}.a3m > $ResDir/${target}.a3m_${i}.gz
	fi
done

## copy predicted distance information
srcdir=${SrcDir}/DistancePred
cp $srcdir/${target}.distanceFeatures.pkl $ResDir/

srcdir=${srcdir}/Dist_Server_EC25CSet10820N11410Atom15Models/
cp $srcdir/${target}*.png $ResDir/
cp $srcdir/${target}.gcnn $ResDir/
cp $srcdir/${target}.CASP.rr $ResDir/${target}.rr
cp $srcdir/${target}.CASP.rr $ResDir/contactmap.txt
cp $srcdir/${target}.CASP.rr $ResDir/${target}.contactmap.txt
cp $srcdir/${target}.bound.txt $ResDir/
cp $srcdir/${target}.correctedDistMatrix412C.pkl $ResDir/

## add code here to copy predicted properties

## copy predicted 3D models
srcdir=${SrcDir}/FoldingResult-Dist_Server_EC25CSet10820N11410Atom15Models-PhiPsiSS8_Server_SeqSet10820Models/${target}_All_t3_c15_A5.0/
if [ -d $srcdir ]; then
	#cp $srcdir/${target}_model1.pdb $ResDir/${target}.pdb
	cp /dev/null $ResDir/${target}.pdb
	modelDir=$ResDir/${target}_models
	mkdir -p $modelDir
	for i in 1 2 3 4 5
	do
		cp $srcdir/${target}_model${i}.pdb $modelDir/${target}_model_${i}.pdb
		cp $srcdir/${target}_model${i}.localQuality.txt $modelDir/
		echo "MODEL $i" >> $ResDir/${target}.pdb
		grep ^ATOM $srcdir/${target}_model${i}.pdb >> $ResDir/${target}.pdb
		echo "ENDMDL" >>  $ResDir/${target}.pdb
	done
	sed s/model/model_/g < $srcdir/${target}.model_summary.txt > $modelDir/${target}.model_summary
	cp $DistanceFoldingHome/0README $modelDir/
	currDir=`pwd`
	cd $ResDir
	zip -r ${target}_models.zip ${target}_models/*
	cd $currDir

	for i in 1 2 3 4 5
	do
		cp $srcdir/${target}_model${i}.pdb.gz $modelDir/${target}_model_${i}.pdb.gz
	done
fi

## generate the all_in_one.zip
##mkdir -p ${target}.all_in_one
cp -r $ResDir ${target}.all_in_one
find ${target}.all_in_one/ -name "*.gz" -type f -delete
find ${target}.all_in_one/ -name "*.zip" -type f -delete

zip -r ${target}.all_in_one.zip ${target}.all_in_one/*
mv ${target}.all_in_one.zip $ResDir/
rm -rf ${target}.all_in_one/
