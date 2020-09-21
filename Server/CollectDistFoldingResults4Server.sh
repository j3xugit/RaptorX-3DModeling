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
if [ ! -d $SrcDir ]; then
	echo "ERROR: invalid $SrcDir"
	exit 1
fi

ResDir=${target}_Results
mkdir -p $ResDir

cp $DistanceFoldingHome/0README_overall $ResDir/

## copy sequence files
cp $SrcDir/${target}.seq $ResDir/
cp $SrcDir/${target}.seq $ResDir/${target}.fasta

## copy MSAs
srcdir=${SrcDir}/${target}_contact/feat_${target}_user
if [ ! -d $srcdir ]; then
	srcdir=$SrcDir/${target}_contact/feat_${target}_uce3_meta
	if [ ! -d $srcdir ]; then
		srcdir=$SrcDir/${target}_contact/feat_${target}_uce3
	fi
fi
if [ ! -d $srcdir ]; then
	echo "ERROR: cannot find $srcdir for MSAs"
	exit 1
fi
cp $srcdir/${target}.ss3 $ResDir/
cp $srcdir/${target}.inputFeatures.pkl $ResDir/

## for display purpose. Since drawing a large MSA is slow, here we display up to 5000 sequences
head -2  $srcdir/${target}.a3m > ${target}.a3m.head
seqLen=`tail -1 ${target}.a3m.head | wc -c | cut -f1 -d' ' `
tail -n +3 $srcdir/${target}.a3m | paste -sd '*\n' | shuf -n 5000 | sed 's/*/\n/g' >> ${target}.a3m.head
cut -c1-$seqLen ${target}.a3m.head > $ResDir/${target}.a3m
rm -f ${target}.a3m.head

for srcdir in ${SrcDir}/${target}_contact/feat_${target}_*
do
	if [ ! -d $srcdir ]; then
		continue
	fi
	method=`basename $srcdir | cut -f2- -d'_'`
	gzip < $srcdir/${target}.a3m > $ResDir/${target}.a3m_${method}.gz
done

## copy predicted distance information
srcdir=${SrcDir}/DistancePred
if [ ! -d $srcdir ]; then
	echo "ERROR: cannot find the folder for predicted distance info: $srcdir"
	exit 1
fi
cp $srcdir/${target}*.png $ResDir/
cp $srcdir/${target}.CM.txt $ResDir/${target}.gcnn
cp $srcdir/${target}.CASP.rr $ResDir/${target}.rr
#cp $srcdir/${target}.CASP.rr $ResDir/contactmap.txt
cp $srcdir/${target}.CASP.rr $ResDir/${target}.contactmap.txt
cp $srcdir/${target}.bound.txt $ResDir/
cp $srcdir/${target}.predictedDistMatrix.pkl $ResDir/
cp $srcdir/${target}.pairPotential.DFIRE16.pkl $ResDir/
cp $srcdir/${target}.epad_prob $ResDir/

## copy predicted properties
srcdir=${SrcDir}/PropertyPred
if [ ! -d $srcdir ]; then
	echo "ERROR: cannot find the folder for predicted properties: $srcdir"
	exit 1
fi
cp $srcdir/${target}.predictedProperties.pkl $ResDir/
cp $srcdir/${target}.propertyFeatures.pkl $ResDir/

## copy predicted 3D models
srcdir=${SrcDir}/${target}-SpickerResults/
if [ -d $srcdir ]; then
	cp /dev/null $ResDir/${target}.pdb
	modelDir=$ResDir/${target}_models
	mkdir -p $modelDir
	cp $ResDir/${target}.contactmap.txt $modelDir/

	cp /dev/null $modelDir/${target}.model_summary

	numModels=0
	for i in 0 1 2 3 4
	do
		if [ ! -f $srcdir/${target}_center${i}.quality.pdb ]; then
			continue
		fi
		numModels=`expr $numModels + 1 `
		j=`expr $i + 1 `

		cp $srcdir/${target}_center${i}.quality.pdb $modelDir/${target}_model_${j}.pdb
		gzip < $modelDir/${target}_model_${j}.pdb > $modelDir/${target}_model_${j}.pdb.gz

		## generate quality file
		qualityFile=$srcdir/Quality-${target}_center${i}.txt
		CaRMSD=`grep ^REMARK $qualityFile | grep -w Ca | cut -f10 -d' ' | cut -c1-6`
		echo ${target}_model_${j}.pdb $CaRMSD >> $modelDir/${target}.model_summary

		qualityFileBaseName=`basename $qualityFile`
		grep -w CA $qualityFile | cut -c1 > $qualityFileBaseName.col1
		grep -w CA $qualityFile | cut -f6 > $qualityFileBaseName.col6
		paste $qualityFileBaseName.col1 $qualityFileBaseName.col6 -d' ' > $qualityFileBaseName.col1-6

		finalQualityFile=$modelDir/${target}_model${j}.localQuality.txt
		cp /dev/null $finalQualityFile
		echo "REMARK: estimated Ca RMSD for ${target}_model_${j}.pdb" >> $finalQualityFile
		echo "REMARK: the three columns are residue index, amino acid and distance deviation" >> $finalQualityFile
		nl -n ln -w 1 -s ' ' $qualityFileBaseName.col1-6 >> $finalQualityFile
		rm -f $qualityFileBaseName.col*

		echo "MODEL $j" >> $ResDir/${target}.pdb
		grep ^ATOM $srcdir/${target}_center${i}.quality.pdb >> $ResDir/${target}.pdb
		echo "ENDMDL" >>  $ResDir/${target}.pdb
	done
	if [ $numModels -eq 0 ]; then
		echo "ERROR: cannot find any 3D models in $srcdir "
		exit 1
	fi
	cp $DistanceFoldingHome/0README $modelDir/
	currDir=`pwd`
	cd $ResDir
	zip -r ${target}_models.zip ${target}_models/*
	cd $currDir
fi

## generate the all_in_one.zip
##mkdir -p ${target}.all_in_one
cp -r $ResDir ${target}.all_in_one
find ${target}.all_in_one/ -name "*.gz" -type f -delete
find ${target}.all_in_one/ -name "*.zip" -type f -delete

zip -r ${target}.all_in_one.zip ${target}.all_in_one/*
mv ${target}.all_in_one.zip $ResDir/
rm -rf ${target}.all_in_one/
rm -f $ResDir/${target}.predictedDistMatrix.pkl
rm -f $ResDir/${target}.pairPotential.DFIRE16.pkl
rm -f $ResDir/${target}.bound.txt
rm -f $ResDir/${target}.epad_prob
