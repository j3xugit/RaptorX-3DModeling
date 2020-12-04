#!/bin/bash

function Usage 
{
	echo "$0 [ -M | -s MSAmethod | -d savefolder ] proteinListFile MetaDir"
	echo "	This script links the input feature files (*.a2m, *.extraCCM.pkl, *.inputFeatures.pkl) dispersed in multiple XXX_OUT/XXX_contact folders into a single folder"
	echo "	It is mainly used for batch contact/distance/orientation prediction of multiple proteins but not for multiple MSAs"
	echo "	To do batch contact/distance/orientation prediction of multiple MSAs, you may use PredictPairwiseRelation4MultiInputs.sh"
	echo "	proteinListFile: a file for a list of protein names, each in one row"
	echo "	MetaDir: a folder contains a set of folders with name XXX_OUT where XXX is the protein names"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-M: when specified, no meta genome data used, default use metagenome data"
	echo "  -s: the MSA generation method, 0 for hhblits, 1 for jackhmmer, 2 for both, 3 for user-provided MSA and 4 for all three, default 4"
}

MSAmethod=4
UseMetaGenomeData=true
ResDir=`pwd`

while getopts ":Ms:d:" opt; do
        case ${opt} in
                M )
                  UseMetaGenomeData=false
                  ;;
		s )
		  MSAmethod=$OPTARG
		  ;;
                d )
                  ResDir=$OPTARG
                  ;;
                \? )
                  echo "Invalid Option: -$OPTARG" 1>&2
                  exit 1
                  ;;
                : )
                  echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                  exit 1
                  ;;
        esac
done
shift $((OPTIND -1))

if [ $# -lt 2 ]; then
	Usage
	exit 1
fi

targetListFile=$1
if [ ! -f $targetListFile ]; then
	echo "ERROR: invalid protein list file $targets"
	exit 1
fi
targets=`cat $targetListFile `

MetaDir=`readlink -f $2`
if [ ! -d $MetaDir ]; then
	echo "ERROR: invalid folder for input features $MetaDir"
	exit 1
fi

if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi
ResDir=`readlink -f $ResDir`

currDir=`pwd`

if [ $MSAmethod -eq 4 ]; then
        methods="uce3 uce5 ure3 ure5 user"
elif [ $MSAmethod -eq 2 ]; then
        methods="uce3 uce5 ure3 ure5"
elif [ $MSAmethod -eq 0 ]; then
        methods="uce3 uce5"
elif [ $MSAmethod -eq 1 ]; then
        methods="ure3 ure5"
elif [ $MSAmethod -eq 3 ]; then
	methods="user"
else
        echo "ERROR: please specify which MSAs shall be used, hhblits MSA or jackhmmer MSA ?"
        exit 1
fi

for method in $methods
do
	if $UseMetaGenomeData; then
		FeatureDir=$ResDir/features-${method}_meta
	else
		FeatureDir=$ResDir/features-${method}
	fi

	mkdir -p $FeatureDir
	cd $FeatureDir

	for target in $targets
	do
		if $UseMetaGenomeData; then
			srcDir=$MetaDir/${target}_OUT/${target}_contact/feat_${target}_${method}_meta
			if [ ! -d $srcDir ]; then	
				srcDir=$MetaDir/${target}_OUT/${target}_contact/feat_${target}_${method}
			fi
		else
			srcDir=$MetaDir/${target}_OUT/${target}_contact/feat_${target}_${method}
		fi

		if [ ! -d $srcDir ]; then
			echo "WARNING: invalid input folder $srcDir"
			continue
		fi

		ln -s $srcDir/${target}.a2m
		for i in $srcDir/${target}.*.pkl
		do
			if [ ! -f $i ]; then
				continue
			fi
			ln -s $i
		done
		#ln -s $srcDir/${target}.extraCCM.pkl
		#ln -s $srcDir/${target}.inputFeatures.pkl
	done
	cd $currDir

	## remove empty folder
	[ "$(ls -A $FeatureDir)" ] && echo "$FeatureDir Not Empty" || rmdir $FeatureDir
done
