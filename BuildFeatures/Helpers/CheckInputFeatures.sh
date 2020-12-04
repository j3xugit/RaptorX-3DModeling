#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 proteinList featureFolder
	exit 1
fi

listFile=$1
featureFolder=$2

for protein in `cat $listFile`
do
	if [ ! -d $featureFolder/feat_${protein}_contact/ ]; then
		echo feat_${protein}_contact
		continue
	fi
	if [ ! -f $featureFolder/feat_${protein}_contact/${protein}.a2m ]; then
		echo feat_${protein}_contact
		continue
	fi
	if [ ! -f $featureFolder/feat_${protein}_contact/${protein}.inputFeatures.pkl ]; then
		echo feat_${protein}_contact
		continue
	fi
	if [ ! -f $featureFolder/feat_${protein}_contact/${protein}.extraCCM.pkl ]; then
                echo feat_${protein}_contact
		continue
        fi
done
