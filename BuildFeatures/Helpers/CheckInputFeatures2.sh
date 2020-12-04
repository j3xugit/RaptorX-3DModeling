#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 proteinList featureFolder
	exit 1
fi

listFile=$1
featureFolder=$2

for protein in `cat $listFile`
do
	if [ ! -f $featureFolder/${protein}.a2m ]; then
		echo $protein.a2m
	fi
	if [ ! -f $featureFolder/${protein}.inputFeatures.pkl ]; then
		echo $protein.inputFeatures.pkl
	fi
	if [ ! -f $featureFolder/${protein}.extraCCM.pkl ]; then
                echo $protein.extraCCM.pkl
        fi
done
