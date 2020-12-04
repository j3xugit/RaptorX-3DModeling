#!/bin/bash

CCMPRED=/home/jinbo/CCMpred/bin/ccmpred
A2MSrcDir=/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/10820_New_PDB25_Training_Data/10820_New_PDB25_Training_Feature/
#PDB25ListFile=/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/10820_New_PDB25_Training_Data/10820_new_training_fasta_list
PDB25ListFile=xai

for chain in `cat $PDB25ListFile`
do
	A2MFile=${A2MSrcDir}/feat_${chain}_contact/${chain}.a2m
	$CCMPRED -b $chain.ccmpred_raw.mpk $A2MFile ${chain}.ccmpred
done
