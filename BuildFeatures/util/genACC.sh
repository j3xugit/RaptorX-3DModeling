#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./genACC.sh <tgt_name> <tmp_root> "
        exit
fi

# ---- get arguments ----#
tgt_name=$1
tmp_root=$2

# ---- process -----#
if [[ -z "${DistFeatureHome}" ]]; then
	echo "Please set environmental variable DistFeatureHome, e.g., $HOME/3DModeling/BuildFeatures"
	exit 1
fi

ACCPred=$DistFeatureHome/util/ACC_Predict/acc_pred
$ACCPred $tmp_root/$tgt_name.hhm $tmp_root/$tgt_name.ss2 $tmp_root/$tgt_name.ss8 $DistFeatureHome/util/ACC_Predict/model.accpred $tmp_root/$tgt_name.acc

