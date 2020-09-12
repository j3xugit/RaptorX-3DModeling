#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./genSS8 <tgt_name> <tmp_root>"
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

SS8_Pred=$DistFeatureHome/util/SS8_Predict/bin/run_raptorx-ss8.pl
$SS8_Pred $tmp_root/$tgt_name.seq -pssm $tmp_root/$tgt_name.psp -outdir $tmp_root/

