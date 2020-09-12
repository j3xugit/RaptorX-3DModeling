#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./genACC <tgt_name> <tmp_root> "
        exit
fi

# ---- get arguments ----#
tgt_name=$1
tmp_root=$2

# ---- process -----#
RaptorX_HOME=~/GitBucket/casp13_pipeline/TGT_Package
ACCPred=$RaptorX_HOME/util/ACC_Predict/acc_pred
$ACCPred $RaptorX_HOME/$tmp_root/$tgt_name.hhm $RaptorX_HOME/$tmp_root/$tgt_name.ss2 $RaptorX_HOME/$tmp_root/$tgt_name.ss8 $RaptorX_HOME/util/ACC_Predict/model.accpred $RaptorX_HOME/$tmp_root/$tgt_name.acc

