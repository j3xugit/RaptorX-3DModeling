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
RaptorX_HOME=~/GitBucket/casp13_pipeline/TGT_Package
SS8_Pred=$RaptorX_HOME/util/SS8_Predict/bin/run_raptorx-ss8.pl
$SS8_Pred $RaptorX_HOME/$tmp_root/$tgt_name.seq -pssm $RaptorX_HOME/$tmp_root/$tgt_name.psp -outdir $RaptorX_HOME/$tmp_root/

