#!/bin/bash

if [ $# -lt 4 ]
then
	echo "Usage: ./contact_feature_make.sh <input_tgt> <input_a3m> <cpu_num> <out_root> "
	echo "[note]: all input/output files should be absolute directory "
	exit
fi
curdir="$(pwd)"


# ---- input files ---- #
cpu_num=$3
out_root=$4

#-> get job id:
fulnam=`basename $1`
relnam=${fulnam%.*}
#tmp=feat"_"$relnam"_contact"

#-> create temporary folder
tmp=$out_root
mkdir -p $tmp/

#-> bin folder
if [[ -z "${DistFeatureHome}" ]]; then
	echo "Please set the environmental variable DistFeatureHome, e.g., $HOME/3DModeling/BuildFeatures"
	exit 1
fi

bin="$DistFeatureHome/bin/"

# ---- copy TGT and A3M to tmp ---------- #
tgt_file=$relnam.tgt
if [ ! -f "$tmp/$tgt_file" ]
then
	cp $1 $tmp/$tgt_file
fi
a3m_file=$relnam.a3m
if [ ! -f "$tmp/$a3m_file" ]
then
	cp $2 $tmp/$a3m_file
fi

# ---- from TGT file to SEQ file -------- #
seq_file=$relnam.seq
if [ ! -f "$tmp/$seq_file" ]
then
	echo ">$relnam" > $seq_file
	head -n 4 $tmp/$tgt_file | tail -n1 | awk '{print $3}' >> $seq_file
	mv $seq_file $tmp/
fi

# ---- from A3M file to A2M file -------- #
a2m_file=$relnam.a2m
if [ ! -f "$tmp/$a2m_file" ]
then
	$DistFeatureHome/util/A3M_To_PSI $tmp/$a3m_file $a2m_file.tmp 
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "failed in util/A3M_To_PSI $tmp/$a3m_file $a2m_file.tmp"
		exit 1
	fi
	grep -v "ss_pred\|ss_conf" $a2m_file.tmp | awk '{print substr($0,34,length($0)-32) }' > $a2m_file
	rm -f $a2m_file.tmp
	mv $a2m_file $tmp/
fi

# ---- from A2M file, generate ccmpred_file -------- #
ccmpred_file_pre=$relnam.ccmpred
if [ ! -f "$tmp/$ccmpred_file_pre" ]
then
	echo "ERROR: Cannot find $tmp/$ccmpred_file_pre"
        exit 1

	echo "ccmpred gpu start"
	$bin/ccmpred_gpu -R -d 0 $tmp/$a2m_file $ccmpred_file_pre > $relnam.ccmpred_log       #-> using GPU version
	if [ ! -f "$ccmpred_file_pre" ]
	then
		echo "ccmpred cpu start"
		$bin/ccmpred_cpu -R $tmp/$a2m_file $ccmpred_file_pre > $relnam.ccmpred_log    #-> using CPU version
	fi
	echo "ccmpred done"
	#--> file check
	if [ ! -f "$ccmpred_file_pre" ]
	then
		echo "failed in $bin/ccmpred -R $tmp/$a2m_file $ccmpred_file_pre"
		exit 1
	fi
fi
ccmpred_file=$relnam.ccmpred_zscore
if [ ! -f "$tmp/$ccmpred_file" ]
then
	python $bin/normalize_ccmpred_sep.py $tmp/$ccmpred_file_pre > $ccmpred_file
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "failed in python $bin/normalize_ccmpred_sep.py $tmp/$ccmpred_file_pre > $ccmpred_file"
		exit 1
	fi
	mv $ccmpred_file $tmp/
fi

# ---- from A2M file, generate potential file ------- #
pot_file=$relnam.pot
if [ ! -f "$tmp/$pot_file" ]
then
	echo "alnstats_omp start"
	$bin/alnstats_omp $tmp/$a2m_file $relnam.ws1 $pot_file
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "failed in $bin/alnstats_omp $tmp/$a2m_file $relnam.ws1 $pot_file"
		exit 1
	fi
	rm -f $relnam.ws1
	echo "alnstats_omp done"
	mv $pot_file $tmp/
fi

# ---- generate SS3 and ACC---- #
ss3_file=$relnam.ss3
if [ ! -f "$tmp/$ss3_file" ]
then
	$bin/DeepCNF_SS_Con -t $tmp/$tgt_file -s 1 > $ss3_file
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "failed in $bin/DeepCNF_SS_Con -t $tmp/$tgt_file -s 1 > $ss3_file"
		exit 1
	fi
	mv $ss3_file $tmp/
fi
acc_file=$relnam.acc
if [ ! -f "$tmp/$acc_file" ]
then
	$bin/AcconPred $tmp/$tgt_file 1 > $acc_file
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "failed in $bin/AcconPred $tmp/$tgt_file 1 > $acc_file"
		exit 1
	fi
	mv $acc_file $tmp/
fi

# ---- printTGT ----#
profile_file=$relnam.profile
if [ ! -f "$tmp/$profile_file" ]
then
	$bin/printTGT $tmp/$tgt_file > $profile_file
	mv $profile_file $tmp/
fi

# ---- make a pseudo DISO file -----#
diso_file=$relnam.diso
if [ ! -f "$tmp/$relnam.diso" ]
then
	tail -n+2 $tmp/$acc_file | awk '{if(NF==8){print $1" "$2" . 0"}else{print $0}}'  > $diso_file
	mv $diso_file $tmp/
fi


# ================ exit 0 ================= #
exit 0

