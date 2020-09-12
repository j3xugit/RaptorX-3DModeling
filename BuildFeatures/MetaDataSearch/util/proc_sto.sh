#!/bin/bash

#------- usage ------------#
if [ $# -lt 3 ]
then
	echo "Usage: ./proc_sto.sh <input_sto> <seq_nam> <output_sto> "
	exit 1
fi

#------ the arguments ----------#
input_sto=$1
seq_nam=$2
output_sto=$3

#------ get name ------#
fulnam=`basename $input_sto`
relnam=${fulnam%.*}

#--- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="TMP_STO_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root

#--- process ----#
grep -v "^#" $input_sto | awk '{if(NF==2){print $0}}' > $tmp_root/$relnam.sto
line=`grep -n "^$seq_nam" $tmp_root/$relnam.sto | cut -d ':' -f 1`
prev=$((line-1))
next=$((line+1))

#--- get new file ---#
head -n $line $tmp_root/$relnam.sto | tail -n1 > $tmp_root/$relnam.first
head -n $prev $tmp_root/$relnam.sto > $tmp_root/$relnam.prev
tail -n+$next $tmp_root/$relnam.sto > $tmp_root/$relnam.next
cat $tmp_root/$relnam.first $tmp_root/$relnam.prev $tmp_root/$relnam.next > $output_sto

#--- remove -----#
rm -rf $tmp_root
exit 0


