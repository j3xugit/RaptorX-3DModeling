#!/bin/bash
if [ $# -ne 4 ]
then
        echo "Usage: ./combine_file <ss3> <ss8> <acc> <tgt> "
        exit
fi


# get ss3
a=`wc $1 | awk '{print $1}'`
b=$((a-2))
tail -n $b $1 > $1".tmp"
awk '{$1="";$2="";$3="";a=$4;$4=$5;$5=$6;$6=a;print $0 }' $1".tmp" > $1".tmp2"
awk '{print $2" "$1}' $1".tmp" > $1".tmp3"

# get ss8
c=`wc $2 | awk '{print $1}'`
d=$((c-3))
tail -n $d $2 > $2".tmp"
awk '{$1="";$2="";$3=""; print $0 }' $2".tmp" > $2".tmp2"

# check length
e=`wc $3 | awk '{print $1}'`
if [ "$b" != "$d" ]
then
	echo "ss3 length not equal to ss8 "
	exit
fi
if [ "$b" != "$e" ]
then
        echo "ss3 length not equal to acc "
        exit
fi

# combine file
paste $1".tmp2" $2".tmp2" $3 $1".tmp3" > $1".tmp4"
echo "#" > $4".tmp2"
cat $4 $1".tmp4" $4".tmp2" > $4".tmp"
mv $4".tmp" $4

# delete file
rm -f $1".tmp"
rm -f $1".tmp2"
rm -f $1".tmp3"
rm -f $1".tmp4"
rm -f $2".tmp"
rm -f $2".tmp2"
rm -f $4".tmp2"

