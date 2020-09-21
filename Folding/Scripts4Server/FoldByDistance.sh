#!/bin/bash

# ----- usage ------ #
usage 
{
	echo "FoldByDistance v0.01 [Mar-14-2017] "
	echo "    Fold proteins by predicted inter-residue distance and secondary structure"
	echo ""
	echo "USAGE:  ./FoldByDistance.sh <-i input_fasta> <-s input_ss> < -B distance_bound_PKL> "
	echo "         [-o out_root] [-m model_num] [-c CPU_num] [-k topk] [-r rand_seed] "
	echo "         [-n NOE_scale] [-d Angle_scale] [ -b dist_cutoff] [-g sigma_value] [ -t bound_type]"
	echo ""
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta  :  input protein sequence file in FASTA format"
	echo "-s input_ss     :  input predicted secondary structure in FASTA format"
	echo "-B distance_bound_PKL    :  estimated distance and lower and upper deviations in symmetrical matrix format in a PKL file"
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_root     :  default output directory is tmp/ at the current directory"
	echo "-m model_num    :  the number of models per CPU [default 20]"
	echo "-c CPU_num      :  the number of CPUs to be used [default 1]"
	echo "-k topk_num     :  output TopK minimal NOE models [default 5]"
	echo "-r rand_seed    :  random seed for 'dgsa.inp' file [default -1]"
	echo "-n NOE_scale    :  scale factor for distance restraints [default 5.0]"
	echo "-d Angle_scale  :  scale factor for dihedral restraints [default 10.0]"
	echo "-b dist_cutoff  :  distance cutoff (default 12.0). Only distance less than this cutoff will be used for model building"
	echo "-g sigma_value  :  distance bound with dist +/- sigma*std [default 1.0]"
	echo "-t bound_type   :  the bound type to be used, ranging from 0 to 4 (default 1). See ConvertDistBounds2CNSTBL.py and EstimateDistances.py for explanation"
	exit 1
}

if [ $# -lt 6 ];
then
        usage
fi
##curdir="$(pwd)"

# ----- main directory ---#
##confold_server_home=$curdir
confold_server_home=/home/jinbo/3DModeling/Folding/
#-> check directory
if [ ! -f "$confold_server_home/FoldByDistance.sh" ]
then
	echo "The program file $confold_server_home/FoldByDistance.sh not exist."
	##echo "please run './setup.pl' to configure the package."
	echo "please set confold_server_home in this script program to its installation directory."
	exit 1
fi

# ----- vice directory ---#
cns_solve=$confold_server_home/util/cns_solve
util=$confold_server_home/util
bin=$confold_server_home/bin

# ----- get arguments ----- #
#-> required arguments
input_fasta=""
input_ss=""
dist_bound=""
output=""
#-> optional arguments
CPU_num=1
model_num=20
topk=5
rand_seed=-1
NOE_scale=5.0
Angle_scale=10.0
sigma_value=1
dist_cutoff=12
bound_type=1

#-> parse arguments
while getopts ":i:s:B:o:m:c:k:r:n:d:b:g:t:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input_fasta=$OPTARG
		;;
	s)
		input_ss=$OPTARG
		;;
	B)
		dist_bound=$OPTARG
		;;
	#-> optional arguments
	o)
		output=$OPTARG
		;;
	m)
		model_num=$OPTARG
		;;
	c)
		CPU_num=$OPTARG
		;;
	k)
		topk=$OPTARG
		;;
	r)
		rand_seed=$OPTARG
		;;
	n)
		NOE_scale=$OPTARG
		;;
	d)
		Angle_scale=$OPTARG
		;;
	b)      
		dist_cutoff=$OPTARG
		;;
	t) 
		bound_type=$OPTARG
		;;
	g)
		sigma_value=$OPTARG
		;;
	#-> default
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done

# ------ get basefilename -------- #
filename=`basename $input_fasta`
extension=${filename##*.}
basefilename=${filename%.*}

# ----- setup process root ----#
tmp=tmp/$basefilename
if [ "$output" != "" ]
then
	tmp=$output
fi
mkdir -p $tmp

# ------ get ld_lib for cns ----#
#ld_lib=$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=~/anaconda2/lib/:$LD_LIBRARY_PATH
source $cns_solve/cns_solve_env.sh
export CNS_CUSTOMMODULE=$util
export KMP_AFFINITY=none

#================ Step 1. generate extended initial structure ==============#
$bin/SEQ_To_Input $input_fasta $tmp/input.seq
cd $tmp
	$cns_solve/intel-x86_64bit-linux/bin/cns_solve < $util/gseq.inp > gseq.log
	$cns_solve/intel-x86_64bit-linux/bin/cns_solve < $util/extn.inp > extn.log
cd $confold_server_home

#================ Step 2. generate all relevant TBL files ==============#
#-> 2.1 distance relevant TBL file
#$bin/EstiDist_To_TBL $input_fasta $esti_dist $esti_dist $lower_std $upper_std $tmp/contact.tbl $sigma_value
python $bin/ConvertDistBounds2CNSTBL.py -c $dist_cutoff -t $bound_type -g $sigma_value $input_fasta $dist_bound > $tmp/contact.tbl

#-> 2.2 secondary structure relevant TBL files
$bin/SSE_To_TBL $input_fasta $input_ss $tmp/ssnoe.tbl $tmp/dihedral.tbl $tmp/hbond.tbl

#================ Step 3. run CNS_Solve for 3D structure folding =======#
#-> 3.1 modify 'dgsa.inp' file
for ((i=1;i<=$CPU_num;i++))
do
	#-> modify 'dgsa.inp' file with given parameters
#	if [ $rand_seed -lt 0 ]
#	then
#		rand_seed=`od -An -N4 -i < /dev/urandom | awk '{if($0<0){print -$0}else {print $0}}'`
#	fi
	$bin/DGSA_File_Mod $input_fasta $util/dgsa.inp_raw $tmp/dgsa_${i}.inp $model_num $rand_seed ${i}_${basefilename} $NOE_scale $Angle_scale
	#-> touch iam.running
	touch $tmp/iam.running_${i}
done

#-> 3.2 running cns_solve
iamfailed=0
cd $tmp
	#-> run cns_solve
	for ((i=1;i<=$CPU_num;i++))
	do
		($cns_solve/intel-x86_64bit-linux/bin/cns_solve < dgsa_${i}.inp > dgsa_${i}.log)&
	done
	wait
	#-> check results
	for ((i=1;i<=$CPU_num;i++))
	do
		if [ -f "${i}_${basefilename}_${model_num}.pdb" ]
		then
			rm iam.running_${i}
			rm ${i}_${basefilename}_*embed*.pdb
		else
			tail -n 30 dgsa_${i}.log
			echo "ERROR! Final structures not found for batch $i!"
			mv iam.running_${i} iam.failed_${i}
			iamfailed=1
		fi
	done
cd $confold_server_home
##export LD_LIBRARY_PATH=$ld_lib

#-> 3.3 check failed
if [ $iamfailed -eq 1 ]
then
	echo "CNS FAILED!"
	exit 2
fi


#================ Step 4. collect all models, and return final topk sorted by NOE =======#
#-> 4.1 collect noe in all generated models
rm -f $tmp/${basefilename}.noe_statis
for ((i=1;i<=$CPU_num;i++))
do
	for ((k=1;k<=${model_num};k++))
	do
		a=`grep "REMARK noe     =" $tmp/${i}_${basefilename}_${k}.pdb | awk '{print $NF}'`
		echo "${i}_${basefilename}_${k}.pdb $a" >> $tmp/${basefilename}.noe_statis

	done
done

#-> 4.2 sort these models by NOE in ascending order, and return topk =====#
sort -n -k2 $tmp/${basefilename}.noe_statis > $tmp/${basefilename}.noe_statis_sort
for ((i=1;i<=$topk;i++))
do
	file=`head -n $i $tmp/${basefilename}.noe_statis_sort | tail -n1 | awk '{print $1}'`
	cp $tmp/$file $tmp/${basefilename}_model$i.pdb
done


#====== remove temporary files ==========#
##rm -f $tmp/*.log

# ---- exit ----- #
exit 0


