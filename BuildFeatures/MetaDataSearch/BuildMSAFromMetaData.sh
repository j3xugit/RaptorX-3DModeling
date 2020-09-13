#!/bin/bash

out_root=""          #-> output to current directory
cpu=4                #-> use 4 CPUs
merge=0              #-> default: 0 for NOT merge
e_value=0.001        #-> default: 0.001
seqid=0.90           #-> default: 0.90
#-> others
kill_tmp=1           #-> default: kill temporary root

cmd=`readlink -f $0`
home=`dirname $cmd`    #-> home directory
#data_db=$home/databases/metaclust50.fasta  #-> can be uniref90, or other sequence databases, with suffix ".fasta"
data_db=$MetaDB

function usage
{
	echo "$0 [-d database | -o out_root | -c CPU_num | -m merge | -e evalue | -s seqid | -K remove_tmp] input_fasta input_a3m"
	echo "    This script generates an extra MSA by searching through a metagenome database with a given A3M file and a sequence file"
	echo "***** required arguments *****"
	echo "-i input_fasta  : Query protein sequence in FASTA format. "
	echo ""
	echo "-I input_a3m    : Initial multiple sequence alignment in A3M format. "
	echo ""
	echo "Options:"
	echo ""
	echo "-d database     : The metagenome database for homology search, default $data_db"
	echo ""
	echo "-o out_root     : the folder for output, if not provided, proteinName_AddiA3M/ will be created in current work directory"
	echo ""
	echo "-c CPU_num      : number of CPUs to be used, default $cpu"
	echo ""
	echo "-m merge        : 1 for merge the additional A3M with input A3M, 0 for not merge, default $merge"
	echo "                  if -1, then realign the input A3M with the additional A3M. "
	echo ""
	echo "-e evalue       : E-value for homology search, default ${e_value}"
	echo ""
	echo "-s seqid        : the seqID cutoff used to remove redundancy in the additional A3M, default $seqid"
	echo "                  If set, this value MUST > 0.65. Redundancy not removed if set to 0"
	echo ""
	echo "-K remove_tmp   : 1 for remove temporary folder, and 0 not; default ${kill_tmp}"
}


#-> parse arguments
while getopts ":d:o:c:m:e:s:K:H:" opt;
do
	case $opt in
	d)
		data_db=$OPTARG
		;;
	o)
		out_root=$OPTARG
		;;
	c)
		cpu=$OPTARG
		;;
	m)
		merge=$OPTARG
		;;
	e)
		e_value=$OPTARG
		;;
	s)
		seqid=$OPTARG
		;;
	K)
		kill_tmp=$OPTARG
		;;
	H)
		home=$OPTARG
		;;
	#-> others
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
shift $((OPTIND -1))

if [ $# -ne 2 ]; then
	usage
	exit 1
fi

in_seq=$1
in_a3m=$2

if [ ! -f $data_db ]; then
	echo "ERROR: invalid metagenome database file $data_db"
	echo "	You may fix this by setting an environmental variable MetaDB to the metagenome database file"
	echo "	or providing the metagenome database file as an input argument of this script"
	exit 1
fi
 
if [ ! -f "$in_seq" ]; then
	echo "ERROR: invalid input seq file $in_seq"
	exit 1
fi

in_seq=`readlink -f $in_seq`
if [ ! -f "$in_a3m" ]; then
	echo "ERROR: invalid A3M file $in_a3m"
	exit 1
fi
in_a3m=`readlink -f $in_a3m`

fulnam=`basename $in_seq`
relnam=${fulnam%.*}

if [ "$out_root" == "" ]; then
	out_root=${relnam}_AddiA3M
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`

thres=0.65
if [ "$seqid" != "0" ]; then
	retv=`echo $seqid | awk '{if($0<a){print 0}else{print 1}}' a=$thres`
	if [ $retv -eq 0 ]; then
		echo "ERROR: current seqid value $seqid is smaller than 0.65">&2
		exit 1
	fi
fi

#---- step 0: create temporary folder ---#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="${out_root}/TMP_AddiA3M_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root

#---- step 1: convert A3M to A2M -----#
$home/util/A3M_To_A2M $in_a3m $tmp_root/$relnam.a2m 0

#---- step 2: build HMM ----#
$home/bin/hmmbuild --symfrac 0 $tmp_root/$relnam.hmm $tmp_root/$relnam.a2m

#---- step 3: search against sequence database  ---#
database=${data_db}
$home/bin/hmmsearch -A $tmp_root/$relnam.sto -o $tmp_root/$relnam.out --cpu $cpu --noali --notextw \
	-E $e_value --domE $e_value --incE $e_value --incdomE $e_value \
	--tblout $tmp_root/$relnam.tblout --domtblout $tmp_root/$relnam.domtblout $tmp_root/$relnam.hmm $database

#---- step 4: process aligned sequences ---#
grep -v "^#" $tmp_root/$relnam.sto | awk '{if(NF==2){print $0}}' > $tmp_root/$relnam.sto2
$home/util/reformat.pl sto a2m $tmp_root/$relnam.sto2 $tmp_root/$relnam.a2m
$home/util/MSA_To_SEQ $tmp_root/$relnam.a2m $tmp_root/$relnam.msa_addi
if [ $merge -eq -1 ]; then
	$home/util/MSA_To_SEQ $in_a3m $tmp_root/$relnam.msa_orig
	cat $tmp_root/$relnam.msa_orig $tmp_root/$relnam.msa_addi > $tmp_root/$relnam.msa
else
	cp $tmp_root/$relnam.msa_addi $tmp_root/$relnam.msa
fi
cat $in_seq $tmp_root/$relnam.msa > $tmp_root/$relnam.msa_seq

#---- step 4.5: remove redundancy ----#
if [ "$seqid" != "0" ]; then
	$home/bin/cd-hit -i $tmp_root/$relnam.msa_seq -o $tmp_root/$relnam.msa_nored -c $seqid -s 0.9
else
	cp $tmp_root/$relnam.msa_seq $tmp_root/$relnam.msa_nored
fi

#---- step 5: re-run hmmsearch ----#
$home/bin/hmmsearch -A $tmp_root/$relnam.sto_ii -o $tmp_root/$relnam.out_ii --cpu $cpu --noali --notextw \
	-E $e_value --domE $e_value --incE $e_value --incdomE $e_value \
	--tblout $tmp_root/$relnam.tblout_ii --domtblout $tmp_root/$relnam.domtblout_ii $tmp_root/$relnam.hmm $tmp_root/$relnam.msa_nored

#---- step 6: put self to first ---#
$home/util/proc_sto.sh $tmp_root/$relnam.sto_ii $relnam $tmp_root/$relnam.sto_fin

#---- step 7: generate final A3M ---#
$home/util/reformat.pl sto a3m $tmp_root/$relnam.sto_fin $tmp_root/$relnam.a3m -M first
$home/util/A3M_Seq_Refine $tmp_root/$relnam.a3m $in_seq $tmp_root/$relnam.a3m_fin

#---- step 8: merge or not ----#
if [ $merge -eq 1 ]; then
	sed '1,2d' $tmp_root/$relnam.a3m_fin > $tmp_root/$relnam.a3m_fin_noquery
	cat $in_a3m $tmp_root/$relnam.a3m_fin_noquery > $out_root/$relnam.a3m
else
	cp $tmp_root/$relnam.a3m_fin $out_root/$relnam.a3m
fi


#------- remove and exit -----#
if [ $kill_tmp -eq 1 ]; then
	rm -rf $tmp_root
else
	rm -rf $out_root/"TMP_AddiA3M_"${relnam}
	mv $tmp_root $out_root/"TMP_AddiA3M_"${relnam}
fi
cp $in_seq $out_root/$relnam.fasta_raw
cp $in_a3m $out_root/$relnam.a3m_orig
