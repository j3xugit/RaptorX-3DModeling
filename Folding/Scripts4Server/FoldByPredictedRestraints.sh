#!/bin/bash

# ----- usage ------ #
Usage 
{
	echo "This script folds one protein by predicted inter-atom contacts/distance and predicted secondary structure/torsion angles"
	echo ""
	echo "USAGE: " $0 " <-i seqFile> < -d predDistMatrix> <-p predProperties> -o resultDir"
	echo "		-a atomPairType -t boundType -c cutoff -g sigma -r minRatio:maxRatio -s seqSeparation"
	echo "		-e PhiPsiRestraintConstant -w DistanceScale:AngleScale"
	echo "		-m numModels:numCPUs -k topK -b randSeed"
	echo ""
	echo "***** input and output information *****"
	echo "-i seqFile:		input protein sequence file in FASTA format"
	echo "-d predDistMatrix:	predicted distance and contact matrix file in PKL format"
	echo "-p predProperties:	predicted property file in PKL format"
	echo "-o resultDir:		output directory, tmp/targetName/ by default"
	echo ""
	echo "*********atom pair type, bound type, and some relevant parameters*********"
	echo "-a atomPairTypes:	atom pair types to be used, e.g., CbCb (default), CaCa, CaCg, CgCg, CaCa+CbCb, CbCb+CgCg, CaCa+CbCb+CgCg or All"
	echo "-t boundType:	distance bound type, can be 0 to 4 (default 3). If negative, then contact instead of distance is used"
	echo "-c cutoff:	distance cutoff (default 15.0). Only distance less than this cutoff is used for building distance restraints"
	echo "				when boundType is contact, it shall be interpreted as Zscore_cutoff of contact probability (default 3.0)"
	echo "-g sigma:		amplication factor for distance deviation, i.e., we construct distance bound by estimated_dist +/- sigma*std [default 1.0]"
	echo ""
	echo "*********usually you may use the default values for the below parameters********"
	echo "-r minRatio:maxRatio:The min and max number of top contacts to be selected divided by sequence length (default 1.0:3.0)"
	echo "-s seqSeparation:    Only two residues whose sequence separation is at least this value will appear in a distance restraint (default 2 for distance and 6 for contacts)"
	echo ""
	echo "-w DistanceScale:AngleScale:	scale factor for distance restraints (default 5.0) and Phi/Psi restraints (default 10.0)"
	echo "					When contacts are used as restraints, the distance and dihedral factors are set to 10 and 5, respectively"
	echo "-e PhiPsiRestraintConstant:	coeffcient for the Phi/Psi angle restraint (default 5.0)"
	echo ""
	echo "-m numModels:numCPUs: The number of models per CPU (default 20) and the number of CPUs to be used (default 3)"
	echo "-k topK:		    Output topK 3D models (default 5)"
	echo ""
	echo "-b randSeed:	random seed used in the 'dgsa.inp' file (default -1, in this case computers will pick up one random seed by itself)"
	echo "		   	By manually setting this parameter to a given value, you can check if each run can generate the same 3D models"
}

if [ $# -lt 6 ];
then
        Usage
	exit -1
fi

# ----- path information ---#
FoldingHome=$DistanceFoldingHome

if [ -z "${FoldingHome}" ]; then
	echo "Please set up the environmental variable DistanceFoldingHome to the correct path, e.g., $ModelingHome/Folding/"
	exit 1
fi


util=$FoldingHome/util
bin=$FoldingHome/bin
cns_solve=$FoldingHome/util/cns_solve

# ----- set default values for arguments ----- #
#-> input and output files and folders
seqFile=""
distFile=""
angleFile=""
output=""

## distance bound type, atom pair type and cutoff
boundType=3
atomPairTypes="CbCb"
distCutoff=15
Zcutoff=3
minRatio=1.0
maxRatio=3.0
sigma_value=1.0
seqSeparation=2

#-> energy coefficients
NOEScale=5
AngleScale=10
##NOEScaleMore=1
dihedralEnergyConstant=5.0

## run folding on 3 CPUs and on each CPU generate 20 3D models. Finally pick up 5 top models.
numCPUs=6
numModels=20
topk=5

## when randSeed is negative, computer will pick up a random seed
randSeed=-1

#-> parse arguments
while getopts ":i:d:p:o:a:t:c:g:r:s:e:w:m:k:b:" opt;
do
	case $opt in
	#-> required arguments
	i)
		seqFile=$OPTARG
		;;
	d)
		distFile=$OPTARG
		;;
	p)
		angleFile=$OPTARG
		;;
	o)
		output=$OPTARG
		;;
	a) 
		atomPairTypes=$OPTARG
		;;
	t) 
		boundType=$OPTARG
		;;
	c)      
		distCutoff=$OPTARG
		Zcutoff=$OPTARG
		;;
	g)
		sigma_value=$OPTARG
		;;
	r)
		minRatio="$(echo $OPTARG | cut -d':' -f1)"
		maxRatio="$(echo $OPTARG | cut -d':' -f2)"
		;;
	s)
		seqSeparation=$OPTARG
		;;
	e)
		dihedralEnergyConstant=$OPTARG
		;;
	w)
		NOEScale="$(echo $OPTARG | cut -d':' -f1)"
		AngleScale="$(echo $OPTARG | cut -d':' -f2)"
		;;
	m)
		numModels="$(echo $OPTARG | cut -d':' -f1)"
		numCPUs="$(echo $OPTARG | cut -d':' -f2)"
		;;
	k)
		topk=$OPTARG
		;;
	b)
		randSeed=$OPTARG
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

##numAtomPairTypes=$(grep -o '+' <<< "$atomPairTypes" | wc -l)
##numAtomPairTypes=`expr $numAtomPairTypes + 1`

##lowerAtomPairTypes=` echo "$atomPairTypes" | tr '[:upper:]' '[:lower:]' `
##if [ $lowerAtomPairTypes == 'all' ]; then
##	numAtomPairTypes=5
##fi

#echo '#atom pair types:' $numAtomPairTypes

if [ ${boundType} -lt 0 ]; then
	NOEScale=10
	AngleScale=5
	##NOEScaleMore=0
fi

initNOEScale=$NOEScale

##if [ ${NOEScaleMore} -eq 1 ]; then
##	NOEScale=` echo "scale=3;" ${NOEScale} \* 15 / ${distCutoff} | bc `
##fi

##NOEScale=` echo "scale=1;" ${NOEScale} / ${numAtomPairTypes} | bc `

#echo "initNOEScale=${initNOEScale}, NOEScale=${NOEScale}, AngleScale=${AngleScale}"

if [[ $seqFile == "" || ! -f $seqFile ]]; then
	echo $seqFile is an invalid input sequence file
	exit -1
fi

# ------ get targetName -------- #
filename=`basename $seqFile`
extension=${filename##*.}
targetName=${filename%.*}

# ----- create a folder for the results ----#
tmp=tmp/${targetName}_${atomPairTypes}_${boundType}
if [ "$output" != "" ]; then
	tmp=$output
fi
mkdir -p $tmp

# ------ get ld_lib for CNS ----#
#ld_lib=$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=~/anaconda2/lib/:$LD_LIBRARY_PATH
source $cns_solve/cns_solve_env.sh
export CNS_CUSTOMMODULE=$util
export KMP_AFFINITY=none

#================ Step 1. generate extended initial structure ==============#
$bin/SEQ_To_Input $seqFile $tmp/input.seq

currDir=`pwd`

cd $tmp
	$cns_solve/intel-x86_64bit-linux/bin/cns_solve < $util/gseq.inp > gseq.log
	$cns_solve/intel-x86_64bit-linux/bin/cns_solve < $util/extn.inp > extn.log

#cd $FoldingHome
cd $currDir

#================ Step 2. generate all relevant TBL files ==============#
#-> 2.1 distance relevant TBL file

if [[ $distFile = "" || ! -f $distFile ]]; then
	echo "The predicted distance information file does not exist: $distFile"
	exit 1
fi

if [ ${boundType} -lt 0 ]; then
	python $bin/GenDistRestraints4CNS.py -t $boundType -c $Zcutoff -a $atomPairTypes $seqFile $distFile > $tmp/contact.tbl

else
	python $bin/GenDistRestraints4CNS.py -t $boundType -c $distCutoff -g $sigma_value -a $atomPairTypes $seqFile $distFile > $tmp/contact.tbl
fi

if [ $? != 0 ]; then
	echo "Error in generating distance restraints from predicted information"
	exit 1
fi


#-> 2.2 secondary structure relevant TBL files, we need to modify this part later to use predicted angles
##$bin/SSE_To_TBL $seqFile $ssFile $tmp/ssnoe.tbl $tmp/dihedral.tbl $tmp/hbond.tbl

##if [ $? != 0 ]; then
##	echo "Error in generating restraints from secondary structure"
##	exit 1
##fi

#-> dihedral angle restraints and restraints derived from secondary structure including hydrogen-bonding from helix and some near-range distance restraints
## angleFile shall contain predicted Phi/Psi and secondary structure information
##echo 'angleFile=' $angleFile
if [[ $angleFile != "" && -f $angleFile ]]; then
	if [[ $angleFile == *".ss3" ]]; then
		echo "Using predicted secondary structure file: $angleFile"
		$bin/SSE_To_TBL $seqFile $angleFile $tmp/ssnoe.tbl $tmp/dihedral.tbl $tmp/hbond.tbl
		if [ $? != 0 ]; then
			echo "Error in generating restraints from predicted secondary structure"
			exit 1
		fi
	else
		## use predicted local structure property information
		python $bin/GenPropertyRestraints4CNS.py -k $dihedralEnergyConstant -d $tmp/ssnoe.tbl -a $tmp/dihedral.tbl -h $tmp/hbond.tbl $seqFile $angleFile
		if [ $? != 0 ]; then
			echo "Error in generating local restraints from predicted Phi/Psi angles and secondary structures"
			exit 1
		fi
	fi
fi

#================ Step 3. run CNS_Solve for 3D structure folding =======#
#-> 3.1 modify 'dgsa.inp' file
for ((i=1;i<=$numCPUs;i++))
do
	#-> modify 'dgsa.inp' file with given parameters
	$bin/DGSA_File_Mod $seqFile $util/dgsa.inp_raw $tmp/dgsa_${i}.inp $numModels $randSeed ${i}_${targetName} $initNOEScale $AngleScale ${initNOEScale}

	if [ $? != 0 ]; then
		echo "Error in generating dgsa.inp"
		exit 1
	fi

	touch $tmp/iam.running_${i}
done

#-> 3.2 running cns_solve
iamfailed=0

currDir=`pwd`
cd $tmp
	#-> run cns_solve
	for ((i=1;i<=$numCPUs;i++))
	do
		($cns_solve/intel-x86_64bit-linux/bin/cns_solve < dgsa_${i}.inp > dgsa_${i}.log)&
	done
	wait
	#-> check results
	for ((i=1;i<=$numCPUs;i++))
	do
		if [ -f "${i}_${targetName}_${numModels}.pdb" ]
		then
			rm iam.running_${i}
			rm ${i}_${targetName}_*embed*.pdb
		else
			tail -n 30 dgsa_${i}.log
			echo "ERROR: Final 3D models not found for batch $i!"
			mv iam.running_${i} iam.failed_${i}
			iamfailed=1
		fi
	done

#cd $FoldingHome
cd $currDir

#-> 3.3 check failed
if [ $iamfailed -eq 1 ]
then
	echo "CNS FAILED!"
	exit 1
fi


#================ Step 4. collect all models, and return final topk sorted by NOE =======#
#-> 4.1 collect noe in all generated models
rm -f $tmp/${targetName}.noe_statis
numRestraints=`wc -l $tmp/contact.tbl | cut -f1 -d' ' `
for ((i=1;i<=$numCPUs;i++))
do
	for ((k=1;k<=${numModels};k++))
	do
		NOE=`grep "REMARK noe     =" $tmp/${i}_${targetName}_${k}.pdb | awk '{print $NF}'`
		avgNOE=`echo "scale=4; $NOE / $numRestraints" | bc `
		echo "${i}_${targetName}_${k}.pdb $avgNOE" >> $tmp/${targetName}.noe_statis

	done
done

#-> 4.2 sort these models by average NOE in ascending order, and return topk =====#
sort -g -k2 $tmp/${targetName}.noe_statis > $tmp/${targetName}.noe_statis_sort

cp /dev/null $tmp/${targetName}.model_summary.txt
for ((i=1;i<=$topk;i++))
do
	file=`head -n $i $tmp/${targetName}.noe_statis_sort | tail -n1 | awk '{print $1}'`
	avgNOE=`head -n $i $tmp/${targetName}.noe_statis_sort | tail -n1 | awk '{print $NF}'`
	sed -i -e '1iREMARK Generated by RaptorX-Contact. Consult jinboxu@gmail.com for method details.\' $tmp/$file
	cp $tmp/$file $tmp/${targetName}_model$i.pdb
	gzip < $tmp/$file > $tmp/${targetName}_model$i.pdb.gz
	echo ${targetName}_model$i.pdb $avgNOE >> $tmp/${targetName}.model_summary.txt
done

## for each top model, we calculate the NOE violation at each residue

for ((i=1;i<=$topk;i++))
do
	python $FoldingHome/scripts/CalcNOEViolationPerResidue.py $seqFile $tmp/${targetName}_model$i.pdb $tmp/contact.tbl
	mv ${targetName}_model$i.localQuality.txt $tmp/
done

#====== remove temporary files ==========#
rm -f $tmp/*.log

# ---- exit ----- #
#exit 0
