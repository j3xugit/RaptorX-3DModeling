#!/bin/sh

GPU=-1

if [ $# -lt 2 ]; then
	echo $0 "seqFile inputFeature_PKL [gpu]"
	echo "   this script takes an input sequence (FASTA format) and its feature file in PKL format and predicts inter-residue distance"
	echo "   the results will be saved into subfolders in the same directory as inputFeature_PKL"
	echo "   gpu can be -1 (default), cuda0, cuda1, cuda2 and cuda3. If set to -1, automatically select a GPU"
	exit 1
fi

seqFile=$1
inputFeature=$2

if [ $# -ge 3 ]; then
	GPU=$3
fi

if [ $GPU != "-1" ]; then
        GPU=`echo -n $GPU | tail -c 1 `
fi


if [[ -z "${DL4DistancePredHome}" ]]; then
        echo "Please set the environmental variable DL4DistancePredHome, e.g., $HOME/3DModeling/DL4DistancePrediction2"
        exit 1
fi

if [ ! -f $seqFile ]; then
	echo "Cannot find the input seq file for distance prediction: " $seqFile
	exit 1
fi

if [ ! -f $inputFeature ]; then
        echo "Cannot find the input feature file for distance prediction: " $inputFeature
        exit 1
fi

if [ $GPU == "-1" ]; then
        neededRAM=`$DL4DistancePredHome/Scripts/EstimateGPURAM4DistPred.sh $seqFile`
        GPU=`$ModelingHome/Common/FindOneGPUByMemory.sh $neededRAM 120`
        #bestGPU=`nvidia-smi --query-gpu=index,memory.free --format=csv | tail -n +2 | cut -f1,2 -d' ' | sort -k2,2 -rn | head -1 | cut -f1 -d',' `
        #GPU=cuda$bestGPU
fi

if [ $GPU == "-1" ]; then
        echo "ERROR: cannot find an appropriate GPU to run distance prediction for $seqFile !"
	## we may run distance prediction on CPU ??
        exit 1
else
        GPU=cuda$GPU
fi

$DL4DistancePredHome/Scripts/RunServerDistancePred4OneProtein.sh EC25CSet10820N11410Atom15Models $inputFeature $GPU
if [ $? -ne 0 ]; then
        echo "Failed to predict distance using EC25CSet10820N11410Atom15Models for $inputFeature"
        exit 1
fi

#(../RunDistancePred4OneProtein.sh EC52CSet10820N11410MEAtom17Models $inputFeature ${GPU0}; ../RunDistancePred4OneProtein.sh EC25CSet10820Atom8Models $inputFeature ${GPU0}; ../RunDistancePred4OneProtein.sh EC25CSet11410Atom7Models $inputFeature ${GPU0}; ../RunDistancePred4OneProtein.sh EC25CSet10820N11410Atom15Models $inputFeature ${GPU0} ) &
#(../RunDistancePred4OneProtein.sh EC52CSet10820MEAtom8Models $inputFeature ${GPU1};  ../RunDistancePred4OneProtein.sh EC52CSet10820Atom8Models $inputFeature ${GPU1}; ../RunDistancePred4OneProtein.sh EC52CSet10820Atom10Models $inputFeature ${GPU1} )&
#(../RunDistancePred4OneProtein.sh EC52CSet11410Atom10Models $inputFeature ${GPU2}; ../RunDistancePred4OneProtein.sh EC52CSet11410Atom14Models $inputFeature ${GPU2}; ../RunDistancePred4OneProtein.sh EC52CSet11410MEAtom9Models $inputFeature ${GPU2} )&
#(../RunDistancePred4OneProtein.sh EC52CSet10820N11410Atom18Models $inputFeature ${GPU3}; ../RunDistancePred4OneProtein.sh EC52CSet7952AdamSGNACbCbModels $inputFeature ${GPU3}; ../RunDistancePred4OneProtein.sh EC52CSet7952DilatedCbCbModels $inputFeature ${GPU3}; ../RunDistancePred4OneProtein.sh EC52CSet7952AdamWCbCbModels $inputFeature ${GPU3} ) &
#wait

