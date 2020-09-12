#!/bin/bash

## This script installs evcouplings for jackhmmer. Looks like that evcouplings works with only python 3.5, but the latest Biopython only supports python 3.6.
## Need to figure out how to fix this issue later.

## We assume the installation folder of anaconda or miniconda is $HOME/anaconda2 or $HOME/miniconda2, respectively.
## if not, please change the below line manually
AnacondaDir=$HOME/anaconda2/
if [ ! -d $AnacondaDir ]; then
	AnacondaDir=$HOME/miniconda2/
fi

if [ ! -d $AnacondaDir ]; then
	echo "ERROR: cannot find a valid folder for anaconda or miniconda. Please make sure one of them is installed at $HOME"
	exit 1
fi

# remove environment
conda remove --name evfold --all -y

# create environment
conda create --name evfold python=3.5 -y
source activate evfold

# pip install evcouplings
pip install evcouplings

# copy modified file
if [ ! -f util/protocol.py ]; then
	echo "Failed in the last step of installation. Cannot find util/protocol.py "
	exit 1
fi

cp util/protocol.py $AnacondaDir/envs/evfold/lib/python3.5/site-packages/evcouplings/align/



