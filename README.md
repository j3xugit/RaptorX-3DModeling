# README #


### What is this repository for? ###

RaptorX predicts protein contact/distance/orientation and local structure properties (e.g, secondary structure and phi/psi angles) by deep convolutional residual networks.
It also builds the 3D model of a protein sequence using predicted distance/orientation and phi/psi angles.
It is mainly tested on the Linux distribution CentOS (>6.0) with bash shell and Python 2.7, but you may also run it if you already have Python 3 installed. 
A version directly supporting Python 3 and tensorflow will be available in a few months. 
This package is also incorporated into our protein structure prediction web server at http://raptorx.uchicago.edu/, which is publicly available for both academia and industry.
If you only want to predict structures for several protein sequences, it is more convenient to use our web server instead of installing this package.

* Version

4.0

* There is no restriction on the usage of this package, but without explicit permission, it shall not be used for commercial purpose and to set up a similar web server for protein folding.

* contact: jinboxu@gmail.com

### How to set up? ###

Download this package by "git clone https://github.com/j3xugit/RaptorX-3DModeling.git" and save it anywhere in your own account, e.g., $HOME/RaptorX-3DModeling/.
It contains the following files and subfolders (and a few others):

BuildFeatures/

DL4DistancePrediction4/

DL4PropertyPrediction/

Folding/

params/

raptorx-path.sh
raptorx-external.sh

README.md

Server/

This RaptorX package consists of 4 major modules: BuildFeatures/ for generating multiple sequence alignment (MSAs) and input features for angle/contact/distance/orientation prediction,
DL4DistancePrediction4/ for contact/distance/orientation prediction, DL4PropertyPrediction/ for local structure property prediction, and Folding/ for building 3D models.

To predict contact/distance/orientation and fold a protein, you may simply run RaptorX-3DModeling/Server/RaptorXFolder.sh, but before doing this some external packages and databases shall be installed and configured.

## Required external packages for all modules ##

1) anaconda or miniconda for Python 2.7

If you have not installed any anaconda or miniconda, you may directly install anaconda or miniconda for Python 2.7.
If you have already installed anaconda/miniconda for Python 3, you may create a virtual enviroment RaptorX by running "conda create --name RaptorX python=2". 
Afterwards, switch to this virtual environment by command "conda activate RaptorX" to install the required packages and run RaptorX. 
See https://docs.python.org/3/tutorial/venv.html for an introduction of virtual environment.

Install numpy by "conda install numpy".

Install msgpack-python by "conda install -c anaconda msgpack-python"; it may not work if you intall it through pip.

2) Biopython (https://biopython.org/)

Needed for both contact/distance/orientation predicton and 3D structure model building.

Install by running "pip install biopython==1.76". Note that a newer version may not work with Python 2.7.

## Required packages for contact/distance/orientation/angle prediction ##

1) Pillow

Needed for visualizing predicted contact and distance; install by running "pip install Pillow"

2) pygpu and Theano 1.0 (http://deeplearning.net/software/theano/install.html)

Needed for train and run deep learning models; install by running "conda install numpy scipy mkl" and then "conda install theano pygpu" .

Please make sure that the CUDA toolkits and CUDNN library have been installed on your machine with GPUs.
Set the environment variable CUDA_ROOT to where cuda is installed, e.g., export CUDA_ROOT=/usr/local/cuda. 
Make sure that the header and lib64 files of CUDNN are in CUDA_ROOT/include and CUDA_ROOT/lib64, respectively. 
We have tested Theano 1.04, CUDA 8 to 10.1 and CUDNN 7 to 7.6.5 . Other versions of CUDA and CUDNN may also work.

In principle, Theano can also run without GPUs, but some minor revisions are needed for some prediction script progtrams in this package.

3) shared_ndarray (https://github.com/crowsonkb/shared_ndarray.git)

Needed for train and run deep learning models for distance/orientation prediction.

Download by "git clone https://github.com/crowsonkb/shared_ndarray.git";
cd to shared_ndarray/ and then run "python setup.py install".

## Required tools and sequence databases for MSA generation ##

1) Install HHblits for MSA generation (https://github.com/soedinglab/hh-suite)

In addition to the HHsuite package itself, please download a sequence database specific to HHsuite and unpack it into a folder, 
e.g. UniRef30_2020_03_hhsuite.tar.gz at http://wwwuser.gwdg.de/~compbiol/uniclust/2020_03/.

2) Install EVcouplings for generating MSAs by jackhmmer (optional, but not recommended since it is too slow)

It is available at https://github.com/debbiemarkslab/EVcouplings . Although the whole package will be installed, only the MSA generation module will be used.

Step 1: Download the package by running "git clone https://github.com/debbiemarkslab/EVcouplings.git". Suppose that it is located at $HOME/EVcouplings. 

Step 2: Run "conda create -n evfold anaconda python=3" to create a virtual environment and then switch to this environment by running "conda activate evfold".

Step 3: cd to $HOME/EVcouplings and run "python setup.py install" to install the whole package.

Note that EVcouplings runs on python 3 while this version of RaptorX runs on python 2.
Without jackhmmer, you may still run RaptorX by generating MSAs using HHblits only.
The sequence database for jackhmmer is uniref90.fasta, which can be downloaded from UniProt.

3) Metagenome data (optional)

Download the data file metaclust_50.fasta at https://metaclust.mmseqs.org/current_release/ and install it somewhere.
In the future, we will install much more metagenome data.

4) (IMPORTANT) Revise the file RaptorX-3DModeling/raptorx-external.sh to setup the path information for the above MSA building tools and databases.

## Install deep learning models for contact/distance/orientation/angle/SS/ACC prediction ##

The deep learning model files for contact/distance/orientation prediction are big (each 100-200M). 
You may download them at https://doi.org/10.5281/zenodo.4710337 or http://raptorx.uchicago.edu/download/ .

1) The package RXDeepModels4DistOri-FM.tar.gz has 6 models for contact/distance/orientation/ prediction. Unpack it and place all the deep model files (ending with .pkl) at $DL4DistancePredHome/models/

2) The package RXDeepModels4Property.tar.gz has 7 deep models for Phi/Psi angle, Secondary Structure (SS) and Solvent Accessibility (ACC) prediction. Unpack it and place all the deep model files (ending with .pkl) at $DL4PropertyPredHome/models/ . By default the package will just predict Phi/Psi angles. If you also want to predict SS and ACC, please use "-m AllSeqSet10820Models" in running the script programs in DL4PropertyPrediction/Scripts/

## Required packages for building protein 3D models from predicted information ##

The below two packages are not needed if you do not want to build 3D models of a protein sequence under prediction. 
That is, they are not needed if you just want to predict SS/ACC/angle and contact/distance/orientation. 

1) PyRosetta (http://www.pyrosetta.org/dow)

It is needed to build 3D models of a protein sequence from predicted distance/orientation and phi/psi angles.
Please download the version supporting Python 2.7 and unpack it. Supposing it is located at PyRosetta4.Release.python27.linux.release-224/,
cd to PyRosetta4.Release.python27.linux.release-224/setup/ and run "python setup.py install".

2) GNU parallel (optional, but recommended for running folding jobs on a Linux workstation)

Some scripts in RaptorX-3DModeling/Folding/ (e.g., ParallelFoldNRelaxOneTarget.sh and SRunFoldNRelaxOneTarget.sh) use GNU parallel to run multiple folding jobs on one or multiple computers.
Run "which parallel" to see if GNU parallel is available or not.
If GNU parallel is not installed, you may still run folding jobs using other scripts.

## Basic Usage

You may run shell script RaptorXFolder.sh in RaptorX-3DModeling/Server/ to predict angle/contact/distance/orientation and/or fold a protein. 
The input can be a protein sequence in FASTA format (ending with .fasta or .seq) or an MSA file in a3m format (ending with .a3m).
In the input file, an amino acid shall be represented by a capital letter instead of a 3-letter code.

NOTE that to run RaptorXFolder.sh, you just need to make sure that the environmental variable CUDA_ROOT and those in RaptorX-3DModeling/raptorx-external.sh are correctly set. 
However, to explicitly run any other shell scripts or code, you will need to set and export all needed environmental variables described in the below section "setup and export of environment variables".

Running RaptorXFolder.sh (and other shell scripts in this package) without any arguments will show its help information. Below are some scenarios.

1) When you already have a multiple sequence alignment in a3m format, please use option "-m 0";

2) When you only want to predict angle/contact/distance/orientation but not 3D models, please use option "-n 0";

3) When you do not want to generate MSAs using jackhmmer, you may use option "-m 9" (without using metagenome data) or "-m 25" (using metagenome data). 
Note that jackhmmer usually is slow, so it is not recommended for MSA generation.

Ideally, RaptorXFolder.sh shall run on a computer with GPUs and a reasonable number of CPUs.
It takes minutes (rarely a couple of hours) to generate MSAs for a protein and several (or at most dozens of) minutes to predict distance/orientation on a single GPU.
The running time for building 3D models from predicted angle/distance/orientation depends on the protein sequence length and whether or not to run relaxation.
When relaxation is not applied, on a single CPU it takes <1 hour to build one 3D model for a protein of 300 residues and 2-3 hours for a protein of 1000 residues. 
However, when relaxation is applied, it may increase the running time by 3 or 4 times. 
Usually one hundred 3D models shall be generated, but this number may be reduced if the predicted distance/orientation is of high quality.
When the protein under prediction is big, GPU and CPU memory may be an issue. 
For a protein of >1000 residues, it may need >12G GPU memory to predict distance/orientation and 10G CPU memory to fold (and relax) one 3D model.


All the result files are saved to a folder target_OUT/ where target is the protein name containing the following subfolders:

1) target_contact/ contains MSAs and input features for distance/orientation prediction;

2) target_thread contains files for Phi/Psi prediction (which are also used for threading);

3) DistancePred/ contains predicted distance/orientation/contact and their visualization;

4) PropertyPred/ contains predicted Phi/Psi angles and secondary structure;

5) target-RelaxResults/ contains all generated decoys;

6) target-SpickerResults/ contains the clustering results of the decoys.

In example/, there is the result of one example protein 1pazA generated by running "Server/RaptorXFolder.sh -o example/ -n 40 -r 1 example/1pazA.fasta".
Note that to save disk space, some resultant files are compressed or deleted.

It is possible to run RaptorXFolder.sh on several machines, each in charge of one major module (MSA generation, feature generation and distance/orientation prediction, 3D model building),
without requiring you to manually copy data among different machines. This will be explained in section "Advanced Usage".

## References

1. improved protein structure prediction by deep learning irrespective of co-evolution information. Nature Machine Intelligence, 2021. Also available at https://www.biorxiv.org/content/10.1101/2020.10.12.336859v1

2. Distance-based protein folding powered by deep learning. PNAS, August 2019. A 2-page abstract also appeared at RECOMB2019 in April 2019.

3. Analysis of distance-based protein structure prediction by deep learning in CASP13. PROTEINS, 2019.

4. Accurate De Novo Prediction of Protein Contact Map by Ultra-Deep Learning Model. PLoS CB, Jan 2017

5. Folding Membrane Proteins by Deep Transfer Learning. Cell Systems, September 2017.

## Setup and export of environment variables ##

1) Set the environemntal variable ModelingHome to the install folder of the whole package, e.g., $HOME/RaptorX-3DModeling.
2) Add ". $ModelingHome/raptorx-path.sh" to the .bashrc file in your own Linux account.
3) Revise the setting in $ModelingHome/raptorx-external.sh and add ". $ModelingHome/raptorx-external.sh " to the .bashrc file in your own Linux account to set environmental variables related to external databases and tools.

Supposing that the RaptorX-3DModeling package is located at $HOME/RaptorX-3DModeling/,
below is an example configuration that can be pasted to the .bashrc file (when your account is using the bash shell).

export CUDA_ROOT=/usr/local/cuda/

export ModelingHome=$HOME/RaptorX-3DModeling/

. $ModelingHome/raptorx-path.sh

. $ModelingHome/raptorx-external.sh

If your account uses csh shell, you may add a similar setting to the file .cshrc in your home directory. 
The settings in raptorx-external.sh and raptorx-path.sh shall also be revised according to the grammar of csh.
Most shell scripts in this package have been tested with bash shell.

## Detailed Usage

* How to generate multiple sequence alignments (MSA)

HHblits and Jackhmmer are two popular tools for protein sequence homology search.
Enclosed in this package (located in folder BuildFeatures/) there are some scripts that call HHblits and Jackhmmer to build MSAs.

1) To generate MSA for phi/psi prediction only, run "BuildMSAs.sh -d ResultDir -m 1 SeqFile" where ResultDir is the folder for result saving.
Helpers/BuildMSA4Threading.sh is another script for this purpose.

2) To generate a single MSA from a protein sequence for contact/distance prediction, you may use BuildFeatures/HHblitsWrapper/BuildMSA4DistPred.sh and/or BuildFeatures/EVAlign/BuildMSAByJack.sh .
Note that BuildFeatures/EVAlign/BuildMSAByJack.sh is usually much slower than BuildMSA4DistPred.sh.

3) To generate MSAs for a single protein for contact/distance/orientation prediction, you may use BuildFeatures/BuildMSAs.sh. 
This script may generate multiple MSA files for a single protein depending on your input options.
For example, to generate MSA for phi/psi and contact/distance prediction, you may run "BuildMSAs.sh -d ResultDir -m 9 SeqFile" or "BuildMSAs.sh -d ResultDir -m 25 SeqFile"
The ResultDir contains two subfolders XXX_contact and XXX_thread where XXX is the protein name. XXX_contact has MSA files for contact/distance prediction and XXX_thread has an MSA file for phi/psi prediction.
By default BuildMSAs.sh will not use jackhmmer to generate MSAs since it is slow, but you may enable it by adding 4 to the option value of "-m".

4) To generate MSAs for multiple proteins, you may use BatchBuildMSAs.sh. A set of MSAs will be generated for an individual protein.
By default, BatchBuildA3M.sh will not use jackhmmer to generate MSAs since it is slow. 

5) To generate input features from one MSA for distance/orientation prediction, you may run BuildFeatures/GenDistFeaturesFromMSA.sh
At least one GPU is needed to run CCMpred efficiently. Otherwise it may take a long time to generate features.

6) To directly generate input features for contact/distance/orientation prediction from protein sequences,
you may run BuildFeatures.sh for a single protein and BatchBuildFeatures.sh for multiple proteins, respectively.
These two scripts call BuildMSAs.sh to build MSAs and then derive input features from MSAs.
At least one GPUs is needed to efficiently run CCMpred. When there are no GPUs, it may take a long time on CPUs.

* How to predict contact/distance/orientation from MSAs or input features

Some scripts such as PredictPairRelationRemote.sh and PredictPairRelation4Server.sh may run on a computer without GPUs as long as you may ssh/scp/rsync (without password) to a remote computer with GPUs with RaptorX installed.

1) To predict contact/distance/orientation from a single MSA file in a3m format, you may use DL4DistancePrediction4/Scripts/PredictPairwiseRelationFromMSA.sh
Note that the 1st sequence in the A3M file shall be the query sequence without any gaps.
One result file XXX.predictedDistMatrix.pkl will be generated that can be opened by cPickle.

2) To predict contact/distance/orientation from a list of MSA files (each protein has one MSA file in this list), first use BuildFeatures/BatchGenDistFeaturesFromMSAs.sh to generate input feature files from all MSA files
and then run DL4DistancePrediction4/Scripts/PredictPairwiseRelation4Inputs.sh. This is much faster than running PredictPairwiseRelationFromMSA.sh on each MSA file separately.

3) To predict contact/distance/orientation from several MSAs generated by BuildMSAs.sh for a single protein, first run BuildFeatures/GenDistFeatures4OneProtein.sh and then DL4DistancePrediction4/Scripts/PredictPairwiseRelation4OneProtein.sh

4) To predict contact/distance/orientation for multiple proteins, each of which has a set of MSAs generated by BuildMSAs.sh, first run BuildFeatures/GenDistFeatures4OneProtein.sh or GenDistFeatures4MultiProteins.sh 
and then DL4DistancePrediction4/Scripts/PredictPairwiseRelation4Proteins.sh 

5) To print out contact matrix from predicted distance/orientation files, use PrintContactPrediction.sh or BatchPrintContactPrediction.sh in DL4DistancePrediction4/Scripts/

Note that the input feature files generated by GenDistFeaturesFromMSAs.sh and BatchGenDistFeaturesFromMSAs.sh and the predicted distance/orientation file may be very large for a large protein.
To save disk space, please avoid generating input feature files or predicting distance/orientation for too many proteins in a single batch.

* How to predict protein local structure properties such as phi/psi angles and secondary structure

Several scripts in /home/jinbo/RaptorX-3DModeling/DL4PropertyPrediction/Scripts/ can be used, e.g., PredictPropertyFromMSA.sh, PredictPropertyFromHHMs.sh, PredictProperty4OneProtein.sh and PredictProperty4Proteins.sh. Some scripts such as PredictPropertyRemote.sh and PredictProperty4Server.sh may run on a remote computer with GPUs and to which you may ssh/scp without password.

* How to build 3D models from predicted angle/distance/orientation files

1) To fold a protein from its primary sequence or an MSA without manually running the intermediate steps, run RaptorX-3DModeling/Server/RaptorXFolder.sh

2) When you already have predicted distance/orientation files (ending with .predictedDistMatrix.pkl) and phi/psi file (ending with .predictedProperties.pkl), you may use them to fold a protein by several scripts. For example, in RaptorX-3DModeling/Folding/, there are LocalFoldNRelaxOneTarget.sh, ParallelFoldNRelaxOneTarget.sh, SRunFoldNRelaxOneTarget.sh and SlurmFoldNRelaxOneTarget.sh, developed for different machine types (e.g., Linux workstation and slurm cluster). In aptorX-3DModeling/Folding/Scripts4Rosetta/, there are FoldNRelaxOneTarget.sh, FoldNRelaxTargets.sh, and RelaxOneTarget.sh, which mainly run on a Linux workstation.

## Advanced Usage 

* Run scripts on several machines without manually copying data

Suppose that you have access to three machines: the 1st one has a small number of CPUs but not any GPUs, the 2nd one has GPUs but very few CPUs, and the third one has many CPUs but not GPUs. 
You may start RaptorXFolder.sh on the 1st machine, which will then automatically run GPU tasks on the 2nd machine and the folding tasks on the 3rd machine. 
During this process, you do not need to manually copy files and results among machines. 
To fullfil this, you shall install and configure this RaptorX package (or a portion of it) on the three machines so that on the 1st one you may run MSA generation,
on the 2nd one you can run GPU tasks (e.g., CCMpred and distance/orientation prediction) and on the 3rd one you may run 3D model building. All related environmental variables shall be correctly set on all machines.

To run GPU tasks at a remote machine, please create one file (e.g., GPUMachines.txt) to specify the remote machines with GPUs and to which you may ssh/scp/rsync without password. 
See an example file in RaptorX-3DModeling/params/. A line in this file looks like "raptorx9.uchicago.edu LargeRAM on" or "jinbo@aptorx7.uchicago.edu SmallRAM off" 
where the three fields are the computer name (and account name if needed), GPUs of a small RAM (<=12G) or a large RAM, and enabled/disabled, respectively.
You may save this file at a default location (i.e., RaptorX-3DModeling/params/GPUMachines.txt). It is also possible to save this file at other places, but you will need to slightly modify RaptorXFolder.sh . 

Note that once RaptorX-3DModeling/params/GPUMachines.txt exists, in order to use your local GPUs, you shall add one line to this file for your local GPUs, e.g., "name-of-your-local-machine SmallRAM on" or "name-of-your-local-machine LargeRAM on".

To run folding jobs at a remote machine, you just need to specify a remote account through the -R option of RaptorXFolder.sh. Again please make sure that you may ssh/scp to this remote account without password.  

## Test

1) Generate input feature files for contact/distance prediction from an MSA file in a3m format:

BuildFeatures/GenDistFeaturesFromMSA.sh -o Test_Feat BuildFeatures/example/1pazA.a3m

Three feature files shall be generated in Test_Feat/: 1pazA.inputsFeatures.pkl, 1pazA.extraCCM.pkl and 1pazA.a2m. 
Meanwhile, the first two .pkl files are needed for contact/distance/orientation prediction and 1pazA.a2m may be needed by very few deep models. 

2) Predict distance/orientation from feature files:

DL4DistancePrediction4/Scripts/PredictPairwiseRelation4OneInput.sh -d ./Test_Dist Test_Feat/1pazA.inputsFeatures.pkl

The result file Test_Dist/1pazA.predictedDistMatrix.pkl should be generated.

3) Generate predicted contact matrix in text format:

To print out predicted contact matrix in text format, run "DL4DistancePrediction4/Scripts/PrintContactPrediction.sh Test_Dist/1pazA.predictedDistMatrix.pkl".
This will generate two text files 1pazA.CASP.rr and 1pazA.CM.txt.


## How to train a new model (to be updated) ##

To train a contact/distance prediction deep network, please go to $DL4DistancePredHome/Work/ and run the shell script TrainCbCbEC25CL51-Adam.sh.
Before running the shell script, please make sure that you set the path of the training and test data correctly.
On raptorx4, raptorx5 and raptorx6 machines, you may simply link the data path to your own work directory. Please read TrainCbCbEC25CL51-Adam.sh for instructions. 
Please do not use cuda0 on thee three machines to run the training algorithm. They are reserved for our web server.


### Develop your own deep network (to be updated) ###


To incorporate your own deep network architeture into the contact/distance prediction module, please follow the below procedure:


1) In Model4PairwisePrediction.py, find the sentence starting with "matrixConv=ResNet". 
Let XXYY denote your own network, then you may change this sentence to the following code.


if modelSpecs['network'] == 'XXYY':


        matrixConv=XXYY(....)


else:


        matrixConv=ResNet(rng, input=input_2d, n_in=n_input2d, n_hiddens=n_hiddens_matrix, n_repeats=matrix_repeats, halfWinSize=hwsz_matrix, mask=mask_matrix, activation=modelSpecs['activation'], batchNorm=modelSpecs['batchNorm'])


Note that your implementation of XXYY shall have the same input and output format as our ResNet implementation.
XXYY shall also contain the following variables: output, n_out, params, paramL2 and paramL1. Meanwhile, params, paramL2 and paramL1 are the list of model parameters and their norms.
output and n_out in XXYY shall be consistent with output and n_out in ResNet.
In addition, it is better to use the BatchNormLayer class in ResNet4Distance.py and the convolution layer we implemented.
If you want to use other batch normalization implementation, please make sure 1) calculate mean and standard deviation for each protein instead of each minibatch; 2) remove the impact of the zero-padding.
Please read our batch normalization code for details.


2) In config.py, please add 'XXYY' to allNetworks, which is a list of all allowed network architectures.


3) If needed, in InitializeModelSpecs() (of config.py), please add some code to set the default architecture parameters for your own network.


4) If needed, please revise ParseCommandLine.py to read in your own network architecture information from the shell script Work/TrainCbCbEC25L51-Adam.sh


5) Revise Work/TrainCbCbEC25L51-Adam.sh accordingly.


Currently two minibatches may have different protein lengths. In the same minibatch, all the proteins are aligned at the right bottom corner and zero padded to have the same length.
We use a mask matrix to indicate the zero-padding pattern at the left and top of a contact matrix.
If you use pooling layers, please make sure that you change the mask matrix correctly after every pooling.
The mask matrix has a smaller shape than the contact matrix. For the 2D convolution, the mask matrix has shape (batchSize, maxProteinLen, maxProteinLen-minProteinLen) where maxProteinLen and minProteinLen are the maximum and minimum protein lengths in a minibatch. For the 1D convolution, the mask matrix has shape (batchSize, maxProteinLen - minProteinLen).


### Who do I talk to? ###


Jinbo Xu at jinboxu@gmail.com
