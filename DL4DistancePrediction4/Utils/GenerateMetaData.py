import os
import sys
import json
import glob
import numpy as np
import random

from Common import SequenceUtils

##this script generates a meta data file in json format for training/validation/test data
## groupFile: a list of groups, each line represents one group
## each group has multiple related proteins (e.g., proteins in the same family and/or the same query protein with different template-based information) associated with weights
## All the proteins share SeqDir, NativeDir, a list of inputFeature folders and optionally AliDir and TemplateDir. 

## groupFile has the following format:
## g1 p1 p2 p3 ...
## g2 p4 p5 p6 ...
## here g represents group name, and p represents one query protein or a query-template pair. p can also be associated with a weight, e.g., p1+w1. By default, the weight is set to 1
## when one group does not have any p, the group name is used as the protein name

def Usage():
	print 'python GenerateMetaData.py groupFile inputInfoFile [flag]'
	print '	groupFile: see the top comment of this script for format'
	print '	inputInfoFile: a specification file containing all necessary file path information'
	print '	the resultant file is saved as XXX.json where XXX is the basename of groupFile'
	print '	flag (default 0) is used to indicate if templates shall be used and/or if the data is for train/validation'
	print '		the last/rightmost bit is set to 1 for validation data and 0 for train data'
	print '		the last second bit is set to 1 to indcate that template information shall be used'
	print '		3: template is used and the data is for validation'
	print '		2: template is used and the data is for train'
	print '		1: template is not used and the data is for validation'
	print '		0: template is not used and the data is for train'

def ParseSpecFile(specFile, forValidation=False):
	with open(specFile, 'r') as fh:
		content = [ c.strip() for c in list(fh)]

	metaData = dict()

	for row in content:
		if len(row)<2:
			continue
		if row.startswith('#'):
			continue
		fields = row.split('=')
		if len(fields) < 2:
			continue

		k = fields[0]
		if k.upper() == 'SeqDir'.upper():
			seqDir = fields[1]
			metaData['seqDir'] = os.path.realpath(seqDir)

		elif k.upper() == 'NativeDir'.upper():
			nativeDir = fields[1]
			metaData['nativeDir'] = os.path.realpath(nativeDir)

		elif k.upper() == 'DistDir'.upper():
			distDir = fields[1]
			metaData['distDir'] = os.path.realpath(distDir)

		elif k.upper() == 'oriDir'.upper():
			oriDir = fields[1]
			metaData['oriDir'] = os.path.realpath(oriDir)

		elif k.upper() == 'inDirs'.upper():
			inDirs = fields[1].split()
			if forValidation:
				inDirs = [random.choice(inDirs)]
			metaData['featureFolders'] = [ os.path.realpath(p) for p in inDirs ]

		elif k.upper() == 'aliDir'.upper():
			aliDir = fields[1]
			metaData['aliDir'] = os.path.realpath(aliDir)

		elif k.upper() == 'tplDir'.upper():
			tplDir = fields[1]
			metaData['tplDir'] = os.path.realpath(tplDir)

		elif k.upper() == 'pdbDir'.upper():
			pdbDir = fields[1]
			metaData['pdbDir'] = os.path.realpath(pdbDir)
		else:
			print 'WARNING: unrecognized keyword: ', k
			continue

	return metaData


## groupFile has the following format:
## g1 p1 p2 p3 ...
## g2 p4 p5 p6 ...
## here g represents group name, and p represents one query protein or a query-template pair. p can also be associated with a weight, e.g., p1+w1. By default, the weight is set to 1
## when one group does not have any p, the group name is used as the protein name

## MetaData is a dict() with the following keys: groups, SeqDir, DistDir, OriDir and a set of infeatureDir
## MetaData['groups'] is a list of groups, each group is a dict() with keys: group name, proteins, sequences, weights and templates
## when not None, templates are used to indicate if template-based information shall be used for this protein or not

def ParseGroupFile(groupFile, metaData, forValidation=False):
	SeqDir = metaData['seqDir']

	groups = []
	with open(groupFile) as fh:
		content = [ line.strip() for line in list(fh) ]

	## parse content
	for row in content:
		if row.startswith('#'):
			continue

		fields = row.split()
		if len(fields)<1:
			continue

		## when the groupFile is for validation and there are at least 2 entries in this group, then we randomly choose one
		## by doing so, we fix the validation data so that validation loss is comparable when multiple models are trained
		if forValidation and len(fields)>2:
			fields = [ fields[0], random.choice(fields[1:]) ]

		group = dict()
		group['name'] = fields[0]

		if len(fields) == 1:
			group['proteins'] = fields
			group['weights'] =[1]
			group['templates'] = [None]

		else:
			proteinPairs = []
			proteins = []

			for f in fields[1:]:
				tmp = f.split('+')
				if len(tmp) == 2:
					weight = np.float32(tmp[1])
				else:
					weight = 1
				seqtemplate = tmp[0].split('-')
				seq = seqtemplate[0]
				if seq not in proteins:
					proteins.append(seq)
				if len(seqtemplate) == 2:
					template = seqtemplate[1]
				else:
					template = None
				proteinPairs.append( (seq, template, weight) )

			## recalculate weight for each protein by summing the weights of all protein pairs belonging to one protein
			## reorganize all templates of each protein
			weights = []
			templates = []
			for p in proteins:
				weight = [ w for seq, temp, w in proteinPairs if seq==p ]
				weights.append( sum(weight) )

				template = [ temp for seq, temp, w in proteinPairs if seq==p and temp is not None ]
				templates.append(template)

			group['templates'] = templates
			group['proteins'] = proteins
			group['weights'] = weights

	
		group['sequences'] =[]
		for p in group['proteins']:
			seqFile = os.path.join(SeqDir, p + '.fasta')
			if not os.path.isfile(seqFile):
				seqFile = os.path.join(SeqDir, p + '.seq')
	
			if not os.path.isfile(seqFile):
				print 'ERROR: seq file does not exist: ', seqFile
				exit(1)

			seq = SequenceUtils.LoadFASTAFile(seqFile)
			group['sequences'].append(seq)

		groups.append(group)
	return groups

def CheckFiles4Templates(metaData):
	tplDir = metaData['tplDir']
	aliDir = metaData['aliDir']

	assert tplDir is not None
	assert aliDir is not None

	if not os.path.isdir(tplDir):
		print 'ERROR: invalid tplDir: ', tplDir
		exit(1)

	if not os.path.isdir(aliDir):
		print 'ERROR: invalid aliDir: ', aliDir
		exit(1)

	allalifiles = glob.glob(aliDir + '/*.fasta')
	allalifiles = [ os.path.basename(f) for f in allalifiles if os.path.getsize(f)>0 ]
	allalifiles = set(allalifiles)

	alltplfiles = glob.glob(tplDir + '/*.tpl.pkl')
	alltplfiles = [ os.path.basename(f) for f in alltplfiles if os.path.getsize(f)>0 ]
	alltplfiles = set(alltplfiles)

	missingPool = []
	for group in metaData['groups']:
		for protein, templateList in zip(group['proteins'], group['templates']):
			if not bool(templateList):
				continue

			for template in templateList:
				tplFile = template + '.tpl.pkl'
				if tplFile not in alltplfiles:
					missingPool.append( tplFile )

				alnfile1 = protein + '-' + template + '.fasta'
				alnfile2 = template + '-' + protein + '.fasta'
				if (alnfile1 not in allalifiles) and (alnfile2 not in allalifiles):
					missingPool.append( (alnfile1, alnfile2) )

	if len(missingPool) < 1:
		print 'Great! all files needed for templates and alignments exist.'
	else:
		print 'ERROR: some files for templates and alignments are missing. The first 10 files are: '
		print missingPool[:10]


## check file exists
def CheckFiles(proteins, fileSuffix, fileDir):
	allfiles = glob.glob(fileDir + '/*.' + fileSuffix)
	allfiles = [ os.path.basename(f) for f in allfiles if os.path.getsize(f)>0 ]
	allfiles = set(allfiles)

	missingPool = []

	for protein in proteins:
		filename = protein + '.' + fileSuffix
		if filename not in allfiles:
			missingPool.append(os.path.join(fileDir, filename))

	return missingPool
			
## check if the needed files exist or not
def CheckData(metaData):
	groups = metaData['groups']
	proteins = []
	for g in groups:
		proteins.extend(g['proteins'])

	if not os.path.isdir(metaData['nativeDir']):
		print 'ERROR: invalid native dir: ', metaData['nativeDir']
		exit(1)

	missing = []
	missing.extend( CheckFiles(proteins, 'native.pkl', metaData['nativeDir']) )
	for featureDir in metaData['featureFolders']:
		if not os.path.isdir(featureDir):
			print 'ERROR: invalid feature dir: ', featureDir
			exit(1)

		missing.extend( CheckFiles(proteins, 'inputFeatures.pkl', featureDir) )
		missing.extend( CheckFiles(proteins, 'extraCCM.pkl', featureDir) )
		missing.extend( CheckFiles(proteins, 'a2m', featureDir) )

	if len(missing) < 1:
		print 'Great: all files for basic input features and ground truth exist'
		return

	missingFile = 'missing.txt'
	print 'ERROR: some files are missing. Please see ', missingFile

	allMissing='\n'.join(missing)
	with open(missingFile, 'a') as fh:
		fh.write('*************************************************************************\n')
		fh.write(allMissing)

if len(sys.argv) < 3:
	Usage()
	exit(1)

groupFile = sys.argv[1]
specFile = sys.argv[2]

templateUsed = False
forValidation = False

if len(sys.argv) >=4:
	flag = np.int32(sys.argv[3])
	templateUsed = (flag & 2)>0
	forValidation = (flag & 1)>0

metaData = ParseSpecFile(specFile, forValidation=forValidation)
print metaData

metaData['groups'] = ParseGroupFile(groupFile, metaData, forValidation=forValidation)

print 'Checking the existence of files for basic features and ground truth...'
CheckData(metaData)

if templateUsed:
	print 'Checking the existence of files for templates and alignments...'
	CheckFiles4Templates(metaData)

savefile = os.path.basename(groupFile) + '.metaData.json'
with open(savefile, 'w') as fh:
	json.dump(metaData, fh)
	#fh.writelines(json.dumps(metaData, indent=4, sort_keys=True))
