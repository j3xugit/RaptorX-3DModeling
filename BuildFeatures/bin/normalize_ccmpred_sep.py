#!/usr/bin/env python

import os
import sys
import math


def save_contact(outfile,contact_map):
	outfile_handle = open(outfile,'w')
	for (u,v) in contact_map.items():
		outstring =  str(u[0]) +' '+ str(u[1]) +' '+ str(v)
		outfile_handle.writelines(outstring+'\n')
	outfile_handle.close()

def cal_mean(sorted_list):
	if not sorted_list:
		return 0

        mean = 0
        for (u,v) in sorted_list.items():
                mean += v
        return mean/len(sorted_list)

def cal_sd(sorted_list,mean):
	if not sorted_list:
		return 0
	if len(sorted_list)<2:
		return 0
        sd = 0
        for (u,v) in sorted_list.items():
                sd += (v-mean) * (v-mean)
        return math.sqrt(sd/(len(sorted_list)-1))

def cal_zscore(sorted_list):
        mean = cal_mean(sorted_list)
        sd = cal_sd(sorted_list,mean)
        zscore_map = {}
        for (u,v) in sorted_list.items():
		if sd < 0.000001:
			zscore_map[u] = 0
		else:
                	zscore_map[u] = (v-mean)/sd
        return zscore_map

def read_ccmpred_feature(ccmpred_file):
        file_handle = open(ccmpred_file)
        short_ccmpred_map = {}
        medium_ccmpred_map = {}
        long_ccmpred_map = {}
        index = 0
        protein_length = -1

        for f in file_handle:
                line = f.rstrip().split()
                dim = len(line)
                protein_length = dim
                for i in range(index+6,dim):
                        if index > i:
                                index,i = i,index
                        single_score = (index,i)

                        if math.fabs(index-i) >= 6 and math.fabs(index-i) < 12:
                                short_ccmpred_map[single_score] = float(line[i])
                        if math.fabs(index-i) >= 12 and math.fabs(index-i) < 24:
                                medium_ccmpred_map[single_score] = float(line[i])
                        if math.fabs(index-i) >= 24:
                                long_ccmpred_map[single_score] = float(line[i])

                index = index + 1

        file_handle.close()

        return short_ccmpred_map,medium_ccmpred_map,long_ccmpred_map,protein_length


# --------- program begin -------------- #
if len(sys.argv) != 2:
	print "command is wrong"
	sys.exit(-1)

ccmpred_file = sys.argv[1]

ccmpred_map_short,ccmpred_map_medium,ccmpred_map_long,protein_length = read_ccmpred_feature(ccmpred_file)

short_zscore_map = cal_zscore(ccmpred_map_short)
medium_zscore_map = cal_zscore(ccmpred_map_medium)
long_zscore_map = cal_zscore(ccmpred_map_long)

output_matrix = [[0 for col in range(protein_length)] for row in range(protein_length)]

for (u,v) in short_zscore_map.items():
	output_matrix[u[0]][u[1]] = v
	output_matrix[u[1]][u[0]] = v

for (u,v) in medium_zscore_map.items():
        output_matrix[u[0]][u[1]] = v
        output_matrix[u[1]][u[0]] = v

for (u,v) in long_zscore_map.items():
        output_matrix[u[0]][u[1]] = v
        output_matrix[u[1]][u[0]] = v

for i in range(protein_length):
	for j in range(protein_length):
		print output_matrix[i][j],
	print ''

#print short_zscore_map
