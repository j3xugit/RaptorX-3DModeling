#!/usr/bin/env python

import numpy as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import sys
import os


if len(sys.argv) != 2:
    print 'Usage: PlotContactMapComparison.py protein_name contactMap1 contactMap2 groundTruth'
    print '	contactMap1 or contactMap2: a text file containing a predicted contact map. Each entry is a prob value or confidence score'
    print '	groundTruth: a text file containing the native distance matrix. Each entry is a Cb-Cb distance calculated from native structure'
    sys.exit(1)

protein_name = sys.argv[1]
predFile1 = sys.argv[2]
predFile2 = sys.argv[3]
truth = sys.argv[4]

"""
our_prediction = np.loadtxt("./" + protein_name + ".gcnn")
our_prediction_tmp = our_prediction
meta_prediction = np.loadtxt("./" + protein_name + ".ccmpred")
meta_prediction_tmp = meta_prediction
ground_truth = np.loadtxt("./" + protein_name + ".distcb")
"""

our_prediction_tmp = np.loadtxt(predFile1)
meta_prediction_tmp = np.loadtxt(predFile2)
ground_truth = np.loadtxt(truth)

L = our_prediction.shape[0]

ground_truth1 = (ground_truth > 0).astype(int)
ground_truth2 = (ground_truth <=8).astype(int)
ground_truth3 = ground_truth1 * ground_truth2
ground_axis = np.argwhere(ground_truth3 == 1)


M1s = np.ones_like(ground_truth3, dtype = np.int8)
mask_LR = np.triu(M1s, 6)
mask_LR1 = np.tril(M1s, -6)
mask_MLR = np.triu(M1s, 12)
mask_SMLR = np.triu(M1s, 6)

#topk=L/2
topk=L

#-> our method
our_prediction_tmp1 = our_prediction_tmp * mask_LR
our_prediction = np.triu(our_prediction_tmp1, 1)   # pick up the upper triangle
our_prediction_flatten = our_prediction.flatten() # flatten the upper triangle to sequence
our_prediction_flatten_sort = np.argsort(-our_prediction_flatten) # sort the sequence from max to min and obtain the index
our_prediction_flatten_mask = np.zeros_like(our_prediction_flatten) #
our_prediction_top = our_prediction_flatten_sort[0:topk]               # pick up the top L as 1, other as 0, maybe add a mask
our_prediction_flatten_mask[our_prediction_top] = 1                   # here for long-media-short range prediction
our_prediction_mask = our_prediction_flatten_mask.reshape((L, L))
our_prediction_mask_correct = our_prediction_mask * ground_truth3
our_prediction_mask_wrong = np.triu(np.abs(our_prediction_mask - our_prediction_mask_correct), 1)
our_prediction_mask_correct1 = our_prediction_mask_correct * mask_LR
our_prediction_mask_wrong1 = our_prediction_mask_wrong * mask_LR
our_predict_axis_correct = np.argwhere(our_prediction_mask_correct1 == 1)
our_predict_axis_wrong = np.argwhere(our_prediction_mask_wrong1 == 1)

#-> metaPSICOV
meta_prediction_tmp1 = meta_prediction_tmp * mask_LR1
meta_prediction = np.tril(meta_prediction_tmp1, -1)
meta_prediction_flatten = meta_prediction.flatten()
meta_prediction_flatten_sort = np.argsort(-meta_prediction_flatten)
meta_prediction_flatten_mask = np.zeros_like(meta_prediction_flatten)
meta_prediction_top = meta_prediction_flatten_sort[0:topk]
meta_prediction_flatten_mask[meta_prediction_top] = 1
meta_prediction_mask = meta_prediction_flatten_mask.reshape((L, L))
meta_prediction_mask_correct = meta_prediction_mask * ground_truth3
meta_prediction_mask_wrong = np.tril(np.abs(meta_prediction_mask - meta_prediction_mask_correct), -1)
meta_prediction_mask_correct1 = meta_prediction_mask_correct * mask_LR1
meta_prediction_mask_wrong1 = meta_prediction_mask_wrong * mask_LR1
meta_predict_axis_correct = np.argwhere(meta_prediction_mask_correct1 == 1)
meta_predict_axis_wrong = np.argwhere(meta_prediction_mask_wrong1 == 1)


# ------- draw PNG ---------- #
plt.figure(1, figsize=(8, 8))
plt.axis([0, L, 0, L])
plt.scatter(ground_axis[:,0], ground_axis[:,1], marker = 'o', color = 'grey', s =3)

plt.hold(True)
plt.scatter(our_predict_axis_correct[:,0], our_predict_axis_correct[:,1], marker = '*', color = 'r', s = 10)
plt.hold(True)
plt.scatter(our_predict_axis_wrong[:,0], our_predict_axis_wrong[:,1], marker = 'x', color = 'g', s = 10)
plt.annotate('Our', xy=(20, 1.9*L), xycoords='axes points', horizontalalignment='left', verticalalignment='upper', fontsize=20)

plt.hold(True)
plt.scatter(meta_predict_axis_correct[:,0], meta_predict_axis_correct[:,1], marker = '*', color = 'r', s = 10)
plt.hold(True)
plt.scatter(meta_predict_axis_wrong[:,0], meta_predict_axis_wrong[:,1], marker = 'x', color = 'g', s = 10)
plt.annotate('CCMpred', xy=(2*L-20, 0.1*L), xycoords='axes points', horizontalalignment='right', verticalalignment='bottom', fontsize=20)

plt.savefig(protein_name + "_ccmpred.png",dpi=300)

