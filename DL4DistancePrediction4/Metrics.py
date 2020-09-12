import numpy as np

def F1(TP, FP, TN, FN):
	epsilon = 0.000001
	precision = TP*1./(TP+FP + epsilon)
        recall = TP*1./(TP+FN + epsilon)
        F1 = 2.*precision*recall/(precision + recall + epsilon)
	return F1, precision, recall


def MCC(TP, FP, TN, FN):
	epsilon = 0.000001
	MCC = (TP*TN - FP*FN)/np.sqrt( epsilon + 1.0*(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) )
	return MCC
