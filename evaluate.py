#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Evaluate CAPP contact predictions by comparing them with the gold-standard
structures from PDB files.

Created on Mon May 11 16:59:20 2015

@author: nitin
"""

from __future__ import division
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import classification_report
#import pconpy

#PDB_DIR = './pconpy/tests/pdb_files/'

def plot_roc(y_true_s, y_score, thresholds):
    # Plot of a ROC curve for a range of threshold values
    plt.figure()
    for t in thresholds:
        y_true = y_true_s.copy()
        y_true[y_true <= t] = 1
        y_true[y_true != 1] = 0

        # Compute ROC curve and ROC area 
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        
        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_true, y_score, pos_label=1)
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    
        plt.plot(fpr['micro'], tpr['micro'], label='t = %0.2f A, area = %0.2f)' % (t, roc_auc['micro']))
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic')
        plt.legend(loc="lower right")

    plt.show()
    
def plot_pr(y_true_s, y_score, thresholds):
    # Plot of a ROC curve for a range of threshold values
    plt.figure()
    for t in thresholds:
        y_true = y_true_s.copy()
        y_true[y_true <= t] = 1
        y_true[y_true != 1] = 0

        # Compute ROC curve and ROC area 
        recall = dict()
        precision = dict()
        pr_auc = dict()
        
        # Compute micro-average ROC curve and ROC area
        # Compute micro-average ROC curve and ROC area
        precision["micro"], recall["micro"], _ = precision_recall_curve( 
                                    y_true, y_score, pos_label=1)
        pr_auc["micro"] = average_precision_score(y_true, y_score,
                                                     average="micro")    
        plt.plot(recall['micro'], precision['micro'], label='t = %0.2f A, area = %0.2f)' % (t, pr_auc['micro']))
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])        
        plt.title('Precision Recall Curve')
        plt.legend(loc="lower right")

    plt.show()
    
if __name__ == '__main__':

    # Load gold-standard PDB file
    pdb_file = sys.argv[1]
    # pconpy options
#    metric='CA'
#    residues = pconpy.get_residues(pdb_file)
#
#    gs_mat = pconpy.calc_dist_matrix(residues, metric)
    gs_df = pd.read_table(pdb_file, sep=' ', header=None)
    # Load CAPP prediction
    capp_file = sys.argv[2]
    capp_df = pd.read_table(capp_file, sep=' ', header=None)

    # Check both contact maps have same size
    assert capp_df.shape == gs_df.shape, 'Shape mismatch'
    
    thresholds = range(2, 11) # Angstrom
    
    y_true_s = gs_df.stack()
    y_score_s = capp_df.stack()
    
    plot_roc(y_true_s, y_score_s, thresholds)
    plot_pr(y_true_s, y_score_s, thresholds)

