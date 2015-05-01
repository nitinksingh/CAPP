# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:56:33 2015

@author: nitin
"""

from __future__ import division
import os
import sys
import numpy as np
import pandas as pd
import itertools
from Bio import AlignIO
from Bio.Data.IUPACData import protein_letters

# Define commonly used variables
alphabets = np.array(list(protein_letters + '-'))
I = np.diag(np.ones(len(alphabets)))

def phimix_x_given_y(x, y):
    """ Compute phi-mixing coefficient between amino acid sequence x and y"""
    # A little hack to ensure 21 by 21 cross-tabulation happens
    m = len(x)
    x = np.append(x, alphabets)
    y = np.append(y, alphabets)
    # Crosstab to get joint dist and remedy the appended alphabet
    theta = (pd.crosstab(x, y).values - I)/m
    alpha = np.sum(theta, axis=1)
    beta = np.sum(theta, axis=0)

    # Remove these later for performance?
    assert abs(np.sum(alpha) - 1) < 10**-4, "x marginal sums to %f != 1" %(np.sum(alpha))
    assert abs(np.sum(beta) - 1) < 10**-4, "y marginal sums to %f != 1" %(np.sum(beta))

    product = np.outer(alpha, beta)
    gamma_scaled = np.sum(np.fabs(theta - product), axis=0)/beta
    # Handle 0s in the y marginal
    gamma_scaled[np.isnan(gamma_scaled)] = 0
    phimix = 0.5*np.max(gamma_scaled)

    return phimix


def capp(aa_array, bs=10, threshold=99):
    """ Coupling analysis of co-evolving protein residues with phi-mixing
    coefficient. A method to predict the co-evolutionary interaction between
    pairwise amino acid residues of a protein.

    Parameters
    ---------------------
    aa_array: A m-by-n numpy array of amino acid residues representing
            multiple sequence alignment (MSA) of 'm' proteins with 'n' residues

    bs: Number of bootstrap runs. default 10
    threshold: The percentile threshold to drop low weight interactions

    Return
    ----------------------
    coupling_matrix: An n-by-n numpy array representing the predicted coupling
            strength predicted by the CAPP method.
    """
    m, n = aa_array.shape
    np.random.seed(2015)
    coupling_matrix = np.zeros((n, n))

    for b in xrange(bs):
        rand_idx = np.random.choice(m,size=m)
        coupling_mat_bs = np.diag(np.ones(n))
        # Compute pairwise phi-mixing coefficient
        for i,j in itertools.permutations(range(n), 2):
            x = aa_array[rand_idx, i]
            y = aa_array[rand_idx, j]
            coupling_mat_bs[j, i] = phimix_x_given_y(x, y)

        coupling_mat_bs_prunned = np.copy(coupling_mat_bs)
        # Examine each triplet for DPI prunning
        for i,j,k in itertools.permutations(range(n), 3):
            if coupling_mat_bs[k,i] < min(coupling_mat_bs[k,j], coupling_mat_bs[j,i]):
                coupling_mat_bs_prunned[k,i] = 0

        # Prunning finished. Take the average
        coupling_matrix += coupling_mat_bs_prunned/bs

    # Thresholding step
    thres = np.percentile(coupling_matrix, threshold)
    coupling_matrix[coupling_matrix < thres] = 0

    return coupling_matrix

if __name__=='__main__':
#    if len(sys.argv) < 2:
#        print('%s fasta_file' %(sys.argv[0]))
#        sys.exit(1)
#    fasta_file = sys.argv[1]
    # Read the fasta file -- Alexey

    fasta_file = "/Users/nitin/Dropbox/LMB-Alexey/data/dirty_ArgRS_alignment/P11875.fas"
    alignment = AlignIO.read(fasta_file, "fasta")
    align_array = np.array([list(rec) for rec in alignment], np.character)

    cmat = capp(align_array[0:20, 0:5], 10, 0.5)
    print cmat