# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:56:33 2015

@author: nitin
"""

from __future__ import division
import os
import sys
import numpy as np
import itertools
from Bio import AlignIO
from Bio.Data.IUPACData import protein_letters
from profiler import do_profile, timefunc


def phimix_x_given_y(x, y, m):
    """ Compute phi-mixing coefficient between amino acid sequence x and y"""
    alphabets = list(protein_letters + '-')
    len_alphabets = len(alphabets)
    keys = [a+b for a in alphabets for b in alphabets]
    joint_dict = dict.fromkeys(keys, 0)
    
    # Crosstab to get joint dist 
    for i in range(m):
        joint_dict[x[i]+y[i]] += 1

    theta = np.array([joint_dict[k] for k in keys]).reshape((len_alphabets, len_alphabets))

    theta = theta/m
    alpha = np.sum(theta, axis=1)
    beta = np.sum(theta, axis=0)
    product = np.outer(alpha, beta)
    # Remove these later for performance?
    assert abs(np.sum(alpha) - 1) < 10**-4, "x marginal sums to %f != 1" %(np.sum(alpha))
    assert abs(np.sum(beta) - 1) < 10**-4, "y marginal sums to %f != 1" %(np.sum(beta))
    
    gamma_scaled = np.sum(np.fabs(theta - product), axis=0)/beta
    # Handle 0s in the y marginal
    gamma_scaled[np.isnan(gamma_scaled)] = 0
    phimix = 0.5*np.max(gamma_scaled)

    return phimix


#@timefunc
#@do_profile(follow=[phimix_x_given_y])
def compute_pairwise_phimix(aa_array, rand_idx, m, n):
    coupling_mat = np.diag(np.ones(n))
    # Compute pairwise phi-mixing coefficient
    for i,j in itertools.permutations(range(n), 2):
        x = aa_array[rand_idx, i]
        y = aa_array[rand_idx, j]
        coupling_mat[j, i] = phimix_x_given_y(x, y, m)
    
    return coupling_mat
 
#@timefunc   
def prune_with_dpi(coupling_mat, n):
    """ Data Processing Inequality based prunning """
    coupling_mat_prunned = np.copy(coupling_mat)
    # Examine each triplet for DPI prunning
    for i,j,k in itertools.permutations(range(n), 3):
        if coupling_mat[k,i] < min(coupling_mat[k,j], coupling_mat[j,i]):
            coupling_mat_prunned[k,i] = 0

    return coupling_mat_prunned

@timefunc
def capp(aa_array, bs=10, threshold=0):
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
        if bs > 1:
            rand_idx = np.random.choice(m,size=m)
        else:
            rand_idx = np.arange(m)
            
        # Compute phi-mixing coeff
        coupling_mat_bs = compute_pairwise_phimix(aa_array, rand_idx, m, n)
        # Prunning step
        coupling_mat_bs_prunned = prune_with_dpi(coupling_mat_bs, n)
        # Prunning finished. Take the average
        coupling_matrix += coupling_mat_bs_prunned/bs

    # Thresholding step
    thres = np.percentile(coupling_matrix, threshold)
    coupling_matrix[coupling_matrix < thres] = 0

    return coupling_matrix

if __name__=='__main__':
    if len(sys.argv) < 2:
        print('%s fasta_file' %(sys.argv[0]))
        sys.exit(1)
    fasta_file = os.path.abspath(sys.argv[1])

    # Read the fasta file - Alexey
    alignment = AlignIO.read(fasta_file, "fasta")
    align_array = np.array([list(rec) for rec in alignment], np.character)

    cmat = capp(align_array[:, 0:2], 1, 0)
#    print cmat
