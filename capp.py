#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coupling analysis of co-evolving protein residues with phi-mixing
coefficient. A method to predict the co-evolutionary interaction between
pairwise amino acid residues of a protein.

Created on Thu Apr 30 11:56:33 2015

@author: nitin
"""

from __future__ import division
import os
import sys
import numpy as np
import itertools
import time
import multiprocessing
import argparse
from Bio import AlignIO
from Bio.Data.IUPACData import extended_protein_letters



def compute_copuling_matrix(aa_array):
    """
    Parameters
    ---------------------
    aa_array: A m-by-n numpy array of amino acid residues representing
            multiple sequence alignment (MSA) of 'm' proteins with 'n' residues

    Return
    ----------------------
    coupling_mat_prunned: An n-by-n numpy array representing the prunned
            coupling matrix.
    """
    m, n = aa_array.shape
    rand_idx = np.random.choice(m,size=m)

    # Compute phi-mixing coeff
    coupling_mat = compute_pairwise_phimix(aa_array, rand_idx, m, n)
    # Prunning step
    coupling_mat_prunned = prune_with_dpi(coupling_mat, n)

    return coupling_mat_prunned


def phimix_x_given_y(x, y, m):
    """ Compute phi-mixing coefficient between amino acid sequence x and y"""
    alphabets = list(extended_protein_letters + '-')
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

    # Handle 0s in the y marginal
    beta[beta == 0] = 10**-4
    gamma_scaled = np.sum(np.fabs(theta - product), axis=0)/beta

#    gamma_scaled[np.isnan(gamma_scaled)] = 0
    phimix = 0.5*np.max(gamma_scaled)

    return phimix


def compute_pairwise_phimix(aa_array, rand_idx, m, n):
    """ Compute pairwise phi-mixing coefficient among all residues"""
    coupling_mat = np.diag(np.ones(n))
    # Compute pairwise phi-mixing coefficient
    for i,j in itertools.permutations(range(n), 2):
        x = aa_array[rand_idx, i]
        y = aa_array[rand_idx, j]
        coupling_mat[j, i] = phimix_x_given_y(x, y, m)

    return coupling_mat

def prune_with_dpi(coupling_mat, n):
    """ Data Processing Inequality based prunning """
    coupling_mat_prunned = np.copy(coupling_mat)
    # Examine each triplet for DPI prunning
    for i,j,k in itertools.permutations(range(n), 3):
        if coupling_mat[k,i] < min(coupling_mat[k,j], coupling_mat[j,i]):
            coupling_mat_prunned[k,i] = 0

    return coupling_mat_prunned


def parse_args(arg_list):
    description = """ Coupling analysis of co-evolving protein residues with
        phi-mixing coefficient. This method predicts the co-evolutionary
        interaction between pairwise amino acid residues of a protein given
        multi-alignment sequences in FASTA format.

        """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('fasta_file', action='store', help='Input fasta file')

    parser.add_argument('-o', action='store', default=None, dest='outfile',
                    help='Output file name. default: fasta_file.txt')

    parser.add_argument('-b', action='store', default=100, type=int,
                    dest='bs',
                    help='Number of bootstraps to run. default: 100')

    parser.add_argument('-t', action='store', default=0, type=float,
                    dest='threshold',
                    help='Threshold to drop low confidence interactions. default: 0')

    parser.add_argument('-v', action='version', version='%(prog)s 0.1')

    results = parser.parse_args(arg_list)

    fasta_file = results.fasta_file
    bs = results.bs
    threshold = results.threshold
    outfile = results.outfile

    return (fasta_file, bs, threshold, outfile)


if __name__=='__main__':

    fasta_file, bs, threshold, outfile = parse_args(sys.argv[1:])

    fasta_file = os.path.abspath(fasta_file)

    # Read the fasta file - Alexey
    alignment = AlignIO.read(fasta_file, "fasta")
    align_array = np.array([list(rec) for rec in alignment], np.character)
    print('Input FASTA file read: %d sequences of %d residues' % align_array.shape)


    align_array_list = [align_array[:,:76]]*bs
    np.random.seed(2015)
    print('Number of boostraps: %d, threshold: %f' %(bs, threshold))

    start = time.time()
    # Run paralle workers for each bootstrap
    pool = multiprocessing.Pool()
    coupling_mat_list = pool.map(compute_copuling_matrix, align_array_list)
    end = time.time()


    # Averaging over all of the bootstrap runs
    coupling_matrix = np.mean(np.array(coupling_mat_list), axis=0)

    # Thresholding step
    thres = np.percentile(coupling_matrix, threshold*100)
    coupling_matrix[coupling_matrix < thres] = 0

    if not outfile:
        outfile = fasta_file + '.txt'
    np.savetxt(outfile, coupling_matrix, fmt='%.3f')

    print('Finished. Time elapsed: %f seconds' %(end - start))
