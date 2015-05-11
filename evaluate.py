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
import pconpy

#PDB_DIR = './pconpy/tests/pdb_files/'

if __name__ == '__main__':

    # Load gold-standard PDB file
    pdb_file = sys.argv[1]
    # pconpy options
    metric='CA'
    residues = pconpy.get_residues(pdb_file)

    gs_mat = pconpy.calc_dist_matrix(residues, metric)

    # Load CAPP prediction
    capp_file = sys.argv[2]
    capp_df = pd.read_table(capp_file, sep=' ', header=None)


    # Check both contact maps have same size
    assert capp_df.shape == gs_mat.shape, 'Shape mismatch'
    print gs_mat[:5, :5]

    print capp_df.iloc[:5, :5]




