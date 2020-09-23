#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from os.path import join
import scanpy as sc


def load_WH2anndata(_dir):
    #_dir : the path of WEDGE results(W.csv, H.csv,cellName.csv, and geneNanme.csv)
    H_df = np.loadtxt(join(_dir, 'H.csv'), delimiter=",")
    W_df = np.loadtxt(join(_dir, 'W.csv'), delimiter=",")
    print('data loaded!')
    cellnames = open(join(_dir, 'cellName.csv'), 'r').readlines()
    cellnames = [item.strip('\n') for item in cellnames]
    assert len(cellnames) == H_df.shape[1],'cell number is not equal'

    gene_names = pd.read_csv(join(_dir, 'geneName.csv'), header=None, index_col=0)
    assert len(genenames) == W_df.shape[0],'gene number is not equal'

    mtx = np.dot(W_df, H_df)
    return sc.AnnData(mtx.T, obs=pd.DataFrame(index=cellnames), var=gene_names)
