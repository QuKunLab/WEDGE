#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/7/14 20:26
# @Author  : YinLei Hu
import scanpy as sc
from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scanpy import AnnData
from typing import Optional, Union
from scipy.sparse import issparse
from sklearn import metrics
import scipy.sparse as ss
import warnings
import itertools
warnings.filterwarnings('ignore')
dataName = 'Baron'
path = './'+ dataName
data_name = ['Reference.csv','Observed.csv','WEDGE_recovery.csv']

labels_true = pd.read_csv(path+'/cellType.csv',header=0, index_col=None)
labels_true = np.squeeze(labels_true.values, axis=1)
labels_true = pd.Categorical(labels_true)
labels_true0 = labels_true
for data_ID in  range(3):
    baron_ref = pd.read_csv(path+'/' +data_name[data_ID], header = 0, index_col=0)
    adata = sc.AnnData(np.transpose(baron_ref.values), obs=pd.DataFrame(baron_ref.columns), var=pd.DataFrame(index=baron_ref.index))

    if data_ID<2:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
        sc.pp.log1p(adata)
    else:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)

    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, resolution=0.3)
    scanpy_label = adata.obs['louvain'].values


    if data_ID==0:
        ARI_value = metrics.adjusted_rand_score(labels_true, scanpy_label)
        NMI_value = metrics.normalized_mutual_info_score(labels_true, scanpy_label)
        print(data_name[data_ID], ARI_value, NMI_value)
    else:
        ARI_value = metrics.adjusted_rand_score(labels_true, scanpy_label)
        NMI_value = metrics.normalized_mutual_info_score(labels_true, scanpy_label)
        print(data_name[data_ID], ARI_value, NMI_value)
