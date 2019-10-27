#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/7/14 20:26
# @Author  : YinLei Hu
import scanpy as sc
from os.path import join
import pandas as pd
import numpy as np
from scanpy import AnnData
from sklearn import metrics
import warnings
warnings.filterwarnings('ignore')
path = './SF2/'
data_name = ['truecounts','Observed','WEDGE_recovery']
# data_ID = 3
geneinfo = pd.read_csv(path+'geneinfo.csv',header=0, index_col=0)

DE_gene = []
for gene_id in range(200):
    if geneinfo['DEFacGroup1'][gene_id]!= 1 or geneinfo['DEFacGroup2'][gene_id]!= 1:
        DE_gene.append(gene_id)

labels_true = pd.read_csv(path+"cellinfo.csv",header=0, index_col=0)
labels_true = labels_true['Group']
labels_true = pd.Categorical(labels_true)

for data_ID in  range(3):
    Data_= pd.read_csv(path+data_name[data_ID]+".csv", header = 0, index_col=0)
    Data_ = Data_.iloc[DE_gene ]
    adata = sc.AnnData(np.transpose(Data_.values), obs=pd.DataFrame(Data_.columns), var=pd.DataFrame(index=Data_.index))
    sc.pp.filter_genes(adata, min_cells=3)
    # # adata = sc.read_csv("./Data/"+data_name[data_ID], first_column_names=True, dtype='int')
    # # print(adata.X.shape)
    if data_ID<2:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
        sc.pp.log1p(adata)
    else:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)

    # sc.pp.scale(adata)
    # adata[np.isnan(adata)] = 0.0
    sc.tl.pca(adata, n_comps=2, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, resolution=0.5)

    scanpy_label = adata.obs['louvain'].values

    if data_ID==0:
        ARI_value = metrics.adjusted_rand_score(labels_true, scanpy_label)
        NMI_value = metrics.normalized_mutual_info_score(labels_true, scanpy_label)
        print(data_name[data_ID], ARI_value, NMI_value)
    else:
        ARI_value = metrics.adjusted_rand_score(labels_true, scanpy_label)
        NMI_value = metrics.normalized_mutual_info_score(labels_true, scanpy_label)
        print(data_name[data_ID], ARI_value, NMI_value)
