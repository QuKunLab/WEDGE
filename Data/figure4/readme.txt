For raw Guo et al. 's dataset, you can obtain the scRNA-seq data of the PBMCs from 2 severe COVID-19 patients from the Gene Expression Comprehensive (GEO) database (GSE150861), and download a scRNA-seq data of PBMCs from 2 healthy donors from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_NGSC3_aggr.
. Then, you can filter genes and cells, the same method as Guo et al., to obtain an expression matrix with 23,324 genes and 68,190 cells.

The test code of other imputation methods can refer to Figure3

Cluster.R: cluster the imputed data
Ref_Label.csv: The cell label reported by Guo et al.
WEDGE_Label.csv: The cell label generated from WEDGE imputation data
cellNames.csv: Cell names in Guo's data
geneName.csv: Gene names in Guo's data
plot_F4abc.m: The code used to draw the subgraphs a, b, c in Figure 4
plot_F4d.R: The code used to plot subfigure d in Figure 4
