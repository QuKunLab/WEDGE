library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
library(readr)
setwd("path_of_imputed_or_raw_data")

cellNames = read.csv("./WEDGE/Data/figure4/cellNames.csv",row.names = 1)
geneNames = read.csv("./WEDGE/Data/figure4/geneNames.csv",row.names = 1)

method = "WEDGE"
W<-read.csv('./WEDGE/W.csv',header = F)
H<-read.csv('./WEDGE/H.csv',header = F)
A = as.matrix(W)%*%as.matrix(H)
rownames(A) = geneNames$geneNames 
colnames(A) = cellNames$cellNames

pbmc.integrated <- CreateSeuratObject(counts = A, project = "WEDGE", min.cells = 3, min.features = 2)
rm(A)
ref_Label = read.csv("./WEDGE/Data/figure4/Ref_Label.csv")

NC_DEg = c('CCL4','CXCL2','CXCL3','IL6','IL10','ATF3','TNF', 'HIVEP2')
VlnPlot(pbmc.integrated, features = NC_DEg, pt.size = F, ncol = 4, idents = c(1,8,12,15))