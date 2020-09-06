library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
library(readr)
setwd("path_of_imputed_or_raw_data")

cellNames = read.csv("./WEDGE/Data/figure4/cellNames.csv",row.names = 1)
geneNames = read.csv("./WEDGE/Data/figure4/geneNames.csv",row.names = 1)

method = "DCA"
A<-read.csv("./DCA_res/mean.tsv",row.names =1,sep = "\t")
A<-as.matrix(A)

method = "SAVERX"
A = readRDS("./SAVERX_Raw.rds")
A = A$estimate
A[is.na(A)] = 0

method = "MAGIC"
W<-read.csv('./MAGIC//U.csv',header = F)
H<-read.csv('./MAGIC/pc_imputed.csv',header = F)
A = as.matrix(W)%*%t(as.matrix(H))
rownames(A) = geneNames$geneNames 
colnames(A) = cellNames$cellNames

method = "WEDGE"
W<-read.csv('./WEDGE//W.csv',header = F)
H<-read.csv('./WEDGE/H.csv',header = F)
A = as.matrix(W)%*%as.matrix(H)
rownames(A) = geneNames$geneNames 
colnames(A) = cellNames$cellNames

method = "ALRA"
A = readRDS("./ALRA_Raw_recovery.rds")
A = t(A)
rownames(A) = geneNames$geneNames 
colnames(A) = cellNames$cellNames


method = "ENHANCE"
A<-read.csv("./ENHANCE_recovery.csv",row.names =1)
A<-as.matrix(A)
rownames(A) = geneNames$geneNames 
colnames(A) = cellNames$cellNames


pbmc.integrated <- CreateSeuratObject(counts = A, project = "WEDGE", min.cells = 3, min.features = 2)
rm(A)

pbmc.integrated <- NormalizeData(pbmc.integrated)#Only SAVERX and DCA need to perform this step

all.genes <- rownames(pbmc.integrated)
pbmc.integrated <- ScaleData(pbmc.integrated, features = all.genes)
pbmc.integrated <- RunPCA(pbmc.integrated,features = all.genes)
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:50)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.5)
write.csv(pbmc.integrated@active.ident,paste0("./",method,"_Label.csv"),row.names = F)



