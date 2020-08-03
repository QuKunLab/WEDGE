library(Seurat)
library(dplyr)
library(stringr)
library(readr)
library(Matrix)
library(SparseM)
setwd("~/Imputation/Mouse_5W/")

###################################################################################################################
#SAVERX Cluster
A = readRDS("./1594017289.34807/denoised.rds")$estimate
A[is.na(A)]<-0
A<-as(A, "dgCMatrix")

tiss <- CreateSeuratObject(raw.data = A)
rm(A)
gc()
#Normalize the data, then regress out correlation with total reads
tiss <- NormalizeData(object = tiss)
tiss@raw.data<-c()
tiss <- ScaleData(object = tiss)
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)

#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
tiss@scale.data<-c()

n.pcs = 30
# Set resolution
res.used <- 1.1
tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE)
#To visualize
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, seed.use = 10, perplexity=30, dim.embed = 2)

save(tiss, file="./SAVERX.Robj")

###################################################################################################################
#ALRA Cluster
A <- readRDS("./ALRA_recovery.rds")

tiss <- CreateSeuratObject(raw.data = A)
rm(A)
gc()
#Normalize the data, then regress out correlation with total reads
#tiss <- NormalizeData(object = tiss)
tiss@raw.data<-c()

tiss <- ScaleData(object = tiss, check.for.norm = F)

tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)


#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
tiss@scale.data<-c()
n.pcs = 30
# Set resolution
res.used <- 1.1

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE)

#To visualize
# If cells are too spread out, you can raise the perplexity. If you have few cells, try a lower perplexity (but never less than 10).
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, seed.use = 10, perplexity=30, dim.embed = 2)

save(tiss, file="./ALRA.Robj")


###################################################################################################################
#ENHANCE Cluster

A = read.csv("./ENHANCE_recovery.csv", header = TRUE,row.names =1)
A<-as.matrix(A)
x = rowMeans(A)
y = apply(A,1,var)
Id = x==0 & y==0
A = A[!Id,]
tiss <- CreateSeuratObject(raw.data = A)
rm(A)
gc()
#Normalize the data, then regress out correlation with total reads
tiss <- NormalizeData(object = tiss)
tiss@raw.data<-c()
tiss <- ScaleData(object = tiss)
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)

#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
tiss@scale.data<-c()
gc()
n.pcs = 30
# Set resolution
res.used <- 1.1

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc =TRUE)

#To visualize
# If cells are too spread out, you can raise the perplexity. If you have few cells, try a lower perplexity (but never less than 10).
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs,
                seed.use = 10, perplexity=30, dim.embed = 2,
                check_duplicates= F)

save(tiss, file="./ENHANCE.Robj")

###################################################################################################################
#DCA Cluster

A = read.csv("./DCA_recovery.csv", header = TRUE,row.names =1)
A<-as.matrix(A)
tiss <- CreateSeuratObject(raw.data = A)
rm(A)
gc()
#Normalize the data, then regress out correlation with total reads
tiss <- NormalizeData(object = tiss)
tiss@raw.data<-c()
tiss <- ScaleData(object = tiss)
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)

#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
tiss@scale.data<-c()
gc()
n.pcs = 30
# Set resolution
res.used <- 1.1

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc =TRUE)

#To visualize
# If cells are too spread out, you can raise the perplexity. If you have few cells, try a lower perplexity (but never less than 10).
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs,
                seed.use = 10, perplexity=30, dim.embed = 2,
                check_duplicates= F)

save(tiss, file="./DCA.Robj")

###################################################################################################################
#ENHANCE Cluster

A = read.csv("./MAGIC_recovery.csv", header = TRUE,row.names =1)
A<-as.matrix(A)

tiss <- CreateSeuratObject(raw.data = A)
rm(A)
gc()
#Normalize the data, then regress out correlation with total reads
#tiss <- NormalizeData(object = tiss)
tiss@raw.data<-c()
tiss <- ScaleData(object = tiss)
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)

#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
tiss@scale.data<-c()
gc()
n.pcs = 30
# Set resolution
res.used <- 0.5

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc =TRUE)

#To visualize
# If cells are too spread out, you can raise the perplexity. If you have few cells, try a lower perplexity (but never less than 10).
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs,
                seed.use = 10, perplexity=30, dim.embed = 2,
                check_duplicates= F)

save(tiss, file="./MAGIC_res0_5.Robj")





