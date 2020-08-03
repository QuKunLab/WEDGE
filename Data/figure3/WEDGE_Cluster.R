library(Seurat)
library(dplyr)
library(stringr)
library(readr)
library(Matrix)


load("./Raw.Robj")#Reduce the memory consumption of creating Seurat objects

W<-read.csv('./auto_rank/W.csv',header = F)
H<-read.csv('./auto_rank/H.csv',header = F)

A <- as.matrix(W)%*%as.matrix(H)
A<-A%*%as(diag(10000/colSums(A)), "dgCMatrix")
A<-as.matrix(A)
attr(A,"dimnames")<-tiss@data@Dimnames

tiss@data<-A

rm(A)
#tiss@data<-exp(tiss@data)-1
tiss@scale.data<-c()
tiss@raw.data<-c()
gc()
#tiss <- NormalizeData(object = tiss)
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)
tiss <- ScaleData(object = tiss,check.for.norm = FALSE)

#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)

# Set number of principal components.
n.pcs = 30

# Set resolution
res.used <- 1.1
tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
tiss@data<-c()
tiss@scale.data<-c()
gc()
#To visualize
# If cells are too spread out, you can raise the perplexity. If you have few cells, try a lower perplexity (but never less than 10).
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, seed.use = 10, perplexity=30, dim.embed = 2)

save(tiss, file="./WEDGE.Robj")