.libPaths(.libPaths()[2])
library(Seurat)
library(dplyr)
library(stringr)
library(readr)
library(Matrix)
library(ggplot2)
library(SparseM)
setwd("G:/program1/10x_data/") # the path of WEDGE.Robj and Raw.Robj

#18 B Cell
load("./WEDGE.Robj")
tiss_WEDGE<-tiss

load("./Raw.Robj")
tiss_Raw<-tiss 

#####################################################################################################


TSNEPlot(object = tiss_WEDGE, do.label = TRUE,pt.size = 0.1)

cell_ID<-tiss_Raw@ident==0&(tiss_WEDGE@ident==1|tiss_WEDGE@ident==3|tiss_WEDGE@ident==50|tiss_WEDGE@ident==53)
#cell_ID<-tiss_Raw@ident==0&(tiss_WEDGE@ident==1|tiss_WEDGE@ident==3|tiss_WEDGE@ident==53)

W<-read.csv('./rank_auto/W.csv',header = F)
H_T<-read.csv('./rank_auto/H_new.csv',header = F)
H<-t(as.matrix(H_T))

WEDGE<- as.matrix(W)%*%H
attr(WEDGE,"dimnames")<-tiss_Raw@data@Dimnames

WEDGE<-WEDGE[,cell_ID]
dim<-attr(WEDGE,"dimnames")
WEDGE<-WEDGE%*%as(diag(10000/colSums(WEDGE)), "dgCMatrix")

WEDGE<-as.matrix(WEDGE)
attr(WEDGE,"dimnames")<-dim

tiss_WEDGE_new =tiss_WEDGE
tiss_WEDGE_new@data <- WEDGE
tiss_WEDGE_new@cell.names<-tiss_Raw@cell.names[cell_ID]
Label_new<-factor(tiss_WEDGE@ident[cell_ID],levels = unique(tiss_WEDGE@ident[cell_ID]))
tiss_WEDGE_new@ident<-Label_new
TSNEPlot(object = tiss_WEDGE_new, do.label = TRUE,pt.size = 0.1)

tiss_Raw_new = tiss_Raw
tiss_Raw_new@data<-tiss_Raw@data[,cell_ID]
tiss_Raw_new@cell.names<-tiss_Raw@cell.names[cell_ID]
Label_new<-factor(tiss_WEDGE@ident[cell_ID],levels = unique(tiss_WEDGE@ident[cell_ID]))
tiss_Raw_new@ident<-Label_new
#tiss_Raw_new@dr = tiss_WEDGE_new@dr

TSNEPlot(object = tiss_Raw_new, do.label = TRUE,pt.size = 0.1)
TSNEPlot(object = tiss_WEDGE_new, do.label = TRUE,pt.size = 0.1)

tiss_Raw_new@meta.data$channel = as.character(tiss_Raw_new@meta.data$channel)
ID = tiss_Raw_new@meta.data$channel %in% c("10X_P7_6","10X_P4_7")
tiss_Raw_new@meta.data$channel[!ID] = "others"
color= c("#e6b34a", "#8494ff","#ff61cc")

tiss_WEDGE_new@meta.data$channel = tiss_Raw_new@meta.data$channel

TSNEPlot(object = tiss_Raw_new, pt.size = 0.1, group.by = "channel",colors.use = color)
TSNEPlot(object = tiss_WEDGE_new, pt.size = 0.1, group.by = "channel",colors.use = color)

TSNEPlot(object = tiss_WEDGE_new, do.label = TRUE,pt.size = 0.1, group.by = "tissue")

tiss_Raw_ing = tiss_Raw_new
tiss_Raw_ing@data = as.matrix(tiss_Raw_new@data)

geneName = tiss_Raw_new@data@Dimnames
geneName[[1]][1] = "FO_Raw"
geneName[[1]][2] = "MZ_Raw"
geneName[[1]][3] = "FO_WEDGE"
geneName[[1]][4] = "MZ_WEDGE"
attr(tiss_Raw_ing@data,"dimnames") = geneName

#other gene list of Figure 3 (Kleiman, E. et al.)
########################################################################################################################
# reading gene list (FO)
FO<-read.csv(file="./genelist/Fig3/FO_Up.txt", header = TRUE, sep = "\t")
FO = as.character(FO$X)
#remove undetected gene
gene_id <- FO %in% dim[[1]]
FO_new = c()
p = 1
for(i in gene_id)
{
  if (i) {
    FO_new = c(FO_new, FO[p])
  }
  p = p + 1
}
FO_new = as.character(FO_new)
#####################################
# reading gene list (MZ)
MZ<-read.csv(file="./genelist/Fig3/MZ_Up.txt", header = TRUE, sep = "\t")
MZ = as.character(MZ$X)
#remove undetected gene
gene_id <- MZ %in% dim[[1]]
sum(gene_id)
MZ_new = c()
p = 1
for(i in gene_id)
{
  if (i) {
    MZ_new = c(MZ_new, MZ[p])
  }
  p = p + 1
}
MZ_new = as.character(MZ_new)

#####################################
num_gene_use = 20

#calculating avarage gene expression(FO)
gene_id = dim[[1]] %in% FO_new[1:num_gene_use]
sum(gene_id)
ep_raw = tiss_Raw_new@data[gene_id,]
ep_WEDGE = tiss_WEDGE_new@data[gene_id,]
for (i in c(1:sum(gene_id))) {
  ep_raw[i,] = (ep_raw[i,] - min(ep_raw[i,]))/(max(ep_raw[i,])-min(ep_raw[i,]) + 1e-15) 
  ep_WEDGE[i,] = (ep_WEDGE[i,] - min(ep_WEDGE[i,]))/(max(ep_WEDGE[i,])- min(ep_WEDGE[i,])+ 1e-15)
}
tiss_Raw_ing@data[1,] <-colSums(ep_raw)/sum(gene_id)
tiss_Raw_ing@data[3,] <-colSums(ep_WEDGE)/sum(gene_id)

######################################
#calculating avarage gene expression(FO)
gene_id = dim[[1]] %in% MZ_new[1:num_gene_use]
sum(gene_id)
ep_raw = tiss_Raw_new@data[gene_id,]
ep_WEDGE = tiss_WEDGE_new@data[gene_id,]
for (i in c(1:sum(gene_id))) {
  ep_raw[i,] = (ep_raw[i,] - min(ep_raw[i,]))/(max(ep_raw[i,])-min(ep_raw[i,])+ 1e-15) 
  ep_WEDGE[i,] = (ep_WEDGE[i,] - min(ep_WEDGE[i,]))/(max(ep_WEDGE[i,])- min(ep_WEDGE[i,])+ 1e-15)
}

tiss_Raw_ing@data[2,] <-colSums(ep_raw)/sum(gene_id)
tiss_Raw_ing@data[4,] <-colSums(ep_WEDGE)/sum(gene_id)

######################################
#plot vlnplot
VlnPlot(tiss_Raw_ing, features = c('FO_Raw','FO_WEDGE','MZ_Raw','MZ_WEDGE'),ident.include=c(1,3),
              point.size.use = NA,same.y.lims=TRUE,nCol = 2)

#########################
#t-test
x = tiss_Raw_ing@data[1,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[1,tiss_Raw_ing@ident==3]
p_Fo_Raw = t.test(x,y,paired = FALSE)$p.value

x = tiss_Raw_ing@data[3,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[3,tiss_Raw_ing@ident==3]
p_Fo_WEDGE = t.test(x,y,paired = FALSE)$p.value

x = tiss_Raw_ing@data[2,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[2,tiss_Raw_ing@ident==3]
p_MZ_Raw = t.test(x,y,paired = FALSE)$p.value

x = tiss_Raw_ing@data[4,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[4,tiss_Raw_ing@ident==3]
p_MZ_WEDGE = t.test(x,y,paired = FALSE)$p.value

cat(" p_Fo_Raw =",p_Fo_Raw,"\n p_Fo_WEDGE =",p_Fo_WEDGE,"\n p_MZ_Raw =",p_MZ_Raw,"\n p_MZ_WEDGE =",p_MZ_WEDGE)
###########################################################################################################

##gene list of Supplementary Figure 8 a-b (Newman, R. et al.)
##################################################################################################################
# reading gene list (FO)
FO<-read.csv(file="./genelist/Supplementary_Fig8/WT_FO_cell_SYMBOL.txt", header = TRUE, sep = ",")
FO = as.character(FO$SYMBOL)
#remove undetected gene
gene_id <- FO %in% dim[[1]]
FO_new = c()
p = 1
for(i in gene_id)
{
  if (i) {
    FO_new = c(FO_new, FO[p])
  }
  p = p + 1
}
FO_new = as.character(FO_new)
#####################################
# reading gene list (MZ)
MZ<-read.csv(file="./genelist/Supplementary_Fig8//WT_MZ_cell_SYMBOL.txt", header = TRUE)
MZ = as.character(MZ$SYMBOL)
#remove undetected gene
gene_id <- MZ %in% dim[[1]]
sum(gene_id)
MZ_new = c()
p = 1
for(i in gene_id)
{
  if (i) {
    MZ_new = c(MZ_new, MZ[p])
  }
  p = p + 1
}
MZ_new = as.character(MZ_new)

#####################################
num_gene_use = 20

#calculating avarage gene expression(FO)
gene_id = dim[[1]] %in% FO_new[1:num_gene_use]
sum(gene_id)
ep_raw = tiss_Raw_new@data[gene_id,]
ep_WEDGE = tiss_WEDGE_new@data[gene_id,]
for (i in c(1:sum(gene_id))) {
  ep_raw[i,] = (ep_raw[i,] - min(ep_raw[i,]))/(max(ep_raw[i,])-min(ep_raw[i,]) + 1e-15) 
  ep_WEDGE[i,] = (ep_WEDGE[i,] - min(ep_WEDGE[i,]))/(max(ep_WEDGE[i,])- min(ep_WEDGE[i,])+ 1e-15)
}
tiss_Raw_ing@data[1,] <-colSums(ep_raw)/sum(gene_id)
tiss_Raw_ing@data[3,] <-colSums(ep_WEDGE)/sum(gene_id)

######################################
#calculating avarage gene expression(FO)
gene_id = dim[[1]] %in% MZ_new[1:num_gene_use]
sum(gene_id)
ep_raw = tiss_Raw_new@data[gene_id,]
ep_WEDGE = tiss_WEDGE_new@data[gene_id,]
for (i in c(1:sum(gene_id))) {
  ep_raw[i,] = (ep_raw[i,] - min(ep_raw[i,]))/(max(ep_raw[i,])-min(ep_raw[i,])+ 1e-15) 
  ep_WEDGE[i,] = (ep_WEDGE[i,] - min(ep_WEDGE[i,]))/(max(ep_WEDGE[i,])- min(ep_WEDGE[i,])+ 1e-15)
}
tiss_Raw_ing@data[2,] <-colSums(ep_raw)/sum(gene_id)
tiss_Raw_ing@data[4,] <-colSums(ep_WEDGE)/sum(gene_id)
######################################
#plot vlnplot

VlnPlot(tiss_Raw_ing, features = c('FO_Raw','FO_WEDGE','MZ_Raw','MZ_WEDGE'),ident.include=c(1,3),
              point.size.use = NA,same.y.lims=TRUE,nCol = 2)

#########################
#t-test
x = tiss_Raw_ing@data[1,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[1,tiss_Raw_ing@ident==3]
p_Fo_Raw = t.test(x,y,paired = FALSE)$p.value

x = tiss_Raw_ing@data[3,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[3,tiss_Raw_ing@ident==3]
p_Fo_WEDGE = t.test(x,y,paired = FALSE)$p.value

x = tiss_Raw_ing@data[2,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[2,tiss_Raw_ing@ident==3]
p_MZ_Raw = t.test(x,y,paired = FALSE)$p.value

x = tiss_Raw_ing@data[4,tiss_Raw_ing@ident==1]
y = tiss_Raw_ing@data[4,tiss_Raw_ing@ident==3]
p_MZ_WEDGE = t.test(x,y,paired = FALSE)$p.value

cat(" p_Fo_Raw =",p_Fo_Raw,"\n p_Fo_WEDGE =",p_Fo_WEDGE,"\n p_MZ_Raw =",p_MZ_Raw,"\n p_MZ_WEDGE =",p_MZ_WEDGE)

########################################################################################################################


A_raw = tiss_Raw_new@data[c("Cr2" ,"Fcer2a","Cd93",'Xkr4'),]
geneName = A_raw@Dimnames
geneName[[1]][4] = "nGene"
A_raw@Dimnames = geneName
A_raw['nGene',]=tiss_Raw_new@meta.data$nGene[cell_ID]

A_WEDGE = tiss_WEDGE_new@data[c("Cr2" ,"Fcer2a","Cd93",'Xkr4'),]
attr(A_WEDGE,"dimnames") = geneName
A_WEDGE['nGene',]=tiss_WEDGE_new@meta.data$nGene[cell_ID]

write.csv(as.matrix(A_raw),'./Raw_fig3_Data.csv')
write.csv(tiss_Raw_new@dr$tsne@cell.embeddings[cell_ID,],'./Raw_fig3_TSNE.csv')
write.csv(A_WEDGE,'./WEDGE_fig3_Data.csv')
write.csv(tiss_WEDGE_new@dr$tsne@cell.embeddings[cell_ID,],'./WEDGE_fig3_TSNE.csv')
#we will use plot_GeneExpression.py to plot the gene expression




