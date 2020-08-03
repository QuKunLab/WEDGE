#nead seuat 2.0
.libPaths(.libPaths()[2])#package location
library(Seurat)
library(dplyr)
library(stringr)
library(readr)
library(Matrix)
library(ggplot2)
library(SparseM)
setwd("G:/program1/10x_data/")

#B Cell
load("./WEDGE.Robj")
tiss_WEDGE<-tiss

load("./Raw.Robj")
tiss_Raw<-tiss 
#####################################################################################################

method_name_set = c("Raw","WEDGE","DCA",'SAVERX','MAGIC_res0_5','ALRA','ENHANCE')
method_id = 7 # change this to plot different method results

method_name = method_name_set[method_id]
load(paste0("./", method_name ,".Robj"))
     
     
TSNEPlot(object = tiss, do.label = TRUE,pt.size = 0.1)

#cell_ID<-tiss_Raw@ident==0&(tiss_WEDGE@ident==1|tiss_WEDGE@ident==3|tiss_WEDGE@ident==50|tiss_WEDGE@ident==53)
cell_ID<-(tiss_Raw@ident==0)&(tiss_WEDGE@ident==1|tiss_WEDGE@ident==3|tiss_WEDGE@ident==53)

Label_WEDGE<-factor(tiss_WEDGE@ident[cell_ID],levels = unique(tiss_WEDGE@ident[cell_ID]))

tiss_Raw_new = tiss_Raw
tiss_Raw_new@data<-tiss_Raw@data[,cell_ID]
tiss_Raw_new@cell.names<-tiss_Raw@cell.names[cell_ID]
tiss_Raw_new@ident<-Label_WEDGE
tiss_Raw_new@dr$tsne@cell.embeddings = tiss@dr$tsne@cell.embeddings
attr(tiss_Raw_new@dr$tsne@cell.embeddings,"dimnames") = attr(tiss_Raw@dr$tsne@cell.embeddings,"dimnames")


tiss_new = tiss_Raw
tiss_new@data<-tiss@data[,cell_ID]

tiss_new@data@Dimnames[2] = tiss_Raw_new@data@Dimnames[2]# If it is MAGIC, please skip this step

tiss_new@cell.names<-tiss_Raw@cell.names[cell_ID]
Label_new<-factor(tiss@ident[cell_ID],levels = unique(tiss@ident[cell_ID]))
tiss_new@ident<-Label_new
tiss_new@dr$tsne@cell.embeddings = tiss@dr$tsne@cell.embeddings
attr(tiss_new@dr$tsne@cell.embeddings,"dimnames") = attr(tiss_Raw@dr$tsne@cell.embeddings,"dimnames")


#Label_Id = tiss@ident %in% c(10,1,46,3)
#tiss_ident = as.character(tiss@ident)
#tiss_ident[!Label_Id] = 'other'
tiss_new@meta.data$res.1 = as.character(tiss@ident)
TSNEPlot(object = tiss_new, do.label = TRUE,pt.size = 0.1, group.by = "res.1")


tiss_new@meta.data$channel = as.character(tiss_Raw_new@meta.data$channel)
ID = tiss_new@meta.data$channel %in% c("10X_P7_6","10X_P4_7")
tiss_new@meta.data$channel[!ID] = "others"
color= c("#e6b34a", "#8494ff","#ff61cc")
TSNEPlot(object = tiss_new, pt.size = 0.1, group.by = "channel",colors.use = color)

path = "G:/program1/10x_data/WEDGE_res/WEDGE_paper/code/figure3/" #'./'
A = tiss_new@data[c("Cr2" ,"Fcer2a"),]
write.csv(as.matrix(A),paste0(path,method_name,'_Data.csv'))
write.csv(tiss_new@dr$tsne@cell.embeddings[cell_ID,],paste0(path,method_name,'_TSNE.csv'))
#please use plot_GeneExpression.py to plot the gene expression


FeaturePlot(tiss_Raw_new, features=c('Cr2','Fcer2a'),no.legend = F)#Raw wxpression
FeaturePlot(tiss_new, features=c('Cr2','Fcer2a'),no.legend = F)#Imputation expression
                         