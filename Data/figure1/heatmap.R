library(dplyr)
library(Seurat)
library(ggplot2)
setwd("G:/MySoftware/Imputation/WEDGE/Data/figure1")

dropout_name = 'truecounts';
dropout_name = 'count';
dropout_name = 'WEDGE_recovery';


ref_data <- read.csv(paste0('./' ,dropout_name,'.csv'),row.names= 1)
ref_pbmc <- CreateSeuratObject(counts = ref_data, project = "ref", min.cells = 1, min.features =1)
if (dropout_name != 'WEDGE_recovery')
{
  ref_pbmc <- NormalizeData(object = ref_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
}
all.genes <- rownames(x = ref_pbmc)
ref_pbmc <- ScaleData(object = ref_pbmc, features = all.genes)


ref_label <- read.csv('./trueLable.csv')
Idents(ref_pbmc) <- ref_label

if (dropout_name == 'truecounts')
{
ref_pbmc.markers <- FindAllMarkers(object = ref_pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.15)

top10 <- ref_pbmc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
length(top10$gene)
}

DoHeatmap(object = ref_pbmc , features = top10$gene)

path2 = paste0('./','Mark_gene.csv')
mark_id =c()
mark_bool = all.genes %in% top10$gene
for(i in 1:length(all.genes))
{
 if (mark_bool[i])
{
    mark_id<-c(mark_id,1)
}
else
{
 mark_id<-c(mark_id,0)
}

}
write.csv(mark_id,file=path2)
