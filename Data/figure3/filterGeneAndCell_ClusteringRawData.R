.libPaths("G:/10x_data/conda_seurat/lib/R/library")
library(Seurat)
library(dplyr)
library(stringr)
library(readr)
library(Matrix)

setwd("G:/10x_data/Mouse_Altas")

library(here)

channel_folders = list.dirs(here("00_data_ingest","01_droplet_raw_data","droplet"), recursive = FALSE)

n = length(strsplit(channel_folders[1],"[/]")[[1]])

raw.data.list = list()
channel.list = list()
for (channel_folder in channel_folders){
  raw.data <- Read10X(channel_folder)
  channel = str_split(str_split(channel_folder,"/", simplify = TRUE)[1,n], "-", simplify = TRUE)[1,2]
  colnames(raw.data) <-  lapply(colnames(raw.data), function(x) paste0(channel, '_', x))
  raw.data.list <- append(raw.data.list, raw.data)
  channel.list <- append(channel.list, rep(channel, length(colnames(raw.data))))
}

raw.data <- do.call(cbind, raw.data.list)
cell.channels <- unlist(channel.list)

ordered_cell_names = order(colnames(raw.data))
raw.data = raw.data[,ordered_cell_names]

meta.data <- read.csv(here("00_data_ingest","01_droplet_raw_data", "metadata_droplet.csv"))
rownames(meta.data) <- meta.data$channel

channel_regex = "(.*?_.*?_.*?)_"
cell.channels <- str_match(colnames(raw.data), channel_regex)[,2]

cell.meta.data <- meta.data[cell.channels,]
rownames(cell.meta.data) <- colnames(raw.data)

# Find ERCC's, compute the percent ERCC, and drop them from the raw data.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,]

# Create the Seurat object with all the data
tiss <- CreateSeuratObject(raw.data = raw.data)

tiss <- AddMetaData(object = tiss, cell.meta.data)
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")

ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")

percent.Rn45s <- tiss@raw.data[c('Rn45s'), ]/Matrix::colSums(tiss@raw.data)
tiss <- AddMetaData(object = tiss, metadata = percent.Rn45s, col.name = "percent.Rn45s")

#A sanity check: genes per cell vs reads per cell.
GenePlot(object = tiss, gene1 = "nUMI", gene2 = "nGene", use.raw=T)

#Filter out cells with few reads and few genes.
tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nUMI"),
                    low.thresholds = c(500, 1000))
tiss@raw.data<-raw.data

#save the data
writeMM(tiss@data,"./Mouse_5W/matrix.mtx")
write.table(tiss@data@Dimnames[2],"./Mouse_5W/barcodes.tsv", row.names = F,
            col.names = F)
write.table(data.frame(tiss@data@Dimnames[[1]],tiss@data@Dimnames[[1]]),"./Mouse_5W/genes.tsv",row.names = F,
            col.names = F)

tiss@raw.data = tiss@data
tiss <- NormalizeData(object = tiss)

tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)

tiss <- ScaleData(object = tiss)

dim(tiss@data)

#Run Principal Component Analysis.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 100)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)
n.pcs = 30
tiss@scale.data <-c()
gc()
#The clustering is performed based on a nearest neighbors graph. Cells that have similar expression will be joined together. The Louvain algorithm looks for groups of cells with high modularity--more connections within the group than between groups. The resolution parameter determines the scale...higher resolution will give more clusters, lower resolution will give fewer.
#For the top-level clustering, aim to under-cluster instead of over-cluster. It will be easy to subset groups and further analyze them below.

# Set resolution
res.used <- 1.1

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
                     resolution = res.used, print.output = 0, save.SNN = TRUE)
#To visualize
# If cells are too spread out, you can raise the perplexity. If you have few cells, try a lower perplexity (but never less than 10).
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, seed.use = 10, perplexity=30, dim.embed = 2)


tiss@meta.data['cluster'] <- tiss@ident
tiss@meta.data['cell'] <- rownames(tiss@meta.data)

anno = read_csv(here("00_data_ingest", "18_global_annotation_csv", "annotations_droplet.csv"))
anno %>% group_by(tissue, cell_ontology_class) %>% summarize(count = n())

tissue_colors = read_csv(here("00_data_ingest", "15_color_palette", "tissue_colors.csv"))
tissue_colors <- rename(tissue_colors, tissue = X1)

tiss@meta.data <- tiss@meta.data %>%
  left_join(anno %>% select(cell_ontology_class,cell_ontology_id,free_annotation, cell), by = 'cell') %>%
  left_join(tissue_colors, by = 'tissue')

rownames(tiss@meta.data) <- tiss@meta.data$cell

TSNEPlot(object = tiss, group.by = 'cell_ontology_class', pt.size = 0.1,do.label = TRUE)+ theme(legend.position="none")
TSNEPlot(object = tiss, pt.size = 0.1,do.label = TRUE)
tiss@raw.data = c()
save(tiss,file = "./Raw.Robj")

            