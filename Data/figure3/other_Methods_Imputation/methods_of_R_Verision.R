############################################################################
#SAVERX

Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/SAVERX/bin/python3")
library(Matrix)
library(SAVERX)

setwd("./Mouse_5W/")

ptm <- proc.time()

file <- saverx("./Mouse_5W.rds", data.species = "Mouse",
               use.pretrain = F,ncores =72, is.large.data = T, model.species = "Mouse")
A <- proc.time() - ptm

################################################################################
#VIPER
setwd("./Mouse_5W/")
library(VIPER)
library(Seurat)
data <- Read10X(data.dir = "../all_data_5W/")

VIPER(as.data.frame(data), percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
             report = T, outdir = './VIPER_res/',
             prefix = "Mouse_VIPER_")
###################################################################################
#ALRA
setwd("/home/math/hyl2016/Imputation/Mouse_5W/")
library(Seurat)

data <- Read10X(data.dir = "../all_data_5W/")
ref_pbmc <- CreateSeuratObject(counts = data, project = "ALRA", min.cells = 1, min.features =1)

ptm <- proc.time()

ref_pbmc <- RunALRA(ref_pbmc)

A <- proc.time() - ptm
###################################################################################