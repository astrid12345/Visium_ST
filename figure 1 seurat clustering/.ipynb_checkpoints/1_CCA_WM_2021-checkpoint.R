########################################################################
# Author    : A. Alsema
# Date      : May-July 2021
# Dataset   : Visium Spatial Transcriptomics for MS lesions
# Purpose   : integrate samples - reduce donor effect
# Required inputs    : 1.datasets.rds containing a list of samples
#########################################################################

rm(list =ls())
library(Seurat)
library(hdf5r)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = 30000 * 1024^2)

##################### load data #####################

spot.list <- readRDS('./RData/seurat/1.datasets.rds')

# 
# want to use the RNA (Utility) or Spatial (Seurat) assay for integration
for (i in 1:length(spot.list)){
    DefaultAssay(spot.list[[i]]) <- "Spatial"
}

# normalize and identify variable features for each dataset independently
spot.list <- lapply(X = spot.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = spot.list, nfeatures = 4000)
print("number of features for integration anchors")
print(length(features))
# integration
gc()
st.anchors <- FindIntegrationAnchors(object.list = spot.list, anchor.features = features, reduction = "cca", dims = 1:30)
saveRDS(st.anchors, file = "./RData/seurat/tmp.WM.anchors.rds")
gc()
st.integrated <- IntegrateData(anchorset = st.anchors)
saveRDS(st.integrated, file = "./RData/seurat/2.WM.integrated.rds")

# set default assay
st.integrated <- readRDS(file = "./RData/seurat/2.WM.integrated.rds")
DefaultAssay(st.integrated) <- "integrated"

require(future)
plan("multiprocess", workers = 60)
# normalization for dimensionality reduction
gc()
st.integrated <- ScaleData(st.integrated, vars.to.regress = c("nFeature_Spatial")) # tested but NOT useful vars.to.regress = c("nCount_Spatial"). mito percent could be biological in ST data so I didn't regress it out as unwanted variation.
st.integrated <- RunPCA(st.integrated)
saveRDS(st.integrated, file = "./RData/seurat/2.WM.integrated_scaled.rds")
