########################################################################
# Author    : A. Alsema
# Date      : Feb 2023
# Dataset   : Visium Spatial Transcriptomics for MS lesions, 15 WM sections, 11 GM sections.
# Purpose   : perform spatial clustering using the BayesSpace workflow. Analysis is done per sample.

# Required Inputs: 
# - sampleID_sce.rds, a SingleCellExperiment object class for each sample.
# - optimized q-value per sample
# - indir, directory where sce ojects are saved
# - outdir, directory where sce_enhanced are saved
# Output    :
# - sampleID_sce_qX.rds: sce object with clustering results
# - sampleID_enhanced.rds: sce object with enhanced clustering results


# Note: Multiple q-values were tested before running this script.
# Regular clustering was used for the publication, because find cluster markers only runs for regular spot-level clustering. 
# Enhanced_sce with subspots is used to run plot_enhanced_gene_expr

#########################################################################

rm(list = ls())
# load libraries
suppressMessages(library(BayesSpace))
suppressMessages(library(RColorBrewer))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(BiocSingular))
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(doMC))

# register parallel backend
registerDoMC(cores = 31) 

indir <- "<yourindir>"
outdir <- "<youroutdir>"

# vector of samples to be processed
sampleIDs <- c("ST31", "ST32", "ST33", "ST34", 
"ST37", "ST38", "ST67", "ST68", 
"ST69", "ST70", "ST71", "ST72",
"ST73", "ST74", "ST79", 
"ST55","ST56", "ST57", "ST58", 
"ST59", "ST60", "ST61", "ST62", 
"ST63", "ST64", "ST65")
# Vector of corresponding cluster numbers for each sample. A range of q-values have been tested and optimized beforehand.
sampleQs <- c(2, 2, 10, 2, 
2, 3, 9, 8, 
6, 8, 8, 10,
6, 7, 5,
4, 5, 10, 10,
8, 9, 10, 4,
8, 8, 8) 
# Number of PCs dimensions for clustering
d <- 15  	
# set seed
set.seed(123)

# Parallelize with foreach method, one sample per CPU
foreach (i = 1:length(sampleIDs)) %dopar% {
    # Extract the current sample ID
    sampleID <- sampleIDs[i]
    print(sampleID)
    
    # Print the start time
    print(Sys.time())
    
    # Read the preprocessed SingleCellExperiment object from file
    sce <- readRDS(paste0(indir, sampleID, "_sce.rds"))
    
    # Extract the number of clusters for the current sample
    q <- sampleQs[i]
    
    ### Run BayesSpace clustering ####
    # Perform spatial clustering using BayesSpace
    sce <- spatialCluster(sce, q = q, d = d, platform = 'Visium',
                          nrep = 50000, burn.in = 1000, gamma = 3, save.chain = TRUE)
    
    # Save the clustered SingleCellExperiment object to file
    saveRDS(sce, file = paste0(outdir, sampleID, "_sce_q", q, ".rds"))
    print("saved clustering results")
    
    ### Run enhanced BayesSpace clustering for subspots ####
    # Perform enhanced clustering to create subspots and refine clustering
    print("creating subspots and clustering...")
    sce.enhanced <- spatialEnhance(sce, q = q, d = d, platform = "Visium",
                                   nrep = 200000, burn.in = 10000, model = "t", jitter_scale = 5, jitter_prior = 0.3, 
                                   gamma = 3, verbose = TRUE, save.chain = FALSE)
    print("finished enhanced clustering...")
    
    # Save the enhanced clustered SingleCellExperiment object to file
    saveRDS(sce.enhanced, file = paste0(outdir, sampleID, "_enhanced.rds"))
    
    # Print the end time and a message indicating completion for the current sample
    print(Sys.time())
    print(paste("saved enhanced clustering for sample", i))
}