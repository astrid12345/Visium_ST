########################################################################
# Author    : A. Alsema
# Date      : July 2021
# Dataset   : Visium Spatial Transcriptomics for MS lesions
# Purpose   : run Spotlight deconvolution of spots per sample. The reference snRNAseq data originated from Schirmer et al.
              # 2019, Nature, doi: 10.1038/s41586-019-1404-z.

# Required inputs:
# - reference_SPOTlight/markers_cellclass_2.csv: the results of FindAllMarkers function on the re-analysed Schirmer data
# can be found under data_files
# - reference_SPOTlight/5.Schirmer_downsampled_cellclass2.rds: the downsampled snRNAseq dataset can be found under data_files
# - 1.datasets.rds: the Visium data of MS lesions, a list of samples directed imported after spaceranger
# - path_to_spaceranger: path to outputs of spaceranger

# Outputs:
# - first version of SPOTlight plots stored in Routput/SPOTlight/
# - cell type proportions per sample stored in "RData/SPOTlight/", ST_ID, "_celltype_prop.csv"
########################################################################


rm(list = ls())
library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(ggcorrplot)

# load single-cell cluster markers and down-sampled single-cell data
cluster_markers_all <- read.csv(
        file = "reference_SPOTlight/markers_cellclass_2.csv", row.names = 1)
ms_down <- readRDS("reference_SPOTlight/5.Schirmer_downsampled_cellclass2.rds")


# load ST data in list, 
ST_ls <- readRDS("Rdata/seurat/1.datasets.rds") 

# loop over WM samples
for (i in 1:length(ST_ls)){
	ST <- ST_ls[[i]]
	ST_ID <- unique(ST$orig.ident)
	print(paste("analysing", ST_ID))
	# filter low-content spots
	print("filtering spots with less than 100 counts or features. From:")
	print(ncol(ST))
	ST <-subset(ST, subset = nCount_Spatial>100)
	ST <-subset(ST, subset = nFeature_Spatial>100)
	print("to")
	print(ncol(ST))

	# takes the longest
	start_time <- Sys.time()
	nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers_all, 
							se_sc = ms_down, 
							mtrx_spatial = ST@assays$Spatial@counts,
							clust_vr = "cellclass_2",
							ntop = NULL,
							hvg = 3000, # try later with 3000
							transf = "uv",
							method = "nsNMF")
	end_time <- Sys.time()

	nmf_mod <- nmf_mod_ls[[1]]
	decon_mtrx <- nmf_mod_ls[[2]]

	# h temporary for plots
	h <- NMF::coef(nmf_mod)
	rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")

	topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
	  h = h,
	  train_cell_clust = decon_mtrx)

	# QC of the trained NMF topics
	pdf(paste0("Routput/SPOTlight/", ST_ID, "_NMF_QC_topic.pdf"), width = 12, height = 20)
	print(topic_profile_plts[[1]] + theme(axis.text.x = element_text(angle = 90), 
									axis.text = element_text(size = 12)))
	print(topic_profile_plts[[2]] + ggplot2::theme(
	  axis.text.x = ggplot2::element_text(angle = 90), 
	  axis.text = ggplot2::element_text(size = 12)))
	dev.off()

	# get basis matrix W. Genes x topic
	# the basis matrix (W) with unique cell  type  marker  genes  and  weights
	w <- basis(nmf_mod)
	dim(w)

	# get coefficient matrix H. topics x cells
	# #  coefficient  matrix  (H)  specifying the corresponding relationship of a cell to a topic 
	h <- coef(nmf_mod)
	dim(h)

	# Extract count matrix
	spot_counts <- ST@assays$Spatial@counts

	# Subset to genes used to train the model
	spot_counts <- spot_counts[rownames(spot_counts) %in% rownames(w), ]

	# find celltype specific profiles 
	ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
								  train_cell_clust = decon_mtrx)
	# regression to get the proportions per spot - this regression is fast
	decon_mtrx <- mixture_deconvolution_nmf(nmf_mod = nmf_mod,
							  mixture_transcriptome = spot_counts,
							  transf = "uv",
							  reference_profiles = ct_topic_profiles, 
							  min_cont = 0.09)

	identical(length(colnames(spot_counts)), nrow(decon_mtrx))
	# should print TRUE, else your new deconvolution matrix is wrong and you are missing spots somewhere

	decon_df <- decon_mtrx %>% data.frame() 
	decon_df$barcodes <- colnames(spot_counts)
	row.names(decon_df) <-  colnames(spot_counts)
    head(decon_df)
    write.csv(decon_df, paste0("RData/SPOTlight/", ST_ID, "_celltype_prop.csv"))
    # add deconvolution results to the meta-data of seurat.
	ST@meta.data <- ST@meta.data %>%
	  tibble::rownames_to_column("barcodes") %>%
	  dplyr::left_join(decon_df, by = "barcodes") %>%
	  tibble::column_to_rownames("barcodes")

	## plots - stored in Routput/SPOTlight/

	# plot 1
	cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

	pdf(paste0("Routput/SPOTlight/", ST_ID, "_piechart.pdf"), width = 10, height = 10)
	print(SPOTlight::spatial_scatterpie(se_obj = ST,
								  cell_types_all = cell_types_all,
								  img_path = paste0(path_to_spaceranger, ST_ID, "/spatial/tissue_lowres_image.png"),
								  pie_scale = 0.4))
	dev.off()

	# plot 2
	pdf(paste0("Routput/SPOTlight/", ST_ID, "_spatialplot_proportions.pdf"))
	print(Seurat::SpatialFeaturePlot(
	  object = ST,
	  features = c("Oligodendrocytes", "OPCs", "GM.astrocytes", "WM.astrocytes"),
	  alpha = c(0.1, 1)))
	print(Seurat::SpatialFeaturePlot(
	  object = ST,
	  features = c("ExNeurons", "InhNeurons", "microglia.macrophages", "Vasculature"),
	  alpha = c(0.1, 1)))
	dev.off()

	# plot 3 - correlation of the percentages
	# Remove cell types not predicted to be on the tissue
	decon_mtrx_sub <- decon_mtrx[, cell_types_all]
	decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]

	# Compute correlation
	decon_cor <- cor(decon_mtrx_sub)

	# Compute correlation P-value
	p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)

	pdf(paste0("Routput/SPOTlight/correlation_perc/", ST_ID, "_Correlation of cell-cell types.pdf"), width = 9, height = 12)
	print(ggcorrplot::ggcorrplot(
	  corr = decon_cor,
	  p.mat = p.mat[[1]],
	  hc.order = TRUE,
	  type = "full",
	  insig = "blank",
	  lab = TRUE,
	  outline.col = "lightgrey",
	  method = "square",
	  colors = c("#6D9EC1", "white", "#E46726"),
	  title = "Predicted cell-cell proportion correlation",
	  legend.title = "Correlation\n(Pearson)") +
	  ggplot2::theme(
		plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
		legend.text = ggplot2::element_text(size = 12),
		legend.title = ggplot2::element_text(size = 15),
		axis.text.x = ggplot2::element_text(angle = 90),
		axis.text = ggplot2::element_text(size = 18, vjust = 0.5)))
	dev.off()
    
}


