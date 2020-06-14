#load library ####
library(Seurat)
library(RNAransform)
library(dplyr)
library(ggplot2)
library(GEOquery)
library(DropletUtils)
library(devtools)
library(URD)
library(edgeR)
library(reticulate)
library(scales)
library(forcats)
library(cowplot)
library(magrittr)
library(varhandle)

######################################################################################################
# Clustering analysis and non-linear dimension reduction ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for Clustering analysis and non-linear dimension reduction. 
######################################################################################################
#Load the integrated data saved from previous code (02) ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round0_integrated.Rda")
######################################################################################################
# Run the standard workflow for visualization and clustering ####
Round0@assays
DefaultAssay(Round0) <- "integrated"
Round0 <- ScaleData(Round0, verbose = FALSE)
Round0 <- RunPCA(Round0, npcs = 20, verbose = FALSE)
######################################################################################################
# t-SNE and Clustering ####
# Twenty principal components (PCs) were selected for Uniform Manifold Approximation and Projection (UMAP)
# for Dimension Reduction. With a resolution of 0.8, cells were clustered by the FindClusters function.
Round0 <- FindNeighbors(Round0, reduction = "pca", dims = 1:20)
Round0 <- FindClusters(Round0, resolution = 0.8)
Round0 <- RunUMAP(object = Round0, reductiosdfn = "pca", dims = 1:20)

#####################################################################################################
# Visualization by clusters ####
DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("integrated_snn_res.0.8")
cell.num <- table(Round0@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round0_cluster_count.pdf', width=10,height=8)
DimPlot(object = Round0,label = TRUE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

#####################################################################################################
# Dot plot to figure out clusters ####

features.plot <- c(  "MGP","KRT5","KRT14", "KRT1","KRT10","FABP5",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "CD163",
                     "CD14","LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CD40","CIITA","LY75","LAMP3",
                     "CTLA4","FOXP3","IL2RA","TIGIT", 
                     "CD8B","CD8A","GZMK","GZMH",
                     "TRBC1","TRAC",
                     "CD3D",
                     "GNLY", "KLRB1")

Idents(object = Round0) <- fct_rev(Idents(object = Round0))

pdf('Dot plot_initial.pdf', width=18,height=16)
DotPlot(object = Round0, features = features.plot) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=18))
dev.off()


########################################################################
# (Supplementary Figure S2) Dot plot with labels ####

Idents(object = Round0) <- ("integrated_snn_res.0.8")
DefaultAssay(Round0) <- "RNA"


features.plot <- c(  "MGP","KRT5","KRT14", "KRT1","KRT10","FABP5",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "CD163",
                     "CD14","LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CD40","CIITA","LY75","LAMP3",
                     "CTLA4","FOXP3","IL2RA","TIGIT", 
                     "CD8B","CD8A","GZMK","GZMH",
                     "TRBC1","TRAC",
                     "CD3D",
                     "GNLY", "KLRB1")


new.cluster.ids <- c(
  "S.Corneum.01","Treg.01","Mature_DC.01","S.Corneum.02","CD4_T_cell","S.Corneum.03","S.Basale.01","S.Corneum.04","S.Spinosum","CD8_T_cell","Semimature_DC","ECM","S.Corneum.05","S.Granulosum","Melanocyte","Doublet.01","Treg.02","CD161_T_cell","NK_cell","Doublet.02","Mature_DC.02","Mature_DC.03","S.Basale.02","Doublet.03","Macrophage","Doublet.04","Doublet.05","Doublet.06"
          )	

names(new.cluster.ids) <- levels(Round0)
Round0 <- RenameIdents(Round0, new.cluster.ids)
Idents(object = Round0) <- factor(Idents(object = Round0), 
                                  levels = (c(
                                    "NK_cell","CD161_T_cell","CD8_T_cell","CD4_T_cell","Treg.01","Treg.02","Mature_DC.01","Mature_DC.02","Mature_DC.03","Semimature_DC","Macrophage","Melanocyte","S.Corneum.01","S.Corneum.02","S.Corneum.03","S.Corneum.04","S.Corneum.05","S.Granulosum","S.Spinosum","S.Basale.01","S.Basale.02","ECM","Doublet.01","Doublet.02","Doublet.03","Doublet.04","Doublet.05","Doublet.06"
                                                                          )))
Round0[["ClusterNames_0.8"]] <- Idents(object = Round0)
Idents(object = Round0) <- fct_rev(Idents(object = Round0))

pdf('Dot plot_Roun0.pdf', width=15,height=9.5)
DotPlot(object = Round0, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=18))
dev.off()

Idents(object = Round0) <- fct_rev(Idents(object = Round0))

########################################################################
# (Supplementary Figure S2) Visualization with label ####
DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("ClusterNames_0.8")
cell.num <- table(Round0@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round0.pdf', width=14,height=8)
DimPlot(object = Round0,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()
#####################################################################################################
# Save ####
save(Round0, file = "Round0_integrated_analyzed.Rda")
#################################################################################################