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
library(viridis)
library(data.table) 
library(ggalluvial)
library(stringr)
######################################################################################################
# Subclustering   ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load t
######################################################################################################
# (Figure 1e) Dot plot with labels ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

features.plot <- rev(c("CLEC4C", "FCER1A",  "ELANE", "KRT5","KRT14", "KRT1","KRT10","FABP5",
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
                     "GNLY", "KLRB1"))

Idents(object = Round2) <- fct_rev(Idents(object = Round2))

pdf('Dot plot_Round2_extended.pdf', width=14,height=6)
DotPlot(object = Round2, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=18))
dev.off()

Idents(object = Round2) <- fct_rev(Idents(object = Round2))
######################################################################################################



levels(Round1)
