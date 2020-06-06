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
# Subclustering   ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/KruegerLab_github")
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/KruegerLab_github/Round0_integrated_analyzed.Rda")
######################################################################################################
# Subclustering_Round1 ####
# Clusters expressing gene signatures of more than two different inflammatory cell types 
# and/or clusters spatially separated from similar inflammatory cell types have risks of being multiplets. 
# Therefore, Mature_DC.03, Mature_DC.04, Mature_DC.05, S.Basale.02, and S.Corneum.06 clusters are excluded
# for downstream analysis. 

DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("ClusterNames_0.8")

Round1 <- subset(Round0, idents = c("Mature_DC.03" , "Mature_DC.04",  "Mature_DC.05","S.Basale.02" ,"S.Corneum.06"), invert=TRUE)
rm(Round0)

new.cluster.ids <- c(  "NK_cell"   ,    "CD161_T_cell" , "CD8_T_cell"   , "CD4_T_cell.01" ,"CD4_T_cell.02", "Treg.01"   ,    "Treg.02"    ,   "Mature_DC.01" ,
                       "Mature_DC.02",  "Semimature_DC" ,"Macrophage"   , "Melanocyte"   ,
                       "S.Corneum.01" , "S.Corneum.02" , "S.Corneum.03" , "S.Corneum.04" , "S.Corneum.05" , "S.Granulosum" , "S.Spinosum"  ,  "S.Basale" ,  "ECM")
names(new.cluster.ids) <- levels(Round1)
Round1 <- RenameIdents(Round1, new.cluster.ids)
Round1[["ClusterNames_0.8"]] <- Idents(object = Round1)


########################################################################################################
## (Figure 5) Mature DC01 vs Mature DC02 comparison ####
# The average gene expression of psoriasis vs. control cells within a cluster was calculated by the AverageExpression function. 
Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

Mature_DCs <- subset(Round1, idents = c("Mature_DC.01" ,  "Mature_DC.02" ))
Mature_DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(Mature_DCs), Mature_DCs$stim, sep = "_")
Idents(Mature_DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(Mature_DCs) <- fct_relevel(Idents(Mature_DCs), sort)
Mature_DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = Mature_DCs)
Idents(object = Mature_DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")

Mature_DC.cluster.averages <- AverageExpression(Mature_DCs, return.seurat = TRUE)

marker <- c( "LAMP3", "CD274",     "PDCD1LG2"  ,"KYNU"  ,    "IDO1" ,     "IL23A" , "IL36G" )

pdf('Mature_DC.01vs02_Heatmp.pdf', width=6, height=3)
DoHeatmap(object = Mature_DC.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

####################################################################################
# Save ####
save(Round1, file = "Round1_integrated_analyzed.Rda")
