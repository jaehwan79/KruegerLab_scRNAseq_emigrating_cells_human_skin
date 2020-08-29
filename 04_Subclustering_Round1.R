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
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round0_integrated_analyzed.Rda")
######################################################################################################
# Subclustering_Round1 ####
# Clusters expressing gene signatures of more than two different inflammatory cell types 
# and/or clusters spatially separated from similar inflammatory cell types have risks of being multiplets. 
# Therefore, Mature_DC.03, Mature_DC.04, Mature_DC.05, S.Basale.02, and S.Corneum.06 clusters are excluded
# for downstream analysis. 

DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("ClusterNames_0.8")

Round1 <- subset(Round0, idents = c("Doublet.01" , "Doublet.02",  "Doublet.03","Doublet.04" ,"Doublet.05","Doublet.06","Doublet.07"), invert=TRUE)
rm(Round0)

new.cluster.ids <- c("NK_cell"   ,    "CD161_T_cell" , "CD8_T_cell" ,   "CD4_T_cell"  ,  "Treg.01"    ,   "Treg.02"    ,
                     "Mature_DC.01",  "Mature_DC.02",  "Mature_DC.03" , "Semimature_DC", "Macrophage"  ,  "Melanocyte"   , "S.Corneum.01" ,
                      "S.Corneum.02",  "S.Corneum.03",  "S.Corneum.04",  "S.Corneum.05" , "S.Granulosum" , "S.Spinosum"  ,  "S.Basale.01"  , "S.Basale.02"   )
names(new.cluster.ids) <- levels(Round1)
Round1 <- RenameIdents(Round1, new.cluster.ids)
Round1[["ClusterNames_0.8"]] <- Idents(object = Round1)


####################################################################################
# Save ####
save(Round1, file = "Round1_integrated_analyzed.Rda")

