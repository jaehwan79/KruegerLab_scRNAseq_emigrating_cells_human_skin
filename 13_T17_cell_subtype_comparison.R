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
library(ggeasy)
library(farver)
library(DescTools)
#install.packages('pheatmap') # if not installed already
library(pheatmap)
#install.packages("heatmaply")
library(heatmaply)
######################################################################################################
# Number of cells (psoriasis vs. control) in each cluster calculation   ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")

#####################################################################################################################
## T17 cell subset comparison ####
rm(list= ls()[!(ls() %in% c('Round2'))])

Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"
Round2$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(Round2), Round2$stim, sep = "_")
Idents(Round2) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(Round2) <- fct_relevel(Idents(Round2), sort)
Round2[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = Round2)
Idents(object = Round2) <- ("ClusterNames_0.8_PsoriasisvsControl")
levels(Round2)

#IL17A+ IL17F- T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17A.Tcells <- subset(Tcells, subset = IL17A > 1)
IL17A.Tcells <- subset(IL17A.Tcells, subset = IL17F < 1)
new.cluster.ids <- c("IL17A.Tcells","IL17A.Tcells","IL17A.Tcells","IL17A.Tcells")
names(new.cluster.ids) <- levels(IL17A.Tcells)
IL17A.Tcells <- RenameIdents(IL17A.Tcells, new.cluster.ids)
levels(Idents(IL17A.Tcells))

#IL17A- IL17F+ T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17F.Tcells <- subset(Tcells, subset = IL17F > 1)
IL17F.Tcells <- subset(IL17F.Tcells, subset = IL17A < 1)
new.cluster.ids <- c("IL17F.Tcells","IL17F.Tcells","IL17F.Tcells")
names(new.cluster.ids) <- levels(IL17F.Tcells)
IL17F.Tcells <- RenameIdents(IL17F.Tcells, new.cluster.ids)
levels(Idents(IL17F.Tcells))

#IL17A+ IL17F+ T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17AF.Tcells <- subset(Tcells, subset = IL17A > 1)
IL17AF.Tcells <- subset(IL17AF.Tcells, subset = IL17F > 1)
new.cluster.ids <- c("IL17AF.Tcells","IL17AF.Tcells","IL17AF.Tcells")
names(new.cluster.ids) <- levels(IL17AF.Tcells)
IL17AF.Tcells <- RenameIdents(IL17AF.Tcells, new.cluster.ids)
levels(Idents(IL17AF.Tcells))


#IL17A- IL17F+ IL10+ T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17F.Tcells <- subset(Tcells, subset = IL17F > 1)
IL17F.Tcells <- subset(IL17F.Tcells, subset = IL17A < 1)
IL17F.Tcells.IL10.pos <- subset(IL17F.Tcells, subset = IL10 > 1)
levels(IL17F.Tcells.IL10.pos) 
new.cluster.ids <- c("IL17F.Tcells.IL10.pos","IL17F.Tcells.IL10.pos")
names(new.cluster.ids) <- levels(IL17F.Tcells.IL10.pos)
IL17F.Tcells.IL10.pos <- RenameIdents(IL17F.Tcells.IL10.pos, new.cluster.ids)

#IL17A- IL17F+ IL10- T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17F.Tcells <- subset(Tcells, subset = IL17F > 1)
IL17F.Tcells <- subset(IL17F.Tcells, subset = IL17A < 1)
IL17F.Tcells.IL10.neg <- subset(IL17F.Tcells, subset = IL10 < 1)
levels(IL17F.Tcells.IL10.neg) 
new.cluster.ids <- c("IL17F.Tcells.IL10.neg","IL17F.Tcells.IL10.neg","IL17F.Tcells.IL10.neg")
names(new.cluster.ids) <- levels(IL17F.Tcells.IL10.neg)
IL17F.Tcells.IL10.neg <- RenameIdents(IL17F.Tcells.IL10.neg, new.cluster.ids)

#IL17A+ IL17F- IFNG+ T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17A.Tcells <- subset(Tcells, subset = IL17A > 1)
IL17A.Tcells <- subset(IL17A.Tcells, subset = IL17F < 1)
IL17A.Tcells.IFNG.pos <- subset(IL17A.Tcells, subset = IFNG > 1)
levels(IL17A.Tcells.IFNG.pos) 
new.cluster.ids <- c("IL17A.Tcells.IFNG.pos","IL17A.Tcells.IFNG.pos")
names(new.cluster.ids) <- levels(IL17A.Tcells.IFNG.pos)
IL17A.Tcells.IFNG.pos <- RenameIdents(IL17A.Tcells.IFNG.pos, new.cluster.ids)
levels(Idents(IL17A.Tcells.IFNG.pos))

#IL17A+ IL17F- IFNG- T-cell subset
Tcells <- subset(Round2, idents = c(  "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis", "Treg_Psoriasis"))
IL17A.Tcells <- subset(Tcells, subset = IL17A > 1)
IL17A.Tcells <- subset(IL17A.Tcells, subset = IL17F < 1)
IL17A.Tcells.IFNG.neg <- subset(IL17A.Tcells, subset = IFNG < 1)
levels(IL17A.Tcells.IFNG.neg) 
new.cluster.ids <- c("IL17A.Tcells.IFNG.neg","IL17A.Tcells.IFNG.neg","IL17A.Tcells.IFNG.neg","IL17A.Tcells.IFNG.neg")
names(new.cluster.ids) <- levels(IL17A.Tcells.IFNG.neg)
IL17A.Tcells.IFNG.neg <- RenameIdents(IL17A.Tcells.IFNG.neg, new.cluster.ids)
levels(Idents(IL17A.Tcells.IFNG.neg))






IL17AF.Tcells.IL10.division <- merge(x= IL17A.Tcells.IFNG.pos, y = c(IL17A.Tcells.IFNG.neg, IL17AF.Tcells,IL17F.Tcells.IL10.pos,IL17F.Tcells.IL10.neg ), project = "IL17AF.Tcells")

levels(Idents(IL17AF.Tcells.IL10.division))
Idents(object = IL17AF.Tcells.IL10.division) <- factor(Idents(object = IL17AF.Tcells.IL10.division), 
                                                       levels = (c(  "IL17A.Tcells.IFNG.pos"  ,  "IL17A.Tcells.IFNG.neg"  ,
                                                                     "IL17AF.Tcells"  , "IL17F.Tcells.IL10.neg","IL17F.Tcells.IL10.pos")))
#####################################################################
# Heatmap for average gene expression ####
IL17AF.Tcells.IL10.division.averages <- AverageExpression(IL17AF.Tcells.IL10.division, return.seurat = FALSE)

marker <- c(  
  "CD4"   , 
  "IFNL1",  "GNLY"  , "PRF1"  ,  
  "CD8B"  ,
  "TNFSF10",  
  "IL36G", "TNF", 
  "TNFSF14", "LTA",  "IL24", "IL34", "EBI3", "CCR6", 
  "RORC", 
  "IFNG"  ,  "IL17A"  , "IL17F" ,  "IL26"  ,  
  "IL10",  
  "CCR10",
  "IL22",   
   "FOXP3" , 
  "TNFRSF4",
  "MAF",
  "CTLA4", 
  "CD86", 
  "CSF2",
  "TGFB1",
  "GZMA" , 
  "IL23R", "IL1R1", "IL2", 
  "CD69",  
  "IL1B", "AHR", 
  "SLC7A8",  
  "FOSB", "BATF", "CCL3", "CCL4", "GZMB",
  "NT5E",
  "CCL5",
  "STAT4",
  "CASP1",
  "LAG3",
  "IL1RN",
  "ITGA2"
)


a <- as.data.frame(IL17AF.Tcells.IL10.division.averages$RNA)
b <- a[marker,]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(b, 1, cal_z_score))
data_subset_norm <- na.omit(data_subset_norm)

# Figure 4b
pdf("Heatmap_gene_expression_IL.17.A.F.T.IL10.cellsAverage_expression09.pdf", width=6,height=10)
p <-   pheatmap(data_subset_norm, display_numbers = F,cluster_cols = F, cluster_rows=T
                ,   colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(512)
                ,number_format = "%.0f"
                , main = "T17 cells:Average expression",border_color=T ,
                clustering_method = "complete",
                scale = "row")
print(p)
dev.off()
