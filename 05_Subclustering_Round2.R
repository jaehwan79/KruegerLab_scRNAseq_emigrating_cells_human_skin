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
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round0_integrated_analyzed.Rda")
######################################################################################################
# Subclustering_Round2 ####
# Clusters expressing gene signatures of more than two different inflammatory cell types 
# and/or clusters spatially separated from similar inflammatory cell types have risks of being multiplets. 
# Therefore, Doublet 01-07 clusters are excluded
# for downstream analysis. 

# In addition, clusters adjacent to each other representing the same inflammatory cell subtypes are combined for 
# downstream analysis (Round2), since the purpose of the analysis is to compare psoriasis vs. control cells 
# within the inflammatory cell subtype clusters. 
# Therefore, CD4_T_cell.01 and CD4_T_cell.02 are combined. Treg.01 and Treg.02 are combined. 
# S.Corneum.01, S.Corneum.02, S.Corneum.03, S.Corneum.04, and S.Corneum.05 are combined.  

DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("ClusterNames_0.8")
Round2 <- subset(Round0, idents = c("Doublet.01" , "Doublet.02",  "Doublet.03","Doublet.04" ,"Doublet.05","Doublet.06", "Doublet.07"), invert=TRUE)
rm(Round0)

new.cluster.ids <- c("NK_cell"   ,    "CD161_T_cell" , "CD8_T_cell" ,   "CD4_T_cell"  ,  "Treg"    ,   "Treg"    ,
                     "Mature_DC",  "Mature_DC",  "Mature_DC" , "Semimature_DC", "Macrophage"  ,  "Melanocyte"   , "S.Corneum" ,
                     "S.Corneum",  "S.Corneum",  "S.Corneum",  "S.Corneum" , "S.Granulosum" , "S.Spinosum"  ,  "S.Basale"  , "S.Basale" )

names(new.cluster.ids) <- levels(Round2)
Round2 <- RenameIdents(Round2, new.cluster.ids)
Round2[["ClusterNames_0.8"]] <- Idents(object = Round2)

######################################################################################################
# Dot plot with labels ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

features.plot <- c(  "KRT5","KRT14", "KRT1","KRT10","FABP5",
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

Idents(object = Round2) <- fct_rev(Idents(object = Round2))

pdf('Dot plot_Round2.pdf', width=14,height=6)
DotPlot(object = Round2, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=18))
dev.off()

Idents(object = Round2) <- fct_rev(Idents(object = Round2))

######################################################################################################
# Visualization with label ####
DefaultAssay(Round2) <- "RNA"
Idents(object = Round2) <- ("stim")
Idents(object = Round2) <- factor(Idents(object = Round2), levels = c("Control","Psoriasis"))         
Idents(object = Round2) <- ("ClusterNames_0.8")

pdf('UMAP_Round2_01.pdf', width=15,height=7)
DimPlot(Round2, reduction = "umap", split.by = "stim", cols=c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','dark slate gray', 'indianred4','deepskyblue','royalblue1','navy','slate gray','gainsboro'))
dev.off()


##cluster count
Idents(object = Round2) <- ("ClusterNames_0.8")
cell.num <- table(Round2@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round2_02.pdf', width=13,height=8)
DimPlot(object = Round2,label = FALSE,  reduction = "umap") +
  scale_colour_manual(values = c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','dark slate gray', 'indianred4','deepskyblue','royalblue1','navy','slate gray','gainsboro'), breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()
######################################################################################################
# Heatmap with labels ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

Combined.markers <- FindAllMarkers(object = Round2, only.pos = TRUE, min.pct = 0.25, 
                                   logfc.threshold = 0.25)
top10 <- Combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
Round2 <- ScaleData(Round2, verbose = FALSE)

pdf('Diff_Heatmp_labeled.pdf', width=24, height=28)
DoHeatmap(object = Round2, features = top10$gene,raster = TRUE, group.bar = TRUE, draw.lines=TRUE) 
dev.off()

###################################################################################################
# Average expression - Dendritic cells ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DCs <- subset(Round2, idents = c("Mature_DC" ,   "Semimature_DC" ))

DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$stim, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(DCs) <- fct_relevel(Idents(DCs), sort)
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")

Idents(object = DCs) <- factor(Idents(object = DCs), levels = (c(  "Mature_DC_Control"  ,  "Mature_DC_Psoriasis"  ,
                                                                   "Semimature_DC_Control"  , "Semimature_DC_Psoriasis")))
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)

DC.cluster.averages <- AverageExpression(DCs, return.seurat = TRUE)

marker <- c( "CD207", "IL23A" , "IL36G" , "HLA-DRB5" , "HLA-DRA" ,  "HLA-DPB1" , "CD86"   , "LAMP3"   ,  "LY75"   ,   "CD40"   ,   "CIITA" ,        "CD274"   ,  "PDCD1LG2" ,
            "CD80"   , "CLEC4C",  "KYNU"   ,  "CD1C" ,
            "ITGAX",  "SIRPA" , "CD14"   ,   "AIF1"    , "IL10",  "THBD"  ,   "CD209"  ,       "TNFRSF10C"  ,           "ITGAM" ,  
 "LILRB2" ,   "LILRB1"  ,    "LILRB4"   
           )

pdf('Dendritic_cells.cluster.averages_Round2_Heatmp.pdf', width=6, height=6)
DoHeatmap(object = DC.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
 # scale_fill_viridis(option="magma")
 # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
   scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# Average expression - Keratinocytes ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

KCs <- subset(Round2, idents = c( c("S.Corneum",   "S.Granulosum" , "S.Spinosum"  ,  "S.Basale"  )))
KCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(KCs), KCs$stim, sep = "_")
Idents(object = KCs) <- ("ClusterNames_0.8_PsoriasisvsControl")

Idents(object = KCs) <- factor(Idents(object = KCs), levels = (c("S.Corneum_Control" ,     "S.Corneum_Psoriasis"  ,  
                                                                 "S.Granulosum_Control"  ,    "S.Granulosum_Psoriasis"  ,
                                                                 "S.Spinosum_Control" , "S.Spinosum_Psoriasis",
                                                                 "S.Basale_Control",   "S.Basale_Psoriasis" )))
KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)

KC.cluster.averages <- AverageExpression(KCs, return.seurat = TRUE)

marker <- c("LCE3D", "CDSN" ,  "FLG" ,"FABP5" ,"IL36G","STAT3", "NFKBIZ",  "KRT1",  "KRT10" ,"KRT5" , "KRT14" ,"CCL27", "KRT15" ,"CD34" )

pdf('KC.cluster.averages_Round2_Heatmp.pdf', width=6, height=3.5)
DoHeatmap(object = KC.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# IL-10 Feature plot and violin plot ####
pdf("IL10_FeaturePlot_Round2.pdf", width=4,height=3)
p <- FeaturePlot(object = Round2, features = "IL10",min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
print(p)
dev.off()

pdf("IL10_VlnPlot_Round2.pdf", width=12,height=4)
p <- VlnPlot(Round2, features="IL10",  split.by = "stim", cols = c("cyan", "red", "green"))
print(p)
dev.off()

######################################################################################################
# Save ####
save(Round2, file = "Round2_integrated_analyzed.Rda")
#################################################################################################
