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
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/KruegerLab_github")
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/KruegerLab_github/Round0_integrated_analyzed.Rda")
######################################################################################################
# Subclustering_Round2 ####
# Clusters expressing gene signatures of more than two different inflammatory cell types 
# and/or clusters spatially separated from similar inflammatory cell types have risks of being multiplets. 
# Therefore, Mature_DC.03, Mature_DC.04, Mature_DC.05, S.Basale.02, and S.Corneum.06 clusters are excluded
# for downstream analysis. 

# In addition, clusters adjacent to each other representing the same inflammatory cell subtypes are combined for 
# downstream analysis (Round2), since the purpose of the analysis is to compare psoriasis vs. control cells 
# within the inflammatory cell subtype clusters. 
# Therefore, CD4_T_cell.01 and CD4_T_cell.02 are combined. Treg.01 and Treg.02 are combined. 
# S.Corneum.01, S.Corneum.02, S.Corneum.03, S.Corneum.04, and S.Corneum.05 are combined.  

DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("ClusterNames_0.8")
Round2 <- subset(Round0, idents = c("Mature_DC.03" , "Mature_DC.04",  "Mature_DC.05","S.Basale.02" ,"S.Corneum.06"), invert=TRUE)
rm(Round0)

new.cluster.ids <- c( "NK_cell"  ,     "CD161_T_cell",  "CD8_T_cell"   , "CD4_T_cell" ,"CD4_T_cell" ,"Treg" ,      "Treg"  ,
                      "Mature_DC" , "Mature_DC" , "Semimature_DC", "Macrophage"   , "Melanocyte"   ,
                       "S.Corneum"  ,"S.Corneum"  ,"S.Corneum" , "S.Corneum" , "S.Corneum" , "S.Granulosum" , "S.Spinosum" ,   "S.Basale"  , "ECM"       )
names(new.cluster.ids) <- levels(Round2)
Round2 <- RenameIdents(Round2, new.cluster.ids)
Round2[["ClusterNames_0.8"]] <- Idents(object = Round2)

######################################################################################################
# (Figure 2) Dot plot with labels ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"


features.plot <- c(  "MGP","KRT5","KRT14", "KRT10",  "KRT1",
                     "FABP5",
                     "LCE1A",
                     "AREG",
                     "KRT6B",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "CD163",
                     "CD68",
                     "FCGR3A","ITGAX", 
                     "AIF1","CD14","LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CD40","CIITA","LY75","LAMP3",
                     "CTLA4","FOXP3","IL2RA","TIGIT", 
                     "CD8B","CD8A","GZMK","GZMH",
                     "TRBC1","TRAC",
                     "CD3D",
                     "KLRB1","GNLY")

Idents(object = Round2) <- fct_rev(Idents(object = Round2))

pdf('Dot plot_Round2.pdf', width=14,height=6)
DotPlot(object = Round2, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=18))
dev.off()

Idents(object = Round2) <- fct_rev(Idents(object = Round2))

######################################################################################################
# (Figure 2) Visualization with label ####
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
DimPlot(object = Round2,label = TRUE,  reduction = "umap") +
  scale_colour_manual(values = c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','dark slate gray', 'indianred4','deepskyblue','royalblue1','navy','slate gray','gainsboro'), breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()
######################################################################################################
# (Supplementary Figure S4) Heatmap with labels ####
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
# (Figure 3) Average expression - Dendritic cells ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DCs <- subset(Round2, idents = c("Mature_DC" ,   "Semimature_DC", "Macrophage" ))

DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$stim, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(DCs) <- fct_relevel(Idents(DCs), sort)
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")

Idents(object = DCs) <- factor(Idents(object = DCs), levels = (c(  "Mature_DC_Control"  ,  "Mature_DC_Psoriasis"  ,
                                                                   "Semimature_DC_Control"  , "Semimature_DC_Psoriasis" , "Macrophage_Psoriasis" )))
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)

DC.cluster.averages <- AverageExpression(DCs, return.seurat = TRUE)

marker <- c("HLA-DRB5" , "HLA-DRA" ,  "HLA-DPB1" , "LAMP3"   ,  "LY75"   ,   "CD40"   ,   "CIITA" ,    "CD86"   ,   "CD80"   ,   "CD274"   ,  "PDCD1LG2" ,
            "KYNU"   ,   "IDO1"    ,  "CD209"  ,   "CD14"   ,   "AIF1"    ,
            "CD68"   ,   "CD1C"   ,   "ITGAX"   ,  "THBD"    ,  "F13A1"    , "FCGR3A"  ,  "CD163"   ,  "IL10"    ,  "IL23A"   ,  "TNFRSF10C" ,"LILRB1"  ,  "ITGAM" ,    "LILRB2" ,   "LILRB4"  )

pdf('Dendritic_cells.cluster.averages_Round2_Heatmp.pdf', width=6, height=6)
DoHeatmap(object = DC.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
 # scale_fill_viridis(option="magma")
 # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
   scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# (Figure 6) Average expression - NKT and T cells ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

NK.t.cells <- subset(Round2, idents = c("NK_cell" ,      "CD161_T_cell", "CD8_T_cell"   , "CD4_T_cell", "Treg"       ))
NK.t.cells$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(NK.t.cells), NK.t.cells$stim, sep = "_")
Idents(NK.t.cells) <- "ClusterNames_0.8_PsoriasisvsControl"
NK.t.cells[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = NK.t.cells)
Idents(object = NK.t.cells) <- ("ClusterNames_0.8_PsoriasisvsControl")
Idents(object = NK.t.cells) <- factor(Idents(object = NK.t.cells), levels = (c( 
  "NK_cell_Control","NK_cell_Psoriasis" ,
  "CD161_T_cell_Control","CD161_T_cell_Psoriasis",
  "CD8_T_cell_Control"   , "CD8_T_cell_Psoriasis",
  "CD4_T_cell_Control", "CD4_T_cell_Psoriasis",
  "Treg_Control","Treg_Psoriasis"      
)))
NK.t.cells[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = NK.t.cells)

NK.Tcell.cluster.averages <- AverageExpression(NK.t.cells, return.seurat = TRUE)

marker <- c("CD3D" ,   "CD4"   ,  "KLRB1"  , "PRF1"  ,  "GNLY"  ,  "CD8A"   , "CD8B"  ,
            "CXCR3" ,  "GZMK"   , "GZMH"   , "IFNG"  ,  "IL17A"  , "IL17F" ,  "IL26"  ,  "IL22"  ,
            "TIGIT"  , "IL2RA"  , "FOXP3" ,  "CTLA4"  , "TNFRSF4")

pdf('NK.Tcell.cluster.averages_Round2_Heatmp.pdf', width=12, height=8)
DoHeatmap(object = NK.Tcell.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# (Figure 7) Average expression - Keratinocytes ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

KCs <- subset(Round2, idents = c( c("S.Corneum",   "S.Granulosum" , "S.Spinosum"  ,  "S.Basale"  )))
KCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(KCs), KCs$stim, sep = "_")
Idents(object = KCs) <- factor(Idents(object = KCs), levels = (c("S.Corneum_Control" ,     "S.Corneum_Psoriasis"  ,  
                                                                 "S.Granulosum_Control"  ,    "S.Granulosum_Psoriasis"  ,
                                                                 "S.Spinosum_Control" , "S.Spinosum_Psoriasis",
                                                                 "S.Basale_Control",   "S.Basale_Psoriasis" )))
KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)

KC.cluster.averages <- AverageExpression(KCs, return.seurat = TRUE)

marker <- c("LCE3D", "CDSN" , "FLG" ,  "FABP5" ,"IL36G", "KRT1",  "KRT10" ,"KRT5" , "KRT14" ,"CCL27", "KRT15" ,"CD34" )

pdf('KC.cluster.averages_Round2_Heatmp.pdf', width=6, height=4)
DoHeatmap(object = KC.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# (Supplementary Figure S3) Flow diagram illustrating individual cell distribution ####
Idents(object = Round2) <- ("ClusterNames_0.8")
Round2$ClusterNames_0.8_Samples <- paste(Round2$number, Idents(Round2), sep = "_")
Idents(Round2) <- "ClusterNames_0.8_Samples"
Idents(Round2) <- fct_relevel(Idents(Round2), sort)

data <- as.data.frame(table(Round2@meta.data$ClusterNames_0.8_Samples))
data[is.na(cellnumbers)] <- 0
colnames(data) <- c("Cluster_in_sample","cell_number")
data$Sample <- str_extract(data$Cluster_in_sample, "[^_]+")

data$Category [grepl("Control",data$Sample)] <- "Control"
data$Category [grepl("Psoriasis",data$Sample)] <- "Psoriasis"

data$Cell_Type [grepl("Basal_KC",data$Cluster_in_sample)] <- "Basal_KC"
data$Cell_Type [grepl("CD161",data$Cluster_in_sample)] <- "CD161_T_cell"
data$Cell_Type [grepl("CD4_T_cell",data$Cluster_in_sample)] <- "CD4_T_cell"
data$Cell_Type [grepl("CD8_T_cell",data$Cluster_in_sample)] <- "CD8_T_cell"
data$Cell_Type [grepl("S.Corneum",data$Cluster_in_sample)] <- "S.Corneum"
data$Cell_Type [grepl("S.Basale",data$Cluster_in_sample)] <- "S.Basale"
data$Cell_Type [grepl("Mature_DC",data$Cluster_in_sample)] <- "Mature_DC"
data$Cell_Type [grepl("Melanocyte",data$Cluster_in_sample)] <- "Melanocyte"
data$Cell_Type [grepl("NK_cell",data$Cluster_in_sample)] <- "NK_cell"
data$Cell_Type [grepl("Semimature_DC",data$Cluster_in_sample)] <- "Semimature_DC"
data$Cell_Type [grepl("Treg",data$Cluster_in_sample)] <- "Treg"
data$Cell_Type [grepl("ECM",data$Cluster_in_sample)] <- "ECM"
data$Cell_Type [grepl("S.Granulosum",data$Cluster_in_sample)] <- "S.Granulosum"
data$Cell_Type [grepl("S.Spinosum",data$Cluster_in_sample)] <- "S.Spinosum"
data$Cell_Type [grepl("Macrophage",data$Cluster_in_sample)] <- "Macrophage"

PASI <- read.delim("~/Dropbox/SIngle_cell/Rcode_JKIM/github_upload/phenotype_data.txt")
data <- merge(data, PASI, by.x  ="Sample", by.y = "ID" ,all.x=TRUE)

data$severity <- ifelse(data$PASI>=12,"Severe psoriasis", "Mild-to-moderate psoriasis" )
data$severity <- ifelse(data$PASI==0, "Normal",data$severity)

data$Sample <- factor(data$Sample, levels = c("Psoriasis01", "Psoriasis02", "Psoriasis03" ,"Psoriasis04" ,"Psoriasis05", "Psoriasis06", "Psoriasis07", "Psoriasis08" ,"Psoriasis09","Psoriasis10", "Psoriasis11", "Psoriasis12" ,"Psoriasis13" ,"Control01" ,  "Control02" ,  "Control03",   "Control04"  , "Control05"  ))
data$Cell_Type <- factor(data$Cell_Type, levels = c("NK_cell", "CD161_T_cell", "CD8_T_cell" ,"CD4_T_cell" ,"Treg",  "Mature_DC" ,"Semimature_DC", "Melanocyte", "S.Corneum" , "S.Granulosum" ,  "S.Spinosum",
                                                    "S.Basale","ECM"))
data$severity <- factor(data$severity, levels = c( "Severe psoriasis","Mild-to-moderate psoriasis", "Normal" ))

A_col <- "firebrick4"
B_col <- "coral1"
C_col <- "cyan2"
alpha <- 0.7

p <- ggplot(data,
            aes(y = cell_number, axis1 = Sample, axis2 = severity, axis3 = Cell_Type)) +
  geom_alluvium(aes(fill = severity, color = severity), 
                width = 1/12, alpha = alpha, knot.pos = 0.4) +
  geom_stratum(width = 1/3, color = "gray") +
  geom_label(stat = "stratum", infer.label = TRUE) +
  scale_x_continuous(breaks = 1:3, labels = c("Sample","severity", "Cell_Type"))     +
  scale_fill_manual(values  = c(A_col, B_col, C_col)) +
  scale_color_manual(values = c(A_col, B_col, C_col)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold")
  )

ggsave(file='Flow_diagram.pdf',p,height=10,width=12)

###################################################################################################
## (Figure 4) UMAP plot - Psoriasis vs. Control ####
p <- DimPlot(Round2, reduction = "umap", group.by = "stim", cols = c("royalblue1","indianred1"))
pdf('UMAP_Psoriasis.vs.Control.pdf', width=10,height=7)
print(p)
dev.off()
###################################################################################################
# (Figure 6) UMAP plot - NK cell and T-cell clusters ####
Idents(object = Round2) <- ("ClusterNames_0.8")
cell.num <- table(Round2@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round2_cluster_Lymphocytes.pdf', width=12,height=8)
DimPlot(object = Round2,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_manual(values = c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'white','white','white', 'white','white','white','white','white','white'), breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

###################################################################################################
# (Figure 7) UMAP plot - keratinocyte clusters ####
Idents(object = Round2) <- ("ClusterNames_0.8")
cell.num <- table(Round2@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round2_cluster_Keratinocytes.pdf', width=12,height=8)
DimPlot(object = Round2,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_manual(values = c('white','white','white','white', 'white', 'white','white','white', 'white','deepskyblue','royalblue1','navy','slate gray','white'), breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()
###################################################################################################
# (Figure 4) IL-10 Feature plot and violin plot ####
pdf("IL10_FeaturePlot_Round2.pdf", width=4,height=3)
p <- FeaturePlot(object = Round2, features = "IL10",min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
print(p)
dev.off()

pdf("IL10_VlnPlot_Round2.pdf", width=12,height=4)
p <- VlnPlot(Round2, features="IL10",  split.by = "stim", cols = c("cyan", "red", "green"))
print(p)
dev.off()
#################################################################################################
# (figure 4) Dendritic cell subcluster feature plots with co-expression ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DCs <- subset(Round2, idents = c("Mature_DC" ,   "Semimature_DC" ))
DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$stim, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(DCs) <- fct_relevel(Idents(DCs), sort)
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")


pdf("THBD_IL10_FeaturePlot_DCs.pdf", width=12,height=4)
p <- FeaturePlot(object = DCs, features = c("THBD", "IL10"), blend = TRUE, label=FALSE, order=TRUE, max.cutoff = "q70") 
print(p)
dev.off()

pdf("CD1C_IL10_FeaturePlot_DCs.pdf", width=12,height=4)
p <- FeaturePlot(object = DCs, features = c("CD1C", "IL10"), blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q70" ) 
print(p)
dev.off()

pdf("LILRB2_IL10_FeaturePlot_DCs.pdf", width=12,height=4)
p <- FeaturePlot(object = DCs, features = c("LILRB2", "IL10"), blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q70" ) 
print(p)
dev.off()

pdf("IL23A_IL10_FeaturePlot_DCs_split.pdf", width=12,height=8)
p <- FeaturePlot(object = DCs, features = c("IL23A", "IL10"), blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q70" , split.by = 'stim') 
print(p)
dev.off()
#################################################################################################
# (Figure 5) IDO1 and KYNU Co-expression plot ####
pdf("IDO1_KYNU_FeaturePlot_Roun2.pdf", width=14,height=4)
p <- FeaturePlot(object = Round2, features = c("IDO1", "KYNU"), blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q70" ) 
print(p)
dev.off()
#################################################################################################
# (figure 6) Regulatory T cell  plots with co-expression ####
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Idents(object = Round2))
Treg <- subset(Round2, idents = "Treg")

pdf("IL17A_FOXP3_FeaturePlot_DCs.pdf", width=12,height=2)
p <- FeaturePlot(object = DCs, features = c("IL17A", "FOXP3"), blend = TRUE, label=FALSE, order=TRUE, max.cutoff = "q70") 
print(p)
dev.off()

pdf("IL17F_FOXP3_FeaturePlot_DCs.pdf", width=12,height=2)
p <- FeaturePlot(object = DCs, features = c("IL17F", "FOXP3"), blend = TRUE, label=FALSE, order=TRUE, max.cutoff = "q70") 
print(p)
dev.off()

#################################################################################################
# Save ####
save(Round2, file = "Round2_integrated_analyzed.Rda")
#################################################################################################