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
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")

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
Tcells <- subset(Round2, idents = c( "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis",  "Treg_Psoriasis"  ))
IL17A.Tcells <- subset(Tcells, subset = IL17A > 1)
IL17A.Tcells <- subset(IL17A.Tcells, subset = IL17F < 1)
levels(IL17A.Tcells)
new.cluster.ids <- c("IL17A_Tcells","IL17A_Tcells","IL17A_Tcells","IL17A_Tcells")
names(new.cluster.ids) <- levels(IL17A.Tcells)
IL17A.Tcells <- RenameIdents(IL17A.Tcells, new.cluster.ids)
levels(Idents(IL17A.Tcells))

#IL17A- IL17F+ T-cell subset
Tcells <- subset(Round2, idents = c( "CD161_T_cell_Psoriasis", "CD8_T_cell_Psoriasis"   , "CD4_T_cell_Psoriasis",  "Treg_Psoriasis"  ))
IL17F.Tcells <- subset(Tcells, subset = IL17F > 1)
IL17F.Tcells <- subset(IL17F.Tcells, subset = IL17A < 1)
levels(IL17F.Tcells)
new.cluster.ids <- c("IL17F_Tcells","IL17F_Tcells","IL17F_Tcells")
names(new.cluster.ids) <- levels(IL17F.Tcells)
IL17F.Tcells <- RenameIdents(IL17F.Tcells, new.cluster.ids)
levels(Idents(IL17F.Tcells))


T17_cells <- merge(x= IL17A.Tcells, y = c(IL17F.Tcells ), project = "T17_cells")

levels(Idents(T17_cells))
#####################################################################
# Heatmap for average gene expression ####
T17_cells.averages <- AverageExpression(T17_cells, return.seurat = FALSE)

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


a <- as.data.frame(T17_cells.averages$RNA)
b <- a[marker,]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(b, 1, cal_z_score))
data_subset_norm <- na.omit(data_subset_norm)


pdf("Heatmap_gene_expression_T17_cells_sAverage_expression11.pdf", width=6,height=10)
p <-   pheatmap(data_subset_norm, display_numbers = F,cluster_cols = F, cluster_rows=T
                ,   colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(512)
                ,number_format = "%.0f"
                , main = "T17 cells:Average expression",border_color=T ,
                clustering_method = "complete",
                scale = "row")
print(p)
dev.off()


#########################################################################
# T17 cell volcano plots
#Idents(object = T17_cells) <- ("ClusterNames_0.8")
#DefaultAssay(T17_cells) <- "RNA"
levels(T17_cells)


# Volcano_plot control vs fresh_tissue_enzyme_digestion
DEG <- FindMarkers(T17_cells, ident.1 =  "IL17A_Tcells" , ident.2 ="IL17F_Tcells", logfc.threshold = 0.1,min.pct =0 )


marker <- c(  
  "CD4"   , 
  "GNLY"  , "PRF1"  ,  
  "CD8B"  ,
   "LTA",  "IL24", "IL34",  "CCR6", 
  "IFNG"  ,  "IL17A"  , "IL17F" ,  
  "CCR10",
  "FOXP3" , 
  "MAF",
  "CTLA4", 
  "CD86", 
  "GZMA" , 
  "IL2", 
  "IL1B", 
  "SLC7A8",  
   "BATF",  "CCL4", "GZMB",
  "NT5E",
  "CCL5",
  "LAG3",
  "IL1RN",
  "ITGA2",
  "IL10"
)



pdf("Volcanoplot_T17_cell_IL17AvsIL17F.pdf", width=14,height=10)
p <- EnhancedVolcano(DEG,
                     lab = rownames(DEG),
                     x = 'avg_log2FC',
                     y = 'p_val',
                     title = 'Control: Emigrating cell vs. Fresh skin enzyme dissociation',
                     selectLab = marker,
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     pointSize = 3.0,
                     labSize = 6.0,
                     cutoffLineWidth = 0.8,
                     colAlpha = 1,
                     legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                    'p-value & Log (base 2) FC'),
                     legendPosition = 'right',
                     legendLabSize = 16,
                     legendIconSize = 5.0,
                     xlab = bquote(~Log[2]~ 'fold change'),
                     drawConnectors = TRUE,
                     widthConnectors = 1,
                     boxedLabels = TRUE 
)


print(p)
dev.off()




