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
## T-cell subset comparison ####
rm(list= ls()[!(ls() %in% c('Round2'))])

Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

#T_cell_except_Treg
T_cell_except_Treg <- subset(Round2, idents = c("NK_cell"  , "CD161_T_cell",  "CD8_T_cell"  ,  "CD4_T_cell"))

#FOXP3- Treg
FOXP3.neg.Treg <- subset(Round2, idents = c("Treg"))
FOXP3.neg.Treg <- subset(FOXP3.neg.Treg, subset = FOXP3 < 1)
levels(FOXP3.neg.Treg)
new.cluster.ids <- c("FOXP3.neg.Treg")
names(new.cluster.ids) <- levels(FOXP3.neg.Treg)
FOXP3.neg.Treg <- RenameIdents(FOXP3.neg.Treg, new.cluster.ids)
levels(Idents(FOXP3.neg.Treg))

#FOXP3+ Treg
FOXP3.pos.Treg <- subset(Round2, idents = c("Treg"))
FOXP3.pos.Treg <- subset(FOXP3.pos.Treg, subset = FOXP3 >1)
levels(FOXP3.pos.Treg)
new.cluster.ids <- c("FOXP3.pos.Treg")
names(new.cluster.ids) <- levels(FOXP3.pos.Treg)
FOXP3.pos.Treg <- RenameIdents(FOXP3.pos.Treg, new.cluster.ids)
levels(Idents(FOXP3.pos.Treg))


T_cell_clusters <- merge(x= T_cell_except_Treg, y = c(FOXP3.neg.Treg,FOXP3.pos.Treg), add.cell.ids=c( "T_cell_except_Treg", "FOXP3.neg.Treg","FOXP3.pos.Treg"), project = "T_cells")


T_cell_clusters$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(T_cell_clusters), T_cell_clusters$stim, sep = "_")
Idents(T_cell_clusters) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(T_cell_clusters) <- fct_relevel(Idents(T_cell_clusters), sort)
T_cell_clusters[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = T_cell_clusters)


Idents(object = T_cell_clusters) <- factor(Idents(object = T_cell_clusters), levels = (c( 
  "NK_cell_Control"     ,     "NK_cell_Psoriasis" ,
  "CD161_T_cell_Control"   ,  "CD161_T_cell_Psoriasis" ,
  "CD8_T_cell_Control"    ,   "CD8_T_cell_Psoriasis"    ,
  "CD4_T_cell_Control"    ,   "CD4_T_cell_Psoriasis"     ,
  "FOXP3.neg.Treg_Control"  , "FOXP3.neg.Treg_Psoriasis",
  "FOXP3.pos.Treg_Control" ,  "FOXP3.pos.Treg_Psoriasis"  
)))

#####################################################################
# Heatmap for average gene expression #### Figure 3
levels(T_cell_clusters)
T_cell_clusters.averages <- AverageExpression(T_cell_clusters, return.seurat = FALSE)


#marker_refined <- marker[1:(which(marker == "END")-1)]
#marker_refined <- as.character(marker_refined)

marker <- c("CD3D" ,  "TRAC",  "TRBC1",  "CD4"   ,  "KLRB1"  ,"IFNL1",  "GNLY"  , "PRF1"  ,   "CD8A"   , "CD8B"  ,
            "GZMK"   , "GZMH"   , "TNFSF10", "TNFSF11", "TNFSF12",  "CCL20" , "IL36G", "TNF", "IL13", 
            "TNFSF14", "LTA",  "IL24", "IL34", "EBI3", "CCR6", "STAT3", "RORC", "IL7R", "IL10", "IL6", "IL23A", "IL33", "CCR10",
            "IFNG"  ,  "IL17A"  , "IL17F" ,  "IL26"  ,  "IL22",    "TIGIT"  , "IL2RA"  , "FOXP3" ,  "CTLA4"  , "TNFRSF4")

a <- as.data.frame(T_cell_clusters.averages$RNA)
b <- a[as.character(marker),]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(b, 1, cal_z_score))
data_subset_norm <- na.omit(data_subset_norm)

pdf("Heatmap_gene_expression_T_cell_subset_Average_expression02.pdf", width=6,height=8)
p <-   pheatmap(data_subset_norm, display_numbers = F,cluster_cols = F, cluster_rows=T
                , colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256)
                ,number_format = "%.1f", main = "Average gene expression",border_color=T ,
                scale = "row")
print(p)
dev.off()


#####################################################################
# Heatmap for cell proportion #### Supplemenatry Figure S4
list <- as.data.frame(table(T_cell_clusters@active.ident))
list[is.na(list)] <- 0
colnames(list) <- c("Cluster","total_cell_number")
total_number <- sum(list[2])

count(marker)

for (i in 1:length(marker)){
  string <- paste("list2 <- subset(T_cell_clusters, subset =`",marker[i],"`> 1)", sep="")
  eval(parse( text= string ))
  list2 <- as.data.frame(table(list2@active.ident))
  list2[is.na(list2)] <- 0
  colnames(list2) <- c("Cluster", as.character(marker[i]))
  list <- merge(list, list2, by  ="Cluster", all.x = TRUE)
  list[is.na(list)] <- 0
}

#list <- list[c(1,3,2,4,5),]
list[is.na(list)] <- 0

list_proportion_cluster <- list
list_proportion_cluster[3:ncol(list_proportion_cluster)] <- list_proportion_cluster[,3:ncol(list_proportion_cluster)] / list_proportion_cluster[,2] *100
list_proportion_cluster[is.na(list_proportion_cluster)] <- 0
#list_proportion_cluster <- as.data.frame(t(list_proportion_cluster))
rownames(list_proportion_cluster) <- list_proportion_cluster[,1]
list_proportion_cluster <- list_proportion_cluster[,-c(1,2)]

list_proportion_cluster <- as.data.frame(t(list_proportion_cluster))



pdf("Heatmap_gene_expression_T_cell_subset_proportion_expression.pdf", width=6,height=8)
p <-   pheatmap(list_proportion_cluster, display_numbers = T,cluster_cols = F, cluster_rows=T
                ,colorRampPalette(c('white','red'))(20)
                , number_format = "%.1f", main = "Proportion of gene expressing cells in Cluster (%)" )

print(p)
dev.off()
