#load library#####
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
library(iTALK)
library(dplyr)
library(data.table)
library(igraph)
library(circlize)
######################################################################################################
# Receptor-ligand interaction analysis ####
# We used iTALK R package (version 0.1.0, https://github.com/Coolgenome/iTALK/)
# for the receptor-ligand interaction analysis. 
# Differentially expressed genes between psoriasis vs. control 
# in each cluster were found by the FindAllMarkers function in Seurat R package. 
# Then, the fold changes of differentially expressed genes were entered into the FindLR function 
# in iTALK R package to calculate receptor-ligand interaction between different cell clusters. 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/KruegerLab_github")
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/KruegerLab_github/Round2_integrated_analyzed.Rda")
######################################################################################################
## Identify differential expressed genes between psoriasis vs. control in each cluster ####
# Differentially expressed genes between psoriasis vs. control 
# in each cluster were found by the FindAllMarkers function in Seurat R package. 
DefaultAssay(Round2) <- "RNA"
Idents(object = Round2) <- ("ClusterNames_0.8")

marker <- levels(Idents(object=Round2))
Round2$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(Round2), Round2$stim, sep = "_")
Idents(Round2) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(Round2) <- fct_relevel(Idents(Round2), sort)
Round2[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = Round2)
Idents(object = Round2) <- ("ClusterNames_0.8_PsoriasisvsControl")

for (i in c(1,2,3,4,5,6,7,9,10,11,12,13,14 
            #differential expressed genes of macrophage (8) between psoriasis vs. control is not calculated, since there is no macrophage in control
            )){table <- FindMarkers(Round2, ident.1 = paste(marker[i],"_Psoriasis",sep=''), ident.2 = paste(marker[i],"_Control",sep=''), verbose = FALSE)
    assign (marker[i], table)
    rm(table)
}

######################################################################################################
## Combine differentially expressed genes in all clusters ####
setDT(NK_cell, keep.rownames = TRUE)
NK_cell$cell_type <- "NK_cell"

setDT(CD161_T_cell, keep.rownames = TRUE)
CD161_T_cell$cell_type <- "CD161_T_cell"

setDT(CD8_T_cell, keep.rownames = TRUE)
CD8_T_cell$cell_type <- "CD8_T_cell"

setDT(CD4_T_cell, keep.rownames = TRUE)
CD4_T_cell$cell_type <- "CD4_T_cell"

setDT(Treg, keep.rownames = TRUE)
Treg$cell_type <- "Treg"

setDT(Mature_DC, keep.rownames = TRUE)
Mature_DC$cell_type <- "Mature_DC"

setDT(Semimature_DC, keep.rownames = TRUE)
Semimature_DC$cell_type <- "Semimature_DC"

setDT(Melanocyte, keep.rownames = TRUE)
Melanocyte$cell_type <- "Melanocyte"

setDT(S.Corneum, keep.rownames = TRUE)
S.Corneum$cell_type <- "S.Corneum"

setDT(S.Granulosum, keep.rownames = TRUE)
S.Granulosum$cell_type <- "S.Granulosum"

setDT(S.Spinosum, keep.rownames = TRUE)
S.Spinosum$cell_type <- "S.Spinosum"

setDT(S.Basale, keep.rownames = TRUE)
S.Basale$cell_type <- "S.Basale"

setDT(ECM, keep.rownames = TRUE)
ECM$cell_type <- "ECM"

data <- rbind (NK_cell, CD161_T_cell,CD8_T_cell,CD4_T_cell,Treg,Mature_DC,Semimature_DC,Melanocyte,S.Corneum,S.Granulosum,S.Spinosum,S.Basale,ECM)
data <- data[,c(1,7,3,2,6)]
colnames(data) <- c("gene", "cell_type","logFC","p.value","q.value")

data$cell_type <- as.character(data$cell_type)
data$gene <- as.character(data$gene)

######################################################################################################
# Receptor-ligand interaction analysis ####
# We used iTALK R package (version 0.1.0, https://github.com/Coolgenome/iTALK/)
# for the receptor-ligand interaction analysis. 
# The FindLR function in iTALK R package to calculate receptor-ligand interaction between different cell clusters. 

comm_list<-c('growth factor','other','cytokine','checkpoint')
res<-NULL

for (i in c(1,2,3,4)){
    comm_type<- comm_list[i]
    res_cat<-FindLR(data,datatype='DEG',comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]  
    res_cat <- res_cat[res_cat$cell_from_logFC*res_cat$cell_to_logFC>0,]
    res_cat <- res_cat[res_cat$cell_from_logFC>0,]
    res_cat <- res_cat[!(res_cat$cell_from==res_cat$cell_to),]
    res<-rbind(res,res_cat)
}

result_original <- res

######################################################################################################
# (Figure 4) Top 50 putative receptor interaction of mature and semimature dendritic cells 
# with other skin inflammatory cells - DC "receptor" interaction ####

#select DC "receptor" interaction
res <- res[grepl("DC", res$cell_to),]

#Select the expression of both receptor and ligand were increased in psoriasis compared to control
res <- res[res$cell_from_logFC*res$cell_to_logFC>0,]
res <- res[res$cell_from_logFC>0,]
res <- res[!(res$cell_from==res$cell_to),]

#Common interaction involving HLA, CD3, CD4 molecules is excluded for the downstream analysis
res <- res[-grep("HLA", res$receptor),]
res <- res[-grep("CD3", res$receptor),]
res <- res[-grep("CD4", res$receptor),]
res <- res[-grep("HLA", res$ligand),]

res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),]

limit <- nrow(res)
if (limit >50){limit <- 50}

pdf(file=paste("DEG_PsoriasisvsControl_DC_receptor",".pdf",sep=''),width=8,height=8)
cell_col<-structure(c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','indianred4','deepskyblue','royalblue1','blue','navy','black'),names=unique(data$cell_type))
LRPlot(res[1:limit,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:limit],link.arr.width=res$cell_to_logFC[1:limit])
dev.off()

######################################################################################################
# (Figure 4) Top 50 putative receptor interaction of mature and semimature dendritic cells 
# with other skin inflammatory cells - DC "ligand" interaction ####

res <- result_original
#select DC "ligand" interaction
res <- res[grepl("DC", res$cell_from),]

#Select the expression of both receptor and ligand were increased in psoriasis compared to control
res <- res[res$cell_from_logFC*res$cell_to_logFC>0,]
res <- res[res$cell_from_logFC>0,]
res <- res[!(res$cell_from==res$cell_to),]

#Common interaction involving HLA, CD3, CD4 molecules is excluded for the downstream analysis
res <- res[-grep("HLA", res$receptor),]
res <- res[-grep("CD3", res$receptor),]
res <- res[-grep("CD4", res$receptor),]
res <- res[-grep("HLA", res$ligand),]

#Irrelevant interaction to psoriasis immunopathogenesis involving GNAI2 and GNAS is excluded for the downstream analysis.
res <- res[-grep("GNAI2", res$ligand),]
res <- res[-grep("GNAS", res$ligand),]
res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),]

limit <- nrow(res)
if (limit >50){limit <- 50}

pdf(file=paste("DEG_PsoriasisvsControl_DC_ligand",".pdf",sep=''),width=8,height=8)
cell_col<-structure(c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','indianred4','deepskyblue','royalblue1','blue','navy','black'),names=unique(data$cell_type))
LRPlot(res[1:limit,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:limit],link.arr.width=res$cell_to_logFC[1:limit])
dev.off()