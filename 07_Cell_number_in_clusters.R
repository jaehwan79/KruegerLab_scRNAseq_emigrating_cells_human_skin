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
# Number of cells (psoriasis vs. control) in each cluster calculation   ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

####################################################################################
## Number of cells (psoriasis vs. control) in each cluster calculation - Round2 ####
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")

Idents(object = Round2) <- ("ClusterNames_0.8")
Round2$ClusterNames_0.8_stim <- paste(Idents(Round2),Round2$stim, sep = "_")
Idents(Round2) <- "ClusterNames_0.8_stim"
Idents(Round2) <- fct_relevel(Idents(Round2), sort)
levels(Idents(object = Round2))

list <- as.data.frame(table(Round2@meta.data$ClusterNames_0.8_stim))
list[is.na(list)] <- 0
colnames(list) <- c("Cluster_in_sample","cell_number")




list$Cluster <- list$Cluster_in_sample
list$Cluster= gsub("_Control", "", list$Cluster)
list$Cluster= gsub("_Psoriasis", "", list$Cluster)
list$Cluster <- factor(list$Cluster)

list$Cluster <- factor(list$Cluster, 
                       levels = (c( "NK_cell"  ,     "CD161_T_cell",  "CD8_T_cell"   , "CD4_T_cell" ,"Treg" ,    
                                    "Mature_DC" ,  "Semimature_DC", "Macrophage"   , "Melanocyte"   ,
                                    "S.Corneum"  , "S.Granulosum" , "S.Spinosum" ,   "S.Basale"  , "ECM"       
                       )))
list$PsoriasisvsControl <- NULL
list$PsoriasisvsControl <- ifelse(grepl("Psoriasis", list$Cluster_in_sample), "Psoriasis", "Control")
list$PsoriasisvsControl  <- factor(list$PsoriasisvsControl )
list$PsoriasisvsControl <- factor(list$PsoriasisvsControl, levels = c("Psoriasis","Control"))

list <- list[,c(1,3,4,2)]
write.csv(list, file="Cellnumbers_Round2.csv")

## Bar graph for clusters
list_for_ggplot <- melt(list)
ggplot(data = list_for_ggplot, aes(Cluster, value, fill = PsoriasisvsControl)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("PsoriasisvsControl", values = c("indianred1","royalblue1"))

ggsave("cluster_bargraph_Round2.pdf", width = 10, height = 4)

####################################################################################
## Number of cells (psoriasis vs. control) in each cluster calculation - Round1 ####
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round1_integrated_analyzed.Rda")

Idents(object = Round1) <- ("ClusterNames_0.8")
Round1$ClusterNames_0.8_stim <- paste(Idents(Round1),Round1$stim, sep = "_")
Idents(Round1) <- "ClusterNames_0.8_stim"
Idents(Round1) <- fct_relevel(Idents(Round1), sort)
levels(Idents(object = Round1))

list <- as.data.frame(table(Round1@meta.data$ClusterNames_0.8_stim))
list[is.na(list)] <- 0
colnames(list) <- c("Cluster_in_sample","cell_number")

list$Cluster <- list$Cluster_in_sample
list$Cluster= gsub("_Control", "", list$Cluster)
list$Cluster= gsub("_Psoriasis", "", list$Cluster)
list$Cluster <- factor(list$Cluster)

list$Cluster <- factor(list$Cluster, 
                       levels = (c( "NK_cell"  ,     "CD161_T_cell",  "CD8_T_cell"   , "CD4_T_cell" ,"Treg.01" ,"Treg.02",    
                                    "Mature_DC.01" , "Mature_DC.02" , "Mature_DC.03" , "Semimature_DC", "Macrophage"   , "Melanocyte"   ,
                                    "S.Corneum.01", "S.Corneum.02","S.Corneum.03","S.Corneum.04","S.Corneum.05",
                                    "S.Granulosum" , "S.Spinosum" ,   "S.Basale"  , "ECM"       
                       )))
list$PsoriasisvsControl <- NULL
list$PsoriasisvsControl <- ifelse(grepl("Psoriasis", list$Cluster_in_sample), "Psoriasis", "Control")
list$PsoriasisvsControl  <- factor(list$PsoriasisvsControl )
list$PsoriasisvsControl <- factor(list$PsoriasisvsControl, levels = c("Psoriasis","Control"))

list <- list[,c(1,3,4,2)]
write.csv(list, file="Cellnumbers_Round1.csv")

## Bar graph for clusters
list_for_ggplot <- melt(list)
ggplot(data = list_for_ggplot, aes(Cluster, value, fill = PsoriasisvsControl)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("PsoriasisvsControl", values = c("indianred1","royalblue1"))

ggsave("cluster_bargraph_Round1.pdf", width = 10, height = 4)

########################################################################################################################
## Number of cells (psoriasis vs. control) in each cluster calculation - Mature DC01 vs Mature DC02 comparison ####
Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

Mature_DCs <- subset(Round1, idents = c("Mature_DC.01" ,  "Mature_DC.02","Mature_DC.03" ))

new.cluster.ids <- c("Mature_DC.A" ,  "Mature_DC.B","Mature_DC.A")
names(new.cluster.ids) <- levels(Mature_DCs)
Mature_DCs <- RenameIdents(Mature_DCs, new.cluster.ids)

Mature_DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(Mature_DCs), Mature_DCs$stim, sep = "_")
Idents(Mature_DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(Mature_DCs) <- fct_relevel(Idents(Mature_DCs), sort)
Mature_DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = Mature_DCs)
Idents(object = Mature_DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")

list <- as.data.frame(table(Mature_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
list[is.na(list)] <- 0
colnames(list) <- c("Cluster_in_sample","cell_number")

list$Cluster <- list$Cluster_in_sample
list$Cluster= gsub("_Control", "", list$Cluster)
list$Cluster= gsub("_Psoriasis", "", list$Cluster)
list$Cluster <- factor(list$Cluster)

list$PsoriasisvsControl <- NULL
list$PsoriasisvsControl <- ifelse(grepl("Psoriasis", list$Cluster_in_sample), "Psoriasis", "Control")
list$PsoriasisvsControl  <- factor(list$PsoriasisvsControl )
list$PsoriasisvsControl <- factor(list$PsoriasisvsControl, levels = c("Psoriasis","Control"))

list <- list[,c(1,3,4,2)]
write.csv(list, file="Cellnumbers_Mature_DCs.csv")

## Bar graph for clusters
list_for_ggplot <- melt(list)
ggplot(data = list_for_ggplot, aes(Cluster, value, fill = PsoriasisvsControl)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("PsoriasisvsControl", values = c("indianred1","royalblue1"))

ggsave("cluster_bargraph_Mature_DCs.pdf", width = 4, height = 4)

####################################################################################
## Number of cells in each sample - Round0 ####
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round0_integrated_analyzed.Rda")

Idents(object = Round0) <- ("number")

list <- as.data.frame(table(Round0@meta.data$number))
list[is.na(list)] <- 0
colnames(list) <- c("Cluster_in_sample","cell_number")

write.csv(list, file="Cellnumbers_per_samples_Round0.csv")

