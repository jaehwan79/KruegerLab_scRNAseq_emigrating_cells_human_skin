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
######################################################################################################
# Number of cells (psoriasis vs. control) in each cluster calculation   ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

#####################################################################################################################
## DC Co-expression
## Upload DC cluster data ####
load("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")

Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"

DCs <- subset(Round2, idents = c("Mature_DC",  "Semimature_DC"      ))
DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$stim, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")
Idents(object = DCs) <- factor(Idents(object = DCs), levels = (c( 
  "Mature_DC_Control","Mature_DC_Psoriasis" ,
  "Semimature_DC_Control","Semimature_DC_Psoriasis"
)))
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)



## IL23A and IL36G co-expression ####
IL23A_DCs <- subset(DCs, subset = IL23A > 1)
IL23A_DCs <- subset(IL23A_DCs, subset = IL36G <1)
IL23A_list <- as.data.frame(table(IL23A_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL23A_list[is.na(IL23A_list)] <- 0
colnames(IL23A_list) <- c("Cluster","IL23A")

IL36G_DCs <- subset(DCs, subset = IL36G > 1)
IL36G_DCs <- subset(IL36G_DCs, subset = IL23A <1 )
IL36G_list <- as.data.frame(table(IL36G_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL36G_list[is.na(IL36G_list)] <- 0
colnames(IL36G_list) <- c("Cluster","IL36G")


IL23A_IL36G_DCs <- subset(DCs, subset = IL23A > 1)
IL23A_IL36G_DCs <- subset(IL23A_IL36G_DCs, subset = IL36G > 1)
IL23A_IL36G_list <- as.data.frame(table(IL23A_IL36G_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL23A_IL36G_list[is.na(IL23A_IL36G_list)] <- 0
colnames(IL23A_IL36G_list) <- c("Cluster","IL23A & IL36G")

list <- merge(IL23A_list, IL36G_list, by  ="Cluster")
list <- merge(list, IL23A_IL36G_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

#Bar graph for clusters (count) - Figure 6c
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL36G","IL23A & IL36G", "IL23A" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of IL23A and IL36G expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_IL23A_IL36G_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 6c
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL36G","IL23A & IL36G", "IL23A" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of IL23A and IL36G expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_IL23A_IL36G_proportion.pdf", width = 5, height = 4)



## IL23A and KYNU co-expression ####
IL23A_DCs <- subset(DCs, subset = IL23A > 1)
IL23A_DCs <- subset(IL23A_DCs, subset = KYNU <2)
IL23A_list <- as.data.frame(table(IL23A_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL23A_list[is.na(IL23A_list)] <- 0
colnames(IL23A_list) <- c("Cluster","IL23A")

KYNU_DCs <- subset(DCs, subset = KYNU > 2)
KYNU_DCs <- subset(KYNU_DCs, subset = IL23A <1 )
KYNU_list <- as.data.frame(table(KYNU_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
KYNU_list[is.na(KYNU_list)] <- 0
colnames(KYNU_list) <- c("Cluster","KYNU")


IL23A_KYNU_DCs <- subset(DCs, subset = IL23A > 1)
IL23A_KYNU_DCs <- subset(IL23A_KYNU_DCs, subset = KYNU > 2)
IL23A_KYNU_list <- as.data.frame(table(IL23A_KYNU_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL23A_KYNU_list[is.na(IL23A_KYNU_list)] <- 0
colnames(IL23A_KYNU_list) <- c("Cluster","IL23A & KYNU")

list <- merge(IL23A_list, KYNU_list, by  ="Cluster")
list <- merge(list, IL23A_KYNU_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#Bar graph for clusters (count) - Figure 6d
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("KYNU","IL23A & KYNU", "IL23A" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of IL23A and KYNU expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_IL23A_KYNU_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 6d
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("KYNU","IL23A & KYNU", "IL23A" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of IL23A and KYNU expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_IL23A_KYNU_proportion.pdf", width = 5, height = 4)


## IL36G and KYNU co-expression ####
IL36G_DCs <- subset(DCs, subset = IL36G > 1)
IL36G_DCs <- subset(IL36G_DCs, subset = KYNU <2)
IL36G_list <- as.data.frame(table(IL36G_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL36G_list[is.na(IL36G_list)] <- 0
colnames(IL36G_list) <- c("Cluster","IL36G")

KYNU_DCs <- subset(DCs, subset = KYNU > 2)
KYNU_DCs <- subset(KYNU_DCs, subset = IL36G <1 )
KYNU_list <- as.data.frame(table(KYNU_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
KYNU_list[is.na(KYNU_list)] <- 0
colnames(KYNU_list) <- c("Cluster","KYNU")


IL36G_KYNU_DCs <- subset(DCs, subset = IL36G > 1)
IL36G_KYNU_DCs <- subset(IL36G_KYNU_DCs, subset = KYNU > 2)
IL36G_KYNU_list <- as.data.frame(table(IL36G_KYNU_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL36G_KYNU_list[is.na(IL36G_KYNU_list)] <- 0
colnames(IL36G_KYNU_list) <- c("Cluster","IL36G & KYNU")

list <- merge(IL36G_list, KYNU_list, by  ="Cluster")
list <- merge(list, IL36G_KYNU_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

#Bar graph for clusters (count) - Supplementary Figure S5c
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("KYNU","IL36G & KYNU", "IL36G" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of IL36G and KYNU expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_IL36G_KYNU_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Supplementary Figure S5c
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("KYNU","IL36G & KYNU", "IL36G" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of IL36G and KYNU expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_IL36G_KYNU_proportion.pdf", width = 5, height = 4)




## THBD and IL10 co-expression ####
THBD_DCs <- subset(DCs, subset = THBD > 1)
THBD_DCs <- subset(THBD_DCs, subset = IL10 <1)
THBD_list <- as.data.frame(table(THBD_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
THBD_list[is.na(THBD_list)] <- 0
colnames(THBD_list) <- c("Cluster","THBD")

IL10_DCs <- subset(DCs, subset = IL10 > 1)
IL10_DCs <- subset(IL10_DCs, subset = THBD <1 )
IL10_list <- as.data.frame(table(IL10_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL10_list[is.na(IL10_list)] <- 0
colnames(IL10_list) <- c("Cluster","IL10")


THBD_IL10_DCs <- subset(DCs, subset = THBD > 1)
THBD_IL10_DCs <- subset(THBD_IL10_DCs, subset = IL10 > 1)
THBD_IL10_list <- as.data.frame(table(THBD_IL10_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
THBD_IL10_list[is.na(THBD_IL10_list)] <- 0
colnames(THBD_IL10_list) <- c("Cluster","THBD & IL10")

list <- merge(THBD_list, IL10_list, by  ="Cluster")
list <- merge(list, THBD_IL10_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#Bar graph for clusters (count) - Figure 7f
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL10","THBD & IL10", "THBD" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of THBD and IL10 expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_THBD_IL10_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 7f
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL10","THBD & IL10", "THBD" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of THBD and IL10 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_THBD_IL10_proportion.pdf", width = 5, height = 4)











## THBD and CD1C co-expression ####
THBD_DCs <- subset(DCs, subset = THBD > 1)
THBD_DCs <- subset(THBD_DCs, subset = CD1C <1)
THBD_list <- as.data.frame(table(THBD_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
THBD_list[is.na(THBD_list)] <- 0
colnames(THBD_list) <- c("Cluster","THBD")

CD1C_DCs <- subset(DCs, subset = CD1C > 1)
CD1C_DCs <- subset(CD1C_DCs, subset = THBD <1 )
CD1C_list <- as.data.frame(table(CD1C_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
CD1C_list[is.na(CD1C_list)] <- 0
colnames(CD1C_list) <- c("Cluster","CD1C")


THBD_CD1C_DCs <- subset(DCs, subset = THBD > 1)
THBD_CD1C_DCs <- subset(THBD_CD1C_DCs, subset = CD1C > 1)
THBD_CD1C_list <- as.data.frame(table(THBD_CD1C_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
THBD_CD1C_list[is.na(THBD_CD1C_list)] <- 0
colnames(THBD_CD1C_list) <- c("Cluster","THBD & CD1C")

list <- merge(THBD_list, CD1C_list, by  ="Cluster")
list <- merge(list, THBD_CD1C_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#Bar graph for clusters (count) - Figure 7e
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CD1C","THBD & CD1C", "THBD" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of THBD and CD1C expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_THBD_CD1C_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 7e
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CD1C","THBD & CD1C", "THBD" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of THBD and CD1C expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_THBD_CD1C_proportion.pdf", width = 5, height = 4)








## THBD and LILRB2 co-expression ####
THBD_DCs <- subset(DCs, subset = THBD > 1)
THBD_DCs <- subset(THBD_DCs, subset = LILRB2 <1)
THBD_list <- as.data.frame(table(THBD_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
THBD_list[is.na(THBD_list)] <- 0
colnames(THBD_list) <- c("Cluster","THBD")

LILRB2_DCs <- subset(DCs, subset = LILRB2 > 1)
LILRB2_DCs <- subset(LILRB2_DCs, subset = THBD <1 )
LILRB2_list <- as.data.frame(table(LILRB2_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
LILRB2_list[is.na(LILRB2_list)] <- 0
colnames(LILRB2_list) <- c("Cluster","LILRB2")


THBD_LILRB2_DCs <- subset(DCs, subset = THBD > 1)
THBD_LILRB2_DCs <- subset(THBD_LILRB2_DCs, subset = LILRB2 > 1)
THBD_LILRB2_list <- as.data.frame(table(THBD_LILRB2_DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
THBD_LILRB2_list[is.na(THBD_LILRB2_list)] <- 0
colnames(THBD_LILRB2_list) <- c("Cluster","THBD & LILRB2")

list <- merge(THBD_list, LILRB2_list, by  ="Cluster")
list <- merge(list, THBD_LILRB2_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(DCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#Bar graph for clusters (count) - Supplementary Figure S5d
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("LILRB2","THBD & LILRB2", "THBD" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of THBD and LILRB2 expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_THBD_LILRB2_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion)  - Supplementary Figure S5d
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("LILRB2","THBD & LILRB2", "THBD" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of THBD and LILRB2 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_THBD_LILRB2_proportion.pdf", width = 5, height = 4)



## DC violin plots - Supplementary Figure S5b
Idents(object = DCs) <- ("ClusterNames_0.8")

pdf("CD1C_VlnPlot_Round2.pdf", width=4,height=4)
p <- VlnPlot(DCs, features="CD1C",split.by = "stim", cols = c("cyan", "red", "cyan","red"))
print(p)
dev.off()

pdf("LY75_VlnPlot_Round2.pdf", width=4,height=4)
p <- VlnPlot(DCs, features="LY75",split.by = "stim", cols = c("cyan", "red", "cyan","red"))
print(p)
dev.off()


pdf("HLA-DPB1_VlnPlot_Round2.pdf", width=4,height=4)
p <- VlnPlot(DCs, features="HLA-DPB1",split.by = "stim", cols = c("cyan", "red"))
print(p)
dev.off()

pdf("CD40_VlnPlot_Round2.pdf", width=4,height=4)
p <- VlnPlot(DCs, features="CD40",split.by = "stim", cols = c("cyan", "red"))
print(p)
dev.off()

pdf("THBD_VlnPlot_Round2.pdf", width=4,height=4)
p <- VlnPlot(DCs, features="THBD",split.by = "stim", cols = c("cyan", "red"))
print(p)
dev.off()

pdf("LILRB2_VlnPlot_Round2.pdf", width=4,height=4)
p <- VlnPlot(DCs, features="LILRB2",split.by = "stim", cols = c("cyan", "red"))
print(p)
dev.off()

#################################################################################################
# KC co-expression  ####
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)
KCs <- subset(Round2, idents = c("S.Corneum"  , "S.Granulosum",  "S.Spinosum"    ,"S.Basale"      ))
KCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(KCs), KCs$stim, sep = "_")
Idents(KCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(KCs) <- fct_relevel(Idents(KCs), sort)
KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)
Idents(object = KCs) <- ("ClusterNames_0.8_PsoriasisvsControl")
levels(KCs)


Idents(object = KCs) <- factor(Idents(object = KCs), levels = (c(   "S.Corneum_Control"   ,"S.Corneum_Psoriasis"  ,
                                                                    "S.Granulosum_Control"  , "S.Granulosum_Psoriasis",
                                                                    "S.Spinosum_Control",  "S.Spinosum_Psoriasis" ,
                                                                    "S.Basale_Control"   ,    "S.Basale_Psoriasis"
)))
KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)




## IL36G and NFKBIZ co-expression ####
IL36G_KCs <- subset(KCs, subset = IL36G > 1)
IL36G_KCs <- subset(IL36G_KCs, subset = NFKBIZ <1)
IL36G_list <- as.data.frame(table(IL36G_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL36G_list[is.na(IL36G_list)] <- 0
colnames(IL36G_list) <- c("Cluster","IL36G")

NFKBIZ_KCs <- subset(KCs, subset = NFKBIZ > 1)
NFKBIZ_KCs <- subset(NFKBIZ_KCs, subset = IL36G <1 )
NFKBIZ_list <- as.data.frame(table(NFKBIZ_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
NFKBIZ_list[is.na(NFKBIZ_list)] <- 0
colnames(NFKBIZ_list) <- c("Cluster","NFKBIZ")

IL36G_NFKBIZ_KCs <- subset(KCs, subset = NFKBIZ > 1)
IL36G_NFKBIZ_KCs <- subset(IL36G_NFKBIZ_KCs, subset = IL36G > 1)
IL36G_NFKBIZ_list <- as.data.frame(table(IL36G_NFKBIZ_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
IL36G_NFKBIZ_list[is.na(IL36G_NFKBIZ_list)] <- 0
colnames(IL36G_NFKBIZ_list) <- c("Cluster","IL36G & NFKBIZ")

list <- merge(IL36G_list, NFKBIZ_list, by  ="Cluster")
list <- merge(list, IL36G_NFKBIZ_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

#Bar graph for clusters (count) - Figure 8c
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("NFKBIZ","IL36G & NFKBIZ", "IL36G" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of IL36G and NFKBIZ expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_IL36G_NFKBIZ_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 8c
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("NFKBIZ","IL36G & NFKBIZ", "IL36G" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of IL36G and NFKBIZ expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_IL36G_NFKBIZ_proportion.pdf", width = 5, height = 4)


## FLG -expression ####
levels(KCs)
FLG_KCs <- subset(KCs, subset = FLG > 1)
FLG_list <- as.data.frame(table(FLG_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
FLG_list[is.na(FLG_list)] <- 0
colnames(FLG_list) <- c("Cluster","FLG")


list <- FLG_list

total_cell_count <- as.data.frame(table(KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,4] <- list[,2]/list[,3]*100
colnames(list)[4] <- "FLG_proportion"

#Bar graph for clusters (count) - Figure 8b
list_for_ggplot <- list[c(1,2)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
#list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FLG","SPINK5 & FLG", "SPINK5" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of FLG expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("red")) +
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_FLG_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 8b
list_for_ggplot <- list[c(1,4)]

list_for_ggplot <- reshape2::melt(list_for_ggplot)
#list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FLG","SPINK5 & FLG", "SPINK5" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of FLG expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("red"))+
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_FLG_proportion.pdf", width = 5, height = 4)


## KRT15 and CCL27 co-expression ####
KRT15_KCs <- subset(KCs, subset = KRT15 > 1)
KRT15_KCs <- subset(KRT15_KCs, subset = CCL27 <1)
KRT15_list <- as.data.frame(table(KRT15_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
KRT15_list[is.na(KRT15_list)] <- 0
colnames(KRT15_list) <- c("Cluster","KRT15")

CCL27_KCs <- subset(KCs, subset = CCL27 > 1)
CCL27_KCs <- subset(CCL27_KCs, subset = KRT15 <1 )
CCL27_list <- as.data.frame(table(CCL27_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
CCL27_list[is.na(CCL27_list)] <- 0
colnames(CCL27_list) <- c("Cluster","CCL27")

KRT15_CCL27_KCs <- subset(KCs, subset = CCL27 > 1)
KRT15_CCL27_KCs <- subset(KRT15_CCL27_KCs, subset = KRT15 > 1)
KRT15_CCL27_list <- as.data.frame(table(KRT15_CCL27_KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
KRT15_CCL27_list[is.na(KRT15_CCL27_list)] <- 0
colnames(KRT15_CCL27_list) <- c("Cluster","KRT15 & CCL27")

list <- merge(KRT15_list, CCL27_list, by  ="Cluster")
list <- merge(list, KRT15_CCL27_list, by  ="Cluster")

total_cell_count <- as.data.frame(table(KCs@meta.data$ClusterNames_0.8_PsoriasisvsControl))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#Bar graph for clusters (count) - Figure 8d
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CCL27","KRT15 & CCL27", "KRT15" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of KRT15 and CCL27 expressing cells")+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) +
  theme(axis.title.x = element_blank())+   theme(axis.title.x=element_blank(),         axis.text.x=element_blank(),         axis.ticks.x=element_blank())

ggsave("cluster_bargraph_DC_KRT15_CCL27_count.pdf", width = 5, height = 4)

#Bar graph for clusters (proportion) - Figure 8d
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CCL27","KRT15 & CCL27", "KRT15" )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_light() +   theme(aspect.ratio = 1/1) +   theme(plot.title = element_text(size = 10, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of KRT15 and CCL27 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))+
  theme(axis.title.x = element_blank())

ggsave("cluster_bargraph_DC_KRT15_CCL27_proportion.pdf", width = 5, height = 4)

#################################################################################################
# NK cell and T-cell subcluster feature plots with co-expression Plot ####
Idents(object = Round2) <- ("ClusterNames_0.8")
DefaultAssay(Round2) <- "RNA"
NK.t.cells <- subset(Round2, idents = c("NK_cell" ,      "CD161_T_cell", "CD8_T_cell"   , "CD4_T_cell", "Treg"       ))

pdf("IL17A_IL17F_FeaturePlot_NK_T_cells.pdf", width=12,height=4) #Figure 4a
p <- FeaturePlot(object = NK.t.cells, features = c("IL17A", "IL17F"), blend = TRUE, label=FALSE, order=TRUE, max.cutoff = "q90") 
print(p)
dev.off()

#################################################################################################
# DC co-expression Plot #### Figure 7c and 7d
Idents(object = Round2) <- ("ClusterNames_0.8")
DCs <- subset(Round2, idents = c("Mature_DC" ,   "Semimature_DC" ))
DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$stim, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(DCs) <- fct_relevel(Idents(DCs), sort)
DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")


pdf("THBD_IL10_FeaturePlot_DCs.pdf", width=12,height=8)
p <- FeaturePlot(object = DCs, features = c("THBD", "IL10"), split.by = "stim", blend = TRUE, label=FALSE, order=TRUE, max.cutoff = "q90") 
print(p)
dev.off()

pdf("CD1C_IL10_FeaturePlot_DCs.pdf", width=12,height=8)
p <- FeaturePlot(object = DCs, features = c("CD1C", "IL10"), split.by = "stim",blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q90" ) 
print(p)
dev.off()

pdf("LILRB2_IL10_FeaturePlot_DCs.pdf", width=12,height=8)
p <- FeaturePlot(object = DCs, features = c("LILRB2", "IL10"), split.by = "stim", blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q90" ) 
print(p)
dev.off()

pdf("IL23A_IL10_FeaturePlot_DCs_split.pdf", width=12,height=8)
p <- FeaturePlot(object = DCs, features = c("IL23A", "IL10"), split.by = "stim",blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q90") 
print(p)
dev.off()


pdf("IL23A_IL36G_FeaturePlot_DCs_split.pdf", width=12,height=8)
p <- FeaturePlot(object = DCs, features = c("IL23A", "IL36G"), split.by = "stim",blend = TRUE, label=FALSE, order=TRUE,  max.cutoff = "q90") 
print(p)
dev.off()

