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
#install.packages("Hmisc")
#install.packages("survival")
library(survival)
library(Hmisc)
#install.packages("corrplot")
library(corrplot)
#if(!require(devtools)) install.packages("devtools")
install.packages("broom")
library(broom)
#devtools::install_github("kassambara/ggpubr")
install.packages("ggpubr")
library("ggpubr")
library(dplyr)
######################################################################################################
# Correlation between scRNA-seq of emigrating cells vs. scRNA-seq of frozen tissue by enzyme digestion  ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

##########################################################################################
Psoriasis02_emigrating_cell <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis02.csv")
Psoriasis02_enzyme_digestion <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis02F.csv")
Psoriasis03_emigrating_cell <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis03.csv")
Psoriasis03_enzyme_digestion <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis03F.csv")
Psoriasis04_emigrating_cell <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis04.csv")
Psoriasis04_enzyme_digestion <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis04F.csv")
Psoriasis06_emigrating_cell <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis06.csv")
Psoriasis06_enzyme_digestion <-read.csv("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis06F.csv")



#Psoriasis02
Psoriasis02_emigrating_cell$cell_type[grepl("T_cell", Psoriasis02_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis02_emigrating_cell$cell_type[grepl("Treg", Psoriasis02_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis02_emigrating_cell$cell_type[grepl("DC", Psoriasis02_emigrating_cell$Var1) ] <- "Dendritic_cell"
Psoriasis02_emigrating_cell$cell_type[grepl("Melanocyte", Psoriasis02_emigrating_cell$Var1) ] <- "Melanocyte"
Psoriasis02_emigrating_cell$cell_type[grepl("NK_cell", Psoriasis02_emigrating_cell$Var1) ] <- "NK_cell"
Psoriasis02_emigrating_cell$cell_type[grepl("Basale", Psoriasis02_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis02_emigrating_cell$cell_type[grepl("Corneum", Psoriasis02_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis02_emigrating_cell$cell_type[grepl("Spinosum", Psoriasis02_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis02_emigrating_cell02 <- ddply(Psoriasis02_emigrating_cell,"cell_type",numcolwise(sum))
Psoriasis02_emigrating_cell02$emigrating_cell_proportion <- Psoriasis02_emigrating_cell02$Freq / sum(Psoriasis02_emigrating_cell$Freq) *100
Psoriasis02_emigrating_cell02
Psoriasis02_emigrating_cell02 <- Psoriasis02_emigrating_cell02[,c(1,4)]

Psoriasis02_enzyme_digestion$cell_type[grepl("T_cell", Psoriasis02_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis02_enzyme_digestion$cell_type[grepl("Treg", Psoriasis02_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis02_enzyme_digestion$cell_type[grepl("DC", Psoriasis02_enzyme_digestion$Var1) ] <- "Dendritic_cell"
Psoriasis02_enzyme_digestion$cell_type[grepl("Melanocyte", Psoriasis02_enzyme_digestion$Var1) ] <- "Melanocyte"
Psoriasis02_enzyme_digestion$cell_type[grepl("NK_cell", Psoriasis02_enzyme_digestion$Var1) ] <- "NK_cell"
Psoriasis02_enzyme_digestion$cell_type[grepl("Basale", Psoriasis02_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis02_enzyme_digestion$cell_type[grepl("Corneum", Psoriasis02_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis02_enzyme_digestion$cell_type[grepl("Spinosum", Psoriasis02_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis02_enzyme_digestion02 <- ddply(Psoriasis02_enzyme_digestion,"cell_type",numcolwise(sum))
Psoriasis02_enzyme_digestion02$enzyme_digestion_proportion <- Psoriasis02_enzyme_digestion02$Freq / sum(Psoriasis02_enzyme_digestion$Freq) *100
Psoriasis02_enzyme_digestion02
Psoriasis02_enzyme_digestion02 <- Psoriasis02_enzyme_digestion02[,c(1,4)]

Psoriasis02 <- data.frame("sample" = "Psoriasis02", "cell_type" = c("NK_cell","T_cell","Dendritic_cell","Melanocyte","Keratinocyte"))
Psoriasis02 <- merge(Psoriasis02, Psoriasis02_emigrating_cell02, by='cell_type',all.x=TRUE)
Psoriasis02 <- merge(Psoriasis02, Psoriasis02_enzyme_digestion02, by='cell_type',all.x=TRUE)
Psoriasis02[is.na(Psoriasis02)] <- 0


#Psoriasis03
Psoriasis03_emigrating_cell$cell_type[grepl("T_cell", Psoriasis03_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis03_emigrating_cell$cell_type[grepl("Treg", Psoriasis03_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis03_emigrating_cell$cell_type[grepl("DC", Psoriasis03_emigrating_cell$Var1) ] <- "Dendritic_cell"
Psoriasis03_emigrating_cell$cell_type[grepl("Melanocyte", Psoriasis03_emigrating_cell$Var1) ] <- "Melanocyte"
Psoriasis03_emigrating_cell$cell_type[grepl("NK_cell", Psoriasis03_emigrating_cell$Var1) ] <- "NK_cell"
Psoriasis03_emigrating_cell$cell_type[grepl("Basale", Psoriasis03_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis03_emigrating_cell$cell_type[grepl("Corneum", Psoriasis03_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis03_emigrating_cell$cell_type[grepl("Spinosum", Psoriasis03_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis03_emigrating_cell$cell_type[grepl("Granulosum", Psoriasis03_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis03_emigrating_cell02 <- ddply(Psoriasis03_emigrating_cell,"cell_type",numcolwise(sum))
Psoriasis03_emigrating_cell02$emigrating_cell_proportion <- Psoriasis03_emigrating_cell02$Freq / sum(Psoriasis03_emigrating_cell$Freq) *100
Psoriasis03_emigrating_cell02
Psoriasis03_emigrating_cell02 <- Psoriasis03_emigrating_cell02[,c(1,4)]

Psoriasis03_enzyme_digestion$cell_type[grepl("T_cell", Psoriasis03_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis03_enzyme_digestion$cell_type[grepl("Treg", Psoriasis03_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis03_enzyme_digestion$cell_type[grepl("DC", Psoriasis03_enzyme_digestion$Var1) ] <- "Dendritic_cell"
Psoriasis03_enzyme_digestion$cell_type[grepl("Melanocyte", Psoriasis03_enzyme_digestion$Var1) ] <- "Melanocyte"
Psoriasis03_enzyme_digestion$cell_type[grepl("NK_cell", Psoriasis03_enzyme_digestion$Var1) ] <- "NK_cell"
Psoriasis03_enzyme_digestion$cell_type[grepl("Basale", Psoriasis03_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis03_enzyme_digestion$cell_type[grepl("Corneum", Psoriasis03_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis03_enzyme_digestion$cell_type[grepl("Spinosum", Psoriasis03_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis03_enzyme_digestion$cell_type[grepl("Granulosum", Psoriasis03_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis03_enzyme_digestion02 <- ddply(Psoriasis03_enzyme_digestion,"cell_type",numcolwise(sum))
Psoriasis03_enzyme_digestion02$enzyme_digestion_proportion <- Psoriasis03_enzyme_digestion02$Freq / sum(Psoriasis03_enzyme_digestion$Freq) *100
Psoriasis03_enzyme_digestion02
Psoriasis03_enzyme_digestion02 <- Psoriasis03_enzyme_digestion02[,c(1,4)]

Psoriasis03 <- data.frame("sample" = "Psoriasis03", "cell_type" = c("NK_cell","T_cell","Dendritic_cell","Melanocyte","Keratinocyte"))
Psoriasis03 <- merge(Psoriasis03, Psoriasis03_emigrating_cell02, by='cell_type',all.x=TRUE)
Psoriasis03 <- merge(Psoriasis03, Psoriasis03_enzyme_digestion02, by='cell_type',all.x=TRUE)
Psoriasis03[is.na(Psoriasis03)] <- 0



#Psoriasis04
Psoriasis04_emigrating_cell$cell_type[grepl("T_cell", Psoriasis04_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis04_emigrating_cell$cell_type[grepl("Treg", Psoriasis04_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis04_emigrating_cell$cell_type[grepl("DC", Psoriasis04_emigrating_cell$Var1) ] <- "Dendritic_cell"
Psoriasis04_emigrating_cell$cell_type[grepl("Melanocyte", Psoriasis04_emigrating_cell$Var1) ] <- "Melanocyte"
Psoriasis04_emigrating_cell$cell_type[grepl("NK_cell", Psoriasis04_emigrating_cell$Var1) ] <- "NK_cell"
Psoriasis04_emigrating_cell$cell_type[grepl("Basale", Psoriasis04_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis04_emigrating_cell$cell_type[grepl("Corneum", Psoriasis04_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis04_emigrating_cell$cell_type[grepl("Spinosum", Psoriasis04_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis04_emigrating_cell$cell_type[grepl("Granulosum", Psoriasis04_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis04_emigrating_cell02 <- ddply(Psoriasis04_emigrating_cell,"cell_type",numcolwise(sum))
Psoriasis04_emigrating_cell02$emigrating_cell_proportion <- Psoriasis04_emigrating_cell02$Freq / sum(Psoriasis04_emigrating_cell$Freq) *100
Psoriasis04_emigrating_cell02
Psoriasis04_emigrating_cell02 <- Psoriasis04_emigrating_cell02[,c(1,4)]

Psoriasis04_enzyme_digestion$cell_type[grepl("T_cell", Psoriasis04_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis04_enzyme_digestion$cell_type[grepl("Treg", Psoriasis04_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis04_enzyme_digestion$cell_type[grepl("DC", Psoriasis04_enzyme_digestion$Var1) ] <- "Dendritic_cell"
Psoriasis04_enzyme_digestion$cell_type[grepl("Melanocyte", Psoriasis04_enzyme_digestion$Var1) ] <- "Melanocyte"
Psoriasis04_enzyme_digestion$cell_type[grepl("NK_cell", Psoriasis04_enzyme_digestion$Var1) ] <- "NK_cell"
Psoriasis04_enzyme_digestion$cell_type[grepl("Basale", Psoriasis04_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis04_enzyme_digestion$cell_type[grepl("Corneum", Psoriasis04_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis04_enzyme_digestion$cell_type[grepl("Spinosum", Psoriasis04_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis04_enzyme_digestion$cell_type[grepl("Granulosum", Psoriasis04_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis04_enzyme_digestion02 <- ddply(Psoriasis04_enzyme_digestion,"cell_type",numcolwise(sum))
Psoriasis04_enzyme_digestion02$enzyme_digestion_proportion <- Psoriasis04_enzyme_digestion02$Freq / sum(Psoriasis04_enzyme_digestion$Freq) *100
Psoriasis04_enzyme_digestion02
Psoriasis04_enzyme_digestion02 <- Psoriasis04_enzyme_digestion02[,c(1,4)]

Psoriasis04 <- data.frame("sample" = "Psoriasis04", "cell_type" = c("NK_cell","T_cell","Dendritic_cell","Melanocyte","Keratinocyte"))
Psoriasis04 <- merge(Psoriasis04, Psoriasis04_emigrating_cell02, by='cell_type',all.x=TRUE)
Psoriasis04 <- merge(Psoriasis04, Psoriasis04_enzyme_digestion02, by='cell_type',all.x=TRUE)
Psoriasis04[is.na(Psoriasis04)] <- 0


#Psoriasis06
Psoriasis06_emigrating_cell$cell_type[grepl("T_cell", Psoriasis06_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis06_emigrating_cell$cell_type[grepl("Treg", Psoriasis06_emigrating_cell$Var1) ] <- "T_cell"
Psoriasis06_emigrating_cell$cell_type[grepl("DC", Psoriasis06_emigrating_cell$Var1) ] <- "Dendritic_cell"
Psoriasis06_emigrating_cell$cell_type[grepl("Melanocyte", Psoriasis06_emigrating_cell$Var1) ] <- "Melanocyte"
Psoriasis06_emigrating_cell$cell_type[grepl("NK_cell", Psoriasis06_emigrating_cell$Var1) ] <- "NK_cell"
Psoriasis06_emigrating_cell$cell_type[grepl("Basale", Psoriasis06_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis06_emigrating_cell$cell_type[grepl("Corneum", Psoriasis06_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis06_emigrating_cell$cell_type[grepl("Spinosum", Psoriasis06_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis06_emigrating_cell$cell_type[grepl("Granulosum", Psoriasis06_emigrating_cell$Var1) ] <- "Keratinocyte"
Psoriasis06_emigrating_cell$cell_type[grepl("Macrophage", Psoriasis06_emigrating_cell$Var1) ] <- "Macrophage"
Psoriasis06_emigrating_cell02 <- ddply(Psoriasis06_emigrating_cell,"cell_type",numcolwise(sum))
Psoriasis06_emigrating_cell02$emigrating_cell_proportion <- Psoriasis06_emigrating_cell02$Freq / sum(Psoriasis06_emigrating_cell$Freq) *100
Psoriasis06_emigrating_cell02
Psoriasis06_emigrating_cell02 <- Psoriasis06_emigrating_cell02[,c(1,4)]

Psoriasis06_enzyme_digestion$cell_type[grepl("T_cell", Psoriasis06_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis06_enzyme_digestion$cell_type[grepl("Treg", Psoriasis06_enzyme_digestion$Var1) ] <- "T_cell"
Psoriasis06_enzyme_digestion$cell_type[grepl("DC", Psoriasis06_enzyme_digestion$Var1) ] <- "Dendritic_cell"
Psoriasis06_enzyme_digestion$cell_type[grepl("Melanocyte", Psoriasis06_enzyme_digestion$Var1) ] <- "Melanocyte"
Psoriasis06_enzyme_digestion$cell_type[grepl("NK_cell", Psoriasis06_enzyme_digestion$Var1) ] <- "NK_cell"
Psoriasis06_enzyme_digestion$cell_type[grepl("Basale", Psoriasis06_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis06_enzyme_digestion$cell_type[grepl("Corneum", Psoriasis06_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis06_enzyme_digestion$cell_type[grepl("Spinosum", Psoriasis06_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis06_enzyme_digestion$cell_type[grepl("Granulosum", Psoriasis06_enzyme_digestion$Var1) ] <- "Keratinocyte"
Psoriasis06_enzyme_digestion$cell_type[grepl("Macrophage", Psoriasis06_enzyme_digestion$Var1) ] <- "Macrophage"
Psoriasis06_enzyme_digestion02 <- ddply(Psoriasis06_enzyme_digestion,"cell_type",numcolwise(sum))
Psoriasis06_enzyme_digestion02$enzyme_digestion_proportion <- Psoriasis06_enzyme_digestion02$Freq / sum(Psoriasis06_enzyme_digestion$Freq) *100
Psoriasis06_enzyme_digestion02
Psoriasis06_enzyme_digestion02 <- Psoriasis06_enzyme_digestion02[,c(1,4)]

Psoriasis06 <- data.frame("sample" = "Psoriasis06", "cell_type" = c("NK_cell","T_cell","Dendritic_cell","Macrophage","Melanocyte","Keratinocyte"))
Psoriasis06 <- merge(Psoriasis06, Psoriasis06_emigrating_cell02, by='cell_type',all.x=TRUE)
Psoriasis06 <- merge(Psoriasis06, Psoriasis06_enzyme_digestion02, by='cell_type',all.x=TRUE)
Psoriasis06[is.na(Psoriasis06)] <- 0

#bind dataset
dataset_for_correlation <- rbind (Psoriasis02, Psoriasis03, Psoriasis04, Psoriasis06)

colnames(dataset_for_correlation)

####################################################################################################################
# Correlation plot ####
# Exclude Control03 sample (only kerationcytes)

pdf('Correlation_plot_emigrating_cell_vs_enzyme_digestion.pdf', width=9, height=9)
ggscatter(dataset_for_correlation, y = "emigrating_cell_proportion", x = "enzyme_digestion_proportion", 
          add = "reg.line", conf.int = TRUE,
          color = "cell_type" ,
          label = "cell_type" , repel = TRUE, size = 5,
          cor.coef = TRUE, cor.method = "pearson", use = "complete.obs",
          add.params = list(color = "black",
                            fill = "lightgray"),
          ylab = "Harvesting emigrating cells - Proportion of cells (%)", 
          xlab = "Enzyme digestion - Proportion of cells (%)")
dev.off()
