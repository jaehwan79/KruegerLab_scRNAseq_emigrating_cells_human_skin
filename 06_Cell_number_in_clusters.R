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

##########################################################################################
## Number of cells (psoriasis vs. control) in each sample - Round2 ####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")

Idents(object = Round2) <- ("ClusterNames_0.8")
Round2$ClusterNames_0.8_number <- paste(Idents(Round2),Round2$number, sep = "_")
Idents(Round2) <- "ClusterNames_0.8_number"
Idents(Round2) <- fct_relevel(Idents(Round2), sort)
levels(Idents(object = Round2))

## Make a list
list <- as.data.frame(table(Round2@meta.data$ClusterNames_0.8_number))
colnames(list) <- c("Cluster_in_sample","cell_number")

list$Cluster_in_sample <- as.character(list$Cluster_in_sample)

## Add clusters with 0
list[nrow(list) + 1,1] = c("CD161_T_cell_Control03")
list[nrow(list) + 1,1] = c("CD161_T_cell_Control05")
list[nrow(list) + 1,1] = c("CD161_T_cell_Psoriasis03")
list[nrow(list) + 1,1] = c("CD161_T_cell_Psoriasis05")
list[nrow(list) + 1,1] = c("CD161_T_cell_Psoriasis06")
list[nrow(list) + 1,1] = c("CD161_T_cell_Psoriasis07")
list[nrow(list) + 1,1] = c("CD161_T_cell_Psoriasis08")
list[nrow(list) + 1,1] = c("CD161_T_cell_Psoriasis13")

list[nrow(list) + 1,1] = c("CD4_T_cell_Control03")
list[nrow(list) + 1,1] = c("CD8_T_cell_Control03")
list[nrow(list) + 1,1] = c("CD8_T_cell_Psoriasis06")

list[nrow(list) + 1,1] = c("Macrophage_Control01")
list[nrow(list) + 1,1] = c("Macrophage_Control02")
list[nrow(list) + 1,1] = c("Macrophage_Control03")
list[nrow(list) + 1,1] = c("Macrophage_Control04")
list[nrow(list) + 1,1] = c("Macrophage_Control05")

list[nrow(list) + 1,1] = c("Macrophage_Psoriasis01")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis02")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis03")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis04")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis07")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis09")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis10")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis11")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis12")
list[nrow(list) + 1,1] = c("Macrophage_Psoriasis13")

list[nrow(list) + 1,1] = c("Mature_DC_Control03")

list[nrow(list) + 1,1] = c("Melanocyte_Control05")

list[nrow(list) + 1,1] = c("NK_cell_Control03")

list[nrow(list) + 1,1] = c("NK_cell_Psoriasis05")

list[nrow(list) + 1,1] = c("S.Granulosum_Control02")
list[nrow(list) + 1,1] = c("S.Granulosum_Control03")
list[nrow(list) + 1,1] = c("S.Granulosum_Control05")

list[nrow(list) + 1,1] = c("S.Granulosum_Psoriasis02")
list[nrow(list) + 1,1] = c("S.Granulosum_Psoriasis05")

list[nrow(list) + 1,1] = c("S.Spinosum_Psoriasis04")

list[nrow(list) + 1,1] = c("Semimature_DC_Control03")

list[nrow(list) + 1,1] = c("Treg_Control03")

list[nrow(list) + 1,1] = c("Treg_Psoriasis05")
list[nrow(list) + 1,1] = c("Treg_Psoriasis06")

list <- list[order(list$Cluster_in_sample),]

list[is.na(list)] <- 0


## Add column for sample
list$sample[grepl("Control01", list$Cluster_in_sample)] <- "Control01"
list$sample[grepl("Control02", list$Cluster_in_sample)] <- "Control02"
list$sample[grepl("Control03", list$Cluster_in_sample)] <- "Control03"
list$sample[grepl("Control04", list$Cluster_in_sample)] <- "Control04"
list$sample[grepl("Control05", list$Cluster_in_sample)] <- "Control05"

list$sample[grepl("Psoriasis01", list$Cluster_in_sample)] <- "Psoriasis01"
list$sample[grepl("Psoriasis02", list$Cluster_in_sample)] <- "Psoriasis02"
list$sample[grepl("Psoriasis03", list$Cluster_in_sample)] <- "Psoriasis03"
list$sample[grepl("Psoriasis04", list$Cluster_in_sample)] <- "Psoriasis04"
list$sample[grepl("Psoriasis05", list$Cluster_in_sample)] <- "Psoriasis05"
list$sample[grepl("Psoriasis06", list$Cluster_in_sample)] <- "Psoriasis06"
list$sample[grepl("Psoriasis07", list$Cluster_in_sample)] <- "Psoriasis07"
list$sample[grepl("Psoriasis08", list$Cluster_in_sample)] <- "Psoriasis08"
list$sample[grepl("Psoriasis09", list$Cluster_in_sample)] <- "Psoriasis09"
list$sample[grepl("Psoriasis10", list$Cluster_in_sample)] <- "Psoriasis10"
list$sample[grepl("Psoriasis11", list$Cluster_in_sample)] <- "Psoriasis11"
list$sample[grepl("Psoriasis12", list$Cluster_in_sample)] <- "Psoriasis12"
list$sample[grepl("Psoriasis13", list$Cluster_in_sample)] <- "Psoriasis13"

## Add column for category
list$category[grepl("Psoriasis", list$Cluster_in_sample)] <- "Psoriasis"
list$category[grepl("Control", list$Cluster_in_sample)] <- "Control"

## Add column for cluster
list$Cluster <- list$Cluster_in_sample
list$Cluster= gsub("_Control", "", list$Cluster)
list$Cluster= gsub("_Psoriasis", "", list$Cluster)
list$Cluster= gsub("01", "", list$Cluster)
list$Cluster= gsub("02", "", list$Cluster)
list$Cluster= gsub("03", "", list$Cluster)
list$Cluster= gsub("04", "", list$Cluster)
list$Cluster= gsub("05", "", list$Cluster)
list$Cluster= gsub("06", "", list$Cluster)
list$Cluster= gsub("07", "", list$Cluster)
list$Cluster= gsub("08", "", list$Cluster)
list$Cluster= gsub("09", "", list$Cluster)
list$Cluster= gsub("10", "", list$Cluster)
list$Cluster= gsub("11", "", list$Cluster)
list$Cluster= gsub("12", "", list$Cluster)
list$Cluster= gsub("13", "", list$Cluster)

## Add total number of cells per sample
Idents(object = Round2) <- ("number")
levels(Idents(object = Round2))

list2 <- as.data.frame(table(Round2@meta.data$number))
colnames(list2) <- c("sample","total_number_of_cells_in_sample")

list <- merge(list,list2, by='sample',all.x=TRUE)
list$cell_proportion_in_cluster <- list$cell_number/list$total_number_of_cells_in_sample*100

## Add column for NK cell, T-cells, Dendritic cells, Macrophage, Melanocyte, Keratinocyte
list$Group[grepl("NK_cell", list$Cluster)] <- "NK_cell"
list$Group[grepl("T_cell", list$Cluster)] <- "T_cell"
list$Group[grepl("Treg", list$Cluster)] <- "T_cell"
list$Group[grepl("DC", list$Cluster)] <- "DC"
list$Group[grepl("Melanocyte", list$Cluster)] <- "Melanocyte"
list$Group[grepl("Macrophage", list$Cluster)] <- "Macrophage"

list$Group[grepl("Corneum", list$Cluster)] <- "Keratinocyte"
list$Group[grepl("Granulosum", list$Cluster)] <- "Keratinocyte"
list$Group[grepl("Spinosum", list$Cluster)] <- "Keratinocyte"
list$Group[grepl("Basale", list$Cluster)] <- "Keratinocyte"


list$Group <- factor(list$Group)
list$Group <- factor(list$Group, levels = c("NK_cell", "T_cell","DC","Macrophage", "Melanocyte", "Keratinocyte"))

list$Total_number_of_cells[grepl("Psoriasis", list$category)] <- 18339
list$Total_number_of_cells[grepl("Control", list$category)] <- 4881

list$cell_proportion_in_total_cells <- list$cell_number/list$Total_number_of_cells*100

## Save table
list$Cluster <- factor(list$Cluster)
list$Cluster <- factor(list$Cluster, 
                       levels = (c( "NK_cell"  ,     "CD161_T_cell",  "CD8_T_cell"   , "CD4_T_cell" ,"Treg" ,    
                                    "Mature_DC" ,  "Semimature_DC", "Macrophage"   , "Melanocyte"   ,
                                    "S.Corneum"  , "S.Granulosum" , "S.Spinosum" ,   "S.Basale"        
                       )))



list <- list[,c(2,5,1,4,8,3,7,9,10,6)]

## Mean and SD calculation
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
    length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
    # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#########################################################################################################################
# Convert dose to a factor variable
df2 <- summarySE(list, measurevar = "cell_number", groupvars = c("Cluster", "category"))
df2$cell_number[which(df2$category=="Psoriasis" & df2$Cluster == "CD161_T_cell")]/df2$cell_number[which(df2$category=="Control" & df2$Cluster == "CD161_T_cell")]
df3 <- summarySE(list, measurevar = "cell_proportion_in_cluster", groupvars = c("Cluster", "category"))

list_for_df <- aggregate(list$cell_number, list(sample = list$sample, Group = list$Group, category=list$category),sum)
colnames(list_for_df)[4] <- "cell_number"
df4 <- summarySE(list_for_df, measurevar = "cell_number", groupvars = c("Group", "category"))

list_for_df <- aggregate(list$cell_proportion_in_cluster, list(sample = list$sample, Group = list$Group, category=list$category),sum)
colnames(list_for_df)[4] <- "cell_proportion_in_cluster"
df5 <- summarySE(list_for_df, measurevar = "cell_proportion_in_cluster", groupvars = c("Group", "category"))

## Cell number or cell proportion of each cluster bar graph 
#(Supplementary FIgure S1g)
ggplot(data = df2, aes(x=Cluster, y=cell_number, fill = category)) +
  geom_bar(stat = 'identity', color='black', position=position_dodge()) +
  theme_light() +
  geom_errorbar(aes(ymin=cell_number, ymax=cell_number+se), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Average cell number per sample")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Clusters") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Cellnumbers_per_sample_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1h)
ggplot(data = df3, aes(x=Cluster, y=cell_proportion_in_cluster, fill = category)) +
  geom_bar(stat = 'identity', color='black', position=position_dodge()) +
  theme_light() +
  geom_errorbar(aes(ymin=cell_proportion_in_cluster, ymax=cell_proportion_in_cluster+se), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Average cell proportion per sample")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Cell_proportion_per_sample_in_cluster_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1c)
ggplot(data = df4, aes(x=Group, y=cell_number, fill = category)) +
  geom_bar(stat = 'identity', color='black', position=position_dodge()) +
  theme_light() +
  geom_errorbar(aes(ymin=cell_number, ymax=cell_number+se), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Average cell number per sample")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Groups") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Group_Cellnumbers_per_sample_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1d)
ggplot(data = df5, aes(x=Group, y=cell_proportion_in_cluster, fill = category)) +
  geom_bar(stat = 'identity', color='black', position=position_dodge()) +
  theme_light() +
  geom_errorbar(aes(ymin=cell_proportion_in_cluster, ymax=cell_proportion_in_cluster+se), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Average cell proportion per sample")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Groups") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Group_Cell_proportion_per_sample_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1a)
list_for_ggplot <- aggregate(list$cell_number, list(Group = list$Group, category = list$category),sum)
ggplot(data = list_for_ggplot, aes(x=Group, y=x, fill = category)) +
  geom_bar(stat = 'identity', color="black", position=position_dodge()) +
  theme_light() +
  ggtitle("Cell number")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Total_cell_count_Group_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1b)
list_for_ggplot <- aggregate(list$cell_proportion_in_total_cells, list(Group = list$Group, category = list$category),sum)
ggplot(data = list_for_ggplot, aes(x=Group, y=x, fill = category)) +
  geom_bar(stat = 'identity',color="black",position=position_dodge())+
  theme_light() +
  ggtitle("Cell proportion")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Total_cell_proportion_Group_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1e)
list_for_ggplot <- aggregate(list$cell_number, list(Cluster = list$Cluster, category = list$category),sum)
ggplot(data = list_for_ggplot, aes(x=Cluster, y=x, fill = category)) +
  geom_bar(stat = 'identity', color="black", position=position_dodge()) +
  theme_light() +
  ggtitle("Cell number")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Total_cell_count_Cluster_Round2.pdf", width = 8, height = 4)

#(Supplementary FIgure S1f)
list_for_ggplot <- aggregate(list$cell_proportion_in_total_cells, list(Cluster = list$Cluster, category = list$category),sum)
ggplot(data = list_for_ggplot, aes(x=Cluster, y=x, fill = category)) +
  geom_bar(stat = 'identity',color="black",position=position_dodge())+
  theme_light() +
  ggtitle("Cell proportion")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("royalblue1", "indianred1")) +
  theme(axis.title.x = element_blank())
ggsave("Total_cell_proportion_Cluster_Round2.pdf", width = 8, height = 4)