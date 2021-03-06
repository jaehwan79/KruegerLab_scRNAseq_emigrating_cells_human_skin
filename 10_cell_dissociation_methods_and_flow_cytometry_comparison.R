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
#devtools::install_github("kassambara/ggpubr")
#install.packages("ggpubr")
library("ggpubr")
######################################################################################################
# Number of cells (psoriasis vs. control) in each cluster calculation   ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for subclustering analysis 
######################################################################################################
#Load the integrated data saved from previous code (03) ####
rm(list = ls())
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
load("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")

rm(list= ls()[!(ls() %in% c('Round2'))])

##########################################################################################
## Emigrating cell proportion ####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

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


## Add column for group
list$group[grepl("T_cell", list$Cluster)] <- "T_cell"
list$group[grepl("Macrophage", list$Cluster)] <- "DC"
list$group[grepl("DC", list$Cluster)] <- "DC"
list$group[grepl("Melanocyte", list$Cluster)] <- "Non_leucocyte"
list$group[grepl("NK", list$Cluster)] <- "Other_leucocyte"
list$group[grepl("Treg", list$Cluster)] <- "T_cell"

list$group[grepl("Corneum", list$Cluster)] <- "Non_leucocyte"
list$group[grepl("Granulosum", list$Cluster)] <- "Non_leucocyte"
list$group[grepl("Spinosum", list$Cluster)] <- "Non_leucocyte"
list$group[grepl("Basale", list$Cluster)] <- "Non_leucocyte"


## Aggregate
list <- aggregate(list$cell_number, list(group = list$group, sample = list$sample),sum)

## Add total number of cells per sample
Idents(object = Round2) <- ("number")
levels(Idents(object = Round2))

list2 <- as.data.frame(table(Round2@meta.data$number))
colnames(list2) <- c("sample","total_number_of_cells_in_sample")

list <- merge(list,list2, by='sample',all.x=TRUE)

list$emigrating_cell_proportion <- list$x / list$total_number_of_cells_in_sample * 100

list$category <- paste(list$sample,list$group, sep = "_")

list_emigrating_cell <- list[,c(6,5)]

list_emigrating_cell2 <- list_emigrating_cell
list_emigrating_cell2$category <- paste(list_emigrating_cell2$category,"emigrating_cell", sep = "_")
colnames(list_emigrating_cell2)[2] <- "cell_proportion"

rm(list)
rm(list2)
##########################################################################################
## Enzyme digestion cell proportion ####
# Enzyme digestion sample csv files are saved by code:09_cell_dissociation_methods_comparison.R
Enzyme_Psoriasis02 <- read.csv("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis02F.csv")
Enzyme_Psoriasis02$sample <- "Psoriasis02"

Enzyme_Psoriasis03 <- read.csv("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis03F.csv")
Enzyme_Psoriasis03$sample <- "Psoriasis03"

Enzyme_Psoriasis04 <- read.csv("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis04F.csv")
Enzyme_Psoriasis04$sample <- "Psoriasis04"

Enzyme_Psoriasis06 <- read.csv("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Label_transfer_Psoriasis06F.csv")
Enzyme_Psoriasis06$sample <- "Psoriasis06"


list_enzyme <- rbind (Enzyme_Psoriasis02,Enzyme_Psoriasis03,Enzyme_Psoriasis04,Enzyme_Psoriasis06)

list_enzyme$group[grepl("T_cell", list_enzyme$Var1)] <- "T_cell"
list_enzyme$group[grepl("DC", list_enzyme$Var1)] <- "DC"
list_enzyme$group[grepl("Melanocyte", list_enzyme$Var1)] <- "Non_leucocyte"
list_enzyme$group[grepl("Treg", list_enzyme$Var1)] <- "T_cell"
list_enzyme$group[grepl("NK", list_enzyme$Var1)] <- "Other_leucocyte"

list_enzyme$group[grepl("Corneum", list_enzyme$Var1)] <- "Non_leucocyte"
list_enzyme$group[grepl("Granulosum", list_enzyme$Var1)] <- "Non_leucocyte"
list_enzyme$group[grepl("Spinosum", list_enzyme$Var1)] <- "Non_leucocyte"
list_enzyme$group[grepl("Basale", list_enzyme$Var1)] <- "Non_leucocyte"


## Aggregate
list_enzyme <- aggregate(list_enzyme$Freq, list(group = list_enzyme$group, sample = list_enzyme$sample),sum)

## Add row with 0 value
list_enzyme[nrow(list_enzyme) + 1,] = c("Other_leucocyte",	"Psoriasis02",	0)

list_enzyme[nrow(list_enzyme) + 1,] = c("T_cell",	"Control05",	0)
list_enzyme[nrow(list_enzyme) + 1,] = c("DC",	"Control05",	0)
list_enzyme[nrow(list_enzyme) + 1,] = c("Non_leucocyte",	"Control05",	0)
list_enzyme[nrow(list_enzyme) + 1,] = c("Other_leucocyte",	"Control05",	0)

## Total number cells
list_enzyme$Total_number_of_cells[grepl("Control05", list_enzyme$sample)] <- 62
list_enzyme$Total_number_of_cells[grepl("Psoriasis02", list_enzyme$sample)] <- 556
list_enzyme$Total_number_of_cells[grepl("Psoriasis03", list_enzyme$sample)] <- 200
list_enzyme$Total_number_of_cells[grepl("Psoriasis04", list_enzyme$sample)] <- 1327
list_enzyme$Total_number_of_cells[grepl("Psoriasis06", list_enzyme$sample)] <- 842

list_enzyme$x <- as.numeric(list_enzyme$x)

list_enzyme$enzyme_digestion_proportion <- list_enzyme$x / list_enzyme$Total_number_of_cells*100

list_enzyme$category <- paste(list_enzyme$sample,list_enzyme$group, sep = "_")

list_enzyme <- list_enzyme[,c(6,5)]

list_enzyme2 <- list_enzyme
list_enzyme2$category <- paste(list_enzyme2$category,"enzyme_digestion", sep = "_")
colnames(list_enzyme2)[2] <- "cell_proportion"

rm(Enzyme_Psoriasis02)
rm(Enzyme_Psoriasis03)
rm(Enzyme_Psoriasis04)
rm(Enzyme_Psoriasis06)


##########################################################################################
## Flow cytometry cell proportion ####
#Flow_cytometry_event_count.csv is uploaded on https://github.com/jaehwan79/KruegerLab_scRNAseq_emigrating_cells_human_skin
df <- read.csv("~/Dropbox/SIngle_cell/Rcode_JKIM/github_upload_02/Flow_cytometry_event_count.csv", row.names=1)

df$DC <- df[,5]-df[,6]
df$Other_leucocyte <- df[,6]
df$T_cell <- df[,4]-df[,5]
df$Non_leucocyte <-   df[,3]-df[,4]
  
df[,c(7,8,9,10)] <- df[,c(7,8,9,10)] / df[,3] *100


df_dc <- df[,c("Category","DC")]
df_dc$category <-   paste(df_dc$Category,"DC", sep = "_")
df_dc <- df_dc[,-1]
colnames(df_dc)[1] <- "Flowcytometry_proportion"

colnames(df)

df_Other_leucocyte <- df[,c("Category","Other_leucocyte")]
df_Other_leucocyte$category <-   paste(df_Other_leucocyte$Category,"Other_leucocyte", sep = "_")
df_Other_leucocyte <- df_Other_leucocyte[,-1]
colnames(df_Other_leucocyte)[1] <- "Flowcytometry_proportion"
  
df_T_cell <- df[,c("Category","T_cell")]
df_T_cell$category <-   paste(df_T_cell$Category,"T_cell", sep = "_")
df_T_cell <- df_T_cell[,-1]
colnames(df_T_cell)[1] <- "Flowcytometry_proportion"
  
df_Non_leucocyte <- df[,c("Category","Non_leucocyte")]
df_Non_leucocyte$category <-   paste(df_Non_leucocyte$Category,"Non_leucocyte", sep = "_")
df_Non_leucocyte <- df_Non_leucocyte[,-1]
colnames(df_Non_leucocyte)[1] <- "Flowcytometry_proportion"

list_flowcytometry <- rbind(df_T_cell,df_dc,df_Other_leucocyte,df_Non_leucocyte   )

list_flowcytometry2 <- list_flowcytometry
list_flowcytometry2$category <- paste(list_flowcytometry2$category,"flow_cytometry", sep = "_")
colnames(list_flowcytometry2)[1] <- "cell_proportion"
list_flowcytometry2 <- list_flowcytometry2[,c(2,1)]

rm(df)
rm(df_dc)
rm(df_Non_leucocyte)
rm(df_T_cell)
rm(df_Other_leucocyte)




##########################################################################################
## merge dataset ####
list <- merge(list_emigrating_cell,list_flowcytometry, by='category',all.x=TRUE)
list <- merge(list,list_enzyme, by='category',all.x=TRUE)

## Add column for sample
list$sample[grepl("Control01", list$category)] <- "Control01"
list$sample[grepl("Control02", list$category)] <- "Control02"
list$sample[grepl("Control03", list$category)] <- "Control03"
list$sample[grepl("Control04", list$category)] <- "Control04"
list$sample[grepl("Control05", list$category)] <- "Control05"

list$sample[grepl("Psoriasis01", list$category)] <- "Psoriasis01"
list$sample[grepl("Psoriasis02", list$category)] <- "Psoriasis02"
list$sample[grepl("Psoriasis03", list$category)] <- "Psoriasis03"
list$sample[grepl("Psoriasis04", list$category)] <- "Psoriasis04"
list$sample[grepl("Psoriasis05", list$category)] <- "Psoriasis05"
list$sample[grepl("Psoriasis06", list$category)] <- "Psoriasis06"
list$sample[grepl("Psoriasis07", list$category)] <- "Psoriasis07"
list$sample[grepl("Psoriasis08", list$category)] <- "Psoriasis08"
list$sample[grepl("Psoriasis09", list$category)] <- "Psoriasis09"
list$sample[grepl("Psoriasis10", list$category)] <- "Psoriasis10"
list$sample[grepl("Psoriasis11", list$category)] <- "Psoriasis11"
list$sample[grepl("Psoriasis12", list$category)] <- "Psoriasis12"
list$sample[grepl("Psoriasis13", list$category)] <- "Psoriasis13"

## Add column for category
list$group[grepl("Psoriasis", list$category)] <- "Psoriasis"
list$group[grepl("Control", list$category)] <- "Control"

## add cell type
list$celltype[grepl("T_cell", list$category)] <- "T_cell"
list$celltype[grepl("DC", list$category)] <- "DC"
list$celltype[grepl("Other_leucocyte", list$category)] <- "Other_leucocyte"
list$celltype[grepl("Non_leucocyte", list$category)] <- "Non_leucocyte"

## List for Bar graph
list2 <- rbind(list_emigrating_cell2, list_flowcytometry2, list_enzyme2)

list2$group[grepl("Psoriasis", list2$category)] <- "Psoriasis"
list2$group[grepl("Control", list2$category)] <- "Control"

# add cell type
list2$celltype[grepl("T_cell", list2$category)] <- "T_cell"
list2$celltype[grepl("DC", list2$category)] <- "DC"
list2$celltype[grepl("Other_leucocyte", list2$category)] <- "Other_leucocyte"
list2$celltype[grepl("Non_leucocyte", list2$category)] <- "Non_leucocyte"

# add experiment
list2$experiment[grepl("emigrating_cell", list2$category)] <- "emigrating_cell"
list2$experiment[grepl("flow_cytometry", list2$category)] <- "flow_cytometry"
list2$experiment[grepl("enzyme_digestion", list2$category)] <- "enzyme_digestion"

colnames(list_enzyme2)

##########################################################################################
## Bar graph ####

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


# Convert dose to a factor variable (all)
df2 <- summarySE(list2, measurevar = "cell_proportion", groupvars = c( "celltype", "experiment"))

df2$celltype <- as.factor(df2$celltype)
df2$celltype <- factor(df2$celltype, levels = c("T_cell", "DC", "Other_leucocyte", "Non_leucocyte" ))

df2$experiment <- as.factor(df2$experiment)
df2$experiment <- factor(df2$experiment, levels = c("flow_cytometry", "emigrating_cell", "enzyme_digestion"))

## Figure 2b
ggplot(data = df2, aes(x=celltype, y=cell_proportion, fill = experiment)) +
  geom_bar(stat = 'identity', color='black', position=position_dodge()) +
  theme_light() +
  geom_errorbar(aes(ymin=cell_proportion, ymax=cell_proportion+se), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Composition of human skin cells")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Cell type") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("", values = c("#00BA38", "indianred1","royalblue1")) +
  theme(axis.title.x = element_blank())+
  theme(text = element_text(size=16))
ggsave("Cell_proportion_experiment_comparison_all_round2.pdf", width = 8, height = 6)
