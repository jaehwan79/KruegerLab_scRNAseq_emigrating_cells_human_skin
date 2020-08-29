#install.packages("tibble")
library(tibble)
#install.packages("tidyverse")
library(tidyverse)
library(readxl)

rm(list = ls())

setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

############################################################################
##Proportion heatmap #### Supplementary Figure S5
#Flow_cytometry_data_proportion.csv is uploaded on https://github.com/jaehwan79/KruegerLab_scRNAseq_emigrating_cells_human_skin
df <- read.csv("~/Dropbox/SIngle_cell/Rcode_JKIM/github_upload_02/Flow_cytometry_data_proportion.csv", row.names=1)

df_for_heatmap <- merge(df_mature, df_semimature, by  ="Marker")
rownames(df_for_heatmap) <- df_for_heatmap[,1]
df_for_heatmap<- df_for_heatmap[,-1]

df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "BDCA-2", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "ILT2", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "DCR1", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD207", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD163", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "BDCA-2", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD172", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD11c", ]

dt2 <- df_for_heatmap %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

dt2$value <- as.numeric(as.character(dt2$value))

dt2$rowname <- as.factor(dt2$rowname)
dt2$rowname <- factor(dt2$rowname, levels = rev(c("CD205","CD40","BDCA-3", "ILT4","BDCA-1","CD11c","BDCA-2","CD163","CD172","CD207","DCR1","ILT2")))
dt2$rowname <- factor(dt2$rowname, levels = rev(c( "CD205","CD40", "BDCA-1", "BDCA-3","ILT4")))

ggplot(dt2, aes(x = colname, y = rowname, fill = value)) +
  geom_tile()+
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))+
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

ggsave("Flowcytometry_heatmap_proportion.pdf", width = 6, height = 2.5)
######################################################################################################################################################
##MFI heatmap #### Figure 5b
#Flow_cytometry_data_MFI.csv is uploaded on https://github.com/jaehwan79/KruegerLab_scRNAseq_emigrating_cells_human_skin
rm(list=ls())
df <- read.csv("~/Dropbox/SIngle_cell/Rcode_JKIM/github_upload_02/Flow_cytometry_data_MFI.csv", row.names=1)

df_mature <- df[df$Maturity == "Mature", ]
df_mature <- na.omit(df_mature)
df_mature <- df_mature[,-3]
colnames(df_mature)[c(1,2)] <- c("Mature_Control","Mature_Psoriasis")

df_semimature <- df[df$Maturity == "Semimature", ]
df_semimature <- na.omit(df_semimature)
df_semimature <- df_semimature[,-3]
colnames(df_semimature)[c(1,2)] <- c("Semimature_Control","Semimature_Psoriasis")

#
df_for_heatmap <- merge(df_mature, df_semimature, by  ="Marker")
rownames(df_for_heatmap) <- df_for_heatmap[,1]
df_for_heatmap<- df_for_heatmap[,-1]

df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "BDCA-2", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "ILT2", ]
#df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "DCR1", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD207", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD163", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "BDCA-2", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD172", ]
df_for_heatmap <- df_for_heatmap[!rownames(df_for_heatmap) %in% "CD11c", ]

dt2 <- df_for_heatmap %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

dt2$value <- as.numeric(as.character(dt2$value))
dt2$value[dt2$value>4000] <- 4000

#dt2$value <- scale(dt2$value)
#dt2$value[dt2$value>2] <- 1

dt2$rowname <- factor(dt2$rowname, levels = rev(c("HLA-DR", "CD205","CD40", "BDCA-1", "BDCA-3","DCR1","ILT4")))

ggplot(dt2, aes(x = colname, y = rowname, fill = value)) +
  geom_tile()+
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))+
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

ggsave("Flowcytometry_heatmap_MFI.pdf", width = 6, height = 2.5)

  
