#load library ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for the downstream single-cell RNA sequencing data analy
library(Seurat)
library(ggplot2)
library(forcats)
library(cowplot)
library(dplyr)
library(EnhancedVolcano)
######################################################################################################
#Start ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

#Upload single-cell dataset with emigrating cell method
load("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round4_integrated_analyzed.Rda")
######################################################################################################
#Function to filter out cells with >25% mitochondrial genes 
#and ubiquitously expressed ribosomal protein-coding (RPS and RPL) and MALAT noncoding RNA, miRNA, and snoRNA genes 
remove.genes<-function(sample, BCtoremove){
  sample <- sample[, !colnames(sample) %in% BCtoremove]
  ubiq <- c(rownames(sample[grep("^MT-", rownames(sample)), ]),
            rownames(sample[grep("^RPS", rownames(sample)), ]),
            rownames(sample[grep("^RPL", rownames(sample)), ]),
            rownames(sample[grep("^RBP", rownames(sample)), ]),
            rownames(sample[grep("^MIR", rownames(sample)), ]),
            rownames(sample[grep("^SNOR", rownames(sample)), ]),
            rownames(sample[grep("^MTRNR", rownames(sample)), ]),
            rownames(sample[grep("^MALAT1", rownames(sample)), ]))
  sample <- sample[!rownames(sample) %in% ubiq, ]
  return(sample)
}

######################################################################################################
## Load data Control02Fresh_skin ####
Control02FreshEnzyme <- Read10X(data.dir = "/Volumes/JKIM_20TB/External_data_backup/EC-YE-5298_processed/8/filtered_gene_bc_matrices/GRCh38")
Control02FreshEnzyme <- CreateSeuratObject(counts = Control02FreshEnzyme, project = "Control02FreshEnzyme", min.cells = 3, min.features = 100)
Control02FreshEnzyme[["percent.mt"]] <- PercentageFeatureSet(object = Control02FreshEnzyme, pattern = "^MT-")
Control02FreshEnzyme@meta.data$stim <- "Fresh_tissue_Enzyme_digestion"
Control02FreshEnzyme@meta.data$number <- "Control02"

Control02FreshEnzyme_morethan25percentMT <- subset(Control02FreshEnzyme, subset = percent.mt > 25)
Control02FreshEnzyme_morethan25percentMTbc <- colnames(Control02FreshEnzyme_morethan25percentMT@assays$RNA@data)
rm(Control02FreshEnzyme)
rm(Control02FreshEnzyme_morethan25percentMT)

Control02FreshEnzyme <- Read10X(data.dir = "/Volumes/JKIM_20TB/External_data_backup/EC-YE-5298_processed/8/filtered_gene_bc_matrices/GRCh38")
Control02FreshEnzyme<-remove.genes(Control02FreshEnzyme, Control02FreshEnzyme_morethan25percentMTbc)

Control02FreshEnzyme <- CreateSeuratObject(counts = Control02FreshEnzyme, project = "Control02FreshEnzyme", min.cells = 3, min.features = 100)
Control02FreshEnzyme[["percent.mt"]] <- PercentageFeatureSet(object = Control02FreshEnzyme, pattern = "^MT-")
Control02FreshEnzyme <- subset(Control02FreshEnzyme, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control02FreshEnzyme@meta.data$stim <- "Fresh_tissue_Enzyme_digestion"
Control02FreshEnzyme@meta.data$number <- "Control02Fresh_tissue_Enzyme_digestion"
rm(Control02FreshEnzyme_morethan25percentMTbc)

Control02FreshEnzyme <- NormalizeData(Control02FreshEnzyme)
Control02FreshEnzyme <- FindVariableFeatures(Control02FreshEnzyme, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control02FreshEnzyme)
Control02FreshEnzyme <- ScaleData(Control02FreshEnzyme, features = all.genes)
rm(all.genes)

######################################################################################################
### Data inegration
# To harmonize merged groups into a single dataset without batch effects, 
# correspondences between cells in three merged datasets were identified by the FIndIntegrationAnchors function, 
# and used for data integration with the IntegratedData function as detailed by Butler et al. 
immune.anchors <- FindIntegrationAnchors(object.list = list(Round4, Control02FreshEnzyme), dims = 1:30)
Round3 <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

######################################################################################################
### Save
save(Round3, file = "Round3_combined.Rda")

######################################################################################################
# Run the standard work flow for visualization and clustering ####
Round3@assays
DefaultAssay(Round3) <- "integrated"
Round3 <- ScaleData(Round3, verbose = FALSE)
Round3 <- RunPCA(Round3, npcs = 20, verbose = FALSE)
######################################################################################################
# t-SNE and Clustering ####
# Twenty principal components (PCs) were selected for Uniform Manifold Approximation and Projection (UMAP)
# for Dimension Reduction. With a resolution of 0.8, cells were clustered by the FindClusters function.
Round3 <- FindNeighbors(Round3, reduction = "pca", dims = 1:20)
Round3 <- FindClusters(Round3, resolution = 0.8)
Round3 <- RunUMAP(object = Round3, reductiosdfn = "pca", dims = 1:20)

#####################################################################################################
## Cluster level order change #2 ####
colnames(Round3@meta.data)
Idents(object = Round3) <- ("stim")
levels(Idents(Round3)) 
Idents(object = Round3) <- factor(Idents(object = Round3), levels = c("Psoriasis"        ,             "Control"         ,              "Fresh_tissue_Enzyme_digestion"))      
Round3[["stim"]] <- Idents(object = Round3)

#####################################################################################################
## Visualization by clusters without labels####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

dir.create("./UMAP_clusters")
setwd("./UMAP_clusters")
getwd()

##cluster count 
Idents(object = Round3) <- ("integrated_snn_res.0.8")
colnames(Round3@meta.data)
cell.num <- table(Round3@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round3_cluster_count_integrated_snn_res.0.8.pdf', width=12,height=8)
DimPlot(object = Round3,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

pdf('UMAP_Round3_cluster_count_integrated_snn_res.0.8_noLabel.pdf', width=12,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##DimPlot
p1 <- DimPlot(Round3, reduction = "umap", group.by = "stim")+
  scale_color_manual(values=c('deeppink','green','deep sky blue'))

p2 <- DimPlot(Round3, reduction = "umap", label = FALSE)

pdf('UMAP_DimPlot_integrated_snn_res.0.8.pdf', width=16,height=7)
plot_grid(p1, p2)
dev.off()

##cluster count 
colnames(Round3@meta.data)
cell.num <- table(Round3@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cluster_count_integrated_snn_res.0.8.pdf', width=10,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###
pdf('UMAP_DimPlot_split_integrated_snn_res.0.8.pdf', width=20,height=7)
DimPlot(Round3, reduction = "umap", split.by = "stim")
dev.off()


##Cell number
colnames(Round3@meta.data)
Idents(object = Round3) <- ("number")

cell.num <- table(Round3@meta.data$number)
ClusterLabels = paste(names(cell.num), paste0("(n=", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cellcount_integrated_snn_res.0.8.pdf', width=14,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

##
colnames(Round3@meta.data)
Idents(object = Round3) 

Idents(object = Round3) <- ("stim")
cell.numb <- table(Round3@meta.data$stim)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round3_category_integrated_snn_res.0.8.pdf',width=12,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

#####################################################################################################
## (SKIP) DEG without labeling ####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

Round3@assays
colnames(Round3@meta.data)
DefaultAssay(Round3)
Idents(object = Round3) <- ("integrated_snn_res.0.8")
levels(Idents(object = Round3))

DefaultAssay(Round3) <- "RNA"
Round3 <- ScaleData(Round3, verbose = FALSE)

Combined.markers <- FindAllMarkers(object = Round3, only.pos = TRUE, min.pct = 0.25, 
                                   logfc.threshold = 0.25)

top1000 <- Combined.markers %>% group_by(cluster) 
write.csv(top1000, file = "top1000.csv")

#####################################################################################################
# Dot plot to figure out clusters ####
Idents(object = Round3) <- ("orig.ident")
DefaultAssay(Round3) <- "RNA"

features.plot <- c(  "DCD", "DCN","CD34","LYVE1","ACKR1","CCL21","MGP","KRT5","KRT14", "KRT1","KRT10","FABP5",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "CD163",
                     "CD14","LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CD40","CIITA","LY75","LAMP3",
                     "CTLA4","FOXP3","IL2RA","TIGIT", 
                     "CD8B","CD8A","GZMK","GZMH",
                     "TRBC1","TRAC",
                     "CD3D",
                     "GNLY","NKG7" ,"KLRB1")


features.plot <- rev(features.plot)

Idents(object = Round3) <- fct_rev(Idents(object = Round3))

pdf('Round3_Dot_plot_initial.pdf', width=18,height=12)
DotPlot(object = Round3, features = features.plot) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=18))
dev.off()


########################################################################
# Dot plot with labels ####
Idents(object = Round3) <- ("integrated_snn_res.0.8")
DefaultAssay(Round3) <- "RNA"

features.plot <- c(  "DCN","KRT5","KRT14", "KRT1","KRT10","FABP5",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "CD163",
                     "CD14","LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CD40","CIITA","LY75","LAMP3",
                     "TRBC1","TRAC",
                     "CD3D")
features.plot <- rev(features.plot)

new.cluster.ids <- c(
  "T_cell","S.Corneum","S.Corneum","S.Corneum","T_cell","Mature_DC","S.Spinosum","S.Granulosum","S.Corneum","S.Basale","Semimature_DC","S.Spinosum","S.Basale","Mature_DC","T_cell","S.Corneum","Melanocyte","S.Corneum","T_cell","Doublet","S.Basale","Mature_DC","Mature_DC","Mature_DC","Mature_DC","Doublet","Doublet","Macrophage","Fibroblast"
)	

names(new.cluster.ids) <- levels(Round3)
Round3 <- RenameIdents(Round3, new.cluster.ids)

Round3 <- subset(Round3, idents = "Doublet" ,invert=TRUE)

levels(Round3)

Idents(object = Round3) <- factor(Idents(object = Round3), 
                                  levels = (c(
                                    "T_cell"  ,
                                    "Mature_DC" ,
                                    "Semimature_DC" ,
                                    "Macrophage"  , 
                                    "Melanocyte" ,
                                    "S.Corneum"  ,
                                    "S.Granulosum",
                                     "S.Spinosum"  ,
                                    "S.Basale"  ,
                                     "Fibroblast"  ,
                                    "Doublet"    
                                  )))

Round3[["ClusterNames_0.8"]] <- Idents(object = Round3)
Idents(object = Round3) <- fct_rev(Idents(object = Round3))


pdf('plot_Round3.pdf', width=12,height=4.5)
DotPlot(object = Round3, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()

pdf('plot_Round3_02.pdf', width=12,height=4.5)
DotPlot(object = Round3, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()


Idents(object = Round3) <- fct_rev(Idents(object = Round3))


#######################################################################
# Visualization with label ####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

dir.create("./UMAP_Label")
setwd("./UMAP_Label")

colnames(Round3@meta.data)
DefaultAssay(Round3)
DefaultAssay(Round3) <- "RNA"

Idents(object = Round3) <- ("stim")
levels(Round3)
Idents(object = Round3) <- factor(Idents(object = Round3), levels = c("Psoriasis"   ,                  "Control"      ,                 "Fresh_tissue_Enzyme_digestion"))       
Idents(object = Round3) <- ("ClusterNames_0.8")

#color_code <- c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','indianred4','deepskyblue','royalblue1','navy','deeppink3','deepskyblue4','hotpink4','snow4','snow3')

##cluster count 
colnames(Round3@meta.data)
Idents(object = Round3) <- ("ClusterNames_0.8" )

cell.num <- table(Round3@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round3_cluster_count_integrated_snn_res.0.8.pdf', width=10,height=8)
DimPlot(object = Round3,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

pdf('UMAP_Round3_cluster_count_integrated_snn_res.0.8_noLabel.pdf', width=10,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##DimPlot
p1 <- DimPlot(Round3, reduction = "umap", group.by = "stim")+
  scale_color_manual(values=c('deeppink','green','deep sky blue'))

p2 <- DimPlot(Round3, reduction = "umap", label = FALSE)

pdf('UMAP_DimPlot_integrated_snn_res.0.8.pdf', width=16,height=7)
plot_grid(p1, p2)
dev.off()

##cluster count 
colnames(Round3@meta.data)
cell.num <- table(Round3@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cluster_count_integrated_snn_res.0.8.pdf', width=10,height=8)
DimPlot(object = Round3,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###
pdf('UMAP_DimPlot_split_integrated_snn_res.0.8.pdf', width=20,height=7)
DimPlot(Round3, reduction = "umap", split.by = "stim")
dev.off()


##Cell number
colnames(Round3@meta.data)
Idents(object = Round3) <- ("number")

cell.num <- table(Round3@meta.data$number)
ClusterLabels = paste(names(cell.num), paste0("(n=", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cellcount_integrated_snn_res.0.8.pdf', width=14,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

##
colnames(Round3@meta.data)
Idents(object = Round3) 

Idents(object = Round3) <- ("stim")
cell.numb <- table(Round3@meta.data$stim)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round3_category_integrated_snn_res.0.8.pdf',width=12,height=8)
DimPlot(object = Round3,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


setwd("../")

#####################################################################################################
# Save ####
save(Round3, file = "Round3_integrated_analyzed.Rda")

#######################################################################
# Round4 for control_emigrating cell vs. control_fresh_tissue_enzyme digestion 
colnames(Round3@meta.data)
Idents(object = Round3) <- ("stim")
#Round4 <- subset(Round3, idents = "Psoriasis" ,invert=TRUE)

Round4 <-Round3
colnames(Round4@meta.data)

Idents(object = Round4) <- ("ClusterNames_0.8")
levels(Round4)
Round4 <- subset(Round4, idents = "Macrophage" ,invert=TRUE)

new.cluster.ids <- c(
  "T_cell"  ,      "Dendritic_cell",     "Dendritic_cell",  "Melanocyte" ,   "Keratinocyte"   ,  "Keratinocyte" , "Keratinocyte"  ,  "Keratinocyte"     , "Fibroblast" 
)	

names(new.cluster.ids) <- levels(Round4)
Round4 <- RenameIdents(Round4, new.cluster.ids)

levels(Round4)
Round4[["ClusterNames_0.8"]] <- Idents(object = Round4)

########################################################################
# Round4 Dot plot with labels ####
Idents(object = Round4) <- ("ClusterNames_0.8")
DefaultAssay(Round4) <- "RNA"
levels(Round4)
features.plot <- c(  "DCN","KRT5","KRT14", "KRT1","KRT10",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CIITA","LY75","LAMP3",
                     "TRBC1","TRAC",
                     "CD3D")
features.plot <- rev(features.plot)

Idents(object = Round4) <- fct_rev(Idents(object = Round4))

pdf('plot_Round4.pdf', width=10,height=3)
DotPlot(object = Round4, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()

pdf('plot_Round4_02.pdf', width=10,height=3)
DotPlot(object = Round4, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()


Idents(object = Round4) <- fct_rev(Idents(object = Round4))




########################################################################
# Round4 Dot plot with labels - 02 ####
Idents(object = Round4) <- ("stim")
Round5 <- Round4 
levels(Round5)
Round5 <- subset(Round5, idents = "Psoriasis",invert=TRUE)
levels(Round5)

features.plot <- c(  "DCN","KRT5","KRT14", "KRT1","KRT10",
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "MLANA", "TYRP1", "DCT",
                     "LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CIITA","LY75","LAMP3",
                     "TRBC1","TRAC",
                     "CD3D")
features.plot <- rev(features.plot)

Idents(object = Round5) <- fct_rev(Idents(object = Round5))

pdf('plot_Round5_stim.pdf', width=11,height=1.9)
DotPlot(object = Round5, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()

pdf('plot_Round5_stim_02.pdf', width=11,height=1.9)
DotPlot(object = Round5, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()


Idents(object = Round4) <- fct_rev(Idents(object = Round4))




#######################################################################
# Visualization with label ####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")


Idents(object = Round4) <- ("ClusterNames_0.8")
levels(Idents(object = Round4))


dir.create("./UMAP_Label_Round4")
setwd("./UMAP_Label_Round4")

colnames(Round4@meta.data)
DefaultAssay(Round4)
DefaultAssay(Round4) <- "RNA"


##cluster count 
colnames(Round4@meta.data)
Idents(object = Round4) <- ("ClusterNames_0.8" )

cell.num <- table(Round4@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round4_cluster_count_integrated_snn_res.0.8.pdf', width=10,height=8)
DimPlot(object = Round4,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

pdf('UMAP_Round4_cluster_count_integrated_snn_res.0.8_noLabel.pdf', width=10,height=8)
DimPlot(object = Round4,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##DimPlot
p1 <- DimPlot(Round4, reduction = "umap", group.by = "stim")+
  scale_color_manual(values=c('deeppink','green','deep sky blue'))

p2 <- DimPlot(Round4, reduction = "umap", label = FALSE)

pdf('UMAP_DimPlot_integrated_snn_res.0.8.pdf', width=16,height=7)
plot_grid(p1, p2)
dev.off()

##cluster count 
colnames(Round4@meta.data)
cell.num <- table(Round4@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cluster_count_integrated_snn_res.0.8.pdf', width=10,height=8)
DimPlot(object = Round4,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###
pdf('UMAP_DimPlot_split_integrated_snn_res.0.8.pdf', width=20,height=7)
DimPlot(Round4, reduction = "umap", split.by = "stim")
dev.off()


##Cell number
colnames(Round4@meta.data)
Idents(object = Round4) <- ("number")

cell.num <- table(Round4@meta.data$number)
ClusterLabels = paste(names(cell.num), paste0("(n=", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cellcount_integrated_snn_res.0.8.pdf', width=14,height=8)
DimPlot(object = Round4,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

##
colnames(Round4@meta.data)
Idents(object = Round4) 

Idents(object = Round4) <- ("stim")
cell.numb <- table(Round4@meta.data$stim)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round4_category_integrated_snn_res.0.8.pdf',width=12,height=8)
DimPlot(object = Round4,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


setwd("../")

#####################################################################################################
# Save ####
save(Round4, file = "Round4_integrated_analyzed.Rda")
#################################################################################################
##Cluster gene expression cells by sample ####
setwd("~/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")
setwd("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

Idents(Round4)
colnames(Round4@meta.data)
Idents(object = Round4) <- ("ClusterNames_0.8")
Idents(object = Round4) <- ("stim")
levels(Idents(object = Round4))

Round4$ClusterNames_0.8_stim <- paste(Round4$ClusterNames_0.8, Round4$stim, sep = "_")

Idents(Round4) <- "ClusterNames_0.8"
#Idents(Round4) <- fct_relevel(Idents(Round4), sort)
levels(Idents(object = Round4))

cellnumbers <- as.data.frame(table(Round4@meta.data$ClusterNames_0.8_stim))
cellnumbers[is.na(cellnumbers)] <- 0

colnames(cellnumbers) <- c("number","cell_number")

write.csv(cellnumbers, file="cellnumbers.csv")

rm(list3)
rm(list3_classifier)
rm(list3_data)




## Bar graph for samples
list <- cellnumbers
list <- list[- grep("Psoriasis", list$number),]


list$Cluster_in_sample
list$Cluster <- list$number

list$Cluster= gsub("_Control", "", list$Cluster)
list$Cluster= gsub("_Fresh_tissue_Enzyme_digestion", "", list$Cluster)

list$Cluster <- factor(list$Cluster)
list$Cluster <- factor(list$Cluster, 
                       levels = (c("T_cell","Dendritic_cell", "Melanocyte" ,"Keratinocyte" ,"Fibroblast" 
                       )))

list$EmigratingvsFreshEnzyme <- list$number
list$EmigratingvsFreshEnzyme= gsub("Dendritic_cell_", "", list$EmigratingvsFreshEnzyme)
list$EmigratingvsFreshEnzyme= gsub("Fibroblast_", "", list$EmigratingvsFreshEnzyme)
list$EmigratingvsFreshEnzyme= gsub("Keratinocyte_", "", list$EmigratingvsFreshEnzyme)
list$EmigratingvsFreshEnzyme= gsub("Melanocyte_", "", list$EmigratingvsFreshEnzyme)
list$EmigratingvsFreshEnzyme= gsub("T_cell_", "", list$EmigratingvsFreshEnzyme)
list$EmigratingvsFreshEnzyme  <- factor(list$EmigratingvsFreshEnzyme )


list$total[list$EmigratingvsFreshEnzyme == "Control"] <- 4823
list$total[list$EmigratingvsFreshEnzyme == "Fresh_tissue_Enzyme_digestion"] <- 2799

list$proportion <- list$cell_number/list$total *100

library(reshape2)

list_for_ggplot <- list[,-c(2,6)]
list_for_ggplot <- melt(list_for_ggplot)


#install.packages("RColorBrewer")
library(RColorBrewer)

colourCount = length(unique(list$EmigratingvsFreshEnzyme))
getPalette(colourCount)

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = EmigratingvsFreshEnzyme)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  xlab("Samples") + ylab("Proportion of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) 


ggsave("Samples_in_cluster_bargraph_Round4_proportion.pdf", width = 10, height = 8)

######################################################################################################
load("D:/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round4_integrated_analyzed.Rda")

levels(Round4)
#########################################################################
# T-cell volcano plots
Idents(object = Round4) <- ("ClusterNames_0.8")
DefaultAssay(Round4) <- "RNA"
levels(Round4)

T_cells <- subset(Round4, idents = "T_cell")


Idents(object = T_cells) <- ("stim")
levels(T_cells)


# Volcano_plot control vs fresh_tissue_enzyme_digestion
DEG <- FindMarkers(T_cells, ident.1 =  "Fresh_tissue_Enzyme_digestion" , ident.2 ="Control", logfc.threshold = 0.1,min.pct =0 )

marker <- c('CD3D','TRAC','TRBC1',
            'ATL5E','ATP5L','ATP5G3','ATP5G2','ATP5I','ATP5C1')

pdf("Volcanoplot_Tcells_Control_Emigrating_vs_Fresh_enzyme_digestion-2.pdf", width=10,height=12)
p <- EnhancedVolcano(DEG,
                     lab = rownames(DEG),
                     x = 'avg_log2FC',
                     y = 'p_val',
                     title = 'Control: Emigrating cell vs. Fresh skin enzyme dissociation',
                     selectLab = c('CD3D','TRAC','TRBC1','ATP5E','ATP5G2','ATP5L'),
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

#########################################################################
# DC volcano plots
Idents(object = Round4) <- ("ClusterNames_0.8")
DefaultAssay(Round4) <- "RNA"
levels(Round4)

DC <- subset(Round4, idents = c("Dendritic_cell"  ))
Idents(object = DC) <- ("stim")
levels(DC)

marker <- c('HLA-DQA1','HLA-DQB1','HLA-DRB1','HLA-DRA','LYZ',
            'ATL5E','ATP5L','ATP5G3','ATP5G2','ATP5I','ATP5C1'
)



# Volcano_plot Mature DC control vs fresh_tissue_enzyme_digestion ####
DEG <- FindMarkers(DC, ident.1 =  "Fresh_tissue_Enzyme_digestion" , ident.2 ="Control", logfc.threshold = 0.1,min.pct =0 )


pdf("Volcanoplot_DC__Emigrating_vs_Fresh_enzyme_digestion02.pdf", width=10,height=12)
p <- EnhancedVolcano(DEG,
                     lab = rownames(DEG),
                     selectLab = marker,
                     x = 'avg_log2FC',
                     y = 'p_val',
                     title = 'DC - Control skin: Emigrating cell vs. Fresh skin enzyme dissociation',
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



#########################################################################
# KC volcano plots
Idents(object = Round4) <- ("ClusterNames_0.8")
DefaultAssay(Round4) <- "RNA"
levels(Round4)

KC <- subset(Round4, idents = c("Keratinocyte"  ))
Idents(object = KC) <- ("stim")
levels(KC)

marker <- c('SPRR2G','LCE3D','CDSN','KRT10','KRT1','KRT14','KRT15',
            'ATL5E','ATP5L','ATP5G3','ATP5G2','ATP5I','ATP5C1'
)



# Volcano_plot Mature DC control vs fresh_tissue_enzyme_digestion ####
DEG <- FindMarkers(DC, ident.1 =  "Fresh_tissue_Enzyme_digestion" , ident.2 ="Control", logfc.threshold = 0.1,min.pct =0 )


pdf("Volcanoplot_KC__Emigrating_vs_Fresh_enzyme_digestion02.pdf", width=10,height=12)
p <- EnhancedVolcano(DEG,
                     lab = rownames(DEG),
                     selectLab = marker,
                     x = 'avg_log2FC',
                     y = 'p_val',
                     title = 'KC - Control skin: Emigrating cell vs. Fresh skin enzyme dissociation',
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


