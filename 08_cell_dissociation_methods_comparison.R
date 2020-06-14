#load library ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for the downstream single-cell RNA sequencing data analy
library(Seurat)
######################################################################################################
#Start ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

# Please download all processed data files in one folder. In this code, the folder is ("~/Documents/GEO/processed_data")
setwd("~/Documents/GEO/processed_data")

# Phenotype data is shared in github 
phenotype <- read.delim("~/Dropbox/SIngle_cell/Rcode_JKIM/github_upload/phenotype_data.txt")

#Load the integrated data saved from previous code (05) as an integrated reference data
load("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020/Round2_integrated_analyzed.Rda")
DefaultAssay(Round2) <- "integrated"
######################################################################################################
# The analysis of individual skin biopsy sample's single-cell data ####

# Quality of single-cell data 
# The number of mitochondrially encoded genes and the total number of genes retrieved from each sample are calculated. 

# Cell type classification using an integrated reference.
# To calculate numbers of each cell type retrieved from individual skin biopsy samples, 
# we projected the integrated reference data onto each sample’s single-cell data. 
# We used TransferData function in Seurat R package (version 3.0, https://satijalab.org/seurat/) for the label transfer. 

######################################################################################################
## Load data Control05 ####
files <- Sys.glob('Control05*')
files <- files[c(4:6)]
dir.create("./Control05")
for(file in files) {file.copy(file, "./Control05")}
new_files <- gsub("Control05_","",files)
setwd("./Control05")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Control05 <- Read10X(data.dir = "./Control05")

Control05 <- CreateSeuratObject(counts = Control05, project = "Control05", min.cells = 3, min.features = 100)
Control05[["percent.mt"]] <- PercentageFeatureSet(object = Control05, pattern = "^MT-")
Control05@meta.data$stim <- as.character(phenotype[grep('Control05', phenotype$ID),2])[1]
Control05@meta.data$number <- as.character(phenotype[grep('Control05', phenotype$ID),1])[1]

unlink("./Control05", recursive = TRUE)

Control05 <- NormalizeData(Control05)
Control05 <- FindVariableFeatures(Control05, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control05)
Control05 <- ScaleData(Control05, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Control05,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Control05 <- AddMetaData(Control05, metadata = predictions)
write.csv(table(Control05$predicted.id), file = "Label_transfer_Control05.csv")

######################################################################################################
## Load data Control05F ####
files <- Sys.glob('Control05F*')
dir.create("./Control05F")
for(file in files) {file.copy(file, "./Control05F")}
new_files <- gsub("Control05F_","",files)
setwd("./Control05F")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Control05F <- Read10X(data.dir = "./Control05F")

Control05F <- CreateSeuratObject(counts = Control05F, project = "Control05F", min.cells = 3, min.features = 100)
Control05F[["percent.mt"]] <- PercentageFeatureSet(object = Control05F, pattern = "^MT-")
Control05F@meta.data$stim <- as.character(phenotype[grep('Control05F', phenotype$ID),2])
Control05F@meta.data$number <- as.character(phenotype[grep('Control05F', phenotype$ID),1])

unlink("./Control05F", recursive = TRUE)

Control05F <- NormalizeData(Control05F)
Control05F <- FindVariableFeatures(Control05F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control05F)
Control05F <- ScaleData(Control05F, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Control05F,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Control05F <- AddMetaData(Control05F, metadata = predictions)
write.csv(table(Control05F$predicted.id), file = "Label_transfer_Control05F.csv")

######################################################################################################
## Load data Psoriasis02 ####
files <- Sys.glob('Psoriasis02*')
files <- files[c(4:6)]
dir.create("./Psoriasis02")
for(file in files) {file.copy(file, "./Psoriasis02")}
new_files <- gsub("Psoriasis02_","",files)
setwd("./Psoriasis02")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis02 <- Read10X(data.dir = "./Psoriasis02")

Psoriasis02 <- CreateSeuratObject(counts = Psoriasis02, project = "Psoriasis02", min.cells = 3, min.features = 100)
Psoriasis02[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis02, pattern = "^MT-")
Psoriasis02@meta.data$stim <- as.character(phenotype[grep('Psoriasis02', phenotype$ID),2])[1]
Psoriasis02@meta.data$number <- as.character(phenotype[grep('Psoriasis02', phenotype$ID),1])[1]

unlink("./Psoriasis02", recursive = TRUE)

Psoriasis02 <- NormalizeData(Psoriasis02)
Psoriasis02 <- FindVariableFeatures(Psoriasis02, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis02)
Psoriasis02 <- ScaleData(Psoriasis02, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis02,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis02 <- AddMetaData(Psoriasis02, metadata = predictions)
write.csv(table(Psoriasis02$predicted.id), file = "Label_transfer_Psoriasis02.csv")

######################################################################################################
## Load data Psoriasis02F ####
files <- Sys.glob('Psoriasis02F*')
dir.create("./Psoriasis02F")
for(file in files) {file.copy(file, "./Psoriasis02F")}
new_files <- gsub("Psoriasis02F_","",files)
setwd("./Psoriasis02F")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis02F <- Read10X(data.dir = "./Psoriasis02F")

Psoriasis02F <- CreateSeuratObject(counts = Psoriasis02F, project = "Psoriasis02F", min.cells = 3, min.features = 100)
Psoriasis02F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis02F, pattern = "^MT-")
Psoriasis02F@meta.data$stim <- as.character(phenotype[grep('Psoriasis02F', phenotype$ID),2])
Psoriasis02F@meta.data$number <- as.character(phenotype[grep('Psoriasis02F', phenotype$ID),1])

unlink("./Psoriasis02F", recursive = TRUE)

Psoriasis02F <- NormalizeData(Psoriasis02F)
Psoriasis02F <- FindVariableFeatures(Psoriasis02F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis02F)
Psoriasis02F <- ScaleData(Psoriasis02F, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis02F,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis02F <- AddMetaData(Psoriasis02F, metadata = predictions)
write.csv(table(Psoriasis02F$predicted.id), file = "Label_transfer_Psoriasis02F.csv")

######################################################################################################
## Load data Psoriasis03 ####
files <- Sys.glob('Psoriasis03*')
files <- files[c(4:6)]
dir.create("./Psoriasis03")
for(file in files) {file.copy(file, "./Psoriasis03")}
new_files <- gsub("Psoriasis03_","",files)
setwd("./Psoriasis03")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis03 <- Read10X(data.dir = "./Psoriasis03")

Psoriasis03 <- CreateSeuratObject(counts = Psoriasis03, project = "Psoriasis03", min.cells = 3, min.features = 100)
Psoriasis03[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis03, pattern = "^MT-")
Psoriasis03@meta.data$stim <- as.character(phenotype[grep('Psoriasis03', phenotype$ID),2])[1]
Psoriasis03@meta.data$number <- as.character(phenotype[grep('Psoriasis03', phenotype$ID),1])[1]

unlink("./Psoriasis03", recursive = TRUE)

Psoriasis03 <- NormalizeData(Psoriasis03)
Psoriasis03 <- FindVariableFeatures(Psoriasis03, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis03)
Psoriasis03 <- ScaleData(Psoriasis03, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis03,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis03 <- AddMetaData(Psoriasis03, metadata = predictions)
write.csv(table(Psoriasis03$predicted.id), file = "Label_transfer_Psoriasis03.csv")
######################################################################################################
## Load data Psoriasis03F ####
files <- Sys.glob('Psoriasis03F*')
dir.create("./Psoriasis03F")
for(file in files) {file.copy(file, "./Psoriasis03F")}
new_files <- gsub("Psoriasis03F_","",files)
setwd("./Psoriasis03F")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis03F <- Read10X(data.dir = "./Psoriasis03F")

Psoriasis03F <- CreateSeuratObject(counts = Psoriasis03F, project = "Psoriasis03F", min.cells = 3, min.features = 100)
Psoriasis03F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis03F, pattern = "^MT-")
Psoriasis03F@meta.data$stim <- as.character(phenotype[grep('Psoriasis03F', phenotype$ID),2])
Psoriasis03F@meta.data$number <- as.character(phenotype[grep('Psoriasis03F', phenotype$ID),1])

unlink("./Psoriasis03F", recursive = TRUE)

Psoriasis03F <- NormalizeData(Psoriasis03F)
Psoriasis03F <- FindVariableFeatures(Psoriasis03F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis03F)
Psoriasis03F <- ScaleData(Psoriasis03F, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis03F,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis03F <- AddMetaData(Psoriasis03F, metadata = predictions)
write.csv(table(Psoriasis03F$predicted.id), file = "Label_transfer_Psoriasis03F.csv")
######################################################################################################
## Load data Psoriasis04 ####
files <- Sys.glob('Psoriasis04*')
files <- files[c(4:6)]
dir.create("./Psoriasis04")
for(file in files) {file.copy(file, "./Psoriasis04")}
new_files <- gsub("Psoriasis04_","",files)
setwd("./Psoriasis04")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis04 <- Read10X(data.dir = "./Psoriasis04")

Psoriasis04 <- CreateSeuratObject(counts = Psoriasis04, project = "Psoriasis04", min.cells = 3, min.features = 100)
Psoriasis04[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis04, pattern = "^MT-")
Psoriasis04@meta.data$stim <- as.character(phenotype[grep('Psoriasis04', phenotype$ID),2])[1]
Psoriasis04@meta.data$number <- as.character(phenotype[grep('Psoriasis04', phenotype$ID),1])[1]

unlink("./Psoriasis04", recursive = TRUE)

Psoriasis04 <- NormalizeData(Psoriasis04)
Psoriasis04 <- FindVariableFeatures(Psoriasis04, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis04)
Psoriasis04 <- ScaleData(Psoriasis04, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis04,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis04 <- AddMetaData(Psoriasis04, metadata = predictions)
write.csv(table(Psoriasis04$predicted.id), file = "Label_transfer_Psoriasis04.csv")
######################################################################################################
## Load data Psoriasis04F ####
files <- Sys.glob('Psoriasis04F*')
dir.create("./Psoriasis04F")
for(file in files) {file.copy(file, "./Psoriasis04F")}
new_files <- gsub("Psoriasis04F_","",files)
setwd("./Psoriasis04F")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis04F <- Read10X(data.dir = "./Psoriasis04F")

Psoriasis04F <- CreateSeuratObject(counts = Psoriasis04F, project = "Psoriasis04F", min.cells = 3, min.features = 100)
Psoriasis04F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis04F, pattern = "^MT-")
Psoriasis04F@meta.data$stim <- as.character(phenotype[grep('Psoriasis04F', phenotype$ID),2])
Psoriasis04F@meta.data$number <- as.character(phenotype[grep('Psoriasis04F', phenotype$ID),1])

unlink("./Psoriasis04F", recursive = TRUE)

Psoriasis04F <- NormalizeData(Psoriasis04F)
Psoriasis04F <- FindVariableFeatures(Psoriasis04F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis04F)
Psoriasis04F <- ScaleData(Psoriasis04F, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis04F,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis04F <- AddMetaData(Psoriasis04F, metadata = predictions)
write.csv(table(Psoriasis04F$predicted.id), file = "Label_transfer_Psoriasis04F.csv")
######################################################################################################
## Load data Psoriasis06 ####
files <- Sys.glob('Psoriasis06*')
files <- files[c(4:6)]
dir.create("./Psoriasis06")
for(file in files) {file.copy(file, "./Psoriasis06")}
new_files <- gsub("Psoriasis06_","",files)
setwd("./Psoriasis06")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis06 <- Read10X(data.dir = "./Psoriasis06")

Psoriasis06 <- CreateSeuratObject(counts = Psoriasis06, project = "Psoriasis06", min.cells = 3, min.features = 100)
Psoriasis06[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis06, pattern = "^MT-")
Psoriasis06@meta.data$stim <- as.character(phenotype[grep('Psoriasis06', phenotype$ID),2])[1]
Psoriasis06@meta.data$number <- as.character(phenotype[grep('Psoriasis06', phenotype$ID),1])[1]

unlink("./Psoriasis06", recursive = TRUE)

Psoriasis06 <- NormalizeData(Psoriasis06)
Psoriasis06 <- FindVariableFeatures(Psoriasis06, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis06)
Psoriasis06 <- ScaleData(Psoriasis06, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis06,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis06 <- AddMetaData(Psoriasis06, metadata = predictions)
write.csv(table(Psoriasis06$predicted.id), file = "Label_transfer_Psoriasis06.csv")
######################################################################################################
## Load data Psoriasis06F ####
files <- Sys.glob('Psoriasis06F*')
dir.create("./Psoriasis06F")
for(file in files) {file.copy(file, "./Psoriasis06F")}
new_files <- gsub("Psoriasis06F_","",files)
setwd("./Psoriasis06F")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis06F <- Read10X(data.dir = "./Psoriasis06F")

Psoriasis06F <- CreateSeuratObject(counts = Psoriasis06F, project = "Psoriasis06F", min.cells = 3, min.features = 100)
Psoriasis06F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis06F, pattern = "^MT-")
Psoriasis06F@meta.data$stim <- as.character(phenotype[grep('Psoriasis06F', phenotype$ID),2])
Psoriasis06F@meta.data$number <- as.character(phenotype[grep('Psoriasis06F', phenotype$ID),1])

unlink("./Psoriasis06F", recursive = TRUE)

Psoriasis06F <- NormalizeData(Psoriasis06F)
Psoriasis06F <- FindVariableFeatures(Psoriasis06F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis06F)
Psoriasis06F <- ScaleData(Psoriasis06F, features = all.genes)

# Label Transfer
transfer.anchors <- FindTransferAnchors(reference = Round2, query = Psoriasis06F,  dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = Round2$ClusterNames_0.8, dims = 1:30)
Psoriasis06F <- AddMetaData(Psoriasis06F, metadata = predictions)
write.csv(table(Psoriasis06F$predicted.id), file = "Label_transfer_Psoriasis06F.csv")

######################################################################################################
## (Figure 1) Merge data set and compare cell dissociation methods ####
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020")

Total <- merge(x = Psoriasis02, y = c(Psoriasis02F, Psoriasis03, Psoriasis03F, Psoriasis04,Psoriasis04F,Psoriasis06,Psoriasis06F,
                                       Control05, Control05F), 
add.cell.ids = c( "Psoriasis_02", "Psoriasis_02F", "Psoriasis_03", "Psoriasis_03F", "Psoriasis_04","Psoriasis_04F",
                  "Psoriasis_06","Psoriasis_06F",  "Control_05" ,"Control_05F"), 
project = "Total")


## Vnplot 
pdf('Vnplot_Cell_dissociation_comparison.pdf',width=30,height=5)
VlnPlot(object = Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
dev.off()
######################################################################################################
