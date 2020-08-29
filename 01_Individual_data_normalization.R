#load library ####
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) 
# installed in R (version 3.6.2, https://www.r-project.org/) for the downstream single-cell RNA sequencing data analy
library(Seurat)
######################################################################################################
#Start ####
rm(list = ls())
setwd("/Users/jkim05/Dropbox/10X_JKIM/aggregation/Single_cell_paper_submission_06_12_2020") ##change directory 

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
# Individual single-cell data loading & quality control ####
# Genes expressed in <3 cells, and cells with <100 or > 5,000 genes, 
# and a mitochondrial gene percentage of >25% were filtered out to eliminate partial cells and doublets. 
# Ubiquitously expressed ribosomal protein-coding (RPS and RPL) and MALAT noncoding RNA, miRNA, and snoRNA genes 
# were also excluded from analysis. 
# Seurat objects were created, followed by normalizing data, scaling data, and finding variable 2,000 genes. 

# This strategy has been used for the coauthors' previous skin single-cell RNA sequencing studies
# 1. Der, E. et al. Tubular cell and keratinocyte single-cell transcriptomics applied to lupus nephritis reveal type I IFN and fibrosis relevant pathways. Nature immunology 20, 915-927 (2019)
# 2. He, H. et al. Single-cell transcriptome analysis of human skin identifies novel fibroblast subpopulation and enrichment of immune subsets in atopic dermatitis. The Journal of allergy and clinical immunology (2020).


# Please download all processed data files in one folder. In this code, the folder is ("~/Documents/GEO/processed_data")
setwd("~/Documents/GEO/processed_data")

# Phenotype data is shared in github 
phenotype <- read.delim("~/Dropbox/SIngle_cell/Rcode_JKIM/github_upload/phenotype_data.txt")

######################################################################################################
## Load data Control01 ####
files <- Sys.glob('Control01*')
dir.create("./Control01")
for(file in files) {file.copy(file, "./Control01")}
new_files <- gsub("Control01_","",files)
setwd("./Control01")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Control01 <- Read10X(data.dir = "./Control01")

Control01 <- CreateSeuratObject(counts = Control01, project = "Control01", min.cells = 3, min.features = 100)
Control01[["percent.mt"]] <- PercentageFeatureSet(object = Control01, pattern = "^MT-")
Control01@meta.data$stim <- as.character(phenotype[grep('Control01', phenotype$ID),2])
Control01@meta.data$number <- as.character(phenotype[grep('Control01', phenotype$ID),1])

Control01_morethan25percentMT <- subset(Control01, subset = percent.mt > 25)
Control01_morethan25percentMTbc <- colnames(Control01_morethan25percentMT@assays$RNA@data)
rm(Control01)
rm(Control01_morethan25percentMT)

Control01 <- Read10X(data.dir = "./Control01")
Control01<-remove.genes(Control01, Control01_morethan25percentMTbc)

Control01 <- CreateSeuratObject(counts = Control01, project = "Control01", min.cells = 3, min.features = 100)
Control01[["percent.mt"]] <- PercentageFeatureSet(object = Control01, pattern = "^MT-")
Control01 <- subset(Control01, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control01@meta.data$stim <- as.character(phenotype[grep('Control01', phenotype$ID),2])
Control01@meta.data$number <- as.character(phenotype[grep('Control01', phenotype$ID),1])
rm(Control01_morethan25percentMTbc)

Control01 <- NormalizeData(Control01)
Control01 <- FindVariableFeatures(Control01, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control01)
Control01 <- ScaleData(Control01, features = all.genes)

rm(all.genes)
unlink("./Control01", recursive = TRUE)
######################################################################################################
## Load data Control02 ####
files <- Sys.glob('Control02*')
dir.create("./Control02")
for(file in files) {file.copy(file, "./Control02")}
new_files <- gsub("Control02_","",files)
setwd("./Control02")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Control02 <- Read10X(data.dir = "./Control02")

Control02 <- CreateSeuratObject(counts = Control02, project = "Control02", min.cells = 3, min.features = 100)
Control02[["percent.mt"]] <- PercentageFeatureSet(object = Control02, pattern = "^MT-")
Control02@meta.data$stim <- as.character(phenotype[grep('Control02', phenotype$ID),2])
Control02@meta.data$number <- as.character(phenotype[grep('Control02', phenotype$ID),1])

Control02_morethan25percentMT <- subset(Control02, subset = percent.mt > 25)
Control02_morethan25percentMTbc <- colnames(Control02_morethan25percentMT@assays$RNA@data)
rm(Control02)
rm(Control02_morethan25percentMT)

Control02 <- Read10X(data.dir = "./Control02")
Control02<-remove.genes(Control02, Control02_morethan25percentMTbc)

Control02 <- CreateSeuratObject(counts = Control02, project = "Control02", min.cells = 3, min.features = 100)
Control02[["percent.mt"]] <- PercentageFeatureSet(object = Control02, pattern = "^MT-")
Control02 <- subset(Control02, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control02@meta.data$stim <- as.character(phenotype[grep('Control02', phenotype$ID),2])
Control02@meta.data$number <- as.character(phenotype[grep('Control02', phenotype$ID),1])
rm(Control02_morethan25percentMTbc)

Control02 <- NormalizeData(Control02)
Control02 <- FindVariableFeatures(Control02, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control02)
Control02 <- ScaleData(Control02, features = all.genes)

rm(all.genes)
unlink("./Control02", recursive = TRUE)
######################################################################################################
## Load data Control03 ####
files <- Sys.glob('Control03*')
dir.create("./Control03")
for(file in files) {file.copy(file, "./Control03")}
new_files <- gsub("Control03_","",files)
setwd("./Control03")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Control03 <- Read10X(data.dir = "./Control03")

Control03 <- CreateSeuratObject(counts = Control03, project = "Control03", min.cells = 3, min.features = 100)
Control03[["percent.mt"]] <- PercentageFeatureSet(object = Control03, pattern = "^MT-")
Control03@meta.data$stim <- as.character(phenotype[grep('Control03', phenotype$ID),2])
Control03@meta.data$number <- as.character(phenotype[grep('Control03', phenotype$ID),1])

Control03_morethan25percentMT <- subset(Control03, subset = percent.mt > 25)
Control03_morethan25percentMTbc <- colnames(Control03_morethan25percentMT@assays$RNA@data)
rm(Control03)
rm(Control03_morethan25percentMT)

Control03 <- Read10X(data.dir = "./Control03")
Control03<-remove.genes(Control03, Control03_morethan25percentMTbc)

Control03 <- CreateSeuratObject(counts = Control03, project = "Control03", min.cells = 3, min.features = 100)
Control03[["percent.mt"]] <- PercentageFeatureSet(object = Control03, pattern = "^MT-")
Control03 <- subset(Control03, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control03@meta.data$stim <- as.character(phenotype[grep('Control03', phenotype$ID),2])
Control03@meta.data$number <- as.character(phenotype[grep('Control03', phenotype$ID),1])
rm(Control03_morethan25percentMTbc)

Control03 <- NormalizeData(Control03)
Control03 <- FindVariableFeatures(Control03, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control03)
Control03 <- ScaleData(Control03, features = all.genes)

rm(all.genes)
unlink("./Control03", recursive = TRUE)
######################################################################################################
## Load data Control04 ####
files <- Sys.glob('Control04*')
dir.create("./Control04")
for(file in files) {file.copy(file, "./Control04")}
new_files <- gsub("Control04_","",files)
setwd("./Control04")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Control04 <- Read10X(data.dir = "./Control04")

Control04 <- CreateSeuratObject(counts = Control04, project = "Control04", min.cells = 3, min.features = 100)
Control04[["percent.mt"]] <- PercentageFeatureSet(object = Control04, pattern = "^MT-")
Control04@meta.data$stim <- as.character(phenotype[grep('Control04', phenotype$ID),2])
Control04@meta.data$number <- as.character(phenotype[grep('Control04', phenotype$ID),1])

Control04_morethan25percentMT <- subset(Control04, subset = percent.mt > 25)
Control04_morethan25percentMTbc <- colnames(Control04_morethan25percentMT@assays$RNA@data)
rm(Control04)
rm(Control04_morethan25percentMT)

Control04 <- Read10X(data.dir = "./Control04")
Control04<-remove.genes(Control04, Control04_morethan25percentMTbc)

Control04 <- CreateSeuratObject(counts = Control04, project = "Control04", min.cells = 3, min.features = 100)
Control04[["percent.mt"]] <- PercentageFeatureSet(object = Control04, pattern = "^MT-")
Control04 <- subset(Control04, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control04@meta.data$stim <- as.character(phenotype[grep('Control04', phenotype$ID),2])
Control04@meta.data$number <- as.character(phenotype[grep('Control04', phenotype$ID),1])
rm(Control04_morethan25percentMTbc)

Control04 <- NormalizeData(Control04)
Control04 <- FindVariableFeatures(Control04, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control04)
Control04 <- ScaleData(Control04, features = all.genes)

rm(all.genes)
unlink("./Control04", recursive = TRUE)
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

Control05_morethan25percentMT <- subset(Control05, subset = percent.mt > 25)
Control05_morethan25percentMTbc <- colnames(Control05_morethan25percentMT@assays$RNA@data)
rm(Control05)
rm(Control05_morethan25percentMT)

Control05 <- Read10X(data.dir = "./Control05")
Control05<-remove.genes(Control05, Control05_morethan25percentMTbc)

Control05 <- CreateSeuratObject(counts = Control05, project = "Control05", min.cells = 3, min.features = 100)
Control05[["percent.mt"]] <- PercentageFeatureSet(object = Control05, pattern = "^MT-")
Control05 <- subset(Control05, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control05@meta.data$stim <- as.character(phenotype[grep('Control05', phenotype$ID),2])[1]
Control05@meta.data$number <- as.character(phenotype[grep('Control05', phenotype$ID),1])[1]
rm(Control05_morethan25percentMTbc)

Control05 <- NormalizeData(Control05)
Control05 <- FindVariableFeatures(Control05, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control05)
Control05 <- ScaleData(Control05, features = all.genes)

rm(all.genes)
unlink("./Control05", recursive = TRUE)
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

Control05F_morethan25percentMT <- subset(Control05F, subset = percent.mt > 25)
Control05F_morethan25percentMTbc <- colnames(Control05F_morethan25percentMT@assays$RNA@data)
rm(Control05F)
rm(Control05F_morethan25percentMT)

Control05F <- Read10X(data.dir = "./Control05F")
Control05F<-remove.genes(Control05F, Control05F_morethan25percentMTbc)

Control05F <- CreateSeuratObject(counts = Control05F, project = "Control05F", min.cells = 3, min.features = 100)
Control05F[["percent.mt"]] <- PercentageFeatureSet(object = Control05F, pattern = "^MT-")
Control05F <- subset(Control05F, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control05F@meta.data$stim <- as.character(phenotype[grep('Control05F', phenotype$ID),2])
Control05F@meta.data$number <- as.character(phenotype[grep('Control05F', phenotype$ID),1])
rm(Control05F_morethan25percentMTbc)

Control05F <- NormalizeData(Control05F)
Control05F <- FindVariableFeatures(Control05F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control05F)
Control05F <- ScaleData(Control05F, features = all.genes)

rm(all.genes)
unlink("./Control05F", recursive = TRUE)
######################################################################################################
## Load data Psoriasis01 ####
files <- Sys.glob('Psoriasis01*')
dir.create("./Psoriasis01")
for(file in files) {file.copy(file, "./Psoriasis01")}
new_files <- gsub("Psoriasis01_","",files)
setwd("./Psoriasis01")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis01 <- Read10X(data.dir = "./Psoriasis01")

Psoriasis01 <- CreateSeuratObject(counts = Psoriasis01, project = "Psoriasis01", min.cells = 3, min.features = 100)
Psoriasis01[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis01, pattern = "^MT-")
Psoriasis01@meta.data$stim <- as.character(phenotype[grep('Psoriasis01', phenotype$ID),2])
Psoriasis01@meta.data$number <- as.character(phenotype[grep('Psoriasis01', phenotype$ID),1])

Psoriasis01_morethan25percentMT <- subset(Psoriasis01, subset = percent.mt > 25)
Psoriasis01_morethan25percentMTbc <- colnames(Psoriasis01_morethan25percentMT@assays$RNA@data)
rm(Psoriasis01)
rm(Psoriasis01_morethan25percentMT)

Psoriasis01 <- Read10X(data.dir = "./Psoriasis01")
Psoriasis01<-remove.genes(Psoriasis01, Psoriasis01_morethan25percentMTbc)

Psoriasis01 <- CreateSeuratObject(counts = Psoriasis01, project = "Psoriasis01", min.cells = 3, min.features = 100)
Psoriasis01[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis01, pattern = "^MT-")
Psoriasis01 <- subset(Psoriasis01, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis01@meta.data$stim <- as.character(phenotype[grep('Psoriasis01', phenotype$ID),2])
Psoriasis01@meta.data$number <- as.character(phenotype[grep('Psoriasis01', phenotype$ID),1])
rm(Psoriasis01_morethan25percentMTbc)

Psoriasis01 <- NormalizeData(Psoriasis01)
Psoriasis01 <- FindVariableFeatures(Psoriasis01, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis01)
Psoriasis01 <- ScaleData(Psoriasis01, features = all.genes)

rm(all.genes)
unlink("./Psoriasis01", recursive = TRUE)
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

Psoriasis02_morethan25percentMT <- subset(Psoriasis02, subset = percent.mt > 25)
Psoriasis02_morethan25percentMTbc <- colnames(Psoriasis02_morethan25percentMT@assays$RNA@data)
rm(Psoriasis02)
rm(Psoriasis02_morethan25percentMT)

Psoriasis02 <- Read10X(data.dir = "./Psoriasis02")
Psoriasis02<-remove.genes(Psoriasis02, Psoriasis02_morethan25percentMTbc)

Psoriasis02 <- CreateSeuratObject(counts = Psoriasis02, project = "Psoriasis02", min.cells = 3, min.features = 100)
Psoriasis02[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis02, pattern = "^MT-")
Psoriasis02 <- subset(Psoriasis02, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis02@meta.data$stim <- as.character(phenotype[grep('Psoriasis02', phenotype$ID),2])[1]
Psoriasis02@meta.data$number <- as.character(phenotype[grep('Psoriasis02', phenotype$ID),1])[1]
rm(Psoriasis02_morethan25percentMTbc)

Psoriasis02 <- NormalizeData(Psoriasis02)
Psoriasis02 <- FindVariableFeatures(Psoriasis02, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis02)
Psoriasis02 <- ScaleData(Psoriasis02, features = all.genes)

rm(all.genes)
unlink("./Psoriasis02", recursive = TRUE)
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

Psoriasis02F_morethan25percentMT <- subset(Psoriasis02F, subset = percent.mt > 25)
Psoriasis02F_morethan25percentMTbc <- colnames(Psoriasis02F_morethan25percentMT@assays$RNA@data)
rm(Psoriasis02F)
rm(Psoriasis02F_morethan25percentMT)

Psoriasis02F <- Read10X(data.dir = "./Psoriasis02F")
Psoriasis02F<-remove.genes(Psoriasis02F, Psoriasis02F_morethan25percentMTbc)

Psoriasis02F <- CreateSeuratObject(counts = Psoriasis02F, project = "Psoriasis02F", min.cells = 3, min.features = 100)
Psoriasis02F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis02F, pattern = "^MT-")
Psoriasis02F <- subset(Psoriasis02F, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis02F@meta.data$stim <- as.character(phenotype[grep('Psoriasis02F', phenotype$ID),2])
Psoriasis02F@meta.data$number <- as.character(phenotype[grep('Psoriasis02F', phenotype$ID),1])
rm(Psoriasis02F_morethan25percentMTbc)

Psoriasis02F <- NormalizeData(Psoriasis02F)
Psoriasis02F <- FindVariableFeatures(Psoriasis02F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis02F)
Psoriasis02F <- ScaleData(Psoriasis02F, features = all.genes)

rm(all.genes)
unlink("./Psoriasis02F", recursive = TRUE)
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

Psoriasis03_morethan25percentMT <- subset(Psoriasis03, subset = percent.mt > 25)
Psoriasis03_morethan25percentMTbc <- colnames(Psoriasis03_morethan25percentMT@assays$RNA@data)
rm(Psoriasis03)
rm(Psoriasis03_morethan25percentMT)

Psoriasis03 <- Read10X(data.dir = "./Psoriasis03")
Psoriasis03<-remove.genes(Psoriasis03, Psoriasis03_morethan25percentMTbc)

Psoriasis03 <- CreateSeuratObject(counts = Psoriasis03, project = "Psoriasis03", min.cells = 3, min.features = 100)
Psoriasis03[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis03, pattern = "^MT-")
Psoriasis03 <- subset(Psoriasis03, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis03@meta.data$stim <- as.character(phenotype[grep('Psoriasis03', phenotype$ID),2])[1]
Psoriasis03@meta.data$number <- as.character(phenotype[grep('Psoriasis03', phenotype$ID),1])[1]
rm(Psoriasis03_morethan25percentMTbc)

Psoriasis03 <- NormalizeData(Psoriasis03)
Psoriasis03 <- FindVariableFeatures(Psoriasis03, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis03)
Psoriasis03 <- ScaleData(Psoriasis03, features = all.genes)

rm(all.genes)
unlink("./Psoriasis03", recursive = TRUE)
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

Psoriasis03F_morethan25percentMT <- subset(Psoriasis03F, subset = percent.mt > 25)
Psoriasis03F_morethan25percentMTbc <- colnames(Psoriasis03F_morethan25percentMT@assays$RNA@data)
rm(Psoriasis03F)
rm(Psoriasis03F_morethan25percentMT)

Psoriasis03F <- Read10X(data.dir = "./Psoriasis03F")
Psoriasis03F<-remove.genes(Psoriasis03F, Psoriasis03F_morethan25percentMTbc)

Psoriasis03F <- CreateSeuratObject(counts = Psoriasis03F, project = "Psoriasis03F", min.cells = 3, min.features = 100)
Psoriasis03F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis03F, pattern = "^MT-")
Psoriasis03F <- subset(Psoriasis03F, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis03F@meta.data$stim <- as.character(phenotype[grep('Psoriasis03F', phenotype$ID),2])
Psoriasis03F@meta.data$number <- as.character(phenotype[grep('Psoriasis03F', phenotype$ID),1])
rm(Psoriasis03F_morethan25percentMTbc)

Psoriasis03F <- NormalizeData(Psoriasis03F)
Psoriasis03F <- FindVariableFeatures(Psoriasis03F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis03F)
Psoriasis03F <- ScaleData(Psoriasis03F, features = all.genes)

rm(all.genes)
unlink("./Psoriasis03F", recursive = TRUE)
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

Psoriasis04_morethan25percentMT <- subset(Psoriasis04, subset = percent.mt > 25)
Psoriasis04_morethan25percentMTbc <- colnames(Psoriasis04_morethan25percentMT@assays$RNA@data)
rm(Psoriasis04)
rm(Psoriasis04_morethan25percentMT)

Psoriasis04 <- Read10X(data.dir = "./Psoriasis04")
Psoriasis04<-remove.genes(Psoriasis04, Psoriasis04_morethan25percentMTbc)

Psoriasis04 <- CreateSeuratObject(counts = Psoriasis04, project = "Psoriasis04", min.cells = 3, min.features = 100)
Psoriasis04[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis04, pattern = "^MT-")
Psoriasis04 <- subset(Psoriasis04, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis04@meta.data$stim <- as.character(phenotype[grep('Psoriasis04', phenotype$ID),2])[1]
Psoriasis04@meta.data$number <- as.character(phenotype[grep('Psoriasis04', phenotype$ID),1])[1]
rm(Psoriasis04_morethan25percentMTbc)

Psoriasis04 <- NormalizeData(Psoriasis04)
Psoriasis04 <- FindVariableFeatures(Psoriasis04, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis04)
Psoriasis04 <- ScaleData(Psoriasis04, features = all.genes)

rm(all.genes)
unlink("./Psoriasis04", recursive = TRUE)
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

Psoriasis04F_morethan25percentMT <- subset(Psoriasis04F, subset = percent.mt > 25)
Psoriasis04F_morethan25percentMTbc <- colnames(Psoriasis04F_morethan25percentMT@assays$RNA@data)
rm(Psoriasis04F)
rm(Psoriasis04F_morethan25percentMT)

Psoriasis04F <- Read10X(data.dir = "./Psoriasis04F")
Psoriasis04F<-remove.genes(Psoriasis04F, Psoriasis04F_morethan25percentMTbc)

Psoriasis04F <- CreateSeuratObject(counts = Psoriasis04F, project = "Psoriasis04F", min.cells = 3, min.features = 100)
Psoriasis04F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis04F, pattern = "^MT-")
Psoriasis04F <- subset(Psoriasis04F, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis04F@meta.data$stim <- as.character(phenotype[grep('Psoriasis04F', phenotype$ID),2])
Psoriasis04F@meta.data$number <- as.character(phenotype[grep('Psoriasis04F', phenotype$ID),1])
rm(Psoriasis04F_morethan25percentMTbc)

Psoriasis04F <- NormalizeData(Psoriasis04F)
Psoriasis04F <- FindVariableFeatures(Psoriasis04F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis04F)
Psoriasis04F <- ScaleData(Psoriasis04F, features = all.genes)

rm(all.genes)
unlink("./Psoriasis04F", recursive = TRUE)
######################################################################################################
## Load data Psoriasis05 ####
files <- Sys.glob('Psoriasis05*')
dir.create("./Psoriasis05")
for(file in files) {file.copy(file, "./Psoriasis05")}
new_files <- gsub("Psoriasis05_","",files)
setwd("./Psoriasis05")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis05 <- Read10X(data.dir = "./Psoriasis05")

Psoriasis05 <- CreateSeuratObject(counts = Psoriasis05, project = "Psoriasis05", min.cells = 3, min.features = 100)
Psoriasis05[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis05, pattern = "^MT-")
Psoriasis05@meta.data$stim <- as.character(phenotype[grep('Psoriasis05', phenotype$ID),2])
Psoriasis05@meta.data$number <- as.character(phenotype[grep('Psoriasis05', phenotype$ID),1])

Psoriasis05_morethan25percentMT <- subset(Psoriasis05, subset = percent.mt > 25)
Psoriasis05_morethan25percentMTbc <- colnames(Psoriasis05_morethan25percentMT@assays$RNA@data)
rm(Psoriasis05)
rm(Psoriasis05_morethan25percentMT)

Psoriasis05 <- Read10X(data.dir = "./Psoriasis05")
Psoriasis05<-remove.genes(Psoriasis05, Psoriasis05_morethan25percentMTbc)

Psoriasis05 <- CreateSeuratObject(counts = Psoriasis05, project = "Psoriasis05", min.cells = 3, min.features = 100)
Psoriasis05[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis05, pattern = "^MT-")
Psoriasis05 <- subset(Psoriasis05, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis05@meta.data$stim <- as.character(phenotype[grep('Psoriasis05', phenotype$ID),2])
Psoriasis05@meta.data$number <- as.character(phenotype[grep('Psoriasis05', phenotype$ID),1])
rm(Psoriasis05_morethan25percentMTbc)

Psoriasis05 <- NormalizeData(Psoriasis05)
Psoriasis05 <- FindVariableFeatures(Psoriasis05, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis05)
Psoriasis05 <- ScaleData(Psoriasis05, features = all.genes)

rm(all.genes)
unlink("./Psoriasis05", recursive = TRUE)
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

Psoriasis06_morethan25percentMT <- subset(Psoriasis06, subset = percent.mt > 25)
Psoriasis06_morethan25percentMTbc <- colnames(Psoriasis06_morethan25percentMT@assays$RNA@data)
rm(Psoriasis06)
rm(Psoriasis06_morethan25percentMT)

Psoriasis06 <- Read10X(data.dir = "./Psoriasis06")
Psoriasis06<-remove.genes(Psoriasis06, Psoriasis06_morethan25percentMTbc)

Psoriasis06 <- CreateSeuratObject(counts = Psoriasis06, project = "Psoriasis06", min.cells = 3, min.features = 100)
Psoriasis06[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis06, pattern = "^MT-")
Psoriasis06 <- subset(Psoriasis06, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis06@meta.data$stim <- as.character(phenotype[grep('Psoriasis06', phenotype$ID),2])[1]
Psoriasis06@meta.data$number <- as.character(phenotype[grep('Psoriasis06', phenotype$ID),1])[1]
rm(Psoriasis06_morethan25percentMTbc)

Psoriasis06 <- NormalizeData(Psoriasis06)
Psoriasis06 <- FindVariableFeatures(Psoriasis06, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis06)
Psoriasis06 <- ScaleData(Psoriasis06, features = all.genes)

rm(all.genes)
unlink("./Psoriasis06", recursive = TRUE)
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

Psoriasis06F_morethan25percentMT <- subset(Psoriasis06F, subset = percent.mt > 25)
Psoriasis06F_morethan25percentMTbc <- colnames(Psoriasis06F_morethan25percentMT@assays$RNA@data)
rm(Psoriasis06F)
rm(Psoriasis06F_morethan25percentMT)

Psoriasis06F <- Read10X(data.dir = "./Psoriasis06F")
Psoriasis06F<-remove.genes(Psoriasis06F, Psoriasis06F_morethan25percentMTbc)

Psoriasis06F <- CreateSeuratObject(counts = Psoriasis06F, project = "Psoriasis06F", min.cells = 3, min.features = 100)
Psoriasis06F[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis06F, pattern = "^MT-")
Psoriasis06F <- subset(Psoriasis06F, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis06F@meta.data$stim <- as.character(phenotype[grep('Psoriasis06F', phenotype$ID),2])
Psoriasis06F@meta.data$number <- as.character(phenotype[grep('Psoriasis06F', phenotype$ID),1])
rm(Psoriasis06F_morethan25percentMTbc)

Psoriasis06F <- NormalizeData(Psoriasis06F)
Psoriasis06F <- FindVariableFeatures(Psoriasis06F, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis06F)
Psoriasis06F <- ScaleData(Psoriasis06F, features = all.genes)

rm(all.genes)
unlink("./Psoriasis06F", recursive = TRUE)
######################################################################################################
## Load data Psoriasis07 ####
files <- Sys.glob('Psoriasis07*')
dir.create("./Psoriasis07")
for(file in files) {file.copy(file, "./Psoriasis07")}
new_files <- gsub("Psoriasis07_","",files)
setwd("./Psoriasis07")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis07 <- Read10X(data.dir = "./Psoriasis07")

Psoriasis07 <- CreateSeuratObject(counts = Psoriasis07, project = "Psoriasis07", min.cells = 3, min.features = 100)
Psoriasis07[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis07, pattern = "^MT-")
Psoriasis07@meta.data$stim <- as.character(phenotype[grep('Psoriasis07', phenotype$ID),2])
Psoriasis07@meta.data$number <- as.character(phenotype[grep('Psoriasis07', phenotype$ID),1])

Psoriasis07_morethan25percentMT <- subset(Psoriasis07, subset = percent.mt > 25)
Psoriasis07_morethan25percentMTbc <- colnames(Psoriasis07_morethan25percentMT@assays$RNA@data)
rm(Psoriasis07)
rm(Psoriasis07_morethan25percentMT)

Psoriasis07 <- Read10X(data.dir = "./Psoriasis07")
Psoriasis07<-remove.genes(Psoriasis07, Psoriasis07_morethan25percentMTbc)

Psoriasis07 <- CreateSeuratObject(counts = Psoriasis07, project = "Psoriasis07", min.cells = 3, min.features = 100)
Psoriasis07[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis07, pattern = "^MT-")
Psoriasis07 <- subset(Psoriasis07, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis07@meta.data$stim <- as.character(phenotype[grep('Psoriasis07', phenotype$ID),2])
Psoriasis07@meta.data$number <- as.character(phenotype[grep('Psoriasis07', phenotype$ID),1])
rm(Psoriasis07_morethan25percentMTbc)

Psoriasis07 <- NormalizeData(Psoriasis07)
Psoriasis07 <- FindVariableFeatures(Psoriasis07, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis07)
Psoriasis07 <- ScaleData(Psoriasis07, features = all.genes)

rm(all.genes)
unlink("./Psoriasis07", recursive = TRUE)
######################################################################################################
## Load data Psoriasis08 ####
files <- Sys.glob('Psoriasis08*')
dir.create("./Psoriasis08")
for(file in files) {file.copy(file, "./Psoriasis08")}
new_files <- gsub("Psoriasis08_","",files)
setwd("./Psoriasis08")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis08 <- Read10X(data.dir = "./Psoriasis08")

Psoriasis08 <- CreateSeuratObject(counts = Psoriasis08, project = "Psoriasis08", min.cells = 3, min.features = 100)
Psoriasis08[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis08, pattern = "^MT-")
Psoriasis08@meta.data$stim <- as.character(phenotype[grep('Psoriasis08', phenotype$ID),2])
Psoriasis08@meta.data$number <- as.character(phenotype[grep('Psoriasis08', phenotype$ID),1])

Psoriasis08_morethan25percentMT <- subset(Psoriasis08, subset = percent.mt > 25)
Psoriasis08_morethan25percentMTbc <- colnames(Psoriasis08_morethan25percentMT@assays$RNA@data)
rm(Psoriasis08)
rm(Psoriasis08_morethan25percentMT)

Psoriasis08 <- Read10X(data.dir = "./Psoriasis08")
Psoriasis08<-remove.genes(Psoriasis08, Psoriasis08_morethan25percentMTbc)

Psoriasis08 <- CreateSeuratObject(counts = Psoriasis08, project = "Psoriasis08", min.cells = 3, min.features = 100)
Psoriasis08[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis08, pattern = "^MT-")
Psoriasis08 <- subset(Psoriasis08, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis08@meta.data$stim <- as.character(phenotype[grep('Psoriasis08', phenotype$ID),2])
Psoriasis08@meta.data$number <- as.character(phenotype[grep('Psoriasis08', phenotype$ID),1])
rm(Psoriasis08_morethan25percentMTbc)

Psoriasis08 <- NormalizeData(Psoriasis08)
Psoriasis08 <- FindVariableFeatures(Psoriasis08, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis08)
Psoriasis08 <- ScaleData(Psoriasis08, features = all.genes)

rm(all.genes)
unlink("./Psoriasis08", recursive = TRUE)
######################################################################################################
## Load data Psoriasis09 ####
files <- Sys.glob('Psoriasis09*')
dir.create("./Psoriasis09")
for(file in files) {file.copy(file, "./Psoriasis09")}
new_files <- gsub("Psoriasis09_","",files)
setwd("./Psoriasis09")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis09 <- Read10X(data.dir = "./Psoriasis09")

Psoriasis09 <- CreateSeuratObject(counts = Psoriasis09, project = "Psoriasis09", min.cells = 3, min.features = 100)
Psoriasis09[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis09, pattern = "^MT-")
Psoriasis09@meta.data$stim <- as.character(phenotype[grep('Psoriasis09', phenotype$ID),2])
Psoriasis09@meta.data$number <- as.character(phenotype[grep('Psoriasis09', phenotype$ID),1])

Psoriasis09_morethan25percentMT <- subset(Psoriasis09, subset = percent.mt > 25)
Psoriasis09_morethan25percentMTbc <- colnames(Psoriasis09_morethan25percentMT@assays$RNA@data)
rm(Psoriasis09)
rm(Psoriasis09_morethan25percentMT)

Psoriasis09 <- Read10X(data.dir = "./Psoriasis09")
Psoriasis09<-remove.genes(Psoriasis09, Psoriasis09_morethan25percentMTbc)

Psoriasis09 <- CreateSeuratObject(counts = Psoriasis09, project = "Psoriasis09", min.cells = 3, min.features = 100)
Psoriasis09[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis09, pattern = "^MT-")
Psoriasis09 <- subset(Psoriasis09, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis09@meta.data$stim <- as.character(phenotype[grep('Psoriasis09', phenotype$ID),2])
Psoriasis09@meta.data$number <- as.character(phenotype[grep('Psoriasis09', phenotype$ID),1])
rm(Psoriasis09_morethan25percentMTbc)

Psoriasis09 <- NormalizeData(Psoriasis09)
Psoriasis09 <- FindVariableFeatures(Psoriasis09, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis09)
Psoriasis09 <- ScaleData(Psoriasis09, features = all.genes)

rm(all.genes)
unlink("./Psoriasis09", recursive = TRUE)
######################################################################################################
## Load data Psoriasis10 ####
files <- Sys.glob('Psoriasis10*')
dir.create("./Psoriasis10")
for(file in files) {file.copy(file, "./Psoriasis10")}
new_files <- gsub("Psoriasis10_","",files)
setwd("./Psoriasis10")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis10 <- Read10X(data.dir = "./Psoriasis10")

Psoriasis10 <- CreateSeuratObject(counts = Psoriasis10, project = "Psoriasis10", min.cells = 3, min.features = 100)
Psoriasis10[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis10, pattern = "^MT-")
Psoriasis10@meta.data$stim <- as.character(phenotype[grep('Psoriasis10', phenotype$ID),2])
Psoriasis10@meta.data$number <- as.character(phenotype[grep('Psoriasis10', phenotype$ID),1])

Psoriasis10_morethan25percentMT <- subset(Psoriasis10, subset = percent.mt > 25)
Psoriasis10_morethan25percentMTbc <- colnames(Psoriasis10_morethan25percentMT@assays$RNA@data)
rm(Psoriasis10)
rm(Psoriasis10_morethan25percentMT)

Psoriasis10 <- Read10X(data.dir = "./Psoriasis10")
Psoriasis10<-remove.genes(Psoriasis10, Psoriasis10_morethan25percentMTbc)

Psoriasis10 <- CreateSeuratObject(counts = Psoriasis10, project = "Psoriasis10", min.cells = 3, min.features = 100)
Psoriasis10[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis10, pattern = "^MT-")
Psoriasis10 <- subset(Psoriasis10, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis10@meta.data$stim <- as.character(phenotype[grep('Psoriasis10', phenotype$ID),2])
Psoriasis10@meta.data$number <- as.character(phenotype[grep('Psoriasis10', phenotype$ID),1])
rm(Psoriasis10_morethan25percentMTbc)

Psoriasis10 <- NormalizeData(Psoriasis10)
Psoriasis10 <- FindVariableFeatures(Psoriasis10, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis10)
Psoriasis10 <- ScaleData(Psoriasis10, features = all.genes)

rm(all.genes)
unlink("./Psoriasis10", recursive = TRUE)
######################################################################################################
## Load data Psoriasis11 ####
files <- Sys.glob('Psoriasis11*')
dir.create("./Psoriasis11")
for(file in files) {file.copy(file, "./Psoriasis11")}
new_files <- gsub("Psoriasis11_","",files)
setwd("./Psoriasis11")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis11 <- Read10X(data.dir = "./Psoriasis11")

Psoriasis11 <- CreateSeuratObject(counts = Psoriasis11, project = "Psoriasis11", min.cells = 3, min.features = 100)
Psoriasis11[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis11, pattern = "^MT-")
Psoriasis11@meta.data$stim <- as.character(phenotype[grep('Psoriasis11', phenotype$ID),2])
Psoriasis11@meta.data$number <- as.character(phenotype[grep('Psoriasis11', phenotype$ID),1])

Psoriasis11_morethan25percentMT <- subset(Psoriasis11, subset = percent.mt > 25)
Psoriasis11_morethan25percentMTbc <- colnames(Psoriasis11_morethan25percentMT@assays$RNA@data)
rm(Psoriasis11)
rm(Psoriasis11_morethan25percentMT)

Psoriasis11 <- Read10X(data.dir = "./Psoriasis11")
Psoriasis11<-remove.genes(Psoriasis11, Psoriasis11_morethan25percentMTbc)

Psoriasis11 <- CreateSeuratObject(counts = Psoriasis11, project = "Psoriasis11", min.cells = 3, min.features = 100)
Psoriasis11[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis11, pattern = "^MT-")
Psoriasis11 <- subset(Psoriasis11, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis11@meta.data$stim <- as.character(phenotype[grep('Psoriasis11', phenotype$ID),2])
Psoriasis11@meta.data$number <- as.character(phenotype[grep('Psoriasis11', phenotype$ID),1])
rm(Psoriasis11_morethan25percentMTbc)

Psoriasis11 <- NormalizeData(Psoriasis11)
Psoriasis11 <- FindVariableFeatures(Psoriasis11, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis11)
Psoriasis11 <- ScaleData(Psoriasis11, features = all.genes)

rm(all.genes)
unlink("./Psoriasis11", recursive = TRUE)
######################################################################################################
## Load data Psoriasis12 ####
files <- Sys.glob('Psoriasis12*')
dir.create("./Psoriasis12")
for(file in files) {file.copy(file, "./Psoriasis12")}
new_files <- gsub("Psoriasis12_","",files)
setwd("./Psoriasis12")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis12 <- Read10X(data.dir = "./Psoriasis12")

Psoriasis12 <- CreateSeuratObject(counts = Psoriasis12, project = "Psoriasis12", min.cells = 3, min.features = 100)
Psoriasis12[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis12, pattern = "^MT-")
Psoriasis12@meta.data$stim <- as.character(phenotype[grep('Psoriasis12', phenotype$ID),2])
Psoriasis12@meta.data$number <- as.character(phenotype[grep('Psoriasis12', phenotype$ID),1])

Psoriasis12_morethan25percentMT <- subset(Psoriasis12, subset = percent.mt > 25)
Psoriasis12_morethan25percentMTbc <- colnames(Psoriasis12_morethan25percentMT@assays$RNA@data)
rm(Psoriasis12)
rm(Psoriasis12_morethan25percentMT)

Psoriasis12 <- Read10X(data.dir = "./Psoriasis12")
Psoriasis12<-remove.genes(Psoriasis12, Psoriasis12_morethan25percentMTbc)

Psoriasis12 <- CreateSeuratObject(counts = Psoriasis12, project = "Psoriasis12", min.cells = 3, min.features = 100)
Psoriasis12[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis12, pattern = "^MT-")
Psoriasis12 <- subset(Psoriasis12, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis12@meta.data$stim <- as.character(phenotype[grep('Psoriasis12', phenotype$ID),2])
Psoriasis12@meta.data$number <- as.character(phenotype[grep('Psoriasis12', phenotype$ID),1])
rm(Psoriasis12_morethan25percentMTbc)

Psoriasis12 <- NormalizeData(Psoriasis12)
Psoriasis12 <- FindVariableFeatures(Psoriasis12, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis12)
Psoriasis12 <- ScaleData(Psoriasis12, features = all.genes)

rm(all.genes)
unlink("./Psoriasis12", recursive = TRUE)
######################################################################################################
## Load data Psoriasis13 ####
files <- Sys.glob('Psoriasis13*')
dir.create("./Psoriasis13")
for(file in files) {file.copy(file, "./Psoriasis13")}
new_files <- gsub("Psoriasis13_","",files)
setwd("./Psoriasis13")
file.copy(from = files, to = new_files)
file.remove(files)
setwd("../")
Psoriasis13 <- Read10X(data.dir = "./Psoriasis13")

Psoriasis13 <- CreateSeuratObject(counts = Psoriasis13, project = "Psoriasis13", min.cells = 3, min.features = 100)
Psoriasis13[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis13, pattern = "^MT-")
Psoriasis13@meta.data$stim <- as.character(phenotype[grep('Psoriasis13', phenotype$ID),2])
Psoriasis13@meta.data$number <- as.character(phenotype[grep('Psoriasis13', phenotype$ID),1])

Psoriasis13_morethan25percentMT <- subset(Psoriasis13, subset = percent.mt > 25)
Psoriasis13_morethan25percentMTbc <- colnames(Psoriasis13_morethan25percentMT@assays$RNA@data)
rm(Psoriasis13)
rm(Psoriasis13_morethan25percentMT)

Psoriasis13 <- Read10X(data.dir = "./Psoriasis13")
Psoriasis13<-remove.genes(Psoriasis13, Psoriasis13_morethan25percentMTbc)

Psoriasis13 <- CreateSeuratObject(counts = Psoriasis13, project = "Psoriasis13", min.cells = 3, min.features = 100)
Psoriasis13[["percent.mt"]] <- PercentageFeatureSet(object = Psoriasis13, pattern = "^MT-")
Psoriasis13 <- subset(Psoriasis13, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Psoriasis13@meta.data$stim <- as.character(phenotype[grep('Psoriasis13', phenotype$ID),2])
Psoriasis13@meta.data$number <- as.character(phenotype[grep('Psoriasis13', phenotype$ID),1])
rm(Psoriasis13_morethan25percentMTbc)

Psoriasis13 <- NormalizeData(Psoriasis13)
Psoriasis13 <- FindVariableFeatures(Psoriasis13, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Psoriasis13)
Psoriasis13 <- ScaleData(Psoriasis13, features = all.genes)

rm(all.genes)
unlink("./Psoriasis13", recursive = TRUE)
######################################################################################################
