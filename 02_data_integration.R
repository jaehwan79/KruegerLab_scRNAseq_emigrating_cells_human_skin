#load library
#instRound0.packages('Seurat')
library(Seurat)
library(RNAransform)
library(dplyr)
library(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(GEOquery)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DropletUtils")
library(DropletUtils)
#instRound0.packages("devtools")
library(devtools)
#source("https://raw.githubusercontent.com/farrellja/URD/master/URD-InstRound0.R")
library(URD)
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(reticulate)
library(scales)
library(forcats)
library(cowplot)

######################################################################################################
#Data Integration
# We used Seurat R package (version 3.0, https://satijalab.org/seurat/) for the single-cell data integration. 
# We first merged single-cell data of samples with the common disease status (psoriasis or control) 
# and the identical reagent kit version (V2.0 or V3.0) into three groups: 
#(1) psoriasis samples with reagent kit V2.0 - 10 samples, \
#(2) psoriasis samples with reagent kit V3.0 - 3 samples, and 
# (3) control samples with reagent kit V2.0 - 5 samples. 

Psoriasis_V2.0 <- merge(x = Psoriasis01, y = c(Psoriasis02, Psoriasis03, Psoriasis04,Psoriasis05,Psoriasis06,Psoriasis07,Psoriasis08,
                                                   Psoriasis09, Psoriasis10
), 
add.cell.ids = c("psoriasis01", "psoriasis02", "psoriasis03", "psoriasis04", "psoriasis05","psoriasis06","psoriasis07","psoriasis08",
                 "psoriasis09","psoriasis10"
), 
project = "Psoriasis_V2.0")

Control_V2.0 <- merge(x = Control01, y = c(
  Control02, Control03, Control04,Control05), 
  add.cell.ids = c(
    "cntl01", "cntl02", "cntl03", "cntl04", "cntl05"), 
  project = "Control_V2.0")

Psoriasis_V3.0 <- merge(x= Psoriasis11, y = c(Psoriasis12,Psoriasis13), add.cell.ids=c("psoriasis11","psoriasis12","psoriasis13"), project = "Psoriasis_V3.0")

######################################################################################################
### Data normalization
Psoriasis_V2.0 <- NormalizeData(Psoriasis_V2.0)
Psoriasis_V2.0 <- FindVariableFeatures(Psoriasis_V2.0, selection.method = "vst", nfeatures = 2000)

Psoriasis_V3.0 <- NormalizeData(Psoriasis_V3.0)
Psoriasis_V3.0 <- FindVariableFeatures(Psoriasis_V3.0, selection.method = "vst", nfeatures = 2000)

Control_V2.0 <- NormalizeData(Control_V2.0)
Control_V2.0 <- FindVariableFeatures(Control_V2.0, selection.method = "vst", nfeatures = 2000)

######################################################################################################
### Data inegration
# To harmonize merged groups into a single dataset without batch effects, 
# correspondences between cells in three merged datasets were identified by the FIndIntegrationAnchors function, 
# and used for data integration with the IntegratedData function as detailed by Butler et al. 
immune.anchors <- FindIntegrationAnchors(object.list = list(Psoriasis_V2.0, Psoriasis_V3.0,Control_V2.0), dims = 1:30)
Round0 <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

######################################################################################################
### Perform an integrated analysis
Round0@assays
DefaultAssay(Round0) <- "integrated"

## Vnplot 
pdf('Vnplot_Round0.pdf',width=20,height=5)
VlnPlot(object = Round0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
dev.off()

#FeatureScatter
plot1 <- FeatureScatter(object = Round0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Round0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf('FeatureScatter_Round0.pdf',width=14,height=5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

save(Round0, file = "Round0_integrated.Rda")
