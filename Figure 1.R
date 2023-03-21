library(Seurat)
library(plyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(Matrix)
library(data.table)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(writexl)

#################### PRE-PROCESSING #################### 

#Read scRNAseq expression data

MISC1.data <- Read10X(data.dir = "/mnt/RAID_5/vincent/Project-Kawasaki_Disease/2021_DatafromGEO_GSE_GSE166489-MISC-S/S1(P1.1)/")
rownames(x = MISC1.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", 
                                                      x = rownames(x = MISC1.data[["Antibody Capture"]]))
                                                      
# Create Seurat object. 
MISC1 <- CreateSeuratObject(counts = MISC1.data[["Gene Expression"]], project = "MIS-C1", min.cells = 3, min.features = 200)

#Standard pre-processing workflow
MISC1[["percent.mt"]] <- PercentageFeatureSet(MISC1, pattern = "^MT-") #NKD <- NKD[!grepl("^MT-", rownames(NKD)), ]
MISC1[["percent.hb"]] <- PercentageFeatureSet(MISC1, pattern = "^HB") # For 10% hemoglobin: https://insight.jci.org/articles/view/135678#SEC4

# Visualize QC metrics as a violin plot
VlnPlot(MISC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)
MISC1 <- subset(MISC1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000  & percent.mt < 15 & percent.hb < 5)

#Removal of MALAT1, XIST, RSP4Y1
MISC1 <- MISC1[!grepl("^MT-", rownames(MISC1)), ] #MT genes are already removed but the column of its percent.mt values remain
MISC1 <- MISC1[!grepl("MALAT1", rownames(MISC1)), ]
MISC1 <- MISC1[!grepl("XIST", rownames(MISC1)), ]
MISC1 <- MISC1[!grepl("RSP4Y1", rownames(MISC1)), ]
MISC1 <- MISC1[grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(MISC1), value = T, invert = T), ]
rm(MISC1.data)



















