# Libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RCurl)
library(AnnotationHub)
library(ensembldb)


### IMPORT ###

# Load Data and Convert to Seurat Object
glio.raw <- Read10X_h5(filename = 'Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.h5')
glio.seurat <- CreateSeuratObject(counts = glio.raw, project = "Glioblastoma", min.cells = 3, min.features = 200)


# Add Mitochondrial Info to Metadata
glio.seurat[["percent.mito"]] <- PercentageFeatureSet(glio.seurat, pattern = "^MT-")


### PREPROCESSING ###

# Cellular Filter - Mitochondrial %, nFeature, nCount
glio.seurat <- subset(glio.seurat, subset = 
                        nFeature_RNA > 200 & 
                        nFeature_RNA < 7500 & 
                        nCount_RNA > 50 & 
                        nCount_RNA < 90000 & 
                        percent.mito < 20)


# Genetic Filter - Remove Genes with Zero Counts and Present in <0.1% of Cells
glio.counts <- GetAssayData(object = glio.seurat, slot = "counts")
glio.nonzero <- glio.counts > 0
glio.keep.genes <- Matrix::rowSums(glio.nonzero) >= 0.001 * ncol(glio.counts)
glio.filtered.counts <- glio.counts[glio.keep.genes, ]
glio.seurat <- CreateSeuratObject(counts = glio.filtered.counts, project = "Glioblastoma_Filtered", meta.data = glio.seurat@meta.data)


# Normalization
glio.seurat <- NormalizeData(glio.seurat)


# Variable Feature ID
glio.seurat <- FindVariableFeatures(glio.seurat)


# Scale Data
genes.all <- rownames(glio.seurat)
glio.seurat <- ScaleData(glio.seurat, features = genes.all)


# LDR (PCA)
glio.seurat <- RunPCA(glio.seurat, features = VariableFeatures(object = glio.seurat))


### CELL CYCLE SCORING ###

# Load TinyAtlas Gene List
cc.genes <- read.csv(text = getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv"))

# Retrieve Annotations from AnnotationHub
ah <- AnnotationHub()

ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

edb <- ah[[id]]

annotations <- genes(edb, 
                     return.type = "data.frame")

annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


# Annotate Cell Cycle Genes
cc.markers <- dplyr::left_join(cc.genes, annotations, by = c("geneID" = "gene_id"))

s.genes <- cc.markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

g2m.genes <- cc.markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


# Score for Cell Cycle Phase
glio.seurat.cc <- CellCycleScoring(glio.seurat,
                                   g2m.features = g2m.genes,
                                   s.features = s.genes)

glio.seurat.cc <- RunPCA(glio.seurat.cc)

# Visualize Cells by Phase - See if Regression is Necessary
DimPlot(glio.seurat.cc,
        reduction = "pca",
        group.by= "Phase")
