# Libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(pheatmap)


### IMPORT ###

# Load Data and Convert to Seurat Object
glio.raw <- Read10X_h5(filename = 'Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.h5')
glio.seurat <- CreateSeuratObject(counts = glio.raw, project = "Glioblastoma", min.cells = 3, min.features = 200)


# Add Mitochondrial Info to Metadata
glio.seurat[["percent.mito"]] <- PercentageFeatureSet(glio.seurat, pattern = "^MT-")


### QC ###

# Count Cell Number and Compare to Loaded Cell Number - Extra: Junk Cells Present, Fewer: Low Capture Efficiency
print(glio.seurat)


# Visualize Gene Count, UMI Count + Mitochondrial Percentage Distributions
VlnPlot(glio.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


# Visualize Number of Genes vs Number of UMIs + Mitochondrial Percentage
glio.seurat@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)


# Visualize Gene Complexity Distribution
glio.seurat$log10GenesPerUMI <- log10(glio.seurat$nFeature_RNA) / log10(glio.seurat$nCount_RNA)

glio.seurat@meta.data %>%
  ggplot(aes(x = log10GenesPerUMI)) +
  geom_density(alpha = 0.2, color = "springgreen4", fill = "springgreen4") +
  theme_classic() +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "gray10")


### FILTER + PREPROCESSING ###

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


# Check Cell Count After Filtering
print(glio.seurat)


# Normalization
glio.seurat <- NormalizeData(glio.seurat)


# Variable Feature ID
glio.seurat <- FindVariableFeatures(glio.seurat)


# Store Top 20 Most Variable Genes
t20.var.feat <- head(VariableFeatures(glio.seurat), 20)


# Plot of Top 20 Variable Features
plot.t20 <- VariableFeaturePlot(glio.seurat)
LabelPoints(plot = plot.t20, points = t20.var.feat, repel = TRUE, xnudge = 0, ynudge = 0)


# Scale Data
genes.all <- rownames(glio.seurat)
glio.seurat <- ScaleData(glio.seurat, features = genes.all)


# LDR (PCA) + Visualization
glio.seurat <- RunPCA(glio.seurat, features = VariableFeatures(object = glio.seurat))

print(glio.seurat[["pca"]], dims = 1:5, nfeatures = 10)
DimHeatmap(glio.seurat, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(glio.seurat)


### CLUSTERING ###

# Cluster
glio.seurat <- FindNeighbors(glio.seurat, dims = 1:20)


# Set Cluster Resolution
glio.seurat <- FindClusters(glio.seurat, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
DimPlot(glio.seurat, group.by = "RNA_snn_res.0.1", label = TRUE) # Lower Res = Fewer Clusters, Higher Res = More Clusters. Test and adjust


# Set Cluster Identity to Correct Resolution
Idents(glio.seurat) <- "RNA_snn_res.0.1"
Idents(glio.seurat)


# NL-LDR (UMAP) + Visualization
glio.seurat <- RunUMAP(glio.seurat, dims = 1:20)
DimPlot(glio.seurat, reduction = "umap", label = TRUE)


### DOUBLET FILTERING ###

# Find Optimal pK Value
sweep.res.list_glio <- paramSweep(glio.seurat, PCs = 1:20, sct = FALSE)
sweep.stats_glio <- summarizeSweep(sweep.res.list_glio, GT = FALSE)
bcmvn_glio <- find.pK(sweep.stats_glio)


# Plot pK Results to Locate pkMax
ggplot(bcmvn_glio, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()


# Automatically Store pK Value of pkMax Based on Result
pK <- bcmvn_glio %>%
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Estimate Homotypic Doublet Proportion
glio.annotations <- glio.seurat@meta.data$seurat_clusters
glio.homotypic.prop <- modelHomotypic(glio.annotations)

nExp_poi <- round(0.076*nrow(glio.seurat@meta.data)) # Default, as recovered and loaded information not provided for dataset
nExp_poi.adj <- round(nExp_poi*(1-glio.homotypic.prop))


# Print Doublet Estimates
print(nExp_poi) # Estimate of real doublets in the data
print(nExp_poi.adj) # Adjusted estimate of doublets after homotypic adjustment


# Run DoubletFinder
glio.seurat <- doubletFinder(glio.seurat, 
                             PCs = 1:20, 
                             pN = 0.25, 
                             pK = pK, 
                             nExp = nExp_poi.adj,
                             reuse.pANN = NULL, 
                             sct = FALSE)


# DoubletFinder Visualization
df_column <- grep("DF.classifications", colnames(glio.seurat@meta.data), value = TRUE)
print(df_column) # Use this to replace the "DF.classifications" string with the correct column name manually below
DimPlot(glio.seurat, reduction = 'umap', group.by = df_column)


# Show Singlet and Doublet Statistics
table(glio.seurat@meta.data$"DF.classifications_0.25_0.005_390")


# Doublet Removal
glio.seurat.filtered <- subset(glio.seurat, subset = DF.classifications_0.25_0.005_390 == "Singlet")

# Count Cells After Doublet Removal and Filtering and Compare to Initial Count
table(glio.seurat.filtered@meta.data$"DF.classifications_0.25_0.005_390")


# Visualize Filtered Clusters
DimPlot(glio.seurat.filtered, reduction = "umap", label = TRUE)


# Retrieve Final Cell Number Per Cluster
glio.n.cells <- FetchData(glio.seurat.filtered, 
           vars = c("ident", "orig.ident")) %>%
           dplyr::count(ident, orig.ident) %>%
           tidyr::spread(ident, n)

View(glio.n.cells)


# Visualize Cell Metrics by Cluster
glio.cluster.metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mito")

FeaturePlot(glio.seurat.filtered, 
            reduction = "umap", 
            features = glio.cluster.metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Plot of Top 20 Variable Features after Doublet Finder
plot.t20 <- VariableFeaturePlot(glio.seurat.filtered)
LabelPoints(plot = plot.t20, points = t20.var.feat, repel = TRUE, xnudge = 0, ynudge = 0)

### DEG ANALYSIS ###

# Find All Differentially Expressed Genes per Cluster
DefaultAssay(glio.seurat.filtered) <- 'RNA'

glio.markers.all <- FindAllMarkers(glio.seurat.filtered,
                                   logfc.threshold = 0.25,
                                   min.pct = 0.1,
                                   only.pos = TRUE)


# Store Top DEGs per Cluster + Export Full List
glio.markers.t5 <- glio.markers.all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(glio.markers.all, file = "Glioblastoma_DEG_Results.csv")


# Visualize Top DEGs per Cluster
FeaturePlot(glio.seurat.filtered, features = unique(glio.markers.t5$gene), ncol = 5)
DoHeatmap(glio.seurat.filtered, features = unique(glio.markers.t5$gene)) + NoLegend()
VlnPlot(glio.seurat.filtered, features = unique(glio.markers.t5$gene), ncol = 5)


### CELL ANNOTATION ###

# Obtain Reference Dataset
ref <- celldex::HumanPrimaryCellAtlasData()


# Run SingleR
sr <- SingleR(test = glio.counts, 
              ref = ref, 
              labels = ref$label.main) # Can run label.fine for more specific cell binning

sr


# Visualize SingleR Results
glio.seurat.filtered$singleR.labels <- sr$labels[match(rownames(glio.seurat.filtered@meta.data), rownames(sr))]
DimPlot(glio.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')


# Classification QC Check
plotScoreHeatmap(sr)
plotDeltaDistribution(sr)
