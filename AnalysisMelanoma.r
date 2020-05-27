#remotes::install_github("rnabioco/clustifyr")
library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(usethis)

mat_humanMelanomaDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_human_melanoma_counts_normalized_9315x19445.tsv.gz") %>%
  select(-X1) %>% t()
colnames(mat_humanMelanomaDC) <- seq_len(ncol(mat_humanMelanomaDC))
meta_humanMelanomaDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_human_melanoma_cell_metadata_9315x14.tsv.gz") %>%
  mutate("rn" = row_number())
ref_humanMelanomaDC <- average_clusters(
  mat_humanMelanomaDC,
  meta_humanMelanomaDC,
  "cell_type"
)

usethis::use_data(ref_humanMelanomaDC, compress = "xz", overwrite = TRUE)

#Preprocessing workflow
melanoma <- CreateSeuratObject(counts = mat_humanMelanomaDC, project = "Melanoma", min.cells = 3, min.features = 200) 
melanoma
melanoma@assays$RNA@data <- melanoma@assays$RNA@counts
melanoma[["percent.mt"]] <- PercentageFeatureSet(melanoma, pattern = "^MT")
head(melanoma@meta.data, 5)
VlnPlot(melanoma, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(melanoma, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(melanoma, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
melanoma <- subset(melanoma, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Variable Features
melanoma <- FindVariableFeatures(melanoma, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(melanoma), 10)
plot1 <- VariableFeaturePlot(melanoma)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear dimension reduction
all.genes <- rownames(melanoma)
melanoma <- ScaleData(melanoma, features = all.genes)
melanoma <- RunPCA(melanoma, features = VariableFeatures(object = melanoma))
print(melanoma[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(melanoma, dims = 1:2, reduction = "pca")
DimPlot(melanoma, reduction = "pca")
DimHeatmap(melanoma, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(melanoma, dims = 1:15, cells = 500, balanced = TRUE)

#Dimensionality
melanoma <- JackStraw(melanoma, num.replicate = 100)
melanoma <- ScoreJackStraw(melanoma, dims = 1:20)
JackStrawPlot(melanoma, dims = 1:15)
ElbowPlot(melanoma)

#Cluster cells
melanoma <- FindNeighbors(melanoma, dims = 1:10)
melanoma <- FindClusters(melanoma, resolution = 0.5)
head(Idents(melanoma), 5)

#Non-linear dimensional reduction (UMAP/tSNE)
melanoma <- RunUMAP(melanoma, dims = 1:10)
DimPlot(melanoma, reduction = "umap")

cluster1.markers <- FindMarkers(melanoma, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(melanoma, ident.1 = 5, ident.2 = 0, min.pct = 0.25)
head(cluster5.markers, n = 5)
melanoma.markers <- FindAllMarkers(melanoma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
melanoma.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(melanoma, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(melanoma, features = c("SPP1", "TIMP1"))
VlnPlot(melanoma, features = c("SPP1", "TIMP1"), slot = "counts", log = TRUE)
FeaturePlot(melanoma, features = c("SPP1", "TIMP1", "C1QA", "APOC1", "TYROBP"))
top10 <- melanoma.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(melanoma, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters
new.cluster.ids <- c(meta_humanMelanomaDC)
names(new.cluster.ids) <- levels(melanoma)
melanoma <- RenameIdents(melanoma, new.cluster.ids)
DimPlot(melanoma, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(Idents(melanoma), 5)

