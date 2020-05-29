library(dplyr)
library(Seurat)
library(patchwork)

#Create Seurat objects
retina.data <- Read10X(data.dir = "~/scRNA-seq-Practice/data_Retina/") #Pull data from the barcode tsv
retina <- CreateSeuratObject(counts = retina.data, project = "Retina", min.cells = 3, min.features = 200) #Create Seurat objects
retina #Display active assays 
retina.data[c("FAM138A", "SAMD11", "HES4"), 1:30] #Display top three genes in the first thirthy cells

#Preprocessing work flow
retina[["percent.mt"]] <- PercentageFeatureSet(retina, pattern = "^MT")
head(retina@meta.data, 5)
VlnPlot(retina, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(retina, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(retina, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
retina <- subset(retina, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

retina <- NormalizeData(retina, normalization.method = "LogNormalize", scale.factor = 10000)

retina <- FindVariableFeatures(retina, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(retina), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(retina)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Perform linear dimension reduction
all.genes <- rownames(retina)
retina <- ScaleData(retina, features = all.genes) #This step takes forever
retina <- RunPCA(retina, features = VariableFeatures(object = retina))
print(retina[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(retina, dims = 1:2, reduction = "pca")
DimPlot(retina, reduction = "pca")
DimHeatmap(retina, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(retina, dims = 1:15, cells = 500, balanced = TRUE)

#Determine dimensionality of the data set
retina <- JackStraw(retina, num.replicate = 100)
retina <- ScoreJackStraw(retina, dims = 1:20)
JackStrawPlot(retina, dims = 1:15)
ElbowPlot(retina)

#Cluster cells
retina <- FindNeighbors(retina, dims = 1:10)
retina <- FindClusters(retina, resolution = 0.5)
head(Idents(retina), 5)

#Non-linear dimensional reduction (UMAP/tSNE)
retina <- RunUMAP(retina, dims = 1:10)
DimPlot(retina, reduction = "umap")

cluster1.markers <- FinaMarkers(retina, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(retina, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n = 5)
retina.markers <- FindAllMarkers(retina, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
retina.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(retina, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(retina, features = c("SOX11", "CADPS"))
VlnPlot(retina, features = c("SOX11", "CADPS"), slot = "counts", log = TRUE)
FeaturePlot(retina, features = c("WIF1", "TF", "MEIS2", "TFAP2A", "GNGT1", "RCVRN", "PRDM13", "SOX11", "CADPS", "CA10"))
top10 <- retina.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(retina, features = top10$gene) + NoLegend()
