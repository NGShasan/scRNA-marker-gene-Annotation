
## Load libraries
library(SingleCellExperiment)
library(DESeq2)
library(Biobase)
library(BisqueRNA)
library(dplyr)
library(Seurat)
library(patchwork)
library(convert)
library(Matrix)
library(scater)
library(cellranger)

## check the current directory
getwd()

## set the working directory
setwd("C://WIP/Sample_Datasets")  # windows directory

#setwd("/home/****/Documents/WIP") # linux directory

#### upload scRNA RData ######

load("SRS380524.rdata")
ls()
sm
class(sm) #Matrix
colnames(sm)
rownames(sm)

# Initialize the Seurat object with the raw (non-normalized data).
s_obj <- CreateSeuratObject(counts = sm, project = "bisque", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(s_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(s_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(s_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Normalizing the data ###
sc <- NormalizeData(s_obj, normalization.method = "LogNormalize", scale.factor = 10000)

## STEP: Feature selection ##
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #error: When using repel, set xnudge and ynudge to 0 for optimal results
plot1 + plot2

##STEP:  Scaling the data ##
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)

### Perform linear dimensional reduction ##
sc <- RunPCA(sc, features = VariableFeatures(object = sc))

# Examine and visualize PCA results a few different ways
print(sc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(sc, dims = 1:2, reduction = "pca")

DimPlot(sc, reduction = "pca")

DimHeatmap(sc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc, dims = 1:15, cells = 500, balanced = TRUE)


# The JackStrawPlot() function provides a visualization tool for comparing the distribution
#of p-values for each PC with a uniform distribution (dashed line). 'Significant' PCs will 
#show a strong enrichment of features with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.


## Determine the 'dimensionality' of the dataset##
sc <- JackStraw(sc, num.replicate = 100)
sc <- ScoreJackStraw(sc, dims = 1:20)

JackStrawPlot(sc, dims = 1:15)

# make an elbowplot to find the right k for clustering
ElbowPlot(sc)

##  STEP: Cluster the cells ##
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(sc), 5)


### STEP: Run non-linear dimensional reduction (UMAP/tSNE) ###

# If you haven't installed UMAP, you can do so via 
#reticulate::py_install(packages ='umap-learn')

sc <- RunUMAP(sc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(sc, reduction = "umap")

##  save the object at this point
saveRDS(sc, file = "C:/WIP/Sample/output_sc.rds")


## remove 5s_RNA from seurat object (1)
sc
head(rownames(sc))

# remove 1st, 2nd and 3rd row to remove 5s-rRNA
sc <- sc[-c(1:3),] 

#lets checked 5s-rRNA removed or not
head(rownames(sc))

## create SingleCellExperiment
sce <- as.SingleCellExperiment(sc)
class(sce)
sce

#### split the gene names "-"

sce_split=strsplit(rownames(sce), "-")
sce_split

# sce.1 <- lapply(sce_split, function(x) {x[c(1)]}) ## see the gene names only
sc.1<- lapply(sce_split, function(x) {x[c(1)]}) ## see the gene names only


## unlist the list 
sc.2= print(unlist(sc.1, use.names=FALSE)) ## unlist the genes

class(sc.2)
length(sc.2)

## set the new rownames

rownames(sce) = sc.2
rownames(sce)

#### Marker Genes Annotations

library(SingleR)
library(celldex)

#ref <- HumanPrimaryCellAtlasData() ## This reference did not work well, 

#so lets chek MonacoImmuneData()

ref <- MonacoImmuneData()
ref

## Lets Annotate the cell markers
pred.sce <- SingleR(test = sce, ref = ref, assay.type.test=1,
                    labels = ref$label.main)


# Summarizing the distribution:
table(pred.sce$labels)

## Totals cells in the normalized matrix
ncol(sce)

## Lets assign cells with 

sc$Prediction <- pred.sce$labels
sc$Prediction

## Lets Plot the cluster UMAP with Marker genes
DimPlot(sc, reduction = "umap", group.by = "Prediction", label = T) # group.by will be exact name in sc$name


##  save the object at this point
saveRDS(sc, file = "C:/WIP/Sample_Datasets/output_sc.rds")


sessionInfo()
