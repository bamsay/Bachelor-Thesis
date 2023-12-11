
library(devtools)
devtools::install_github("IMSBCompBio/SpaCo")

### Import Libraries
library(SPACO)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(VennDiagram)
library(tidyverse)
library(ggforce)
library(ggVennDiagram)


### PART ONE
### install anterior1 data from SeuratData
### Treat the anterior1 as a spatial data

InstallData("stxBrain")
brain <- LoadData("stxBrain",  type = "anterior1")
brain
brain <- PercentageFeatureSet(brain, pattern = "^mt-" ,col.name = "percent.mt")
brain <- PercentageFeatureSet(brain, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain <- SCTransform(brain, assay = "Spatial", variable.features.n = 3000, verbose = F)

### Perform downstream processing to obtain the clusters

brainPCA <- RunPCA(brain, assay = "SCT", verbose = FALSE)

brainPCA <- FindNeighbors(brainPCA, reduction = "pca", dims = 1:30, verbose = FALSE)

brainPCA <- FindClusters(brainPCA, resolution = 0.5, verbose = FALSE)

brainPCA <- RunUMAP(brainPCA, reduction = "pca", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(brainPCA, reduction = "umap",  
              label = TRUE, label.size = 5)  + ggtitle('Clusters of Anterior 1 mouse brain spatial dataset')

p1


### find markers for cell annotation
brain_markers <- FindAllMarkers(brainPCA, logfc.threshold = 0.25,
                                only.pos = TRUE, 
                                min.pct = 0.1)
head(brain_markers)

### Compare the marker genes of each cluster to 
### tabula muris and panglao db databases, Identify cells in each cluster.
### Annotate the clusters with the cell names

brainPCA <- RenameIdents(brainPCA, '0' = 'ca1-Neurons'). ### This was done for all clusters

p2 <- DimPlot(brainPCA, reduction = "umap",  
              label = TRUE, label.size = 5)  + ggtitle('Annotated Clusters of Anterior 1 mouse brain spatial dataset')

p2



### Create SPACO object using the default tissue KNN graph

SpaCoObject<- seurat_to_spaco(Seurat = brain, assay = "SCT", n_image= 1, slot = "scale.data")

SpaCoObject <- RunSCA(SpaCoObject, compute_nSpacs = TRUE )

SpaCoObject@nSpacs

#### Compute spatial variable genes (SVG’s).

DE_genes<- SVGTest(SpaCoObject)

head(DE_genes)

DE_genes_sortbrain <- DE_genes[order(DE_genes$score, decreasing = TRUE),]

sigs_brain <- rownames(DE_genes[DE_genes$p.adjust<0.05,])


### Now treat the anterior 1 as a single cell data

brain2 <- LoadData("stxBrain",  type = "anterior1")

brain2 <- PercentageFeatureSet(brain2, pattern = "^mt-" ,col.name = "percent.mt")
brain2 <- PercentageFeatureSet(brain2, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain2 <- NormalizeData(brain2)
brain2 <- FindVariableFeatures(brain2, selection.method = "vst", nfeatures = 3000)
brain2 <- SCTransform(brain2, assay = "Spatial",variable.features.n = 3000, verbose = F)
brain2

### Downstream Processing, dimension reduction, knn graph construction
brain2 <- RunPCA(brain2, features = VariableFeatures(object = brain2), verbose = F)
brain2 <- RunUMAP(brain2,reduction = "pca",dims = 1:30, verbose = F)
brain2 <- FindNeighbors(brain2,reduction = "pca", k.param = 6, dims = 1:30, verbose = F, return.neighbor = T)

### Create a SPACO object using the constructed KNN graph

SpaCoObject2<- create_SpaCoObject_from_KNN(brain2, n = 6)

SpaCoObject2   <- RunSCA(SpaCoObject2,compute_nSpacs = TRUE)

SpaCoObject2@nSpacs

#### Compute spatial variable genes (SVG’s)

DE_genes2<- SVGTest(SpaCoObject2)

head(DE_genes2)

DE_genes_sortbrain2 <- DE_genes2[order(DE_genes2$score, decreasing = TRUE),]

sigs_brain2 <- rownames(DE_genes2[DE_genes2$p.adjust<0.05,])



### Part 2
### install anterior2 data from SeuratData
### Treat the anterior2 as a spatial data

brain3 <- LoadData("stxBrain",  type = "anterior2")
brain3
brain3 <- PercentageFeatureSet(brain3, pattern = "^mt-" ,col.name = "percent.mt")
brain3 <- PercentageFeatureSet(brain3, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain3 <- SCTransform(brain3, assay = "Spatial", variable.features.n = 3000, verbose = F)


### Create SPACO object using default KNN graph

SpaCoObject3<- seurat_to_spaco(Seurat = brain3, assay = "SCT", n_image= 1, slot = "scale.data")

SpaCoObject3 <- RunSCA(SpaCoObject3, compute_nSpacs = TRUE )

SpaCoObject3@nSpacs

### Compute spatial variable genes (SVG’s).

DE_genes3<- SVGTest(SpaCoObject3)

head(DE_genes3)

DE_genes_sortbrain3 <- DE_genes3[order(DE_genes3$score, decreasing = TRUE),]

sigs_brain3 <- rownames(DE_genes3[DE_genes3$p.adjust<0.05,])


### Now treat the anterior 2 as a single cell data

brain4 <- LoadData("stxBrain",  type = "anterior2")

brain4 <- PercentageFeatureSet(brain4, pattern = "^mt-" ,col.name = "percent.mt")
brain4 <- PercentageFeatureSet(brain4, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain4 <- NormalizeData(brain4)
brain4 <- FindVariableFeatures(brain4, selection.method = "vst", nfeatures = 3000)
brain4 <- SCTransform(brain4, assay = "Spatial", variable.features.n = 3000, verbose = F)
brain4

### Downstream Processing, dimension reduction, knn graph construction
brain4 <- RunPCA(brain4, features = VariableFeatures(object = brain4), verbose = F)
brain4 <- RunUMAP(brain4,reduction = "pca",dims = 1:30, verbose = F)
brain4 <- FindNeighbors(brain4,reduction = "pca", k.param = 6, dims = 1:30, verbose = F, return.neighbor = T)

### Create a SPACO object using the constructed KNN graph

SpaCoObject4<- create_SpaCoObject_from_KNN(brain4, n = 6)

SpaCoObject4   <- RunSCA(SpaCoObject4,compute_nSpacs = TRUE)

SpaCoObject4@nSpacs

### Compute spatial variable genes (SVG’s).

DE_genes4<- SVGTest(SpaCoObject4)

head(DE_genes4)

DE_genes_sortbrain4 <- DE_genes4[order(DE_genes4$score, decreasing = TRUE),]

sigs_brain4 <- rownames(DE_genes4[DE_genes4$p.adjust<0.05,])



### PART 3
### Now lets analyse the single cell Hippo dataset

hippo <- readRDS("C:/Users/Bamidele David/Downloads/hippo.rds")

hippo

hippo <- PercentageFeatureSet(hippo, pattern = "^mt-" ,col.name = "percent.mt")
hippo <- PercentageFeatureSet(hippo, pattern = "^Hbb-" ,col.name = "percent.hbb")

hippo <- NormalizeData(hippo)
hippo <- FindVariableFeatures(hippo, selection.method = "vst", nfeatures = 3000)
hippo <- SCTransform(hippo, assay = "RNA", ncells = 3000, variable.features.n = 3000, verbose = F)

### Downstream Processing
hippo <- RunPCA(hippo, features = VariableFeatures(object = hippo), verbose = F)
hippo <- RunUMAP(hippo,reduction = "pca",dims = 1:30, verbose = F)
hippo <- FindNeighbors(hippo,reduction = "pca", k.param = 6, dims = 1:30, return.neighbor = T, verbose = F)
hippo <- FindClusters(hippo, resolution = 0.5, verbose = FALSE)

### Visualize the Louvain clustering in a UMAP plot

p3 <- DimPlot(hippo, reduction = "umap", group.by = 'subclass',
        label = TRUE, label.size = 5)  + ggtitle('Annotated Clusters of the Mouse Hippocampus sc-RNA seq dataset')
p3


### Create a SPACO object using the constructed KNN graph

SpaCoObject5<- create_SpaCoObject_from_KNN(hippo, n = 6)

SpaCoObject5   <- RunSCA(SpaCoObject5,compute_nSpacs = TRUE)

SpaCoObject5@nSpacs


### Compute spatial variable genes (SVG’s).

DE_genes5<- SVGTest(SpaCoObject5)

head(DE_genes5)

DE_genes_sort5 <- DE_genes5[order(DE_genes5$score, decreasing = TRUE),]

sigs_hippo <- rownames(DE_genes5[DE_genes5$p.adjust<0.05,])



### PART 4

### Compare the intersection of SVGs using Venn Diagrams

# Anterior 1 as spatial and scRNA data

Venn_1 <- list(antr1_as_spatial= sigs_brain, antr1_as_sc= sigs_brain2)

ggVennDiagram(Venn_1, color= "blue", lwd = 10, lty =1) +
  scale_fill_gradient(low= '#B9F8D3', high ='#E78EA9') + 
  guides(fill= 'none') +
  ggtitle("Anterior1 as spatial and single cell data sets")

# Anterior 2 as spatial and scRNA data
Venn_2 <- list(antr2_as_spatial= sigs_brain3, antr2_as_sc= sigs_brain4)

ggVennDiagram(Venn_2, color= "blue", lwd = 10, lty =1) +
  scale_fill_gradient(low= '#B9F8D3', high ='#E78EA9') + 
  guides(fill= 'none') +
  ggtitle("Anterior2 as spatial and single cell data sets")

# Anterior 1 and anterior 2 spatial datasets
Venn_3 <- list(anterior1= sigs_brain, anterior2= sigs_brain3)

ggVennDiagram(Venn_3, color= "blue", lwd = 10, lty =1) +
  scale_fill_gradient(low= '#B9F8D3', high ='#E78EA9') + 
  guides(fill= 'none') +
  ggtitle("Anterior1 & 2 spatial datasets")

# Anterior 1 and anterior 2 scRNA datasets
Venn_4 <- list(anterior1= sigs_brain2, anterior2= sigs_brain4)

ggVennDiagram(Venn_4, color= "blue", lwd = 10, lty =1) +
  scale_fill_gradient(low= '#B9F8D3', high ='#E78EA9') + 
  guides(fill= 'none') +
  ggtitle("Anterior1 & 2 as single cell datasets")


# Spatially- treated anterior 1 & 2 and hippocampus scRNA datasets

Venn_5 <- list(antr1= sigs_brain, scRNA_hippo = sigs_hippo, antr2= sigs_brain3)

ggVennDiagram(Venn_5, color= "blue", lwd = 10, lty =1) +
  scale_fill_gradient(low= '#B9F8D3', high ='#E78EA9') + 
  guides(fill= 'none') +
  ggtitle("At sc-RNA seq KNN = 6")

# Single cell - treated anterior 1 & 2 and hippocampus scRNA datasets

Venn_6 <- list(antr1= sigs_brain2, scRNA_hippo = sigs_hippo, antr2= sigs_brain4)

ggVennDiagram(Venn_6, color= "blue", lwd = 10, lty =1) +
  scale_fill_gradient(low= '#B9F8D3', high ='#E78EA9') + 
  guides(fill= 'none') +
  ggtitle("At sc-RNA seq KNN = 6")



### Visualising the SVGs in UMAP and feature plots

### First, we choose 6 SVGs obtained from spatially - treated anterior 1: Adora2a,Nts,3110035E14Rik,Mbp,Lmo3 and Ttr

### We check if these SVGs exists in the single cell treated anterior 1

gene_to_check <- "Ttr" ### This was done for each SVG
gene_exists <- gene_to_check %in%sigs_brain2

print(gene_exists)

### We visualize the 6 SVGs as a spatial feature and UMAP plot

SpatialFeaturePlot(brainPCA, features = c('Adora2a','Nts','3110035E14Rik','Mbp','Lmo3','Ttr'), ncol = 3)  
+ ggtitle('Spatial Layout of Selected SVGs from the anterior 1 dataset')

p3 <- FeaturePlot(brainPCA, features = c('Adora2a','Nts','3110035E14Rik','Mbp','Lmo3','Ttr'), 
                  min.cutoff = 'q10') + ggtitle('UMAP Layout of selected SVGs from the anterior 1 dataset')


p3


### Choose 6 SVGs obtained from the Hippocampus scRNA Data: Calb1, Hpca, Gabrd, Pcp4,Wfs1 and Pvalb

### Visualize the hippo marker genes 
p4 <- FeaturePlot(hippo, features = c('Calb1','Hpca','Gabrd','Pcp4','Wfs1','Pvalb'), 
                  min.cutoff = 'q10')
p4


### Also check if these SVGs are present in the anterior 1 spatial dataset

gene_to_check <- "Calb1" ### This was done for each of the 6 SVG
gene_exists <- gene_to_check %in%sigs_brain

print(gene_exists)


### Visualize the hippo marker genes locations in the spatial dataset

SpatialFeaturePlot(brainPCA, features = c('Calb1','Hpca','Gabrd','Pcp4','Wfs1','Pvalb'), ncol = 3)



























































