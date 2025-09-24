library(Seurat)
library(SeuratObject)
library(readr)
library(haven)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(R.utils)
library(Rcpp)
library(harmony)
library(dplyr)
library(scDblFinder)
library(Matrix)
library(patchwork)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(cowplot)

# input file
seu_data <- readRDS(seurat_file)
seu_data[["percent.mt"]] <- PercentageFeatureSet(seu_data, pattern = "^mt-")

## filter low quality cells
seu_data <- subset(seu_data, subset = percent.mt<20)

## filter doublet
seu_data <- as.SingleCellExperiment(seu_data)
seu_data <- scDblFinder(seu_data, samples = "orig.ident")
table(seu_data$scDblFinder.class)
seu_data <- as.Seurat(seu_data)
DefaultAssay(seu_data) <- "RNA"

seu_data <- subset(seu_data,subset = scDblFinder.class == "singlet")


# PCA and cluster
seu_data <- NormalizeData(seu_data, normalization.method = 'LogNormalize', scale.factor = 10000)
seu_data <- FindVariableFeatures(seu_data, selection.method = "vst", nfeatures = 3000)
seu_data <- ScaleData(seu_data, vars.to.regress = "percent.mt", features = VariableFeatures(seu_data))

# run PCA 
seu_data <- RunPCA(seu_data, features = VariableFeatures(seu_data))

# harmony
seu_data <- RunHarmony(seu_data,"orig.ident")
seu_data <- FindNeighbors(seu_data, reduction = "harmony", dims = 1:35)
seu_data <- FindClusters(seu_data, resolution = 0.8)
seu_data <- RunUMAP(seu_data, reduction = "harmony", dims = 1:35)

# Cell distributions across different groups using Seurat’s DimPlot
DimPlot(seu_data, reduction = "umap",group.by = "Type",label = T,cols = color, 
        raster = FALSE,label.size = 8) +  
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
        axis.line = element_line(colour = "black", size = 0.8),  
        axis.ticks = element_line(colour = "black", size = 0.5),  
        axis.title = element_text(size = 12, face = "bold"),  
        legend.title = element_text(size = 15, face = "bold"),  
        legend.text = element_text(size = 15, face = "bold")   
  ) +
  NoLegend()

# Gene expression levels using Seurat’s FeaturePlot function
FeaturePlot(AML,features = c('YTHDF2'),ncol = 1,raster=FALSE,max.cutoff = 2.5
            ,pt.size = 1) +
  scale_color_gradientn(colours = rainbow)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
  coord_fixed(ratio = 1)