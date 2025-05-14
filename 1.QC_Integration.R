rm(list=ls())
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

#@@@@@@@@@@@@@@@@@@@@@
# STEP-1 QUALUTY CONTROL IN EACH SAMPLE AND MERGING AND DUPLICATE REMOVAL 
#@@@@@@@@@@@@@@@@@@@@@
set.seed(1234)
base_path <- "/cellranger_scran/"
sample_dirs <- c("OX_60912_LC38", "OX_60913_LC52", "OX_60915_LC63",  "OX_60914_LC57",
                 "OX_60916_LC71", "OX_60917_LC104", "OX_60918_LC115", "OX_60919_LC221")

for (sample_dir in sample_dirs) {
  full_path <- file.path(base_path, sample_dir, "outs/filtered_feature_bc_matrix")
  cts <- Read10X(data.dir = full_path)
  seurat_obj_name <- paste0(sample_dir, "_seurat")
  assign(seurat_obj_name, CreateSeuratObject(counts = cts))
}

# Provide proper name and sample names and sample_type for better downstream analysis 
LC38_NR<- OX_60912_LC38_seurat
LC38_NR@meta.data$sample<-"LC38_NR"
LC38_NR@meta.data$sample_type<-"Non_Rec"
LC52_R<- OX_60913_LC52_seurat
LC52_R@meta.data$sample<- "LC52_R"
LC52_R@meta.data$sample_type<-"Rec"
LC57_R<- OX_60914_LC57_seurat
LC57_R@meta.data$sample<-"LC57_R"
LC57_R@meta.data$sample_type<-"Rec"
LC63_NR<- OX_60915_LC63_seurat
LC63_NR@meta.data$sample<-"LC63_NR"
LC63_NR@meta.data$sample_type<-"Non_Rec"
LC71_R<- OX_60916_LC71_seurat
LC71_R@meta.data$sample<- "LC71_R"
LC71_R@meta.data$sample_type<-"Rec"
LC104_R<- OX_60917_LC104_seurat
LC104_R@meta.data$sample<-"LC104_R"
LC104_R@meta.data$sample_type<-"Rec"
LC115_NR<- OX_60918_LC115_seurat
LC115_NR@meta.data$sample<-"LC115_NR"
LC115_NR@meta.data$sample_type<-"Non_Rec"
LC221_NR<- OX_60919_LC221_seurat
LC221_NR@meta.data$sample<-"LC221_NR"
LC221_NR@meta.data$sample_type<-"Non_Rec"

rm(cts)
ls()
# List of sample names to remove UNNECESSORY ONCE ASSIGN  
samples_to_remove <- c("OX_60912_LC38_seurat", "OX_60913_LC52_seurat", "OX_60914_LC57_seurat",
                       "OX_60915_LC63_seurat", "OX_60916_LC71_seurat", "OX_60917_LC104_seurat",
                       "OX_60918_LC115_seurat", "OX_60919_LC221_seurat")
for (sample_name in samples_to_remove) {
  if (exists(sample_name)) {
    rm(list = sample_name)
  }
}

#@@@@@@@@@@@@@@@@@@@      QUALUTY CONTROL IN EACH SAMPLE AND MERGING 
LC38_NR[["percent.mt"]] <- PercentageFeatureSet(LC38_NR, pattern = "^MT-")
LC52_R[["percent.mt"]] <- PercentageFeatureSet(LC52_R, pattern = "^MT-")
LC57_R[["percent.mt"]] <- PercentageFeatureSet(LC57_R, pattern = "^MT-")
LC63_NR[["percent.mt"]] <- PercentageFeatureSet(LC63_NR, pattern = "^MT-")
LC71_R[["percent.mt"]] <- PercentageFeatureSet(LC71_R, pattern = "^MT-")
LC104_R[["percent.mt"]] <- PercentageFeatureSet(LC104_R, pattern = "^MT-")
LC115_NR[["percent.mt"]] <- PercentageFeatureSet(LC115_NR, pattern = "^MT-")
LC221_NR[["percent.mt"]] <- PercentageFeatureSet(LC221_NR, pattern = "^MT-")

# # Because of the sample LC57_R sample has less number of cells, excluded from further analysis
# # Cells with few reads are likely low-quality cells or empty droplets. Cells with an abnormally high number of reads might be doublets or multiplets. 
VlnPlot(LC38_NR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LC52_R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LC63_NR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LC71_R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LC104_R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LC115_NR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LC221_NR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Filtered out "bad" cells # https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/Seurat_QC_to_Clustering/#:~:text=%2C%20percent%20ribosomal).-,QC%20metrics%20are%20stored%20as%20metadata,mapped%20to%20the%20mitochondrial%20genome.&text=This%20entire%20workflow%20is%20exploratory,thresholds%20after%20performing%20downstream%20steps.
LC38_NR <- subset(LC38_NR, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC52_R <- subset(LC52_R, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC57_R <- subset(LC57_R, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC63_NR <- subset(LC63_NR, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC71_R <- subset(LC71_R, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC104_R <- subset(LC104_R, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC115_NR <- subset(LC115_NR, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)
LC221_NR <- subset(LC221_NR, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA >1000 & percent.mt < 20)

# merge 7 Seurat objects together
merged_nsclc<- merge(x=LC38_NR, y=c(LC63_NR, LC115_NR, LC221_NR, LC52_R, LC104_R, LC71_R),
                     project = 'NSCLC')


#@@@@@@@@@@@@@@@   DUPLICATE REMOVAL 
merged_nsclc$orig.ident<-"nsclc_project"
merged_nsclc[["RNA"]] <- JoinLayers(merged_nsclc[["RNA"]])
merged_nsclc <- NormalizeData(object = merged_nsclc)
merged_nsclc <- FindVariableFeatures(object = merged_nsclc)
merged_nsclc <- ScaleData(object = merged_nsclc)
merged_nsclc <- RunPCA(object = merged_nsclc)
merged_nsclc <- FindNeighbors(object = merged_nsclc, dims = 1:20)
merged_nsclc <- FindClusters(object = merged_nsclc) 
merged_nsclc <- RunUMAP(object = merged_nsclc, dims = 1:20)
DimPlot(merged_nsclc, reduction = "umap")

## pK Identification (no ground-truth) OPTIMUM pK values 
sweep.res.list_nsclc <- paramSweep(merged_nsclc, PCs = 1:20, sct = F)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)

ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1)) + # Maximum BCmetric values is the optimum pK values 
  geom_point() +
  geom_line()
pK <- bcmvn_nsclc %>%
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
# select the pK that corresponds to max bcmvn to optimize doublet detection
pK <- as.numeric(as.character(pK[[1]]))
## Homotypic Doublet Proportion Estimate 
annotations <- merged_nsclc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(merged_nsclc@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# run doubletFinder 
merged_nsclc <- doubletFinder(merged_nsclc, 
                                      PCs = 1:20, 
                                      pN = 0.25, 
                                      pK = pK, 
                                      nExp = nExp_poi.adj,
                                      reuse.pANN = FALSE, sct = F)
merged_nsclc<-subset(merged_nsclc, subset = DF.classifications_0.25_0.05_2503 =="Singlet")
View(merged_nsclc@meta.data)

fout1 <- c("filtered_data.rds")
saveRDS(merged_nsclc,file=fout1)

merged_nsclc<-readRDS(fout1)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ INTEGRATE DATA FOR BATCH EFFECT USING CCA METHOD 
sample_list <- SplitObject(merged_nsclc, split.by = "sample")
sample_list <- lapply(sample_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", 
                              nfeatures = 2000, verbose = FALSE)
  return(obj)
})

anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:100) 
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:100)
integrated_data
data<-integrated_data
fout2 <- c("integrated_data.rds")
saveRDS(integrated_data,file=fout2)
