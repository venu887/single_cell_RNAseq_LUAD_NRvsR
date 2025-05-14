#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## *******************************************************************
## Analyze immune cells: using 30PCs and annotate clusters
## *******************************************************************
#https://github.com/hannahsfuchs/OBDS_Course/issues/2
# https://github.com/Zhongqige/HER2_scRNAseq/blob/master/R/cd45_analysis.R
# Python: https://www.sc-best-practices.org/introduction/prior_art.html
# https://broadinstitute.github.io/2019_scWorkshop/pseudotime-cell-trajectories.html
# https://bioconductor.org/books/3.12/OSCA/clustering.html 
# https://link.springer.com/article/10.1007/s12672-024-01357-7

rm(list=ls())
library(SingleCellExperiment)
library(Seurat)
library(dplyr) 
library(ggplot2)
library(NMF)
library(reshape2)
packageVersion("Seurat")
sessionInfo()
source("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/source_fxns4seurat.R")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 1: load Seurat object
set.seed(42)
fin2 <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45pos_nsclc.rds")
cd45pos <-readRDS(fin2)
cd45pos
View(cd45pos@meta.data)
cd45pos@meta.data <- cd45pos@meta.data[, -c(7:14)]
DimPlot(cd45pos, reduction = "umap")
DefaultAssay(cd45pos)<-"integrated"
cd45pos
table(cd45pos@meta.data$sample)
## Step 2: scale and dimentional reduction
ndim=30  #*****
cd45pos <- FindNeighbors(cd45pos, dims = 1:ndim)

# # Find optimal resolution
# library(clustree)
# resolutions <- seq(0, 1.5, 0.1)
# n.dims <- 30
# combined <- FindClusters(cd45pos, reduction.type = "cca.aligned",
#                          dims.use = 1:n.dims, resolution = resolutions)
# 
# clustree(combined)
# rm(combined)

cd45pos <- FindClusters(cd45pos)
cd45pos <- RunUMAP(cd45pos, dims = 1:ndim)


allClust_mks <- FindAllMarkers(cd45pos, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.25)

cd45pos
## export data
write.table(allClust_mks, c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45pos_nsclc_allClust.txt"),
            quote = F, row.names = F)


#@@@@@@@@@@@ Raw UMAP
Imm_plot_1<-DimPlot(cd45pos,  reduction = "umap", label = T)
print(Imm_plot_1)

pdf("/home/u251079/scRNA/Lung_scRNA/3_Imm_UMAP.pdf", width = 6, height =5)
print(Imm_plot_1)
dev.off()
#@@@@@@@@@@@ 







#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 2: Calculate Enrichment scores and average expression using Panglao database and keep major cell types
inp<-"/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/"
fin_db1 <- paste(inp,"PanglaoDB_markers_27_Mar_2020.tsv",sep="")
panglao <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)
listOrgan = c("Lungs","Immune system","Epithelium",
              "Connective tissue","Vasculature")
panglao = panglao[is.element(panglao$organ,listOrgan),]
unique(panglao$cell.type)
immune_cell_types  = c("T cells","B cells","Dendritic cells",
              "Macrophages","Alveolar macrophages",
              "Mast cells","Monocytes","Neutrophils",
              "NK cells","Plasma cells",
              "Plasmacytoid dendritic cells",
              "Gamma delta T cells")

tmpIMM  = grepl("Hs",panglao$species) & 
  is.element(panglao$cell.type,immune_cell_types)
IMM_sub = panglao[tmpIMM,2:3]
colnames(IMM_sub)=c("Symbol","cellType")



## Enrichment approach by using panglao database
enrichedData <- enrichScoreCalc_mt(cd45pos,
                                   allClust_mks,IMM_sub)

nclusters <- nrow(enrichedData$score_Mt)
Ftab_candidate <- matrix(NA,nclusters,4)

for (i in c(1:nclusters)){
  tmpi_hi <- order(enrichedData$score_Mt[i,],decreasing = TRUE)[c(1,2)]
  tmpscore_hi <- enrichedData$score_Mt[i,tmpi_hi]
  tmpCell_hi <- colnames(enrichedData$score_Mt)[tmpi_hi]
  Ftab_candidate[i,] <- c(tmpscore_hi,tmpCell_hi)
}

colnames(Ftab_candidate) <- c("score1","score2","cellType1","cellType2")
rownames(Ftab_candidate) <- rownames(enrichedData$score_Mt)

## Using canonical markers from either literature (FACS markers) 
## or overlap genes based on enrichment analysis
gList = c("CD3D","CD4","LTB","IL7R",  ##CD4 
          "CD8A","GZMB","PRF1", #CD8, NK
          "CD79A",       #Plasma
          "MS4A1",       #B
          "HLA-DQA1",    #DC from Panglao & lit "HLA-DQB1"
          "FCGR3A",       #CD16 macro> mono
          "CD14",         #mono > macro literature/schelker
          "CD68",         #macro > mono 
          "C1QA",         #macro Tirosh   
          "LILRA4",       #pDC
          "CD1C","CADM1", #cDC1 vs. cDC2 both Panglao and literature
          "S100A12",      #neutro Danaher&LM
          "CPA3",         #mast
          "MKI67")        #proliferation,Tgd    


# Get (arithmetic) average of markers: non-log scale
tmp_averages <- AverageExpression(cd45pos,
                                  assays="RNA",
                                  features=gList)$RNA
tmp_averages <- t(tmp_averages)
tmp_zscores <- scale(tmp_averages)
mks_clust_avg <- cbind(rownames(tmp_zscores),tmp_zscores)
mks_clust_avg<-as.data.frame(mks_clust_avg)
mks_clust_avg1<-mks_clust_avg[,-1]

## export output from enrichment and marker expression
Ftab_out <- cbind(mks_clust_avg1,Ftab_candidate)
Ftab_out[, 1:22] <- lapply(Ftab_out[, 1:22], as.numeric)
str(Ftab_out)

# Calculate the score difference
Ftab_out$score_diff <- Ftab_out$score1 - Ftab_out$score2
Ftab_out$cluster_cell<-ifelse(Ftab_out$score_diff >1, Ftab_out$cellType1, NA)
Ftab_out$cluster<-rownames(Ftab_out)
#@@@@@@@@@@@@
Step 3.4: Manually inspect the 3.3 output and assign cell type:
  (1) if difference between 1st and 2nd highest score > 1,
assign cell type associated with 1st highest score
(2) else
  assign based on cell type expressing canonical markers:
  Rules for assigning immune cell type based on enrichment and zscore
NK:  T-cells or NK by enrichment & CD3D<0 & PRF1>0.5
NKT: T-cells or NK by enrichment & CD3D>0 & PRF1>0.5
CD4: T-cells by enrichment & sum(CD3D,LTB,IL7R) >2
CD8: T-cells or NK by enrichment & sum(CD3D,CD8A)
B:  by enrichment & MS4A1>3
plasma: by enrichment & 0<CD79A<2
mast: by enrichment & CPA3>1
pDC: by enrichment & LILRA4>1
Neutrophil: either Macrophage or Monocyte or DC by enrichment,
AND S100A12>1
DC: either Macrophage or Monocyte or DC by enrichment,
AND sum(CD1C,CADM1,HLA-DQA1>2)
Mono: either Macrophage or Monocyte or DC by enrichment,
AND (CD14>2)
Macophage: either Macrophage or Monocyte or DC by enrichment,
AND sum(FCGR3A,CD68,C1QA)
#@@@@@@@@@@@@
str(Ftab_out)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster_cell %in% "Plasmacytoid dendritic cells", "pDC", Ftab_out$cluster_cell)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster1 %in% "Dendritic cells", "DC", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster1 %in% "Gamma delta T cells", "GaDelT", Ftab_out$cluster1)
# manual inspection of expression of marker genes, where if the cluster is not assigned by enrichment method, 
# Cluster 0: Average expression (AE) CD8A, CD3D and IL7R :CD8
# Cluster 1: AE CD8A, CD3D, LOW GZMB, LOW PRF1 : CD8
# Cluster 2: High CD79A, MS4A1 : Plasma cells
# Cluster 3: high expression of CD14 : Monocytes 
# Cluster 5: high expression of S100A2 :Neutrophils 
# Cluster 6:  LOW CD4, HIGH IL7R, LTB, CD3D : CD4
# Cluster 7: high LTB, CD4, Low IL7R, CD3D, CADM1 : CD4 
# Cluster 8: High C1QA, FCGR3A, CD68 : Macophage
# Cluster 10: low IL7R, low CD3D, high CD4 : CD4
# Cluster 11: High C1QA, FCGR3A, CD68, CD14 low, HLA : Macophage/DC
# Cluster 12: Low CD79A, MS4A1: Plasma cells
# Cluster 15: CD8A, IL7R, LTB, HIGH PRF1 : CD8 # T cells from both NATs and the TME were mainly composed of CD8-positive T cells and expressed high levels of cytotoxic markers, such as GZMA/B/H/K, PRF1, NKG7, IFNG, GNLY and CXCL13 https://www.nature.com/articles/s41392-022-01150-4 
# Cluster 16: High C1QA, FCGR3A, CD68, High CD1C : Macophage
# Cluster 17: High CPA3 : MAST cells

#@@@@@@@@@@@@@@@@@@@@@@@@
str(Ftab_out)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g0", "CD8", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g1", "CD8", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g2", "B cells", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g3", "Mono", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g5", "Neutrophil", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g6", "CD4", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g7", "CD4", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g8", "Mac", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g10", "CD4", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g11", "Mac|DC", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g12", "Plasma cells", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g15", "CD8", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g16", "Prolifirating Mac", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g17", "Mast", Ftab_out$cluster1)
Ftab_out$cluster1<-ifelse(Ftab_out$cluster %in% "g20", "Prolif Mac|DC", Ftab_out$cluster1)


fout <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45pos_Ftab_out.txt")
write.table(Ftab_out, fout, sep="\t",
            quote = F, row.names = F)









#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 3: Update cell types to the immune data
fin3 <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45pos_Ftab_out.txt") # Manual assessment
cellAnnotInfo <- read.delim(fin3,sep="\t",header=T,
                            stringsAsFactors=FALSE)
# Assign cell types to meta data
names(Ftab_out)
num <- length(Ftab_out$cluster1)
numbers <- seq(0, (num - 1))
Ftab_out$cluster_assign <- numbers
DefaultAssay(cd45pos)<-"RNA"
# Match cluster numbers from Ftab_out to integrated_snn_res.1 in allSam.integrated2 and assign cluster cell names
cd45pos@meta.data$Cell_type <- ifelse(cd45pos@meta.data$seurat_clusters %in% Ftab_out$cluster_assign, 
                                             Ftab_out$cluster1[match(cd45pos@meta.data$seurat_clusters, Ftab_out$cluster_assign)], 
                                             NA)
table(cd45pos@meta.data$Cell_type)
# Combine cluster number and cluster name into Cell_type
cd45pos@meta.data$Cell_type1 <- paste0( cd45pos@meta.data$Cell_type,"_","C",cd45pos@meta.data$integrated_snn_res.0.8)

#+++++++++++ Plot-2 UMAP
plot1<-DimPlot(cd45pos, group.by = 'Cell_type1', split.by = "sample_type", label = TRUE, label.size = 6, repel = T) +NoLegend()
plot1
pdf("/home/u251079/scRNA/Lung_scRNA/3_Imm_cell_UMAP1.pdf", width = 16, height =8)
print(plot1)
dev.off()

#+++++++++++ Plot-3 UMAP
p2<-DimPlot(cd45pos, reduction = "umap", group.by = 'Cell_type1', label = TRUE, label.size = 6, repel = T) 
p2
pdf("/home/u251079/scRNA/Lung_scRNA/3_Imm_cell_UMAP_2.pdf", width = 15, height =10)
print(p2)
dev.off()

#+++++++++++ Plot-4 All markerr genes 
pdf("/home/u251079/scRNA/Lung_scRNA/3_Marker_gene_Feature_Plot.pdf", width = 15, height =15)
FeaturePlot(cd45pos, 
            features = gList,
            pt.size = 0.5,                          # Point size for "bubble" effect
            order = TRUE) 
dev.off()


#+++++++++++ Plot-5 Violin Plots
pdf("/home/u251079/scRNA/Lung_scRNA/3_Marker_gene_VlnPlot_Plot.pdf", width = 20, height =15)
VlnPlot(cd45pos, 
        features = gList)
dev.off()


#+++++++++++ Plot-this 6 DOT PLOTs https://www.nature.com/articles/s41597-023-02074-6/tables/4 
file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/cd45pos_nsclc_Immune_cells.rds")
cd45pos<-readRDS(file)
cd45pos

markers <- list(
  "Markers: B-cell" = c("CD79A", "MS4A1", "CD22", "BANK1"),
  "CD4+" = c("CD3D","CTLA4", "CD4","LTB"), 
  "CD8+" = c("PRF1","GZMA", "GZMB", "CD8A","CD8B"),
  "DC" = c("CLEC10A","CD1C","CADM1", "HLA-DQA1"),
  "Macrophages" = c("CD68", "MCEMP1", "FCGR3A", "CD63","MS4A7", "IL1B", "IL4I1",  "APOE", "C1QA", "C1QB", "C1QC"),
  "Mast" = c("KIT", "CPA3"),
  "Mono" = c("CD14", "CSF3R"),
  "Neu" = c( "CSF3R", "S100A8"),
  "NK" = c("KLRF1", "KLRD1", "GNLY", "NKG7"),
  "pDCs" = c("IRF7", "IRF8","LILRA4"),
  "Plasma" = c("IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4"), 
  "Proliferating" = c("TOP2A", "MKI67", "NUSAP1")
)
#+++++++++++ Plot-6.1
all_markers <- unique(unlist(markers))
DefaultAssay(cd45pos) <- "RNA"
valid_markers <- all_markers[all_markers %in% rownames(cd45pos@assays$RNA$counts)]
dot_data <- DotPlot(cd45pos, features = valid_markers, group.by = 'Cell_type1')
plot1<-dot_data+RotatedAxis()+ scale_size(range = c(1,9))

pdf("/home/u251079/scRNA/Lung_scRNA/3_Marker_gene_Dot_Plot.pdf", width = 17, height =7)
print(plot1)
dev.off()

#+++++++++++ Plot-6.2
library(viridis)
all_markers <- unique(unlist(markers))
valid_markers <- all_markers[all_markers %in% rownames(cd45pos@assays$RNA$counts)]
plot2<-DotPlot(cd45pos, features =valid_markers, group.by = 'Cell_type1',
               assay = "RNA", dot.scale = 1, 
               cluster.idents = FALSE) +
  scale_size(range = c(0, 5)) +
  scale_size_area(max_size = 10) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) +
  scale_color_gradientn(colours = c("lightblue", "blue", "darkblue"), limits = c(0, 1), 
                        oob = scales::squish, name = 'log2 (count + 1)')

plot2

#@@@@@@@  Plot-6.2 ADDing name on top on the plot.
markers_df <- stack(markers)
colnames(markers_df) <- c("Symbol", "CellType")
valid_markers_df <- markers_df[markers_df$Symbol %in% rownames(cd45pos@assays$RNA$counts), ]
valid_markers_df <- valid_markers_df[order(valid_markers_df$CellType), ]
celltype_labels <- unique(valid_markers_df$CellType) 
gene_list <- split(valid_markers_df$Symbol, valid_markers_df$CellType) 
valid_markers <- unlist(gene_list)  
names(valid_markers) <- rep(celltype_labels, times = lengths(gene_list))  
valid_markers <- valid_markers[!duplicated(valid_markers)]

plot2 <- DotPlot(cd45pos, 
                 features = valid_markers, 
                 group.by = 'Cell_type1', 
                 assay = "RNA", 
                 dot.scale = 1, 
                 cluster.idents = FALSE) +
  scale_size(range = c(0, 5)) +
  scale_size_area(max_size = 10) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) +
  scale_color_gradientn(colours = c("lightblue", "blue", "darkblue"), 
                        limits = c(0, 1), 
                        oob = scales::squish, 
                        name = 'log2 (count + 1)')

plot2



pdf("/home/u251079/scRNA/Lung_scRNA/3_DotPlot_IMM.pdf", width = 18, height =7.5)
print(plot2)
dev.off()








##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@# Save data 
saveRDS(cd45pos,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/cd45pos_nsclc_Immune_cells.rds"))

fin2=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/cd45pos_nsclc_Immune_cells.rds")
cd45pos <-readRDS(fin2)

DefaultAssay(cd45pos)<-"integrated"
unique(cd45pos@meta.data$Cell_type)
table(cd45pos@meta.data$Cell_type)
# save the data according to the cell type of clusters 
Macrophges <- subset(cd45pos,subset = Cell_type == c("Mac|DC","Mac", "Mono", "Neutrophil", "Prolifirating Mac", "Prolif Mac|DC"))
saveRDS(Macrophges,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/Macrophges.rds"))


# Subset the Seurat object for T cells (CD4, CD8, and γδT cells)
T_cells <- subset(cd45pos, subset = Cell_type %in% c("CD4", "CD8", "GaDelT"))
saveRDS(T_cells,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/T_cells.rds"))


NK_NKT<-subset(cd45pos, subset = Cell_type %in% c("NK cells" ,  "NKT"))
saveRDS(NK_NKT,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/NK_NKT.rds"))


B_cells<-subset(cd45pos, subset = Cell_type %in% c("Plasma cells", "B cells" ))
saveRDS(B_cells,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/B_cells.rds"))

Mast_NA_cells<-subset(cd45pos, subset = seurat_clusters %in% c(17,20))
saveRDS(Mast_NA_cells,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/Mast_NA_cells.rds"))

DC<-subset(cd45pos, subset = Cell_type %in% c("DC", "pDC"))
saveRDS(DC,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/DC_cells.rds"))


# Bingo with cell type assessment
# Each cell type analysis can find in 3.1.Sub_clust_anno.R
# Statistic analysis with all cell and proportion of all cells in sample wise and group wise in 3.2.FDR_P-val_CD45pos.R



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@# SingleR Annotation for CD45 cells
fin2=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/Immune_cells_cd45pos.rds")
cd45pos <-readRDS(fin2)
names(cd45pos@meta.data)
View(cd45pos@meta.data)
# cd45pos<-integrated_data
DefaultAssay(cd45pos)<-"integrated"
# Check the object's processing history
cd45pos@commands

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))
# expression values are log counts (log normalized counts)

# run SingleR (default mode) ---------

DefaultAssay(cd45pos) <- "RNA"
cd45pos[["RNA"]]<- JoinLayers(cd45pos[["RNA"]])
pbmc_counts <- GetAssayData(cd45pos, layer = 'counts', assay = "RNA")

# Run SingleR
pred <- SingleR(test = pbmc_counts,
                ref = ref,
                labels = ref$label.main)

pred
View(cd45pos@meta.data)
unique(cd45pos@meta.data$singleR.labels)
DimPlot(cd45pos, reduction = 'umap')

cd45pos$singleR.labels <- pred$labels[match(rownames(cd45pos@meta.data), rownames(pred))]
plot1<-DimPlot(cd45pos, reduction = 'umap', group.by = 'singleR.labels', label=T)
print(plot1)

# Define immune and non-immune cell types
immune_cell_types <- c(
  "NK_cell", "T_cells", "DC", "Monocyte", "B_cell", 
  "Macrophage", "Neutrophils", "Platelets", "Pro-B_cell_CD34+",
  "HSC_CD34+", "HSC_-G-CSF", "GMP", "MEP", "CMP",
  "Pre-B_cell_CD34-", "BM & Prog."
)

non_immune_cell_types <- c(
  "Smooth_muscle_cells", "Epithelial_cells", "Fibroblasts", 
  "Endothelial_cells", "Astrocyte", "Tissue_stem_cells", 
  "Chondrocytes", "Gametocytes", "Neurons", "Osteoblasts",
  "Neuroepithelial_cell", "Hepatocytes", "Embryonic_stem_cells",
  "iPS_cells"
)

# Add immune/non-immune classification to metadata
cd45pos@meta.data$cell_type <- cd45pos$singleR.labels
cd45pos@meta.data$category <- ifelse(cd45pos@meta.data$cell_type %in% immune_cell_types, "Immune", "Non-Immune")


# View metadata
View(cd45pos@meta.data)
unique(cd45pos@meta.data$singleR.labels)

# Plot annotated clusters on UMAP
DimPlot(cd45pos, reduction = "umap", group.by = "cell_type") + ggtitle("Cell Type")
p1<-DimPlot(cd45pos, reduction = "umap", group.by = "category", label = T) + ggtitle("SingleR: Immune and Non-Immune cells")+NoLegend()
pdf("/home/u251079/scRNA/plots_Rec_nsclc/SingleR_Imm_Non_UMAP.pdf", width = 5, height =5)
print(p1)
dev.off()

cd45pos
table(cd45pos@meta.data$category)

# https://bioconductor.org/books/3.17/SingleRBook/annotation-diagnostics.html
pdf("/home/u251079/scRNA/plots_Rec_nsclc/SingleR_Imm_cell_DotPlot.pdf", width = 12, height =7)
print(plot1)
dev.off()

# Filter out non-immune cells
cd45pos_immune <- subset(cd45pos, category == "Immune")

# Annotation diagnostics ----------
# ...Based on the scores within cells -----------
pred
pred$scores

# Extract the annotation column and ensure it is a vector
annotation_col <- cd45pos_immune@meta.data[,"integrated_snn_res.1", drop=FALSE]

# Ensure the annotation_col is a data frame
annotation_col <- as.data.frame(annotation_col)
colnames(annotation_col) <- "Cluster"

# Ensure the annotation_col names match the cell names in pred
rownames(annotation_col) <- rownames(cd45pos_immune@meta.data)

# Check the head of annotation_col to ensure it is correctly formatted
head(annotation_col)

# Order the cells in pred according to the annotation
ordered_cells <- order(annotation_col$Cluster)
pred_ordered <- pred[ordered_cells, , drop = FALSE]

# Order the annotation_col accordingly
annotation_col_ordered <- annotation_col[ordered_cells, , drop = FALSE]

# Create the heatmap plot with the ordered annotation
plot2 <- plotScoreHeatmap(pred_ordered, annotation_col = annotation_col_ordered)

# Display the plot
print(plot2)



plot2<-plotScoreHeatmap(pred, annotation_col=as.data.frame(cd45pos_immune@meta.data[,"integrated_snn_res.1",drop=T]))
print(plot2)


pdf("/home/u251079/scRNA/plots_Rec_nsclc/SingleR_Imm_Non_Heatmap.pdf", width = 9, height =12)
print(plot2)
dev.off()






















#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@. OTHER PLOTS CODE:

#@@@@@@@@. DOT PLOTs.1
# Define markers for each cell type/subtype: https://www.nature.com/articles/s41597-023-02074-6/tables/4 
View(cd45pos@meta.data)
markers <- list(
  "Proliferating T/NK" = c("TOP2A", "MKI67", "NUSAP1"),
  "NK" = c("KLRF1", "KLRD1", "KLRB1", "GNLY", "NKG7"),
  "Naïve T" = c("PTPRC"),
  "CD8+ Tem" = c("GZMA", "GZMM", "CD8A", "CD8B"),
  "CD4+ Treg" = c("CTLA4", "CD4", "CD3D","LTB","IL7R"),
  "Mature naïve B" = c("CD22", "CD53", "CD79A"),
  "Plasma" = c("IGHA2", "IGHM", "TNFRSF17"),
  "pDCs" = c("IRF7", "IRF8"),
  "cDC2/moDCs" = c("CLEC10A"),
  "Mast" = c("KIT", "CPA3", "CD63"),
  "Monocytes" = c("CD14", "CSF3R"),
  "Low quality Mφ" = c("LYZ", "FTL"),  # Plus high number of MT- and RPL/S genes
  "Lipid-associated Mφ" = c("MS4A7", "IL1B", "IL4I1", "FOLR2", "APOE", "C1QA", "C1QB", "C1QC", "CTSB", "CTSD"),
  "Alveolar Mφ" = c("MCEMP1", "PPARG", "MRC1"),
  "Proliferating Mφ" = c("CDCA8", "MKI67", "CENPF", "CD14", "TOP2A"),
  "Neutrophils" = c("FCGR3B", "CSF3R", "S100A12", "S100A8")
)


#@@@@@@@@. PLOT.2
all_markers <- unique(unlist(markers))
DefaultAssay(cd45pos) <- "RNA"

valid_markers <- all_markers[all_markers %in% rownames(cd45pos@assays$RNA$counts)]

dot_data <- DotPlot(cd45pos, features = valid_markers, group.by = 'Cell_type1')
dot_data_df <- dot_data$data

markers_df <- stack(markers)
colnames(markers_df) <- c("Symbol", "cellType")

dot_data_df <- dot_data_df %>%
  left_join(markers_df, by = c("features.plot" = "Symbol"))

dot_plot <- ggplot(dot_data_df, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "blue") +
  labs(
    title = "Cell Markers from Non-Small Cell Lung Cancer",
    x = "Marker Genes",
    y = "Cell Clusters",
    size = "Percent Expressed",
    color = "Avg. Scaled Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold") # Y-axis title
  ) +
  facet_wrap(~ cellType, scales = "free_x", nrow = 1)   # One row facet

print(dot_plot)

pdf("/home/u251079/scRNA/Lung_scRNA/3_1integrating_NSCLC_markgenes.pdf", width = 25, height =6)
print(dot_plot)
dev.off()




file_path= c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45pos_nsclc_allClust.txt")
# Adjust reading the file
allClust_mks <- read.delim(file_path, header = TRUE, sep = " ", stringsAsFactors = FALSE)

top10 = allClust_mks %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#dotplot
plot = DotPlot(object = cd45pos, features = unique(top10$gene) , assay="RNA") # group.by = "Cell_type1"
plot + RotatedAxis()


