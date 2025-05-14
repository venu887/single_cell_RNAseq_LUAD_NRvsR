#@@@@@@@@@@@@@@@@@@@@@@@@@@
# Separation of IMMUNE and NON-Immune cells in all samples
#@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)


source("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/source_fxns4seurat.R")

##@@@@@@@@@@@@@@@@@@@@@@@@@@ Step-1:  Panglo data
inp<-"/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/"
fin_db1 <- paste(inp,"PanglaoDB_markers_27_Mar_2020.tsv",sep="")
panglao <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)
listOrgan = c("Lungs","Immune system","Epithelium",
              "Connective tissue","Vasculature")
panglao = panglao[is.element(panglao$organ,listOrgan),]

unique(panglao$cell.type)

listSelCT = c("T cells","B cells","NK cells","Plasma cells", 
              "Dendritic cells",
              "Macrophages","Alveolar macrophages",
              "Mast cells","Monocytes","Neutrophils",
              "Plasmacytoid dendritic cells",
              "Pulmonary alveolar type I cells",
              "Pulmonary alveolar type II cells",
              "Epithelial cells", 
              "Airway epithelial cells",
              "Fibroblasts", 
              "Stromal cells",
              "Endothelial cells","Ciliated cells","Ionocytes")

tmpsel  = grepl("Hs",panglao$species) & 
  is.element(panglao$cell.type,listSelCT)
db_sub = panglao[tmpsel,2:3]
colnames(db_sub)=c("Symbol","cellType")
unique(db_sub$cellType)
 
##@@@@@@@@@@@@@@@@@@@@@@@@@@ Step-2 Export all data 
fout2 <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/data1_integrated_ndim100.rds")
integrated_data<- readRDS(fout2)

integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 100, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, dims = 1:50) # 20 good

ElbowPlot(integrated_data)
# Find optimal resolution
# library(clustree)
# resolutions <- seq(0, 1.5, 0.1)
# n.dims <- 20
# combined <- FindClusters(integrated_data, reduction.type = "cca.aligned",
#                          dims.use = 1:n.dims, resolution = resolutions)
# 
# clustree(combined)
# rm(combined)


integrated_data <- FindClusters(integrated_data, resolution =1)
integrated_data <- RunUMAP(integrated_data, dims = 1:49)
p1<-DimPlot(integrated_data, reduction = "umap", label = T)
p1
pdf("/home/u251079/scRNA/Lung_scRNA/2_Raw_UMAP.pdf", width = 6, height =5)
print(p1)
dev.off()

allClust_mks <- FindAllMarkers(integrated_data, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
fin_4mks <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/data1_integrated_allClust_markers.txt")
write.table(allClust_mks, fin_4mks, sep = "\t", quote = FALSE, row.names = TRUE)


###@@@@@@@@@@@@@@@@@@@@@@@@@@ Step 3: Calculate scores using Enrichment and Average expression from step-1 and step-2
## Step 3.1: ## Find the top two highest scores and associated cell types for each cluster using ENRICHMENT method
enrichedData <- enrichScoreCalc_mt(integrated_data,
                                   allClust_mks, db_sub)

nclusters <- nrow(enrichedData$score_Mt)
Ftab_candidate <- matrix(NA,nclusters,4) # Fisher exact test

for (i in c(1:nclusters)){
  tmpi_hi <- order(enrichedData$score_Mt[i,],decreasing = TRUE)[c(1,2)]
  tmpscore_hi <- enrichedData$score_Mt[i,tmpi_hi]
  tmpCell_hi <- colnames(enrichedData$score_Mt)[tmpi_hi]
  Ftab_candidate[i,] <- c(tmpscore_hi,tmpCell_hi)
}

colnames(Ftab_candidate) <- c("score1","score2","cellType1","cellType2")
rownames(Ftab_candidate) <- rownames(enrichedData$score_Mt)

## Step 3.2: using canonical markers from literature
gList = c("PTPRC",    #CD45 for immune
          "CD3D",     #T_cells
          "CD79A",    #B and plasma
          "CD68",     # myeloid
          "COL1A1",   # fibroblast
          "PECAM1",   #CD31 for endothelial
          "EPCAM","SCGB1A1")  #lung epithelial marker    


# Get (arithmetic) average of markers: non-log scale
tmp_averages <- AverageExpression(integrated_data,
                                  assays="RNA",
                                  features=gList)$RNA
# Assuming tmp_averages is a transposed matrix with numeric data
# Convert the matrix to a data frame first
tmp_averages_df <- as.data.frame(tmp_averages)
tmp_averages_df$Gene = rownames(tmp_averages)
mks_clust_avg <- tmp_averages_df[, c(ncol(tmp_averages_df), 1:(ncol(tmp_averages_df)-1))]
mks_clust_avg<-t(mks_clust_avg)
mks_clust_avg<-mks_clust_avg[-1,]
mks_clust_avg1<-as.data.frame(mks_clust_avg)
mks_clust_avg1$max_col <- colnames(mks_clust_avg1)[max.col(mks_clust_avg1, ties.method = 'first')]
Ftab_out <- cbind(mks_clust_avg1,Ftab_candidate)
class(Ftab_out)
Ftab_out$score1 <- as.numeric(Ftab_out$score1)
Ftab_out$score2 <- as.numeric(Ftab_out$score2)

# ## Step 4: Manually inspect the Step.3 output and categorize into:
# endothelial --> fibroblast --> epithelial, including secteroty cells
# --> immune based on:
#   (1) if difference between 1st and 2nd highest score > 1,
# assign cell type associated with 1st highest score
# (2) else
#   assigne based on cell type whose marker expression is highest
# Create matrix of clusters X cell_Type for splitting immune vs non-immune
Ftab_out$score_diff <- Ftab_out$score1 - Ftab_out$score2
Ftab_out$cluster_cell<-ifelse(Ftab_out$score_diff >1, Ftab_out$cellType1, Ftab_out$max_col)
table(Ftab_out$cluster_cell)
unique(Ftab_out$cluster_cell)

Ftab_out$cluster_cell1<-ifelse(Ftab_out$cluster_cell %in% "EPCAM","Epithelial",Ftab_out$cluster_cell)
Ftab_out$cluster_cell1<-ifelse(Ftab_out$cluster_cell1 %in% "PTPRC", "CD45_immune",Ftab_out$cluster_cell1)
Ftab_out$cluster_cell1<-ifelse(Ftab_out$cluster_cell1 %in% "COL1A1", "Fibroblast",Ftab_out$cluster_cell1)
unique(Ftab_out$cluster_cell1)
table(Ftab_out$cluster_cell1)

Non_Imm<-c( "Pulmonary alveolar type I cells",
            "Pulmonary alveolar type II cells",
            "Clara cells","Airway epithelial cells",
            "Airway goblet cells","Ciliated cells",
            "Fibroblast","Stromal cells","Ionocytes",
            "Endothelial cells", "Endothelial",
            "lung epithelial", "fibroblast", "Epithelial", "Fibroblast")
Ftab_out$cluster_type_category <- ifelse(Ftab_out$cluster_cell1 %in% Non_Imm,"non-immune" ,"immune")

# Print the updated table
print(table(Ftab_out$cluster_cell1, Ftab_out$cluster_type_category))
Ftab_out$cluster_assign<-rownames(Ftab_out)
Ftab_out$cluster_assign<- gsub("^g", "", rownames(Ftab_out))

library(dplyr)
Ftab_out <- Ftab_out %>%
  mutate(
    broad_classification = case_when(
      cluster_cell1 == "Epithelial" ~ "Epithelial",
      cluster_cell1 == "Endothelial cells" ~ "Endothelial",
      cluster_cell1 == "Fibroblast" ~ "Fibroblasts",
      cluster_cell1 %in% c("T cells", "B cells", "NK cells", "Plasma cells") ~ "Immune",
      cluster_cell1 == "CD45_immune" ~ "Immune",
      TRUE ~ "Other" # Default case for unclassified cells
    )
  )

fout <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/data1_integrated_annot.txt")
write.table(Ftab_out, fout, sep="\t", quote = F, row.names = F)


finAnnot <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/data1_integrated_annot.txt")
tmp <- read.delim(finAnnot,sep="\t",head=T, stringsAsFactors=FALSE)
# Immune and non immune categories
# tmp$cluster_type_category  for Immune and non_immune plot
# tmp$cluster_cell1 for cell type
# tmp$broad_classification for broad classification
tmp<-Ftab_out
cellID_vec <- tmp$cluster_type_category 
num<-length(cellID_vec)
numbers <- seq(0, (num-1))
names(cellID_vec) <- numbers

## check available cell attributions
integrated_data2<-integrated_data
DimPlot(integrated_data,reduction = "umap", label = T)
DefaultAssay(integrated_data2) <- "RNA"
integrated_data2[["RNA"]] <- JoinLayers(integrated_data2[["RNA"]])
## add cell type annotation
tmpClust = Idents(integrated_data2)
integrated_data2 <- AddMetaData(integrated_data2,
                                metadata=tmpClust,col.name='cellID_vec') # cellID_lab or cellID_vec
View(integrated_data2@meta.data)
integrated_data2 <- RenameIdents(integrated_data2,cellID_vec) # cellID_lab or cellID_vec
tmp = Idents(integrated_data2)
integrated_data2 <- AddMetaData(integrated_data2,tmp,
                                col.name = 'cellIDAll')

tmpClust = Idents(integrated_data2)
table(integrated_data2@meta.data$cellIDAll)
table(integrated_data2@meta.data$integrated_snn_res.1, integrated_data2@meta.data$cellIDAll)
integrated_data2@meta.data$Cell_type1 <- paste0(integrated_data2@meta.data$integrated_snn_res.1,"_",integrated_data2@meta.data$cellIDAll)

gList = c("PTPRC",    #CD45 for immune
          "CD3D",     #T_cells
          "CD79A",    #B and plasma
          "CD68",     # myeloid
          "COL1A1",   # fibroblast
          "PECAM1",   #CD31 for endothelial
          "EPCAM","NKX2-1","SFTPC")  #lung epithelial marker    

plot1<-FeaturePlot(integrated_data2, features =gList,order = T)
plot1

pdf("/home/u251079/scRNA/Lung_scRNA/2_1_mark_FeaturePlot.pdf", width = 15, height =10)
print(plot1)
dev.off()


#++++++++++  PLOTs
p1<-DimPlot(integrated_data2, reduction = "umap", group.by = "cellIDAll", label = T)+NoLegend()+ 
  ggtitle("")  
p1

pdf("/home/u251079/scRNA/Lung_scRNA/2_1_subFig1B_Imm_non_imm.pdf", width = 5, height =4)
print(p1)
dev.off()


p3<-DimPlot(integrated_data2, reduction = "umap", group.by = "cellIDAll", split.by = "sample_type", label = T)+ggtitle("Batch wise samples: Non-Recurrence and Recurrence")+NoLegend()
p4<-DimPlot(integrated_data2, reduction = "umap", group.by = "cellIDAll", split.by = "sample", ncol = 4,label = F)+ggtitle("Individual samples")+NoLegend()
combined_plot<-p3 |p4
# Display the combined plot
print(combined_plot)

pdf("/home/u251079/scRNA/Lung_scRNA/2_1Supp_2_Imm_non_imm_batch.pdf", width = 13, height =5)
print(combined_plot)
dev.off()


# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Evaluation of cluster-8,
# # Evaluation of cluster-8, which is associated with Immune markers and non-immune markers 
# mark_list = c("PTPRC",    #CD45 for immune
#           "CD68",     #myeloid
#           "COL1A1",   #fibroblast
#           "PECAM1",   #CD31 for endothelial
#           "EPCAM", "SFTPC")  #lung epithelial marker  
# Idents(integrated_data) <- "integrated_snn_res.1"
# cluster_cells <- WhichCells(integrated_data, idents = "8") 
# p1<-DimPlot(integrated_data, reduction = "umap", cells = cluster_cells)+NoLegend() +
#   ggtitle("Cluster 8: Immune")
# p1
# p2<-FeaturePlot(integrated_data, features = mark_list, cells = cluster_cells)
# p2
# 
# # Subset the data to include only the selected cells
# subset_data <- subset(integrated_data, cells = cluster_cells)
# p3<-DotPlot(subset_data, features = mark_list, cols = c("blue", "red"), dot.scale = 8)+RotatedAxis()+coord_flip()
# p3
# pdf("/home/u251079/scRNA/plots_Rec_nsclc/Supp_2_Imm_Clust_8.pdf", width = 15, height =7)
# print(p3)
# dev.off()
# 
# combine_plot_8_1<-p1|p3
# pdf("/home/u251079/scRNA/plots_Rec_nsclc/Supp_2_Imm_Clust_8_1.pdf", width = 8, height =3.5)
# print(combine_plot_8_1)
# dev.off()

# https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/ 
#install.packages("scCustomize")
library(scCustomize)
scCustomize::Clustered_DotPlot(integrated_data2, features = gList, group.by = "Cell_type1",
                               plot_km_elbow = FALSE)


integrated_data2@meta.data$Dot_plot<-paste0(integrated_data2@meta.data$clustID,"_", integrated_data2@meta.data$cellIDAll)
pdf("/home/u251079/scRNA/Lung_scRNA/2_Supp_Clustered_DotPlot.pdf", width = 11, height =6)
scCustomize::Clustered_DotPlot(integrated_data2, features = gList, group.by = "Cell_type1",
                               plot_km_elbow = FALSE)
dev.off()





saveRDS(integrated_data2,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/Integrated_cd45pos_neg.rds"))


Idents(integrated_data2)
cd45pos_seurat <- subset(integrated_data2,idents='immune')
View(cd45pos_seurat@meta.data)
cd45pos_seurat
saveRDS(cd45pos_seurat,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45pos_nsclc.rds"))

cd45neg_seurat <- subset(integrated_data2,idents='immune',invert=TRUE)
cd45neg_seurat
saveRDS(cd45neg_seurat,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45neg_nsclc.rds"))




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@# SingleR Annotation
# https://bioconductor.org/books/3.17/SingleRBook/annotation-diagnostics.html
## Step 1:Automatic Annotation with SingleR all cells for immune and non-immune cell types 
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
inf <-c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/Integrated_cd45pos_neg.rds") # merged_nsclc_Singlet_integrated_data_Allclust_seurat_data.rds
data <-readRDS(inf)
names(data@meta.data)
View(data@meta.data)
data

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# run SingleR (default mode) ---------
DefaultAssay(data) <- "RNA"
data
data[["RNA"]]<- JoinLayers(data[["RNA"]])
counts <- GetAssayData(data, layer = 'counts', assay = "RNA")

# Run SingleR
pred <- SingleR(test = counts,
                ref = ref,
                labels = ref$label.main)

pred

data$singleR.labels <- pred$labels[match(rownames(data@meta.data), rownames(pred))]
plot1<-DimPlot(data, reduction = 'umap', group.by = 'singleR.labels', label=T)
print(plot1)

# immune and non-immune cell types
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

data@meta.data$cell_type <- data$singleR.labels # Add immune/non-immune classification to metadata
data@meta.data$category <- ifelse(data@meta.data$cell_type %in% immune_cell_types, "Immune", "Non-Immune")
unique(data@meta.data$singleR.labels)
p1<-DimPlot(data, reduction = "umap", group.by = "category", label = T) + ggtitle("SingleR: Immune and Non-Immune cells")+NoLegend()
print(p1)
pdf("/home/u251079/scRNA/Lung_scRNA/2_SingleR_Imm_Non-Immune_plot.pdf", width = 5, height =4)
print(p1)
dev.off()









