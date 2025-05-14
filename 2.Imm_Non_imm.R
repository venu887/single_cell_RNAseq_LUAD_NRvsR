#@@@@@@@@@@@@@@@@@@@@@@@@@@
# Separation of IMMUNE and NON-Immune cells in all samples
#@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

#Download this file from paper :10.1158/0008-5472.CAN-23-0128
#Single-Cell Characterization of Pulmonary Nodules Implicates Suppression of Immunosurveillance across Early Stages of Lung Adenocarcinoma 
source("source_fxns4seurat.R") 

##@@@@@@@@@@@@@@@@@@@@@@@@@@ Step-1:  Panglo data
inp<-"cellranger_scran/NSCLC_output/"
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
 
##@@@@@@@@@@@@@@@@@@@@@@@@@@ Step-2 Export all data 
fout2 <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/data1_integrated_ndim100.rds")
integrated_data<- readRDS(fout2)
integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 100, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, dims = 1:50) 

ElbowPlot(integrated_data)
# Find optimal resolution
# library(clustree)
# resolutions <- seq(0, 1.5, 0.1)
# n.dims <- 20
# combined <- FindClusters(integrated_data, reduction.type = "cca.aligned",
#                          dims.use = 1:n.dims, resolution = resolutions)
# clustree(combined)
# rm(combined)
integrated_data <- FindClusters(integrated_data, resolution =1)
integrated_data <- RunUMAP(integrated_data, dims = 1:49)


#@@@@@@@@@@@@@@@@@@  Supplementary Figure (SF1) @@@@@@@@@@@@@@@@@@@@@@@@@@
gList = c("PTPRC",    #CD45 for immune
          "CD3D",     #T_cells
          "CD79A",    #B and plasma
          "CD68",     # myeloid
          "COL1A1",   # fibroblast
          "PECAM1",   #CD31 for endothelial
          "EPCAM","SCGB1A1")  #lung epithelial marker   
Plot_1<- FeaturePlot(integrated_data, features = gList)
# Save figure

#@@@@@@@@@@@@@@@@@  Supplementary Figure (SF2A) @@@@@@@@@@@@@@@@@@@@@@@@@@
Plot_2<-DimPlot(integrated_data, reduction = "umap", label = T)
# Save figure


allClust_mks <- FindAllMarkers(integrated_data, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
fin_4mks <- c("NSCLC_output/integrated_allmarkClust.txt")
write.table(allClust_mks, fin_4mks, sep = "\t", quote = FALSE, row.names = TRUE)


###@@@@@@@@@@@@@@@@@@@@@@@@@@ Step 3: Calculate scores using Enrichment and Average expression from step-1 and step-2
## Step 3.1: ## Find the top two highest scores and associated cell types for each cluster using ENRICHMENT method
enrichedData <- enrichScoreCalc_mt(integrated_data,
                                   allClust_mks, db_sub)
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

## Step 3.2: using canonical markers from literature
# Get (arithmetic) average of markers: non-log scale
tmp_averages <- AverageExpression(integrated_data,
                                  assays="RNA",
                                  features=gList)$RNA
tmp_averages_df <- as.data.frame(tmp_averages)
tmp_averages_df$Gene = rownames(tmp_averages)
mks_clust_avg <- tmp_averages_df[, c(ncol(tmp_averages_df), 1:(ncol(tmp_averages_df)-1))]
mks_clust_avg<-t(mks_clust_avg)
mks_clust_avg<-mks_clust_avg[-1,]
mks_clust_avg1<-as.data.frame(mks_clust_avg)
mks_clust_avg1$max_col <- colnames(mks_clust_avg1)[max.col(mks_clust_avg1, ties.method = 'first')]
Ftab_out <- cbind(mks_clust_avg1,Ftab_candidate)
Ftab_out$score1 <- as.numeric(Ftab_out$score1)
Ftab_out$score2 <- as.numeric(Ftab_out$score2)

# ## Step 4: Manually inspect the Step.3 output and categorize into:
# Go to here for understanding https://github.com/LinhTranUCLA/scRNA_subsolid
Ftab_out$score_diff <- Ftab_out$score1 - Ftab_out$score2
Ftab_out$cluster_cell<-ifelse(Ftab_out$score_diff >1, Ftab_out$cellType1, Ftab_out$max_col)
Ftab_out$cluster_cell1<-ifelse(Ftab_out$cluster_cell %in% "EPCAM","Epithelial",Ftab_out$cluster_cell)
Ftab_out$cluster_cell1<-ifelse(Ftab_out$cluster_cell1 %in% "PTPRC", "CD45_immune",Ftab_out$cluster_cell1)
Ftab_out$cluster_cell1<-ifelse(Ftab_out$cluster_cell1 %in% "COL1A1", "Fibroblast",Ftab_out$cluster_cell1)
Non_Imm<-c( "Pulmonary alveolar type I cells",
            "Pulmonary alveolar type II cells",
            "Clara cells","Airway epithelial cells",
            "Airway goblet cells","Ciliated cells",
            "Fibroblast","Stromal cells","Ionocytes",
            "Endothelial cells", "Endothelial",
            "lung epithelial", "fibroblast", "Epithelial", "Fibroblast")
Ftab_out$cluster_type_category <- ifelse(Ftab_out$cluster_cell1 %in% Non_Imm,"non-immune" ,"immune")
# Print the updated table
Ftab_out$cluster_assign<-rownames(Ftab_out)
Ftab_out$cluster_assign<- gsub("^g", "", rownames(Ftab_out))
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

fout <- c("cellranger_scran/NSCLC_output/integrated_annot.txt")
write.table(Ftab_out, fout, sep="\t", quote = F, row.names = F)


finAnnot <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/data1_integrated_annot.txt")
tmp <- read.delim(finAnnot,sep="\t",head=T, stringsAsFactors=FALSE)
# OR 
tmp<-Ftab_out
cellID_vec <- tmp$cluster_type_category # broad_classification replace the classification cell types and cluster_type_category cell types from the data  
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



#@@@@@@@@@@@@@@@@@@@@@ Figure 2A @@@@@@@@@@@@@@@@@@@@@@@@@@
Plot_2<-DimPlot(integrated_data2, reduction = "umap", group.by = "cluster_type_category", label = T)+NoLegend()+ 
  ggtitle("")  
Plot_2

#@@@@@@@@@@@@@@@@@@@@@ Figure 2B @@@@@@@@@@@@@@@@@@@@@@@@@@
Plot_3<-DimPlot(integrated_data2, reduction = "umap", group.by = "cellIDAll", label = T)+NoLegend()+ 
  ggtitle("")  
Plot_3

#@@@@@@@@@@@@@@@@@@@@  Supplementary Figure (SF2B) @@@@@@@@@@@@@@@@@@@@@@@@@@
# https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/ 
#install.packages("scCustomize")
#Plot_4
library(scCustomize)
scCustomize::Clustered_DotPlot(integrated_data2, features = gList, group.by = "Cell_type1",
                               plot_km_elbow = FALSE)
integrated_data2@meta.data$Dot_plot<-paste0(integrated_data2@meta.data$clustID,"_", integrated_data2@meta.data$cellIDAll)
pdf("Supp_Clustered_DotPlot.pdf", width = 11, height =6)
scCustomize::Clustered_DotPlot(integrated_data2, features = gList, group.by = "Cell_type1",
                               plot_km_elbow = FALSE)
dev.off()



#@@@@@@@@@@@@@@@@@@@ Supplementary Figure (SF2C) @@@@@@@@@@@@@@@@@@@@@@@@@@
Plot_5<-DimPlot(integrated_data2, reduction = "umap", group.by = "broad_classification", label = T)+NoLegend()+ 
  ggtitle("")  
# Save plot

#@@@@@@@@@@@@@@@@@  Supplementary Figure (SF2D) @@@@@@@@@@@@@@@@@@@@@@@@@@
Plot_6<-DimPlot(integrated_data2, reduction = "umap", group.by = "Cell_type1", label = T)+NoLegend()+ 
  ggtitle("")  
# Save plot

#@@@@@@@@@@@@@@@@@  Supplementary Figure (SF2F) @@@@@@@@@@@@@@@@@@@@@@@@@@
Plot_7A<-DimPlot(integrated_data2, reduction = "umap", group.by = "cellIDAll", split.by = "sample_type", label = T)+ggtitle("Batch wise samples: Non-Recurrence and Recurrence")+NoLegend()
Plot_7B<-DimPlot(integrated_data2, reduction = "umap", group.by = "cellIDAll", split.by = "sample", ncol = 4,label = F)+ggtitle("Individual samples")+NoLegend()
combined_plot<-Plot_7A | Plot_7B
print(combined_plot)

saveRDS(integrated_data2,file=c("NSCLC_output/Integrated_cd45pos_neg.rds"))

cd45pos_seurat <- subset(integrated_data2,idents='immune')
saveRDS(cd45pos_seurat,file=c("NSCLC_output/cd45pos_nsclc.rds"))

cd45neg_seurat <- subset(integrated_data2,idents='immune',invert=TRUE)
saveRDS(cd45neg_seurat,file=c("NSCLC_output/cd45neg_nsclc.rds"))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@# SingleR Annotation  Supplementary Figure (SF2E):
# https://bioconductor.org/books/3.17/SingleRBook/annotation-diagnostics.html
## Step 1:Automatic Annotation with SingleR all cells for immune and non-immune cell types 
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
inf <-c("NSCLC_output/Integrated_cd45pos_neg.rds")
data <-readRDS(inf)
# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))
# run SingleR (default mode) ---------
DefaultAssay(data) <- "RNA"
data[["RNA"]]<- JoinLayers(data[["RNA"]])
counts <- GetAssayData(data, layer = 'counts', assay = "RNA")
# Run SingleR
pred <- SingleR(test = counts,
                ref = ref,
                labels = ref$label.main)

pred
data$singleR.labels <- pred$labels[match(rownames(data@meta.data), rownames(pred))]
Plot_8<-DimPlot(data, reduction = 'umap', group.by = 'singleR.labels', label=T)
print(Plot_8)

# immune and non-immune cell types
# Seperate immune and non-immune cells based on mapped names.
data@meta.data$cell_type <- data$singleR.labels # Add immune/non-immune classification to metadata
data@meta.data$category <- ifelse(data@meta.data$cell_type %in% immune_cell_types, "Immune", "Non-Immune")
unique(data@meta.data$singleR.labels)
Plot_9_singleR<-DimPlot(data, reduction = "umap", group.by = "category", label = T) + ggtitle("SingleR: Immune and Non-Immune cells")+NoLegend()
# Save the plot.
