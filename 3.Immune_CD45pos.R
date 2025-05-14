#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## *******************************************************************
## Analyze immune cells: Annotate clusters
## *******************************************************************
rm(list=ls())
library(SingleCellExperiment)
library(Seurat)
library(dplyr) 
library(ggplot2)
library(NMF)
library(reshape2)
source("NSCLC_output/source_fxns4seurat.R")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 1: load Seurat object
set.seed(42)
fin2 <- c("NSCLC_output/cd45pos_nsclc.rds")
cd45pos <-readRDS(fin2)
cd45pos@meta.data <- cd45pos@meta.data[, -c(7:14)] # Remove excess data from previous step 
DimPlot(cd45pos, reduction = "umap")
DefaultAssay(cd45pos)<-"integrated"

## Step 2: scale and dimentional reduction
ndim=30  #Good number of PC
cd45pos <- FindNeighbors(cd45pos, dims = 1:ndim)
# # Find optimal resolution
# library(clustree)
# resolutions <- seq(0, 1.5, 0.1)
# n.dims <- 30
# combined <- FindClusters(cd45pos, reduction.type = "cca.aligned",
#                          dims.use = 1:n.dims, resolution = resolutions)
# clustree(combined)
# rm(combined)

cd45pos <- FindClusters(cd45pos)
cd45pos <- RunUMAP(cd45pos, dims = 1:ndim)
allClust_mks <- FindAllMarkers(cd45pos, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.25)
cd45pos
write.table(allClust_mks, c("NSCLC_output/cd45pos_nsclc_allClust.txt"),
            quote = F, row.names = F)

#@@@@@@@@@@@ Supplementary Figure 2A Raw UMAP @@@@@@@@@@@ 
Imm_plot_1<-DimPlot(cd45pos,  reduction = "umap", label = T)
print(Imm_plot_1)
# Save the figure

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 2: Calculate Enrichment scores and average expression using Panglao database and keep major cell types
inp<-"NSCLC_output/"
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
fout <- c("NSCLC_output/cd45pos_Ftab_out.txt")
write.table(Ftab_out, fout, sep="\t",
            quote = F, row.names = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 3: Update cell types to the immune data
fin3 <- c("NSCLC_output/cd45pos_Ftab_out.txt") 
cellAnnotInfo <- read.delim(fin3,sep="\t",header=T,
                            stringsAsFactors=FALSE)
num <- length(Ftab_out$cluster1)
numbers <- seq(0, (num - 1))
Ftab_out$cluster_assign <- numbers
DefaultAssay(cd45pos)<-"RNA"
# Match cluster numbers from Ftab_out to integrated_snn_res.1 in allSam.integrated2 and assign cluster cell names
cd45pos@meta.data$Cell_type <- ifelse(cd45pos@meta.data$seurat_clusters %in% Ftab_out$cluster_assign, 
                                             Ftab_out$cluster1[match(cd45pos@meta.data$seurat_clusters, Ftab_out$cluster_assign)], 
                                             NA)
cd45pos@meta.data$Cell_type1 <- paste0( cd45pos@meta.data$Cell_type,"_","C",cd45pos@meta.data$integrated_snn_res.0.8)
file=c("cd45pos_nsclc_Immune_cells.rds")
saveRDS(cd45pos,file) # Save the data in case if you lost

#@@@@@@@@@@@@@@@ Plot-2C @@@@@@@@@@@@@@@@@@@@@@
p2<-DimPlot(cd45pos, reduction = "umap", group.by = 'Cell_type1', label = TRUE, label.size = 6, repel = T) 
p2
# Save the figure 

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure 2B @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
plot1<-DimPlot(cd45pos, group.by = 'Cell_type1', split.by = "sample_type", label = TRUE, label.size = 6, repel = T) +NoLegend()
plot1
# Save the figure 

library(viridis)
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure 2D @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
all_markers <- unique(unlist(markers))
DefaultAssay(cd45pos) <- "RNA"
valid_markers <- all_markers[all_markers %in% rownames(cd45pos@assays$RNA$counts)]
dot_data <- DotPlot(cd45pos, features = valid_markers, group.by = 'Cell_type1')
plot1<-dot_data+RotatedAxis()+ scale_size(range = c(1,9))
#save the plot.

#@@@@@@@@@@@@@@@ Plot-2D @@@@@@@@@@@@@@@@@@@@@@
all_markers <- unique(unlist(markers))
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
# save the plot

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@# Save data 
saveRDS(cd45pos,file=c("Immune_cells_data/cd45pos_nsclc_Immune_cells.rds"))
# save the data according to the cell type of clusters 
Macrophges <- subset(cd45pos,subset = Cell_type == c("Mac|DC","Mac", "Mono", "Neutrophil", "Prolifirating Mac", "Prolif Mac|DC"))
saveRDS(Macrophges,file=c("Immune_cells_data/Macrophges.rds"))


# Subset the Seurat object for T cells (CD4, CD8, and γδT cells)
T_cells <- subset(cd45pos, subset = Cell_type %in% c("CD4", "CD8", "GaDelT"))
saveRDS(T_cells,file=c("Immune_cells_data/T_cells.rds"))


NK_NKT<-subset(cd45pos, subset = Cell_type %in% c("NK cells" ,  "NKT"))
saveRDS(NK_NKT,file=c("Immune_cells_data/NK_NKT.rds"))

B_cells<-subset(cd45pos, subset = Cell_type %in% c("Plasma cells", "B cells" ))
saveRDS(B_cells,file=c("Immune_cells_data/B_cells.rds"))

Mast_NA_cells<-subset(cd45pos, subset = seurat_clusters %in% c(17,20))
saveRDS(Mast_NA_cells,file=c("Immune_cells_data/Mast_NA_cells.rds"))

DC<-subset(cd45pos, subset = Cell_type %in% c("DC", "pDC"))
saveRDS(DC,file=c("Immune_cells_data/DC_cells.rds"))


# Bingo with immune cell type assessment

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@. OTHER PLOTS CODE 
#@@@@@@@@@@@@@@@@@@@@@@ Method-1 FDR and Fisher p-value, calculating the proportion with fisher extract test between recurrence vs non-recurrence samples
# Fold change calculations in each individual clusters
# calculate total cells in Non-Rec and Rec samples
# count the number of cells present in each cluster
# Divide the total cells  by specific cluster types to get proportion 
# Divide the proportion of Non-Rec with Rec will give us Fold Change
library(SingleCellExperiment)
library(Seurat)
library(dplyr) 
library(ggplot2)
library(NMF)
library(reshape2)
library(limma)
library(speckle)
library(speckle)
library(knitr)
library(tidyr)
library(ggpubr)
library(ggrepel)
library(stringr)
file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/cd45pos_nsclc_Immune_cells.rds")
cd45pos<-readRDS(file)
cd45pos

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure 2C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# proportion plots of each samples in each immune cell type
#@@@@@@@@@@@@@@@@@@@@@@@@@@ 
names(cd45pos@meta.data)
cont_table <- table(cd45pos$Cell_type1, cd45pos$sample) # prop1: sample and prop2: sample_type
cont_df <- as.data.frame(cont_table)
colnames(cont_df) <- c("Cell_type", "Sample", "Count")

contingency_df <- cont_df %>%
  group_by(Cell_type) %>%
  mutate(Proportion = Count / sum(Count))
xx<-c("Plasma cells_C12", "B cells_C2","CD4_C10", "GaDelT_C14", "B cells_C13",
      "CD8_C1", "NK cells_C4", "CD4_C6", "pDC_C19", "CD4_C7", "CD8_C0", "Plasma cells_C9",
      "Mast_C17", "CD8_C15", "Mono_C3", "Neutrophil_C5","DC_C18","Mac_C8", "Prolifirating Mac_C16",
      "Mac|DC_C11", "Prolif Mac|DC_C20")
# Reverse the order of 'xx' and set it as levels for 'Cell_type'
contingency_df$Cell_type <- factor(contingency_df$Cell_type, levels = rev(xx))

# Create the plot with the reversed x-axis order
plot_prapo <- ggplot(contingency_df, aes(x = Cell_type, y = Proportion, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = NULL, y = "Proportion", fill = "Cell_type") +
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black", angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(color = "black")) +
  coord_flip()
# Save the plot

# Generate final_data for calculation of proportion for each cell type in group wise
final_data <- cd45pos@meta.data %>%
  count(Cell_type1, sample, name = "Freq") %>%
  pivot_wider(names_from = sample, values_from = Freq, values_fill = 0) %>%
  # Calculate proportions
  mutate(across(-Cell_type1, ~ .x / sum(.x), .names = "Prop_{col}")) %>%
  # Calculate mean proportions for R and NR groups
  rowwise() %>%
  mutate(
    Mean_R_prop = mean(c_across(starts_with("Prop_LC") & ends_with("_R")), na.rm = TRUE),
    Mean_NR_prop = mean(c_across(starts_with("Prop_LC") & ends_with("_NR")), na.rm = TRUE),
    FC_mean_prop = ifelse(Mean_NR_prop > 0, Mean_R_prop / Mean_NR_prop, NA),
    Enriched_in = case_when(
      FC_mean_prop > 1 ~ "R",
      FC_mean_prop < 1 ~ "NR",
      TRUE ~ "None"
    ),
    R_sum_prop = sum(c_across(starts_with("Prop_LC") & ends_with("_R")), na.rm = TRUE),
    NR_sum_prop = sum(c_across(starts_with("Prop_LC") & ends_with("_NR")), na.rm = TRUE),
    FC_sum = ifelse(NR_sum_prop > 0, R_sum_prop / NR_sum_prop, NA),
    Enriched_in_sum = case_when(
      FC_sum > 1 ~ "R",
      FC_sum < 1 ~ "NR",
      TRUE ~ "None"
    )
  ) %>%
  # Summarize counts for R and NR
  mutate(
    R_sum = sum(c_across(starts_with("LC") & ends_with("_R")), na.rm = TRUE),
    NR_sum = sum(c_across(starts_with("LC") & ends_with("_NR")), na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate total counts for R and NR groups
total_R <- sum(final_data$R_sum, na.rm = TRUE)
total_NR <- sum(final_data$NR_sum, na.rm = TRUE)

# Step 2: Perform Fisher's Exact Test for each cell type
final_data <- final_data %>%
  rowwise() %>%
  mutate(
    Fisher_p_value = fisher.test(matrix(
      c(R_sum, NR_sum, total_R - R_sum, total_NR - NR_sum), 
      nrow = 2
    ))$p.value,
    Significant = Fisher_p_value < 0.05 & abs(log2(FC_mean_prop)) > 1
  ) %>%
  ungroup()

final_data <- final_data %>%
  rowwise() %>%
  mutate(
    Fisher_p_value = fisher.test(matrix(
      c(R_sum, NR_sum, total_R - R_sum, total_NR - NR_sum), 
      nrow = 2
    ))$p.value,
    Significant = Fisher_p_value < 0.05 & abs(log2(FC_mean_prop)) > 1
  ) %>%
  ungroup()
# Adjust p-values using BH method
final_data$Fisher_p_value_adj <- p.adjust(final_data$Fisher_p_value, method = "BH")
final_data$Significant_adj <- final_data$Fisher_p_value_adj < 0.05 & abs(log2(final_data$FC_mean_prop)) > 1
fin1 <- c("Immune_cells_data/Immune_Cells_FC_P_val.csv")
write.csv(final_data, fin1)


#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure 2E @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
names(final_data)
final_data1<-final_data
xx<-c("Prop_LC104_R", "Prop_LC115_NR","Prop_LC221_NR", "Prop_LC38_NR", "Prop_LC52_R", "Prop_LC63_NR", "Prop_LC71_R")
proportions <- final_data1[,xx]
proportions$Cluster <- final_data1$Cell_type1
long_data <- melt(proportions)
colnames(long_data) <- c("Cluster", "Sample", "Proportion")
long_data$Group <- ifelse(grepl("_R$", long_data$Sample), "R", "NR")
# Merge Fisher's p-values with long_data for annotation
long_data <- left_join(long_data, final_data %>% select(Cell_type1, Fisher_p_value_adj), 
                       by = c("Cluster" = "Cell_type1"))

# Format Fisher’s p-values for display
long_data <- long_data %>%
  mutate(Fisher_label = case_when(
    Fisher_p_value_adj < 0.001 ~ "***",
    Fisher_p_value_adj < 0.01  ~ "**",
    Fisher_p_value_adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

pdf("Imm_Prop_boxplot.pdf", width = 10, height = 11)
ggplot(long_data, aes(x = Group, y = Proportion, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) + 
  scale_color_manual(values = c("R" = "darkorange", "NR" = "purple")) +
  facet_wrap(~ Cluster, scales = "free_y", ncol = 4) +
  geom_text(data = distinct(long_data, Cluster, Fisher_label), 
            aes(x = 1.5, y = max(long_data$Proportion, na.rm = TRUE) * 0.1, 
                label = Fisher_label), 
            inherit.aes = FALSE, size = 5) +  
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Group", y = "Proportion", title = "Immune Proportions Across Clusters (Fisher extract)")
dev.off()




