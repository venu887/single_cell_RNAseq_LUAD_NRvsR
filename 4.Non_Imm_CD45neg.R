#@@@@@@@@@@@@@@@@@@@@@@ 
# Non_Immune cells analysis 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## *******************************************************************
## Analyze non-immune cells: using 30PCs and annotate clusters
## *******************************************************************
rm(list = ls())
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(ggplot2)
library(NMF)
library(reshape2)
library(tidyr)
## loading Panglao database and keep major cell types
inp<-"/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/"
fin_db1 <- paste(inp,"PanglaoDB_markers_27_Mar_2020.tsv",sep="")
panglao <- read.delim(fin_db1,sep="\t",stringsAsFactors=FALSE)
listOrgan = c("Lungs","Immune system","Epithelium",
              "Connective tissue","Vasculature")
panglao = panglao[is.element(panglao$organ,listOrgan),]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
non_Imm<-c("Pulmonary alveolar type I cells",
           "Pulmonary alveolar type II cells",
           "Clara cells","Airway epithelial cells",
           "Airway goblet cells","Ciliated cells",
           "Ionocytes",
           "Fibroblasts","Stromal cells",
           "Endothelial cells")
tmp_non_IMM  = grepl("Hs",panglao$species) & 
  is.element(panglao$cell.type,non_Imm)
non_IMM_sub = panglao[tmp_non_IMM,2:3]
colnames(non_IMM_sub)=c("Symbol","cellType")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Step 1: load Seurat object
fin1 <- "NSCLC_output/cd45neg_nsclc.rds"
cd45neg <-readRDS(fin1)
DefaultAssay(cd45neg)<-"integrated"
cd45neg@meta.data <- cd45neg@meta.data[, -c(7:14)]
DimPlot(cd45neg, reduction = "umap")

## Step 2: scale and dimentional reduction
ndim=30  #*****
cd45neg <- FindNeighbors(cd45neg, dims = 1:ndim)
cd45neg <- FindClusters(cd45neg, resolution =0.3)
cd45neg <- RunUMAP(cd45neg,dims = 1:ndim)
# input for cluster annotation
target_seu <- cd45neg   
## Step 3: cluster annotation
## Step 3.1: enrichment approach
## Step 3.1.1: find cluster markers
allClust_mks <- FindAllMarkers(target_seu, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
write.table(allClust_mks, c("NSCLC_output/cd45neg_Allmarkers.txt"),
            quote = F, row.names = F)

## Step 3.1.2: enrichment approach by using panglao database
source("NSCLC_output/source_fxns4seurat.R")
enrichedData <- enrichScoreCalc_mt(target_seu,
                                   allClust_mks,non_IMM_sub)
enrichedData$score_Mt
## Find the top two highest scores and associated cell types for each cluster
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
## Step 3.2: using cannonical markers from either literature (FACS markers) 
## or overlap genes based on enrichment analysis
gList = c("COL1A1","COL1A2", #fibro
          "PECAM1","ECSCR",   #endo 
          "EPCAM","NKX2-1","SFTPC","SFTPA1",  ##AT2 or Clara
          "SCGB1A1","SCGB3A1",  #Clara or BASC
          "PDPN",   #AT1 
          "TFCP2L1",  #ionocyte, 
          "FOXJ1")    #ciliated
tmp_averages <- AverageExpression(target_seu,
                                  assays="RNA",
                                  features=gList)$RNA
tmp_zscores <- scale(log(tmp_averages+1))   
tmp_averages_df <- as.data.frame(tmp_zscores)
tmp_averages_df$Gene = rownames(tmp_averages)

mks_clust_avg <- tmp_averages_df[, c(ncol(tmp_averages_df), 1:(ncol(tmp_averages_df)-1))]
mks_clust_avg<-t(mks_clust_avg)
mks_clust_avg<-mks_clust_avg[-1,]
mks_clust_avg1<-as.data.frame(mks_clust_avg)
mks_clust_avg1$max_col <- colnames(mks_clust_avg1)[max.col(mks_clust_avg1, ties.method = 'first')]
## Step 3.3: export output from enrichment and marker expression
Ftab_out <- cbind(mks_clust_avg1,Ftab_candidate)
Ftab_out[, c(1:13, 15, 16)] <- lapply(Ftab_out[, c(1:13, 15, 16)], as.numeric)
Ftab_out$score_diff <- Ftab_out$score1 - Ftab_out$score2
Ftab_out$cluster_cell<-ifelse(Ftab_out$score_diff >1, Ftab_out$cellType1, Ftab_out$max_col)

## Step 3.4: Manually inspect the 3.3 output and assign cell type:
# (1) if difference between 1st and 2nd highest score > 1,
#         assign cell type associated with 1st highest score
# (2) else
#         assigne based on cell type expressing canonical markers
#         AT2 = EPCAM and NKX2-1
#         BASC = SCGB1A1 and SCGB1A1 and PECAM (intermediated level)
#         Clara = EPCAM, NKX2-1, SCGB1A1,SCGB1A1
# Create matrix of clusters X cell_Type for next step

###########Manual assessment 
# Define the cluster-to-cell type mapping
Ftab_out$cluster<-rownames(Ftab_out)
Ftab_out$cluster_name<-ifelse(Ftab_out$cluster %in% names(mapping), mapping, NA)
fout <- c("NSCLC_output/cd45neg_Ftab_out.txt")
write.table(Ftab_out, fout, sep="\t",
            quote = F, row.names = F)
# Assign cell types to meta data
names(Ftab_out)
num <- length(Ftab_out$cluster_name)
numbers <- seq(0, (num - 1))
Ftab_out$cluster_assign <- numbers
# Match cluster numbers from Ftab_out to integrated_snn_res.1 in allSam.integrated2 and assign cluster cell names
target_seu@meta.data$Cell_type <- ifelse(target_seu@meta.data$integrated_snn_res.0.3 %in% Ftab_out$cluster_assign, 
                                         Ftab_out$cluster_name[match(target_seu@meta.data$integrated_snn_res.0.3, Ftab_out$cluster_assign)], 
                                         NA)
# Combine cluster number and cluster name into Cell_type
target_seu@meta.data$Cell_type2 <- ifelse(target_seu@meta.data$integrated_snn_res.0.3 %in% Ftab_out$cluster_assign, 
                                          paste0(Ftab_out$cluster_cell[match(target_seu@meta.data$integrated_snn_res.0.3, Ftab_out$cluster_assign)],"-","C", target_seu@meta.data$integrated_snn_res.0.3), 
                                          NA)
# Combine cluster number and cluster name into Cell_type
target_seu@meta.data$Cell_type1 <- ifelse(target_seu@meta.data$integrated_snn_res.0.3 %in% Ftab_out$cluster_assign, 
                                          paste0(Ftab_out$cluster_name[match(target_seu@meta.data$integrated_snn_res.0.3, Ftab_out$cluster_assign)],"-","C", target_seu@meta.data$integrated_snn_res.0.3 ), 
                                          NA)

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 2E @@@@@@@@@@@@@@@@@@@@@@@@
Plot_1 <- DimPlot(target_seu, reduction = "umap", group.by = 'Cell_type',label = TRUE, label.size = 4)+NoLegend()+ggtitle(label = NULL)



#@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 2F @@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@ DOT PLOT using all Top Marker genes
#@@@@@@@@. Updated
non_immune_markers <- list(
  Fibro = c("COL1A1", "COL1A2"),
  Endo = c("PECAM1", "ECSCR"),
  AT2_Clara = c("EPCAM", "CDH1", "NKX2-1", "SFTPC", "SFTPA1"),
  Clara_BASC = c("SCGB1A1", "SCGB3A1"),
  AT1 = c("PDPN"),
  Ciliated = c("FOXJ1")
)
markers_df <- stack(non_immune_markers)
colnames(markers_df) <- c("Symbol", "CellType")
valid_markers_df <- markers_df[markers_df$Symbol %in% rownames(cd45neg@assays$RNA$counts), ]
valid_markers_df <- valid_markers_df[order(valid_markers_df$CellType), ]
celltype_labels <- unique(valid_markers_df$CellType)  
gene_list <- split(valid_markers_df$Symbol, valid_markers_df$CellType)  
valid_markers <- unlist(gene_list)
names(valid_markers) <- rep(celltype_labels, times = lengths(gene_list)) 
valid_markers <- valid_markers[!duplicated(valid_markers)]
Plot_2 <- DotPlot(cd45neg, 
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
Plot_1

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 3C @@@@@@@@@@@@@@@@@@@@@@@@
Plot_3<-DimPlot(target_seu , reduction = "umap", split.by ='sample_type')+ggtitle(label = NULL)+theme(legend.position = c(0.2, 0.2))


#@@@@@@@@@@@@@@@@@@@@@@@@@@ Figure 3D @@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@ Proportion test 
#install.packages("devtools") # https://github.com/rpolicastro/scProportionTest
#devtools::install_github("rpolicastro/scProportionTest")
#Plot_4
library("scProportionTest")
prop_test <- sc_utils(target_seu)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "Cell_type1",
  sample_1 = "Non_Rec", 
  sample_2 = "Rec",
  sample_identity = "sample_type" # Correct column name here
)
permutation_plot(prop_test)
# Save plot










#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure SF6A @@@@@@@@@@@@@@@@@@@@@@@@
#DimPlot( target_seu , reduction = "umap", group.by ='Cell_type2', label = TRUE, label.size = 6)
Plot_5<-DimPlot( target_seu , reduction = "umap",label = TRUE, label.size = 4)+NoLegend()+ggtitle(label = NULL)
Plot_6 <- DimPlot(target_seu, reduction = "umap", group.by = 'sample_type',label = TRUE, label.size = 4)+NoLegend()+ggtitle(label = NULL)
# Save bothe the plots Plot_5+Plot_6
saveRDS(target_seu,file=c("NSCLC_output/cd45negc_celltype_data.rds"))





#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure SF6B @@@@@@@@@@@@@@@@@@@@@@@@
# Remove all and re-upload to reduce the burden in R environment
library(dplyr)
library(knitr)
library(ggplot2)
library(tidyr)
file=c("NSCLC_output/cd45negc_celltype_data.rds")
cd45neg<-readRDS(file) 
target_seu<-cd45neg
cont_table <- table(target_seu$Cell_type1, cd45neg$sample) 
cont_df <- as.data.frame(cont_table)
colnames(cont_df) <- c("Cell_type", "Sample", "Count")

contingency_df <- cont_df %>%
  group_by(Cell_type) %>%
  mutate(Proportion = Count / sum(Count))
xx<-c("Fibroblasts-C6", "Endothelial cells-C4", "AT2-C2", "Clara-C5", "Fibroblasts-C9", 
      "Ciliated cells-C7", "AT2-C0", "BASC-C10", "AT2-C1", "Clara-C3", "BASC-C8")
contingency_df$Cell_type <- factor(contingency_df$Cell_type, levels = rev(xx))
Plot_7<-ggplot(contingency_df, aes(x = Cell_type, y = Proportion, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = NULL, y = "Proportion", fill = "Cell_type") +
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black",angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text(color = "black"))+coord_flip()
Plot_7
# save this plot

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure SF6C @@@@@@@@@@@@@@@@@@@@@@@@
gList = c("COL1A1","COL1A2", #fibro
          "PECAM1","ECSCR",   #endo 
          "EPCAM", "CDH1","NKX2-1","SFTPC","SFTPA1",  ##AT2 or Clara
          "SCGB1A1","SCGB3A1",  #Clara or BASC
          "PDPN", #AT1 
          "FOXJ1")    #ciliated
all_markers <- unique(unlist(gList))
valid_markers <- all_markers[all_markers %in% rownames(cd45neg@assays$RNA$counts)]
dot_data <- DotPlot(cd45neg, features = valid_markers, group.by = 'Cell_type1')
Plot_8<-dot_data+RotatedAxis()+ scale_size(range = c(1,9))
Plot_8 # save this plot




#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure SF6D @@@@@@@@@@@@@@@@@@@@@@@@
# Create the contingency table
Table_nonImm <- table(target_seu@meta.data$Cell_type1, target_seu@meta.data$sample_type)
Table_nonImm <- as.data.frame(Table_nonImm)
# Extract relevant rows and columns
Rec <- Table_nonImm[12:22, 3]
Table_nonImm <- Table_nonImm[-c(12:22), ]
Table_nonImm <- cbind(Table_nonImm, Rec)
Table_nonImm <- Table_nonImm[, -2]
colnames(Table_nonImm) <- c("Cluster", "NR_cell", "R_cells")
# Calculate total and proportions
n1 <- sum(Table_nonImm$NR_cell)
n2 <- sum(Table_nonImm$R_cells)
Table_nonImm$Total <- Table_nonImm$NR_cell + Table_nonImm$R_cells
Table_nonImm$NR_prapor <- Table_nonImm$NR_cell / n1
Table_nonImm$R_prapor <- Table_nonImm$R_cells / n2
Table_nonImm$FC <- Table_nonImm$NR_prapor / Table_nonImm$R_prapor
# Prepare the data for Fisher's test
data <- Table_nonImm[, 1:3]
total_NR <- sum(data$NR_cell)
total_R <- sum(data$R_cells)
# Function to calculate Fisher's p-value
fisher_test_p_value <- function(nr, r) {
  table <- matrix(c(nr, r, total_NR - nr, total_R - r), nrow = 2)
  test <- fisher.test(table)
  return(test$p.value)
}
# Apply Fisher's test to each row
data$p_value <- mapply(fisher_test_p_value, data$NR_cell, data$R_cells)
data$p_value1 <- format(data$p_value, scientific = TRUE, digits = 4)

# Add p-value and scientific p-value to Table_nonImm
Table_nonImm$p_value <- data$p_value
Table_nonImm$p_value_scientific <- data$p_value1
Table_nonImm$p_value_adj <- p.adjust(Table_nonImm$p_value, method = "BH")

# Add rounded columns for NR proportion, R proportion, and Fold Change (FC)
Tab_round <- round(Table_nonImm[, c("NR_prapor", "R_prapor", "FC")], 5)
names(Tab_round) <- c("Round5_NR_prapor", "Round5_R_prapor", "Round5_FC")
Table_nonImm1 <- cbind(Table_nonImm, Tab_round)

# Summarize counts by Cell_type1 and sample, and calculate proportions
sample_wise_counts <- target_seu@meta.data %>%
  count(Cell_type1, sample, name = "Freq") %>%
  pivot_wider(names_from = sample, values_from = Freq, values_fill = 0) %>%
  mutate(across(-Cell_type1, ~ .x / sum(.x), .names = "Prop_{col}")) %>%
  rowwise() %>%
  mutate(
    Mean_R_prop = mean(c_across(starts_with("Prop_LC") & ends_with("_R")), na.rm = TRUE),
    Mean_NR_prop = mean(c_across(starts_with("Prop_LC") & ends_with("_NR")), na.rm = TRUE),
    FC_mean_prop = ifelse(Mean_NR_prop > 0, Mean_R_prop / Mean_NR_prop, NA),
    Enriched_in = case_when(
      FC_mean_prop > 1 ~ "R",
      FC_mean_prop < 1 ~ "NR",
      TRUE ~ "None"
    )
  ) %>%
  ungroup()

# Merge the sample-wise counts with Table_nonImm1
Table_nonImm_samplewise <- merge(Table_nonImm1, sample_wise_counts, by.x = "Cluster", by.y = "Cell_type1", all = TRUE)

xx<-c("Cluster", "Prop_LC104_R","Prop_LC115_NR","Prop_LC221_NR","Prop_LC38_NR","Prop_LC52_R","Prop_LC63_NR","Prop_LC71_R", "p_value_adj")
nonImm<-Table_nonImm_samplewise[,xx]
long_nonImm <- melt(nonImm, id.vars = c("Cluster", "p_value_adj"),
                    variable.name = "Sample", value.name = "Proportion")
long_nonImm$Group <- ifelse(grepl("_R$", long_nonImm$Sample), "R", "NR")
long_nonImm <- long_nonImm %>%
  mutate(Fisher_label = case_when(
    p_value_adj < 0.001 ~ "***",
    p_value_adj< 0.01  ~ "**",
    p_value_adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

Plot_9<-ggplot(long_nonImm, aes(x = Group, y = Proportion, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  
  scale_color_manual(values = c("R" = "darkorange", "NR" = "purple")) +
  facet_wrap(~ Cluster, scales = "free_y", ncol = 4) +
  geom_text(data = distinct(long_nonImm, Cluster, Fisher_label), 
            aes(x = 1.5, y = max(long_nonImm$Proportion, na.rm = TRUE) * 0.1, 
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
  labs(x = "Group", y = "Proportion", title = "Non-Immune Proportions Across Clusters (Fisher test)")
Plot_9 # Save the figure.

#@@@@@@@@@@@@@@@@@@@@@@@@@@ Supplementary Figure SF6E @@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@ AT2 cells
target_seu_subset <- subset(target_seu, idents = c(0, 1, 2))
Plot_10<-VlnPlot(target_seu_subset, 
            features = c("SFTPC", "NKX2-1"), 
            #      group.by = "sample_type", 
            split.by = "Cell_type1") 

Plot_10
Plot_11<-VlnPlot(target_seu_subset, 
            features = c("SFTPC", "NKX2-1"), 
            group.by = "sample_type" )
Plot_11
Plot_12<-Plot_11|Plot_12
