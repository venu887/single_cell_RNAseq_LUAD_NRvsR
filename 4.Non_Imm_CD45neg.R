#@@@@@@@@@@@@@@@@@@@@@@ 
# Non_Immune cells analysis 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## *******************************************************************
## Analyze non-immune cells: using 30PCs and annotate clusters
## Input from Immune_NonImmune_split.R
## *******************************************************************
# Marker genes: https://www.nature.com/articles/s41586-024-07113-9#Sec2
# https://desailab.stanford.edu/sites/g/files/sbiybj24296/files/media/file/li2.pdf 
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
# non Immune cells that used in previous step
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
fin1 <- "/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45neg_nsclc.rds"
cd45neg <-readRDS(fin1)
dim(cd45neg)
cd45neg
DefaultAssay(cd45neg)<-"integrated"
cd45neg@meta.data <- cd45neg@meta.data[, -c(7:14)]
names(cd45neg@meta.data)
View(cd45neg@meta.data)
DimPlot(cd45neg, reduction = "umap")

## Step 2: scale and dimentional reduction
ndim=30  #*****
cd45neg <- FindNeighbors(cd45neg, dims = 1:ndim)
# Find optimal resolution
# library(clustree)
# resolutions <- seq(0, 1.5, 0.1)
# n.dims <- 20
# combined <- FindClusters(cd45neg, reduction.type = "cca.aligned",
#                          dims.use = 1:n.dims, resolution = resolutions)
# 
# clustree(combined)
# rm(combined)

cd45neg <- FindClusters(cd45neg, resolution =0.3)
cd45neg <- RunUMAP(cd45neg,dims = 1:ndim)

DimPlot(cd45neg, label = TRUE)


# input for cluster annotation
target_seu <- cd45neg   

## Step 3: cluster annotation
## Step 3.1: enrichment approach
## Step 3.1.1: find cluster markers
allClust_mks <- FindAllMarkers(target_seu, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)


## export data
write.table(allClust_mks, c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45neg_Allmarkers.txt"),
            quote = F, row.names = F)

## Step 3.1.2: enrichment approach by using panglao database
source("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/source_fxns4seurat.R")
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
  #        "MUC5AC","MUC16",   #Goblet
          "PDPN",   #AT1 
          "TFCP2L1",  #ionocyte, 
          "FOXJ1")    #ciliated


tmp_averages <- AverageExpression(target_seu,
                                  assays="RNA",
                                  features=gList)$RNA
# tmp_averages <- t(tmp_averages)
tmp_zscores <- scale(log(tmp_averages+1))   ##high dynamic rank --> log transform to norm distribution
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

# Display the result
Ftab_out$cluster_cell<-ifelse(Ftab_out$score_diff >1, Ftab_out$cellType1, Ftab_out$max_col)

unique(Ftab_out$cluster_cell)
c("COL1A1","COL1A2", #fibro
  "PECAM1","ECSCR",   #endo 
  "EPCAM","NKX2-1","SFTPC","SFTPA1",  ##AT2 or Clara
  "SCGB1A1","SCGB3A1",  #Clara or BASC
  "MUC5AC","MUC16",   #Goblet
  "PDPN",   #AT1 
  "TFCP2L1",  #ionocyte, 
  "FOXJ1")    #ciliated

## Step 3.4: Manually inspect the 3.3 output and assign cell type:
# (1) if difference between 1st and 2nd highest score > 1,
#         assign cell type associated with 1st highest score
# (2) else
#         assigne based on cell type expressing canonical markers
#         AT2 = EPCAM and NKX2-1
#         BASC = SCGB1A1 and SCGB1A1 and PECAM (intermediated level)
#         Clara = EPCAM, NKX2-1, SCGB1A1,SCGB1A1
# Create matrix of clusters X cell_Type for next step
Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_cell %in% "TFCP2L1" ,"AT2",Ftab_out$cluster_cell)
Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "SFTPC" ,"AT2",Ftab_out$cluster_name)
Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "COL1A2" ,"Fibroblasts",Ftab_out$cluster_name)
Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "SCGB3A1" ,"BASC",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "EPCAM" ,"Clara",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "PECAM1" ,"Endothelial cells",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "MUC16" ,"Goblet",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "SFTPA1" ,"AT2",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "PDPN" ,"AT1",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "MUC5AC" ,"Goblet",Ftab_out$cluster_name)
# Ftab_out$cluster_name<-ifelse(Ftab_out$cluster_name %in% "Pulmonary alveolar type II cells" ,"AT2",Ftab_out$cluster_name)

unique(Ftab_out$cluster_name)
###########Manual assessment 
# Define the cluster-to-cell type mapping
Ftab_out$cluster<-rownames(Ftab_out)
mapping <- c(
  "g0" = "AT2",
  "g1" = "AT2",
  "g2" = "AT2",
  "g3" = "Clara",
  "g4" = "Endothelial cells",
  "g5" = "Clara",
  "g6" = "Fibroblasts",
  "g7" = "Ciliated cells",
  "g8" = "BASC",
  "g9" = "Fibroblasts",
  "g10" = "BASC"
)
Ftab_out$cluster_name<-ifelse(Ftab_out$cluster %in% names(mapping), mapping, NA)

fout <- c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45neg_Ftab_out.txt")
write.table(Ftab_out, fout, sep="\t",
            quote = F, row.names = F)

# Assign cell types to meta data
names(Ftab_out)
num <- length(Ftab_out$cluster_name)
numbers <- seq(0, (num - 1))
Ftab_out$cluster_assign <- numbers
names(target_seu@meta.data)

View(target_seu@meta.data)

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

View(target_seu@meta.data)
table(cd45neg@meta.data$integrated_snn_res.0.3)
#DimPlot( target_seu , reduction = "umap", group.by ='Cell_type2', label = TRUE, label.size = 6)
plot1<-DimPlot( target_seu , reduction = "umap",label = TRUE, label.size = 4)+NoLegend()+ggtitle(label = NULL)
plot2<-DimPlot( target_seu , reduction = "umap", group.by ='sample_type')+ggtitle(label = NULL)+theme(legend.position = c(0.2, 0.2))
plot3 <- DimPlot(target_seu, reduction = "umap", group.by = 'Cell_type',label = TRUE, label.size = 4)+NoLegend()+ggtitle(label = NULL)
plot1
plot2
plot3
# theme(legend.position = c(0.76, 0.3)) # Adjust the coordinates as needed

#plot3
DefaultAssay(target_seu)<-"RNA"
gList = c("COL1A1","COL1A2", #fibro
          "PECAM1","ECSCR",   #endo 
          "EPCAM", "CDH1","NKX2-1","SFTPC","SFTPA1",  ##AT2 or Clara
          "SCGB1A1","SCGB3A1",  #Clara or BASC
          "PDPN", #AT1 
          "FOXJ1")    #ciliated

library(patchwork)
p<-(plot1+plot3+plot2)
p
pdf("/home/u251079/scRNA/Lung_scRNA/Non_Imm_Umap.pdf", width = 12, height =5)
print(p)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@ DOT PLOT using all Top Marker genes
Symbols_and_types <- non_IMM_sub
DefaultAssay(cd45neg)<-"RNA"
valid_Symbols <- Symbols_and_types$Symbol[Symbols_and_types$Symbol %in% rownames(cd45neg@assays$RNA$counts)]
Symbols_and_types <- Symbols_and_types %>%
  filter(Symbol %in% valid_Symbols)

dot_data <- DotPlot(cd45neg, features = unique(valid_Symbols), group.by = "integrated_snn_res.0.3")
dot_data_df <- dot_data$data
filtered_dot_data <- dot_data_df %>%
  filter(pct.exp > 30 & avg.exp.scaled > 0.2)

filtered_dot_data <- filtered_dot_data %>%
  distinct(features.plot, id, .keep_all = TRUE)

filtered_dot_data <- filtered_dot_data %>%
  left_join(Symbols_and_types, by = c("features.plot" = "Symbol"))

dot_plot <- ggplot(filtered_dot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Top marker genes per.exp >30 & Avg.exp >0.5",
    x = "Symbols",
    y = "Cell Clusters",
    size = "Percent Expressed",
    color = "Avg. Scaled Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ cellType, scales = "free_x", nrow = 1)
dot_plot   

file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45negc_celltype_data.rds")
saveRDS(target_seu,file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45negc_celltype_data.rds"))







file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45negc_celltype_data.rds")
cd45neg<-readRDS(file) 
cd45neg
View(cd45neg@meta.data)
unique(cd45neg$Cell_type1)
unique(cd45neg$Cell_type2)
unique(cd45neg$Cell_type)

#+++++++++++ Plot UMAP
DefaultAssay(cd45neg)<-"integrated"
p1<-DimPlot(cd45neg, group.by = 'Cell_type1', split.by = "sample_type", label = TRUE, label.size = 6, repel = T) +NoLegend()
p1
pdf("/home/u251079/scRNA/Lung_scRNA/4_NonImm_cell_UMAP1.pdf", width = 12, height =5)
print(p1)
dev.off()



gList = c("COL1A1","COL1A2", #fibro
          "PECAM1","ECSCR",   #endo 
          "EPCAM", "CDH1","NKX2-1","SFTPC","SFTPA1",  ##AT2 or Clara
          "SCGB1A1","SCGB3A1",  #Clara or BASC
          "PDPN", #AT1 
          "FOXJ1")    #ciliated
#@@@@@@@@. PLOT.1
all_markers <- unique(unlist(gList))
valid_markers <- all_markers[all_markers %in% rownames(cd45neg@assays$RNA$counts)]
dot_data <- DotPlot(cd45neg, features = valid_markers, group.by = 'Cell_type1')
plot1<-dot_data+RotatedAxis()+ scale_size(range = c(1,9))
plot1
pdf("/home/u251079/scRNA/Lung_scRNA/4_Marker_gene_Dot_Plot.pdf", width = 8, height =4)
print(plot1)
dev.off()

#@@@@@@@@. PLOT.2
library(viridis)
all_markers <- unique(unlist(markers))
valid_markers <- all_markers[all_markers %in% rownames(cd45neg@assays$RNA$counts)]
plot2<-DotPlot(cd45neg, features =valid_markers, group.by = 'Cell_type1',
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
plot2 <- DotPlot(cd45neg, 
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

pdf("/home/u251079/scRNA/Lung_scRNA/4_DotPlot_Non_imm.pdf", width = 10, height =4.5)
print(plot2)
dev.off()




#@@@@@@@@@@@@@@@@@@@@@@@@@@ Proportion plot 
library(dplyr)
library(knitr)
library(ggplot2)
library(tidyr)
target_seu<-cd45neg
View(target_seu@meta.data)
names(target_seu@meta.data)
cont_table <- table(target_seu$Cell_type1, cd45neg$sample_type) # prop1: sample and prop2: sample_type
cont_df <- as.data.frame(cont_table)
colnames(cont_df) <- c("Cell_type", "Sample", "Count")

contingency_df <- cont_df %>%
  group_by(Cell_type) %>%
  mutate(Proportion = Count / sum(Count))
xx<-c("Fibroblasts-C6", "Endothelial cells-C4", "AT2-C2", "Clara-C5", "Fibroblasts-C9", 
      "Ciliated cells-C7", "AT2-C0", "BASC-C10", "AT2-C1", "Clara-C3", "BASC-C8")
contingency_df$Cell_type <- factor(contingency_df$Cell_type, levels = rev(xx))
plot_prapo<-ggplot(contingency_df, aes(x = Cell_type, y = Proportion, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = NULL, y = "Proportion", fill = "Cell_type") +
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black",angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text(color = "black"))+coord_flip()

plot_prapo1<- plot_prapo # Sample
plot_prapo2<-plot_prapo # sample_type
plot_save<-plot_prapo1+plot_prapo2
print(plot_save)
pdf("/home/u251079/scRNA/Lung_scRNA/4_NI_Prap.pdf", width = 15, height =5)
print(plot_save)
dev.off()


#@@@@@@@@@@@@@@@@@@  Compute contingency table
cont_table <- table(cd45neg$Cell_type1, cd45neg$sample)
cont_df <- as.data.frame(cont_table)
colnames(cont_df) <- c("Cell_type", "Sample", "Count")
# Compute total non-immune cells per sample
total_cells_df <- cont_df %>%
  group_by(Sample) %>%
  summarize(Total_NonImmune_Cells = sum(Count))
# Compute overall total for normalization
total_cells <- sum(total_cells_df$Total_NonImmune_Cells)
# Compute fraction for each cell type
fraction_df <- cont_df %>%
  left_join(total_cells_df, by = "Sample") %>%
  mutate(Cell_Fraction = Count / total_cells)
# Assign Sample Groups (R vs. NR)
fraction_df <- fraction_df %>%
  mutate(Sample_Group = ifelse(str_detect(Sample, "_R$"), "R", "NR"))
# Generate plots for each cell type
plot_list <- lapply(unique(fraction_df$Cell_type), function(cell) {
  df_subset <- fraction_df %>% filter(Cell_type == cell)
  # Order samples within each group (R and NR) by fraction
  df_subset <- df_subset %>%
    arrange(Sample_Group, desc(Cell_Fraction))
  # Convert Sample to a factor **for this specific cell type**
  df_subset$Sample <- factor(df_subset$Sample, levels = df_subset$Sample)
  
  ggplot(df_subset, aes(x = Sample, y = Cell_Fraction, fill = Sample_Group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("R" = "#1f78b4", "NR" = "#e31a1c")) +
    labs(x = "Sample", y ="Fraction", fill = "Group",
         title = paste(cell, "in Total Non-Immune Cells")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

final_plot <- wrap_plots(plot_list, ncol = 3)  # Adjust columns for better layout
final_plot
ggsave("/home/u251079/scRNA/Lung_scRNA/Non_Immune_Cells_Fraction.pdf", final_plot, width = 15, height = 10)










#@@@@@@@@@@@@@@@@@@ Proportion test 
#install.packages("devtools") # https://github.com/rpolicastro/scProportionTest
#devtools::install_github("rpolicastro/scProportionTest")
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


pdf("/home/u251079/scRNA/Lung_scRNA/4_NI_scProp.pdf", width = 4.5, height =3.5)
permutation_plot(prop_test)
dev.off()



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

# Adjust p-values using BH method
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
View(Table_nonImm_samplewise)

write.csv(Table_nonImm_samplewise, file = "/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/CD45neg_FC_P_val.csv.csv")


names(Table_nonImm_samplewise)
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

pdf("/home/u251079/scRNA/Lung_scRNA/4_NonImm_Prop_boxplot.pdf", width = 7, height = 7)
ggplot(long_nonImm, aes(x = Group, y = Proportion, color = Group)) +
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
dev.off()





names(Table_nonImm_samplewise)
xx<-c("Cluster","FC_mean_prop", "p_value_adj")
vol_data<-Table_nonImm_samplewise[,xx]
vol_data$log_pval <- -log10(vol_data$p_value)
vol_data$Significance <- ifelse(vol_data$p_value < 0.05 & abs(log2(vol_data$FC_mean_prop)) > 1, "Significant", "Not Significant")
vol_data$log2_FC <- log2(vol_data$FC_mean_prop)

# Assign colors based on logFC values (red if logFC > 1, blue if logFC < -1)
vol_data$Color <- ifelse(vol_data$log2_FC > 1, "#458B74", ifelse(vol_data$log2_FC < -1, "#CD6600", "gray"))

# Create the volcano plot
p<-ggplot(vol_data, aes(x = log2_FC, y = log_pval, color = Color)) +
  geom_point(alpha = 0.8, size = 3) + 
  scale_color_identity() +  # Use the color column directly
  geom_text_repel(aes(label = Cluster), size = 4, max.overlaps = 10) +  # Add labels for clusters
  theme_classic() + 
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", title = "Volcano Plot: Non-Immune Cell Proportions") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +  # Significance line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Vertical line at log2FC = 1
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +  # Vertical line at log2FC = -1
  theme(
    text = element_text(size = 14),
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

pdf("/home/u251079/scRNA/Lung_scRNA/4_NonImm_Prop_volcano.pdf", width = 10, height = 7)
print(p)
dev.off()
#@@@@@@@@@@@@@@@@@@@@







#@@@@@@@@@@@@@@@@ AT2 cells
target_seu_subset <- subset(target_seu, idents = c(0, 1, 2))
p1<-VlnPlot(target_seu_subset, 
            features = c("SFTPC", "NKX2-1"), 
            #      group.by = "sample_type", 
            split.by = "Cell_type1") 

p1
p2<-VlnPlot(target_seu_subset, 
            features = c("SFTPC", "NKX2-1"), 
            group.by = "sample_type" )
#split.by = "Cell_type1") 

p2
p3<-p1|p2
p3
pdf("/home/u251079/scRNA/Lung_scRNA/AT2_Vln.pdf", width = 10, height =4)
print(p3)
dev.off()




#@@@@@@@@@@@@@@@@ Endothelial cells
library(Seurat)
library(ggpubr)
endo <- subset(target_seu, idents = c(4))
endo_filtered <- subset(endo, subset = PECAM1 > 0)
p4 <- VlnPlot(endo_filtered, 
              features = c("PECAM1"), 
              group.by = "sample_type") +
  labs(x = NULL) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)) +
  NoLegend() +
  stat_compare_means(method = "wilcox.test")

# Display plot
p4






pdf("/home/u251079/scRNA/Lung_scRNA/Endo_Vln.pdf", width = 4, height =3.5)
print(p4)
dev.off()
