This is a comprehensive bioinformatics analysis pipeline for lung cancer (specifically adenocarcinoma) data, focusing on recurrence vs non-recurrence comparisons across multiple data types and cell populations. 

## Overview
The code performs integrated analysis of:
- Bulk RNA-seq data (BCM dataset)
- Single-cell RNA-seq data (immune and non-immune cells)
- Survival analysis across multiple datasets
- Multi-omics integration with clinical outcomes

## Main Analysis Sections:

### 1. Bulk RNA-seq Analysis (Recurrence vs Non-Recurrence)
- Differential Expression: Uses limma/edgeR to identify DEGs between recurrent and non-recurrent tumors
- Pathway Analysis: GO and KEGG enrichment of DEGs
- GSEA: Hallmark pathway enrichment using ranked genes
- PCA: Principal component analysis to visualize sample separation
- Volcano Plots: Visualization of significant DEGs

### 2. Single-Cell RNA-seq Analysis
#### NK Cells Analysis
- Marker expression profiling (cytotoxic, inhibitory, activating markers)

#### T Cell Analysis 
- CD4/CD8/Gamma Delta T cell subpopulations
- Exhaustion, memory, and activation markers
- Heatmap and dot plot visualizations

#### Macrophage/Myeloid Analysis
- Multiple macrophage subsets (immunosuppressive, proliferating, etc.)
- Monocyte and neutrophil populations
- Comprehensive marker expression profiling

#### B Cell/Plasma Cell Analysis
- B cell maturation states (naive, memory, plasma cells)
- Immunoglobulin expression patterns
- Detailed characterization of plasma cell subsets

### 3. Non-Immune Cell Analysis
- Fibroblasts, endothelial cells, epithelial subsets
- Cancer cell markers and proliferation signatures
- Tumor microenvironment characterization

### 4. Survival Analysis
- Plasma Cell Signature: C12 plasma cells as prognostic marker
- Multiple Validation Cohorts: BCM and Okayama datasets
- Stratified Analysis: By stage and smoking status
- Kaplan-Meier plots with confidence intervals

## Key Biological Insights:
1. Plasma cells (C12) emerge as a strong prognostic factor for better recurrence-free survival
2. NK cell dysfunction in recurrent tumors
3. Macrophage polarization differences between recurrence states
4. T cell exhaustion patterns in recurrent disease
5. Multi-cellular immune interactions in tumor microenvironment

## Technical Approach:
- Integrated multi-scale data: Bulk + single-cell + clinical
- Comprehensive cell type characterization
- Pathway-centric interpretation
- Clinical validation across independent cohorts
- Rigorous statistical testing with multiple testing correction

This represents a sophisticated systems biology approach to understanding lung cancer recurrence mechanisms through integrated analysis of transcriptional programs across different cellular compartments and data modalities.

# Perform DEG analysis between Rec and Non-Rec samples of BCMLC samples.
# Then perform Pathway analysis.
rm(list = ls())
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
load("/Data/BCM_LC_RNAseq.rda")
output<-"/LUNG_scRNAseq/Images/New_Plots"
myinf1 = "/Data/BCM_prot_gene_TPM_Symbol.csv"
data=read.csv(myinf1, row.names = 1, header = T)
info<-ORI.info

data = log2(data+1)
dim(data)
se = grep("T", colnames(data))
data = data[,se]

se = grep("Adenocarcinoma|adenocarcinoma", info$Histology)
info = info[se,]

colnames(data)
colnames(data) <- gsub("^X", "LC", colnames(data))  
colnames(data) <- gsub("T$", "", colnames(data))    
head(data)
head(info)
colnames(data)
rownames(info)
common_samples <- intersect(colnames(data), rownames(info))
data <- data[, common_samples]
info <- info[common_samples, ]
info$Recurrence
info$Recurrence_Status <- ifelse(info$Recurrence == 1, "Rec", "NonRec")
info$Recurrence_Status <- factor(info$Recurrence_Status, levels = c("NonRec", "Rec"))
table(info$Recurrence_Status)
sample_group <- factor(info$Recurrence_Status, levels = c("NonRec", "Rec"))
design <- model.matrix(~0 + sample_group)
colnames(design) <- levels(sample_group)

contrast.matrix <- makeContrasts(Recurrence_vs_NonRecurrence = Rec - NonRec, 
                                 levels = design)
expr_matrix <- as.matrix(data)
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
deg_results <- topTable(fit3, n = Inf, coef = 1, adjust.method = "BH")
dim(deg_results)
head(deg_results)
deg_filtered <- deg_results[deg_results$adj.P.Val < 0.1, ]

library(ggplot2)
deg_results$Significant <- ifelse(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1, 
                                  "Significant", "Not Significant")

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Recurrence vs. Non-Recurrence",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value")

# --- PCA between Recurrence vs Non-Recurrence ---
top_var_genes <- order(apply(expr_matrix, 1, var), decreasing = TRUE)[1:100]
expr_pca <- t(expr_matrix[top_var_genes, ])
# Run PCA
pca_res <- prcomp(expr_pca, scale. = TRUE)
# Make data frame for plotting
pca_df <- data.frame(Sample = rownames(pca_res$x),
                     PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Recurrence = info$Recurrence_Status)

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Recurrence, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "PCA: Recurrence vs Non-Recurrence",
       x = paste0("PC1 (", round(100*summary(pca_res)$importance[2,1],1), "%)"),
       y = paste0("PC2 (", round(100*summary(pca_res)$importance[2,2],1), "%)"))

# Save PCA results
# --- Pathway Analysis (using clusterProfiler) ---
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
# Get significant DEGs
deg_genes <- rownames(deg_results)
# Convert to ENTREZ IDs
gene_map <- bitr(deg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_genes <- unique(gene_map$ENTREZID)

# Run GO enrichment
ego <- enrichGO(gene = entrez_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Run KEGG enrichment
ekegg <- enrichKEGG(gene = entrez_genes,
                    organism = 'hsa',
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)
# GO dotplot
dotplot(ego, showCategory=20, title="GO Biological Process")
# KEGG barplot
barplot(ekegg, showCategory=20, title="KEGG Pathways")


#@@@@@@@@@@@@@@@@@@@@@@@ NK cells SF3
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(scales)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(grid)
output<-"/LUNG_scRNAseq/Images/New_Plots"
cd45pos <- readRDS("/Data/cd45pos_nsclc_Immune_cells.rds")
table(cd45pos@meta.data$Cell_type)
table(cd45pos@meta.data$sample_type)
# file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/NK_NKT.rds")
# NK= readRDS(file)

NK_clusters <- c("NK cells")
NK <- subset(cd45pos, subset = Cell_type %in% NK_clusters)
table(NK@meta.data$sample)

NK$CellType_Sample <- paste("NK", NK$sample_type, sep = "_")
NK$CellType_Sample <- factor(NK$CellType_Sample, levels = c("NK_Non_Rec", "NK_Rec"))

# -------------------------------
# 2️⃣ Define NK markers by type
# -------------------------------
marker_list <- list(
  Cytotoxic = c("PRF1", "GZMB", "GZMH", "GZMA", "NKG7", "NCAM1", "KLRK1", "KLRD1"),
  Inhibitory = c("KLRC1", "KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "PDCD1", "TIGIT", "LAG3", "HAVCR2"),
  Activating = c("KLRK1", "CD226", "NCR1","NCR3", "FCGR3A"),
  Cytokine = c("IFNG", "TNF", "CSF2", "CCL5", "CCL22"),
  Trafficking = c("CCR7", "CXCR4", "CXCR6", "CCR5", "CX3CR1", "SELL"),
  Tissue_Resident = c("ITGA1", "CD69", "ITGAE")
)
NK_markers <- unique(unlist(marker_list))

# Create gene_map for MarkerType annotation
gene_map <- bind_rows(
  lapply(names(marker_list), function(type) {
    tibble(Gene = marker_list[[type]], MarkerType = type)
  })
) %>% distinct(Gene, .keep_all = TRUE)

# -------------------------------
# 3️⃣ DotPlot
# -------------------------------
dp <- DotPlot(NK,
              features = NK_markers,
              group.by = "CellType_Sample",
              assay = "RNA",
              dot.scale = 1,
              cluster.idents = FALSE)

dp_data <- dp$data %>%
  left_join(gene_map, by = c("features.plot" = "Gene"))

p1 <- ggplot(dp_data, aes(x = id, y = features.plot, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_size(range = c(0, 4)) +
  scale_size_area(max_size = 6) +
  scale_color_gradientn(colours = c("lightblue", "blue", "darkblue"),
                        limits = c(0, 1),
                        oob = squish,
                        name = "log2(count + 1)") +
  facet_wrap(~MarkerType, scales = "free_y", ncol = 1) +
  theme_cowplot() +labs(x=NULL, y=NULL)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = "bold"))
p1
pdf(file = file.path(output, "NK_dot_plot.pdf"),
    width = 5,height = 12)  
print(p1)
dev.off()


# -------------------------------
# 4️⃣ Heatmap (log2(count +1))
# -------------------------------
# -------------------------------
# Prepare numeric matrix from DotPlot summary
# -------------------------------
dp_summary <- dp_data %>%
  group_by(features.plot, id, MarkerType) %>%
  summarise(mean_scaled = mean(avg.exp.scaled), .groups = "drop") %>%
  tidyr::spread(key = id, value = mean_scaled)

mat <- dp_summary %>%
  select(-MarkerType) %>%
  column_to_rownames("features.plot") %>%
  as.matrix()

mat[is.na(mat)] <- 0
row_annotation <- dp_summary %>%
  select(features.plot, MarkerType) %>%
  column_to_rownames("features.plot")
marker_order <- unique(row_annotation$MarkerType)
row_annotation <- row_annotation %>%
  arrange(factor(MarkerType, levels = marker_order))
mat <- mat[rownames(row_annotation), ]

# -------------------------------
# Plot heatmap as grob
# -------------------------------
heatmap_grob <- pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = row_annotation,
  color = colorRampPalette(c("red", "white", "blue"))(50),
  silent = TRUE
)

heatmap_gg <- as.ggplot(heatmap_grob)
pdf(file = file.path(output, "NK_Heat.pdf"),
    width = 5,height = 10)  
print(heatmap_gg)
dev.off()
# -------------------------------
# 5️⃣ DEG analysis: NK_Rec vs NK_Non_Rec
# -------------------------------
# Set identity to CellType_Sample
Idents(NK) <- NK$CellType_Sample

# DEGs: Non-Recurrent vs Recurrent
DEG_nonrec_vs_rec <- FindMarkers(
  NK,
  ident.1 = "NK_Non_Rec",  # now Non-Recurrent is the "up" group
  ident.2 = "NK_Rec",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

DEG_nonrec_vs_rec <- DEG_nonrec_vs_rec %>%
  rownames_to_column("Gene")
# Up-regulated in Non-Recurrent (higher in Non-Rec)
up_genes_nonrec <- DEG_nonrec_vs_rec %>%
  filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
  pull(Gene)

# Down-regulated in Non-Recurrent (higher in Rec)
down_genes_nonrec <- DEG_nonrec_vs_rec %>%
  filter(avg_log2FC < 0 & p_val_adj < 0.05) %>%
  pull(Gene)

library(msigdbr)
library(clusterProfiler)

# Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# Enrichment for up-regulated genes in Non-Recurrent
ego_up_nonrec <- enricher(
  gene = up_genes_nonrec,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Enrichment for down-regulated genes in Non-Recurrent
ego_down_nonrec <- enricher(
  gene = down_genes_nonrec,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

p1<-dotplot(ego_up_nonrec, showCategory = 15) + ggtitle("Up-regulated pathways in NK_Non_Rec")
p2<-dotplot(ego_down_nonrec, showCategory = 15) + ggtitle("Down-regulated pathways in NK_Non_Rec")

pdf(file = file.path(output, "NK_Hallmark.pdf"),
    width = 10,height = 8)  
print(p1/p2)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@. T cells SF5
output<-"/LUNG_scRNAseq/Images/New_Plots/"
cd45pos <- readRDS("/Data/cd45pos_nsclc_Immune_cells.rds")
table(cd45pos@meta.data$Cell_type)
table(cd45pos@meta.data$sample_type)

T_cells_clusters <- c( "CD4", "CD8", "GaDelT")
T_cells <- subset(cd45pos, subset = Cell_type %in% T_cells_clusters)
# file=c("/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/T_cells.rds")
# T_cells<-readRDS(file)
table(T_cells@meta.data$Cell_type1)
"CD4_C10     CD4_C6     CD4_C7     CD8_C0     CD8_C1    CD8_C15 GaDelT_C14 
943       1909       1694       3965       3597        293        412"
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "CD4_C6", "C6_CD4_Tcm", NA)
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "CD4_C7", "C7_CD4_Trag", T_cells@meta.data$Subcells)
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "CD4_C10", "C10_CD4_Ext", T_cells@meta.data$Subcells)
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "CD8_C0", "C0_CD8_Tem", T_cells@meta.data$Subcells)
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "CD8_C1", "C1_CD8_Ext", T_cells@meta.data$Subcells)
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "CD8_C15", "C15_CD8_actTRM", T_cells@meta.data$Subcells)
T_cells@meta.data$Subcells<-ifelse(T_cells@meta.data$Cell_type1 %in% "GaDelT_C14", "C14_GaDelT", T_cells@meta.data$Subcells)
table(T_cells@meta.data$Subcells)

T_cells@meta.data$Subcells<-paste0(T_cells@meta.data$Subcells, "_",T_cells@meta.data$sample_type)


View(T_cells@meta.data)
DimPlot(T_cells)
T_cells

DefaultAssay(T_cells) <- "RNA" 
T_marker<-c("CD3D", "TRAC", "CD4", "CD8A", "CD8B","GZMB", "GZMH", "GZMK", "IFNG",
            "IL7R", "CD40LG","TCF7","SELL", "CCR7", "LEF1", "MAL", "CD28", "CD27",
            "ITGA1", "ITGAE", "PDCD1", "TOX", "HAVCR2", "LAG3", "CXCL13", "ZBED2",
            "ETV1", "LAYN", "ENTPD1", "TIGIT", "TNFRSF9", "TNFRSF18", "CTLA4", "ICOS",
            "MAGEH1", "IL2RA", "FOXP3", "TNFRSF4", "HLA-DRA", "HLA-DRB1", "MKI67", "STMN1")
T_marker1<-DotPlot(T_cells, features = T_marker, group.by = "Subcells")+coord_flip()+RotatedAxis()+labs(x=NULL, y=NULL)+scale_size(range = c(1,7))
print(T_marker1)
tmp_averages <- AverageExpression(T_cells,
                                  assays="RNA",
                                  features=T_marker, group.by = "Subcells")$RNA
# Assuming tmp_averages is already calculated
# Extract the average expression matrix
average_expression_matrix <- as.matrix(tmp_averages)
# Create a heatmap using pheatmap
T_marker2<-pheatmap::pheatmap(average_expression_matrix,
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                              angle_col = 45)
library(cowplot)
pdf(paste0(output,"T_markers.pdf"), width = 15, height =10)
plot_grid(T_marker1, as.ggplot(T_marker2), ncol = 2, labels = c("A", "B"))
dev.off()


library(ggplot2)
library(cowplot)
library(ggplotify)
library(pheatmap)

pdf(paste0(output, "T_markers.pdf"), width = 15, height = 10)

plot_grid(
  T_marker1,
  as.ggplot(T_marker2),  # convert the pheatmap to ggplot
  ncol = 2,
  labels = c("A", "B")
)

dev.off()






#@@@@@@@@@@@@@@@@@@@@@@@@ Macrophages SF6 and Plasma cells 
# https://www.nature.com/articles/s41597-024-03885-x
rm(list = ls())
library(Seurat)
library(dplyr)
library(cowplot)
library(SeuratObject)
library(ggplot2)
library(scales)
library(ggrepel)

output<-"/LUNG_scRNAseq/Images/New_Plots"
cd45pos <- readRDS("/Data/cd45pos_nsclc_Immune_cells.rds")
table(cd45pos@meta.data$Cell_type1,cd45pos@meta.data$sample_type)

macrophage_clusters <- c("Mac_C8", "Mac|DC_C11", "Prolifirating Mac_C16", "Prolif Mac|DC_C20", "Mono_C3", "Neutrophil_C5")
Macrophages <- subset(cd45pos, subset = Cell_type1 %in% macrophage_clusters)
DefaultAssay(Macrophages) <- "RNA"
table(Macrophages@meta.data$sample_type, Macrophages@meta.data$Cell_type1)

M_G <- list("Monocytes" = c("CD14", "CSF3R"),
            "Neutrophils" = c("FCGR3B", "CSF3R", "S100A12", "S100A8"),
  "Low_quality_Macro" = c("LYZ", "FTL"),  
  "Lipid_associated_Macro" = c("MS4A7", "IL1B", "IL4I1", "FOLR2", "APOE", "C1QA", "C1QB", "C1QC", "CTSB", "CTSD"),
  "Alveolar_Macro" = c("MCEMP1", "PPARG", "MRC1"),
  "Proliferating_Macro" = c("CDCA8", "MKI67", "CENPF", "CD14", "TOP2A"),
#)
#M_G <- list(
  "Immunosuppressive_Macrophages" = c("APOE", "FOLR2", "C1QA", "C1QB", "C1QC", "MRC1", "CD163"),
  "Lipid_Associated_Macrophages" = c("FABP4", "APOE", "LPL", "TREM2"),
  "Alveolar_Macrophages" = c("PPARG", "MRC1", "MCEMP1", "MARCO"),
  "Proinflammatory_M1_Macrophages" = c("CD80", "CD86", "IL12B", "NOS2", "TNF", "IL1B"),
  "AntiInflammatory_M2_Macrophages" = c("CD163", "MRC1", "ARG1", "TGFBI", "IL10"),
  "Interstitial_Macrophages" = c("CD11b", "ITGAM", "HLA-DRA", "F13A1"),
  "Proliferating_Macrophages" = c("MKI67", "TOP2A", "CDCA8", "CENPF"),
  "Tumor_Associated_Macrophages" = c("VEGFA", "MMP9", "PDL1 (CD274)", "TGFB1")
)

markers_df <- stack(M_G)
colnames(markers_df) <- c("Symbol", "CellType")

valid_markers_df <- markers_df[markers_df$Symbol %in% rownames(Macrophages[["RNA"]]), ]
valid_markers_df <- valid_markers_df[order(valid_markers_df$CellType), ]

celltype_labels <- unique(valid_markers_df$CellType)
gene_list <- split(valid_markers_df$Symbol, valid_markers_df$CellType)
valid_markers <- unlist(gene_list)
names(valid_markers) <- rep(celltype_labels, times = lengths(gene_list))
valid_markers <- valid_markers[!duplicated(valid_markers)]

# Keep group mapping as a dataframe
gene_map <- valid_markers_df[, c("Symbol", "CellType")]

# Create a combined column
Macrophages$CellType_Sample <- paste(Macrophages$Cell_type1, Macrophages$sample_type, sep = "_")
unique(Macrophages$CellType_Sample)
Macrophages$CellType_Sample <- factor(
  Macrophages$CellType_Sample,
  levels = c(
    "Mac_C8_Non_Rec", "Mac_C8_Rec",
    "Mac|DC_C11_Non_Rec", "Mac|DC_C11_Rec",
    "Prolifirating Mac_C16_Non_Rec", "Prolifirating Mac_C16_Rec",
    "Prolif Mac|DC_C20_Rec",
    "Mono_C3_Non_Rec", "Mono_C3_Rec",
    "Neutrophil_C5_Non_Rec", "Neutrophil_C5_Rec"
  )
)

dp <- DotPlot(Macrophages, 
              features = unique(gene_map$Symbol), 
              group.by = "CellType_Sample", 
              assay = "RNA", 
              dot.scale = 1, 
              cluster.idents = FALSE)

# Add marker group info
dp_data <- dp$data %>%
  left_join(gene_map, by = c("features.plot" = "Symbol"))

# Plot with ggplot2 and facet_wrap for marker groups
p<-ggplot(dp_data, aes(x = id, y = features.plot, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_size(range = c(0, 4)) +
  scale_size_area(max_size = 6) +
  scale_color_gradientn(colours = c("lightblue", "blue", "darkblue"),
                        limits = c(0, 1),
                        oob = scales::squish,
                        name = "log2 (count + 1)") +
  facet_wrap(~CellType, scales = "free_y", ncol = 1) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(face = "bold"))
p

pdf(file = file.path(output, "Macrophages_DotPlot_CellType.pdf"),
    width = 6,height = 20)  
print(p)
dev.off()


# Compute mean scaled expression
dp_summary <- dp_data %>%
  group_by(id, CellType) %>%
  summarise(mean_scaled = mean(avg.exp.scaled), .groups = "drop") %>%
  tidyr::spread(CellType, mean_scaled)

# Convert to matrix for pheatmap
mat <- as.matrix(dp_summary[,-1])
rownames(mat) <- dp_summary$id

# Heatmap
pdf(file = file.path(output, "Macrophages_Heatmap.pdf"),
    width = 7,height = 5)  
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("red","white", "blue"))(50),
         main = "Macrophages marker expression and sample type")
dev.off()



#@@@@@@@@@@@@@@@@@@@. PLASMA B cells 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Subset plasma cells C12 
# plasma cells https://pmc.ncbi.nlm.nih.gov/articles/PMC9633381/#S2
# B cells https://www.cell.com/cell/fulltext/S0092-8674(24)00712-8#sec-2
# https://ashpublications.org/blood/article/144/5/496/515800/Integrative-single-cell-chromatin-and 
# https://ashpublications.org/blood/article/114/25/5173/26511/An-in-vitro-model-of-differentiation-of-memory-B 
rm(list = ls())
output<-"/LUNG_scRNAseq/Images/New_Plots/"
cd45pos <- readRDS("/Data/cd45pos_nsclc_Immune_cells.rds")

unique(cd45pos@meta.data$Cell_type1)
plasma_cells <- subset(cd45pos, Cell_type1 %in% c("B cells_C2","B cells_C13", "Plasma cells_C12", "Plasma cells_C9")) #Plasma cells_C9 or Plasma cells_C12
plasma_cells <- subset(cd45pos, Cell_type1 %in% c("B cells_C2","B cells_C13", "Plasma cells_C12", "Plasma cells_C9")) #Plasma cells_C9 or Plasma cells_C12

#@@@@@@@@@@@@@ # STUDY-1 for mg_1
#@@@@@@@@@@@@@@@@@@@@ STUDY-2  mg_2 for supplement discussion figure 
mg_1 <- c(
  # Naïve B cells
  "MS4A1", "CD19", "BANK1", "SELL", "IGHD", "IL4R", "FCER2", "TCL1A", "BACH2",
  # Memory B cells
  "CD27", "CD24",
  # Plasmablasts and Plasma Cells
  "CD38", "SDC1", "TNFRSF17", "MZB1", "JCHAIN",
  # Immunoglobulin genes (commonly expressed in various B cells, including plasma cells)
  "IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",
  # Proliferation markers (can be associated with plasmablasts, plasma cells, or proliferating B cells)
  "MKI67", "TUBB", "STMN1", "TYMS", "BAG3", "HSPA6", "HSPB1",
  # MHC class II genes
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2",
  "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1",
  "HLA-DRB5"
)
mg_2 <- c(
  "FCRL5", "TCL1A","IL4R", "CD72","BACH2","IGHD", "IGHM", "NR4A1","NR4A2","CREM", "CD83",
  "ISG15", "IF44L","IFIT3","CD27", "TNFRSF13B","TXNIP","GPR183","HSPA1A", "HSPA1B","DNAJB1", "EGR1",
  "FCRL5","CCR1", "CXCR3","PDCD1","HCK", "FCRL3","FGR", "NME1","APEX1","POLD2", "POLE3", "MYC", "BCL6",
  "RGS13", "AICDA", "IL21R", "HMGB2", "CD38", "MZB1", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DQA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4",
  "IGHA1", "IGHG2", "IL10", "IL12A", "EBI3", "TGFB1", "IL35",
  "ITGAX", "TBX21", "CR2", "FCRL4", "PDCD1", 
  "PRDM1", "IRF4", "XBP1", "MS4A1", "LTB", "RGS13", 
  "MKI67", "STMN1", "TOP2A", "IGHG", "IGHA"
)

unique(mg_1, mg_2)
mg_combined <- unique(c(mg_1, mg_2))


plot1 <- DotPlot(plasma_cells, features =mg_combined, group.by = "Cell_type1", split.by = "sample_type") + 
  scale_size(range = c(1, 6)) +  
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_flip()

print(plot1)

pdf(paste0(output,"B_Plasma_Dotplot.pdf"), width = 7, height =17)
print(plot1)
dev.off()

#@@@@@@@@@@@@@@@@@@@ Study-3 

library(Seurat)
library(pheatmap)
library(dplyr)

plasma_cells@meta.data$Cell_type1 <- recode(plasma_cells@meta.data$Cell_type1, 
                                            "Plasma cells_C9"  = "stressed_PC_C9", 
                                            "Plasma cells_C12" = "prePB_C12",
                                            "B cells_C13"      = "MBC_C13",
                                            "B cells_C2"       = "actBC_C2")
mg_3 <- c(
  "TSC22D3", "KLF6", "EZR", "MS4A1", "CD37","JUND", "HLA-DPB1", "CD69","KLF2", "FOXP1",
  "ACTG1", "ACTB", "TUBA1B", "TUBB", "ENO1", "LDHA", "HSPD1", "RAN", "PFN1", "IL2RA",
  "MYDGF", "TP53INP1", "PPIB", "IGHG1", "PDIA4", "ITM2A", "CD44", "UBE2J1", "SEC11C", "RPN2",
  "IFI6", "ISG15", "XAF1", "MX1","IFI44L", "SSR4", "IRF7", "XBP1", "CD79A", "SAMD9L"
)

MBC  <- c("TSC22D3", "KLF6", "EZR", "MS4A1", "CD37","JUND", "HLA-DPB1", "CD69","KLF2", "FOXP1")
prePB <- c("ACTG1","ACTB", "TUBA1B", "TUBB", "ENO1", "LDHA", "HSPD1", "RAN", "PFN1", "IL2RA")
PB   <- c("MYDGF","TP53INP1", "PPIB", "IGHG1", "PDIA4", "ITM2A", "CD44", "UBE2J1", "SEC11C", "RPN2") 
PC   <- c("IFI6", "ISG15", "XAF1", "MX1","IFI44L", "SSR4", "IRF7", "XBP1", "CD79A", "SAMD9L")

plasma_cells$combined_group <- paste(plasma_cells$Cell_type1, plasma_cells$sample_type, sep = "_")

tmp_averages_split <- AverageExpression(
  plasma_cells,
  assays = "RNA",
  features = mg_3,
  group.by = "combined_group"
)$RNA
average_expression_matrix_split <- as.matrix(tmp_averages_split)

marker_annotation <- data.frame(Cell_Type = rep(NA, length(mg_3)), row.names = mg_3)
marker_annotation$Cell_Type[mg_3 %in% MBC]   <- "MBC"
marker_annotation$Cell_Type[mg_3 %in% prePB] <- "prePB"
marker_annotation$Cell_Type[mg_3 %in% PB]    <- "PB"
marker_annotation$Cell_Type[mg_3 %in% PC]    <- "PC"
marker_annotation$Cell_Type <- factor(marker_annotation$Cell_Type, levels = c("MBC", "prePB", "PB", "PC"))

ann_colors <- list(
  Cell_Type = c(
    "MBC" = "#548B54",
    "prePB" = "#CD6839",
    "PB" = "#27408B",
    "PC" = "red"
  )
)

ordered_genes <- rownames(marker_annotation)[order(marker_annotation$Cell_Type)]
average_expression_matrix_split <- average_expression_matrix_split[ordered_genes, ]

BC_PC_marker_split <- pheatmap::pheatmap(
  average_expression_matrix_split,
  cluster_rows = T,
  cluster_cols = TRUE,
  scale = "row",
  color = colorRampPalette(c("#104E8B", "white", "#B22222"))(50),
  angle_col = 45,
  border_color = "#363636",
  annotation_row = marker_annotation,
  fontsize = 9,
  treeheight_row = 10,
  treeheight_col = 10
)

pdf(paste0(output, "B_Plasma_Heatmap_split.pdf"), width = 7, height = 8)
print(BC_PC_marker_split)
dev.off()



#@@@@@@@@@@@@@@@@@@@ BCMLC GSEA
#============================================================
# 0️⃣ Setup
#============================================================
rm(list = ls())
output <- "/LUNG_scRNAseq/Images/New_Plots/"
load("/Data/BCM_LC_RNAseq.rda")
BCMLC_data <- read.csv("/Data/BCM_prot_gene_TPM_Symbol.csv")
info <- ORI.info

#------------------------------------------------------------
# Load packages
#------------------------------------------------------------
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
  library(ggrepel)
})

expr <- as.matrix(BCMLC_data)
colnames(expr) <- gsub("^X", "", colnames(expr))
expr <- expr[, grepl("T$", colnames(expr))]  # only tumor samples

# Metadata
info$Recurrence_Status <- ifelse(info$Recurrence == 1, "Recurrence", "Non-Recurrence")
info$Recurrence_Status <- factor(info$Recurrence_Status, levels = c("Non-Recurrence", "Recurrence"))

# Match sample order
common_samples <- intersect(colnames(expr), info$Sample.id)
expr <- expr[, common_samples]
info <- info[match(common_samples, info$Sample.id), ]

#------------------------------------------------------------
# 1️⃣ Log2-transform TPM for limma
#------------------------------------------------------------
expr_log <- log2(expr + 1)

design <- model.matrix(~ Recurrence_Status, data = info)
colnames(design) <- c("Intercept", "Recurrence")

fit <- lmFit(expr_log, design)
fit <- eBayes(fit)
deg_limma <- topTable(fit, coef = "Recurrence", number = Inf, sort.by = "P")

deg_limma <- deg_limma %>% 
  mutate(gene = rownames(.),
         status = ifelse(P.Value < 0.05 & logFC > 0, "Up",
                         ifelse(P.Value < 0.05 & logFC < 0, "Down", "NS")))

write.table(deg_limma, paste0(output, "DEG_limma_TPM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------------------------
# 2️⃣ Volcano plot
#------------------------------------------------------------
deg_volcano <- deg_limma %>%
  mutate(sig = ifelse(P.Value < 0.01, "Significant", "NS"))

ggplot(deg_volcano, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "#D55E00", "NS" = "grey")) +
  theme_classic(base_size = 14) +
  labs(title = "Volcano Plot of DEGs (TPM)",
       x = "log2 Fold Change (Recurrence vs Non-Recurrence)",
       y = "-log10(P-value)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed")

#------------------------------------------------------------
# 3️⃣ GSEA using ranked genes
#------------------------------------------------------------
deg_limma$ranking_metric <- deg_limma$logFC * -log10(deg_limma$P.Value)
ranks <- deg_limma$ranking_metric
names(ranks) <- deg_limma$gene
ranks <- sort(ranks, decreasing = TRUE)

msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_res <- GSEA(ranks,
                 TERM2GENE = msig_h,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

gsea_df <- as.data.frame(gsea_res@result)
write.table(gsea_df, paste0(output, "GSEA_results_TPM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Top 10 pos/neg pathways
pos_pathways <- gsea_df %>% filter(NES > 0) %>% top_n(10, NES)
neg_pathways <- gsea_df %>% filter(NES < 0) %>% top_n(-10, NES)
top_pathways <- rbind(pos_pathways, neg_pathways)
top_pathways$Description
p <- ggplot(
  top_pathways %>% 
    mutate(Group = ifelse(NES > 0, "Recurrence", "Non-Recurrence")),
  aes(x = reorder(Description, NES), y = NES, fill = Group)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(
    values = c("Non-Recurrence" = "#0072B2",   # blue
               "Recurrence" = "#D55E00")       # orange
  ) +
  labs(
    title = "Top Enriched Pathways \n(Recurrence vs Non-Recurrence)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "Upregulated"
  ) +
  theme_classic(base_size = 12)

p
ggsave(paste0(output, "GSEA_TopPathways_TPM.pdf"), height = 4.5, width = 8)

#------------------------------------------------------------
# 4️⃣ PCA of TPM
#------------------------------------------------------------
expr_log <- log2(expr + 1)
pca_res <- prcomp(t(expr_log), scale. = TRUE)
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Recurrence_Status = info$Recurrence_Status,
  Sample = info$Sample.id
)
library(ggrepel)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Recurrence_Status, label = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3) +
  scale_color_manual(values = c("Non-Recurrence" = "#0072B2", 
                                "Recurrence" = "#D55E00")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +
  ggtitle("Unbiased PCA of Bulk RNA-seq Samples (TPM)")

p1<-ggplot(pca_df, aes(x = PC1, y = PC2, color = Recurrence_Status, label = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = Inf) +  # show all labels
  scale_color_manual(values = c("Non-Recurrence" = "#0072B2", 
                                "Recurrence" = "#D55E00")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +
  ggtitle("Unbiased PCA of Bulk RNA-seq Samples (TPM)")

p1
library(ggplot2)
library(ggrepel)
library(dplyr)
library(matrixStats)

# Select top 50 DEGs based on adjusted p-value and absolute logFC
top_DEGs <- deg_limma %>%
  filter(status %in% c("Up", "Down")) %>%
  arrange(adj.P.Val) %>%
  slice(1:50) %>%
  pull(gene)

# Subset expression matrix for top 50 DEGs
expr_top50 <- expr_log_filtered[top_DEGs, ]

# PCA
pca_res_top50 <- prcomp(t(expr_top50), scale. = TRUE)
percentVar_top50 <- round(100 * (pca_res_top50$sdev^2 / sum(pca_res_top50$sdev^2)), 1)

# PCA dataframe
pca_df_top50 <- data.frame(
  PC1 = pca_res_top50$x[,1],
  PC2 = pca_res_top50$x[,2],
  Recurrence_Status = info$Recurrence_Status,
  Sample = info$Sample.id
)

# Convex hull for each group
hull_data <- pca_df_top50 %>%
  group_by(Recurrence_Status) %>%
  slice(chull(PC1, PC2))

# PCA plot
p2<-ggplot(pca_df_top50, aes(x = PC1, y = PC2, color = Recurrence_Status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 20) +
  geom_polygon(data = hull_data, aes(fill = Recurrence_Status), alpha = 0.15, color = NA) +
  scale_color_manual(values = c("Non-Recurrence" = "#0072B2", 
                                "Recurrence" = "#D55E00")) +
  scale_fill_manual(values = c("Non-Recurrence" = "#0072B2", 
                               "Recurrence" = "#D55E00")) +
  xlab(paste0("PC1: ", percentVar_top50[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_top50[2], "% variance")) +
  theme_classic(base_size = 14) +
  ggtitle("PCA of Bulk RNA-seq Samples (Top 50 DEGs, TPM)") +
  theme(legend.position = "right")
p2

library(patchwork)

# Add labels
p1_labeled <- p1 + ggtitle("A. Unbiased PCA of Bulk RNA-seq Samples (TPM)")
p2_labeled <- p2 + ggtitle("B. PCA of Top 50 DEGs (TPM)")

# Combine
combined_plot <- p1_labeled / p2_labeled  # stacked vertically

# Save
ggsave(filename = paste0(output, "PCA_combined_top50DEGs.pdf"),
       plot = combined_plot, width = 10, height = 12)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Non_Immune cells
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
library(Seurat)
library(pheatmap)
output <- "/LUNG_scRNAseq/Images/New_Plots/"
file=c("/Data/cd45negc_celltype_data.rds")
cd45neg<-readRDS(file) 
cd45neg
unique(cd45neg$Cell_type1)
unique(cd45neg$Cell_type2)
unique(cd45neg$Cell_type)
unique(cd45neg$sample_type)
table(cd45neg$Cell_type1, cd45neg$sample_type)
cd45neg$Cell_type1_sample <- paste0(cd45neg$Cell_type1, "_", cd45neg$sample_type)
head(cd45neg$Cell_type1_sample)
unique(cd45neg$Cell_type1_sample)

non_immune_markers <- list(
  Fibro = c("COL1A1", "COL1A2"),
  Endo = c("PECAM1", "ECSCR"),
  AT2_Clara = c("EPCAM", "CDH1", "NKX2-1", "SFTPC", "SFTPA1"),
  Clara_BASC = c("SCGB1A1", "SCGB3A1"),
  AT1 = c("PDPN"),
  Ciliated = c("FOXJ1")
)

# https://static-content.springer.com/esm/art%3A10.1038%2Fs41597-023-02074-6/MediaObjects/41597_2023_2074_MOESM1_ESM.pdf
non_immune_markers <- list(
  Fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "DCN", "FBLN1"),
  Endothelial = c("PECAM1", "ECSCR", "VWF", "KDR"),
  AT2_Clara = c("EPCAM", "CDH1", "NKX2-1", "SFTPC", "SFTPA1"),
  Clara_BASC_Cili = c("SCGB1A1", "SCGB3A1", "FOXJ1"),
  Epithelial_Cancer = c("SFTPC", "SFTPA1", "SFTPB", "MUC1", "KRT8","CDKN2A", "SOX2", "CXCL1", "LAMC2"),
  Prolif = c("MKI67", "TOP2A")
)

markers_df <- stack(non_immune_markers)
colnames(markers_df) <- c("Gene", "CellType")

valid_markers_df <- markers_df %>% 
  filter(Gene %in% rownames(cd45neg@assays$RNA$counts)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  arrange(CellType)

valid_markers <- valid_markers_df$Gene

avg_exp <- AverageExpression(
  cd45neg,
  features = valid_markers,
  group.by = "Cell_type1_sample",
  assays = "RNA"
)$RNA

avg_exp_log <- log1p(avg_exp)

p <- pheatmap(
  avg_exp_log[valid_markers, ],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  color = colorRampPalette(c("darkblue", "white", "darkred"))(50),
  main = "Average Expression Heatmap of Non-immune & Cancer Markers",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA
)

print(p)

pdf(paste0(output, "Non_immune_marker_expr.pdf"), height = 5, width = 8)
print(p)
dev.off()



######@@@@@@@@@@@@ Dot plot
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
                 group.by = 'Cell_type1_sample', 
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

pdf(paste0(output, "Non_immune_marker_dot_plot.pdf"), height = 8, width = 15)
print(plot2)
dev.off()


############# SUMMARY of Marker genes
DefaultAssay(cd45neg) <- "RNA"

# Marker list
non_immune_markers <- list(
  Fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "DCN", "FBLN1"),
  Endothelial = c("PECAM1", "ECSCR", "VWF", "KDR"),
  AT2_Clara = c("EPCAM", "CDH1", "NKX2-1", "SFTPC", "SFTPA1"),
  Clara_BASC_Cili = c("SCGB1A1", "SCGB3A1", "FOXJ1"),
  Epithelial_Cancer = c("SFTPC", "SFTPA1", "SFTPB", "MUC1", "KRT8","CDKN2A", "SOX2", "CXCL1", "LAMC2"),
  Prolif = c("MKI67", "TOP2A")
)

# Flatten marker list
markers_df <- stack(non_immune_markers)
colnames(markers_df) <- c("Gene", "CellType")

# Filter valid markers present in RNA assay
valid_markers <- markers_df$Gene[markers_df$Gene %in% rownames(cd45neg@assays$RNA$counts)]

# Average expression by sample type (NR vs R)
avg_exp <- AverageExpression(cd45neg, features = valid_markers, group.by = "sample_type", assays = "RNA")$RNA
avg_exp_log <- log1p(avg_exp)

# Fetch expression matrix for valid markers
expr_matrix <- FetchData(cd45neg, vars = valid_markers)
expr_matrix$sample_type <- cd45neg$sample_type

# Calculate proportion of cells expressing each gene
prop_expr <- expr_matrix %>%
  pivot_longer(cols = all_of(valid_markers), names_to = "Gene", values_to = "Expression") %>%
  mutate(Detected = Expression > 0) %>%
  group_by(sample_type, Gene) %>%
  summarise(Proportion = mean(Detected), .groups = "drop") %>%
  pivot_wider(names_from = sample_type, values_from = Proportion, names_prefix = "Prop_")

colnames(avg_exp_log)
summary_df <- data.frame(
  Gene = rownames(avg_exp_log),
  AvgExp_NR = avg_exp_log[, "Non_Rec"],
  AvgExp_R = avg_exp_log[, "Rec"]
) %>%
  left_join(prop_expr, by = "Gene")

summary_df

###########################################
# Plasma Cells (C12) DEG & Pathway Analysis
# Goal: Identify anti-tumor mechanisms associated with plasma cells (NR vs R)
###########################################
rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(msigdbr)
library(cowplot)
library(ggplotify)
library(ggpubr)
library(viridis)

# ----------------------
# Load Seurat object
# ----------------------
file <- "/Data/cd45pos_nsclc_Immune_cells.rds"
cd45pos <- readRDS(file)

# ----------------------
# 1. Subset Plasma Cells (C12)
# ----------------------
plasma_cells <- subset(cd45pos, Cell_type1 == "Plasma cells_C12")
plasma_cells@meta.data <- plasma_cells@meta.data %>%
  mutate(Sample_group = ifelse(grepl("NR", sample), "NR", "R"))
Idents(plasma_cells) <- plasma_cells@meta.data$Sample_group

# ----------------------
# 2. Differential Expression
# ----------------------
deg_wilcox <- FindMarkers(
  object = plasma_cells,
  ident.1 = "NR",
  ident.2 = "R",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)
deg_wilcox$gene <- rownames(deg_wilcox)

# Classify DEGs
logFC_threshold <- 0.25
pval_threshold <- 0.05
deg_wilcox$color_group <- "NS"
deg_wilcox$color_group[deg_wilcox$p_val_adj < pval_threshold & deg_wilcox$avg_log2FC > logFC_threshold] <- "Up"
deg_wilcox$color_group[deg_wilcox$p_val_adj < pval_threshold & deg_wilcox$avg_log2FC < -logFC_threshold] <- "Down"
deg_wilcox$color_group <- factor(deg_wilcox$color_group, levels = c("NS","Down","Up"))

# ----------------------
# 3. Volcano plot
# ----------------------
top_up <- deg_wilcox %>% filter(color_group=="Up") %>% slice_max(avg_log2FC, n=10)
top_down <- deg_wilcox %>% filter(color_group=="Down") %>% slice_min(avg_log2FC, n=10)
top_genes <- rbind(top_up, top_down)

volcano_plot <- ggplot(deg_wilcox, aes(x=avg_log2FC, y=-log10(p_val_adj), color=color_group)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("grey","#00688B","#EE2C2C")) +
  geom_hline(yintercept=-log10(pval_threshold), linetype="dashed", color="blue") +
  geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), linetype="dashed", color="blue") +
  geom_text_repel(data=top_genes, aes(label=gene), size=3, max.overlaps=10) +
  theme_classic() +
  labs(title="DEGs: NR vs R - Plasma Cells",
       x="log2 Fold Change", y="-log10(adj. P-value)", color="Significance")

# ----------------------
# 4. Ranked gene list for GSEA
# ----------------------
ranked_genes <- deg_wilcox %>%
  dplyr::select(gene, avg_log2FC) %>%
  filter(!duplicated(gene)) %>%
  arrange(desc(avg_log2FC))
ranked_vec <- setNames(ranked_genes$avg_log2FC, ranked_genes$gene)

# Map gene symbols to Entrez IDs
ranked_genes$EntrezID <- mapIds(org.Hs.eg.db, keys=ranked_genes$gene, keytype="SYMBOL", column="ENTREZID")
entrez_ids <- na.omit(ranked_genes$EntrezID)

# ----------------------
# 5. GO ORA
# ----------------------
go_res_all <- enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="ALL", pvalueCutoff=0.05)
go_df <- as.data.frame(go_res_all)

# Add Regulation and triangle shape
go_df$Regulation <- sapply(go_df$geneID, function(ids){
  gene_list <- unlist(strsplit(ids,"/"))
  gene_fc <- ranked_genes$avg_log2FC[ranked_genes$EntrezID %in% gene_list]
  if(length(gene_fc)==0) return(NA)
  ifelse(mean(gene_fc) > 0, "Up in NR", "Down in NR")
})
go_df$Shape <- ifelse(go_df$Regulation=="Up in NR", 24, ifelse(go_df$Regulation=="Down in NR", 25, 21))
go_sig <- go_df %>% filter(p.adjust <= 0.01) %>% arrange(p.adjust)
go_top <- go_sig %>% slice_head(n=20)
go_top$fill_padj <- -log10(go_top$p.adjust)

go_facet_plot <- ggplot(go_top, aes(x=reorder(Description, Count),
                                    y=Count,
                                    color=Regulation,
                                    shape=Shape,
                                    size=Count,
                                    fill=fill_padj)) +
  geom_point(stroke=1.2) +
  coord_flip() +
  facet_wrap(~ONTOLOGY, scales="free_y", ncol=1) +
  scale_color_manual(values=c("Up in NR"="#EE2C2C", "Down in NR"="#00688B", "NA"="grey")) +
  scale_shape_identity() +
  scale_fill_viridis_c(option="D", name="-log10(adj P)") +
  scale_size(range=c(3,8)) +
  guides(size = guide_legend(override.aes = list(shape=24, fill="grey"))) +
  theme_bw(base_size=12) +
  labs(title="GO Enrichment (All Ontologies) with NR Regulation",
       x="GO Term", y="Gene Count", color="Regulation in NR", size="Gene Count")
go_facet_plot 
# ----------------------
# 6. KEGG ORA
# ----------------------
kegg_res <- enrichKEGG(gene=entrez_ids, organism="hsa", pvalueCutoff=0.05)
kegg_df <- as.data.frame(kegg_res)

kegg_df$Regulation <- sapply(kegg_df$geneID, function(ids){
  gene_list <- unlist(strsplit(ids,"/"))
  gene_fc <- ranked_genes$avg_log2FC[ranked_genes$EntrezID %in% gene_list]
  if(length(gene_fc)==0) return(NA)
  ifelse(mean(gene_fc) > 0, "Up in NR", "Down in NR")
})
kegg_df$Shape <- ifelse(kegg_df$Regulation=="Up in NR", 24, ifelse(kegg_df$Regulation=="Down in NR", 25, 21))
kegg_sig <- kegg_df %>% filter(p.adjust <= 0.01) %>% arrange(p.adjust)
kegg_top <- kegg_sig %>% slice_head(n=20)
kegg_top$fill_padj <- -log10(kegg_top$p.adjust)

kegg_facet_plot <- ggplot(kegg_top, aes(x=reorder(Description, Count),
                                        y=Count,
                                        color=Regulation,
                                        shape=Shape,
                                        size=Count,
                                        fill=fill_padj)) +
  geom_point(stroke=1.2) +
  coord_flip() +
  facet_wrap(~category, scales="free_y", ncol=1) +
  scale_color_manual(values=c("Up in NR"="#EE2C2C", "Down in NR"="#00688B", "NA"="grey")) +
  scale_shape_identity() +
  scale_fill_viridis_c(option="D", name="-log10(adj P)") +
  scale_size(range=c(3,8)) +
  guides(size = guide_legend(override.aes = list(shape=24, fill="grey"))) +
  theme_bw(base_size=12) +
  labs(title="KEGG Pathway Enrichment (NR vs R)", x="KEGG Pathway", y="Gene Count", color="Regulation in NR", size="Gene Count")

# ----------------------
# 7. GSEA: Hallmark pathways
# ----------------------
msig_h <- msigdbr(species="Homo sapiens", collection="H") %>% dplyr::select(gs_name, gene_symbol)
hallmark_pathways <- split(msig_h$gene_symbol, msig_h$gs_name)
fgsea_res <- fgseaMultilevel(pathways = hallmark_pathways, stats = ranked_vec)
fgsea_sig <- fgsea_res %>% filter(padj <= 0.05) %>% arrange(padj)
fgsea_top20 <- head(fgsea_sig, 20)
fgsea_top20$Regulation <- ifelse(fgsea_top20$NES > 0, "Up in NR", "Down in NR")

# Bar plot for NES
gsea_bar <- ggplot(fgsea_top20, aes(x=reorder(pathway, NES), y=NES, fill=Regulation)) +
  geom_col(width=0.7) +
  coord_flip() +
  scale_fill_manual(values=c("Up in NR"="#0072B2", "Down in NR"="#D55E00")) +
  theme_minimal(base_size=12) +
  labs(title="Top Hallmark Pathways (GSEA) with NR Regulation",
       x="Pathway", y="Normalized Enrichment Score (NES)",
       fill="Regulation in NR")

# ----------------------
# 8. Save all plots
# ----------------------
output <- "/LUNG_scRNAseq/Images/New_Plots/"

pdf(paste0(output,"Plasma_C12_DEG.pdf"), height=5, width=5); print(volcano_plot); dev.off()
pdf(paste0(output,"Plasma_C12_GO.pdf"), height=10, width=8); print(go_facet_plot); dev.off()
pdf(paste0(output,"Plasma_C12_KEGG.pdf"), height=10, width=8); print(kegg_facet_plot); dev.off()
pdf(paste0(output,"Plasma_C12_Hallmark.pdf"), height=5, width=8); print(gsea_bar); dev.off()

output <- "/LUNG_scRNAseq/Images/New_Plots/"
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.3] KM plot for Plasma.cells_C12.ES (BCM dataset, with CI)
rm(list=ls())

library(survival)
output <- "/LUNG_scRNAseq/Images/New_Plots/"
myinf1 <- "/Data/BCM_LC_RNAseq.rda"
myFig <- paste0(output, "Fig3_KM_BCM_PlasmaC12_CI.pdf")

load(myinf1)
data <- ORI.data
info <- ORI.info

comxx <- intersect(row.names(data), row.names(info))
data <- data[comxx,]; info <- info[comxx,]

se <- which(colnames(data)=="Plasma.cells_C12.ES")
mys <- as.numeric(data[,se])
xx <- cbind(mys, info)
xx <- xx[xx[, "t.surv"]>0,]
mydat <- cbind(info, mys)
colnames(mydat)[1:2] <- c("time", "event")

mycat <- ifelse(mydat$mys > median(mydat$mys), 1, 0)
xx <- cbind(mycat, mydat)
mytex <- paste(c("High-Score","Low-Score")," (n=",
               c(sum(xx$mycat==1,na.rm=T), sum(xx$mycat==0,na.rm=T)),")",sep="")

tmp <- summary(coxph(Surv(time,event)~mycat, data=xx))$coefficients
myp <- tmp[5]; myh <- round(tmp[2],3)
myp <- ifelse(myp<0.001,formatC(myp,format="e",digits=0),signif(myp,1))
myp <- paste("P=",myp,sep=""); myh <- paste("HR=",myh,sep="")

fit <- survfit(Surv(time,event)~mycat,data=xx)
mycol <- c('lightcoral','skyblue')

pdf(myFig, width=4, height=3)
par(mfrow=c(1,1), lend=2, tcl=-0.15, mar=c(3,3,1,1)+0.1, mgp=c(1.1,0.15,0))
plot(fit, col=mycol[2:1], lwd=2, conf.int=TRUE, ylim=c(0,1),
     xlab="Survival time (months)", ylab="Probability of survival (RFS)",
     cex.axis=0.9, cex.lab=0.9)
legend("bottomleft", mytex, lwd=2, col=mycol, bty="n", cex=0.9)
text(x=max(xx$time)*0.9, y=0.85, labels=myp, cex=0.9, adj=c(1,0))
text(x=max(xx$time)*0.9, y=0.80, labels=myh, cex=0.9, adj=c(1,0))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.2] KM plot for Plasma.cells_C12.ES (Okayama dataset, with CI)
rm(list=ls())

library(survival)
output <- "/LUNG_scRNAseq/Images/New_Plots/"
myinf1 <- "/Data/Okayama_GSE31210_Data.rda"
myFig <- paste0(output, "Fig4_KM_Okayama_PlasmaC12_CI.pdf")

load(myinf1)
data <- ORI.data; info <- ORI.info
comxx <- intersect(row.names(data), row.names(info))
data <- data[comxx,]; info <- info[comxx,]

se <- which(colnames(data)=="Plasma.cells_C12.ES")
mys <- as.numeric(data[,se])
xx <- cbind(mys, info)
xx <- xx[xx[, "t.surv"]>0,]
mydat <- cbind(info, mys)
colnames(mydat)[1:2] <- c("time","event")

mycat <- ifelse(mydat$mys > median(mydat$mys), 1, 0)
xx <- cbind(mycat, mydat)
mytex <- paste(c("High-Score","Low-Score")," (n=",
               c(sum(xx$mycat==1,na.rm=T), sum(xx$mycat==0,na.rm=T)),")",sep="")

tmp <- summary(coxph(Surv(time,event)~mycat, data=xx))$coefficients
myp <- tmp[5]; myh <- round(tmp[2],3)
myp <- ifelse(myp<0.001,formatC(myp,format="e",digits=0),signif(myp,1))
myp <- paste("P=",myp,sep=""); myh <- paste("HR=",myh,sep="")

fit <- survfit(Surv(time,event)~mycat,data=xx)
mycol <- c('lightcoral','skyblue')

pdf(myFig, width=4, height=3)
par(mfrow=c(1,1), lend=2, tcl=-0.15, mar=c(3,3,1,1)+0.1, mgp=c(1.1,0.15,0))
plot(fit, col=mycol[2:1], lwd=2, conf.int=TRUE, ylim=c(0,1),
     xlab="Survival time (days)", ylab="Probability of survival (RFS)",
     cex.axis=0.9, cex.lab=0.9)
legend("bottomleft", mytex, lwd=2, col=mycol, bty="n", cex=0.9)
text(x=max(xx$time)*0.9, y=0.85, labels=myp, cex=0.9, adj=c(1,0))
text(x=max(xx$time)*0.9, y=0.80, labels=myh, cex=0.9, adj=c(1,0))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.4] KM plot for Plasma.cells_C12.ES (Stage I, with CI)
rm(list=ls())

library(survival)
output <- "/LUNG_scRNAseq/Images/New_Plots/"
myinf1 <- "/Data/Okayama_GSE31210_Data.rda"
myFig <- paste0(output, "Fig4_KM_Okayama_PlasmaC12_StageI_CI.pdf")

load(myinf1)
data <- ORI.data; info <- ORI.info
comxx <- intersect(row.names(data), row.names(info))
data <- data[comxx,]; info <- info[comxx,]

se <- which(colnames(data)=="Plasma.cells_C12.ES")
mys <- as.numeric(data[,se])
mydat <- cbind(info, mys)
colnames(mydat)[1:2] <- c("time","event")
mydat <- mydat[mydat$pstage.iorii==" I" & mydat$time>0,]

mycat <- ifelse(mydat$mys > median(mydat$mys), 1, 0)
xx <- cbind(mycat, mydat)
mytex <- paste(c("High-Score","Low-Score")," (n=",
               c(sum(xx$mycat==1,na.rm=T), sum(xx$mycat==0,na.rm=T)),")",sep="")

tmp <- summary(coxph(Surv(time,event)~mycat, data=xx))$coefficients
myp <- tmp[5]; myh <- round(tmp[2],3)
myp <- ifelse(myp<0.001,formatC(myp,format="e",digits=0),signif(myp,1))
myp <- paste("P=",myp,sep=""); myh <- paste("HR=",myh,sep="")

fit <- survfit(Surv(time,event)~mycat,data=xx)
mycol <- c('lightcoral','skyblue')

pdf(myFig, width=4, height=3)
par(mfrow=c(1,1), lend=2, tcl=-0.15, mar=c(3,3,1,1)+0.1, mgp=c(1.1,0.15,0))
plot(fit, col=mycol[2:1], lwd=2, conf.int=TRUE, ylim=c(0,1),
     xlab="Survival time (days)", ylab="Probability of survival (RFS)",
     cex.axis=0.9, cex.lab=0.9)
legend("bottomleft", mytex, lwd=2, col=mycol, bty="n", cex=0.9)
text(x=max(xx$time)*0.9, y=0.85, labels=myp, cex=0.9, adj=c(1,0))
text(x=max(xx$time)*0.9, y=0.80, labels=myh, cex=0.9, adj=c(1,0))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.5] KM plot for Plasma.cells_C12.ES (Stage II, with CI)
rm(list=ls())

library(survival)
output <- "/LUNG_scRNAseq/Images/New_Plots/"
myinf1 <- "/Data/Okayama_GSE31210_Data.rda"
myFig  <- paste0(output, "Fig4_KM_Okayama_PlasmaC12_StageII_CI.pdf")

load(myinf1)
data <- ORI.data; info <- ORI.info
comxx <- intersect(row.names(data), row.names(info))
data <- data[comxx,]; info <- info[comxx,]

se  <- which(colnames(data)=="Plasma.cells_C12.ES")
mys <- as.numeric(data[,se])
mydat <- cbind(info, mys)
colnames(mydat)[1:2] <- c("time","event")
mydat <- mydat[mydat$pstage.iorii==" II" & mydat$time>0,]

mycat <- ifelse(mydat$mys>median(mydat$mys),1,0)
xx <- cbind(mycat,mydat)
mytex <- paste(c("High-Score","Low-Score")," (n=",
               c(sum(xx$mycat==1,na.rm=T),sum(xx$mycat==0,na.rm=T)),")",sep="")

tmp  <- summary(coxph(Surv(time,event)~mycat,data=xx))$coefficients
myp  <- tmp[5]; myh <- round(tmp[2],3)
myp  <- ifelse(myp<0.001,formatC(myp,format="e",digits=0),signif(myp,1))
myp  <- paste("P=",myp,sep=""); myh <- paste("HR=",myh,sep="")

fit  <- survfit(Surv(time,event)~mycat,data=xx)
mycol <- c('lightcoral','skyblue')

pdf(myFig,width=4,height=3)
par(mfrow=c(1,1),lend=2,tcl=-0.15,mar=c(3,3,1,1)+0.1,mgp=c(1.1,0.15,0))
plot(fit,col=mycol[2:1],lwd=2,conf.int=TRUE,ylim=c(0,1),
     xlab="Survival time (days)",ylab="Probability of survival (RFS)",
     cex.axis=0.9,cex.lab=0.9)
legend("bottomleft",mytex,lwd=2,col=mycol,bty="n",cex=0.9)
text(x=max(xx$time)*0.9,y=0.85,labels=myp,cex=0.9,adj=c(1,0))
text(x=max(xx$time)*0.9,y=0.80,labels=myh,cex=0.9,adj=c(1,0))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.6] KM plot for Plasma.cells_C12.ES -- Smokers (with CI)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list=ls())
library(survival)

#-----------------------------
# Paths
output <- "/LUNG_scRNAseq/Images/New_Plots/"
myinf1 <- "/Data/Okayama_GSE31210_Data.rda"
myFig  <- paste0(output, "Fig4_KM_Okayama_PlasmaC12_Smoker_CI.pdf")

#-----------------------------
# Load data
load(myinf1)
data <- ORI.data
info <- ORI.info

comxx <- intersect(row.names(data), row.names(info))
data <- data[comxx,]
info <- info[comxx,]

se  <- which(colnames(data) == "Plasma.cells_C12.ES")
mys <- as.numeric(data[,se])
mydat <- cbind(info, mys)
colnames(mydat)[1:2] <- c("time", "event")

#-----------------------------
# Filter for smokers
cat("Unique smoking.status categories:\n")
print(unique(as.character(mydat$smoking.status)))

# Adjust this string if needed (check printed categories)
mydat <- mydat[which(as.character(mydat$smoking.status) == " Ever-smoker" & mydat$time > 0), ]

if(nrow(mydat) < 2){
  stop("No or too few smoker samples after filtering. Check 'smoking.status' values.")
}

#-----------------------------
# KM + COX
mycat <- ifelse(mydat$mys > median(mydat$mys, na.rm=TRUE), 1, 0)
xx <- cbind(mycat, mydat)
mytex <- paste(c("High-Score","Low-Score")," (n=",
               c(sum(xx$mycat==1,na.rm=TRUE),sum(xx$mycat==0,na.rm=TRUE)),")",sep="")

fit.cox <- coxph(Surv(time,event) ~ mycat, data=xx)
tmp <- summary(fit.cox)$coefficients
myp <- tmp[5]
myh <- round(tmp[2],3)
myp <- ifelse(myp<0.001,formatC(myp,format="e",digits=0),signif(myp,1))
myp <- paste("P=",myp,sep="")
myh <- paste("HR=",myh,sep="")

fit <- survfit(Surv(time,event) ~ mycat, data=xx)
mycol <- c('lightcoral','skyblue')

#-----------------------------
pdf(myFig, width=4, height=3)
par(mfrow=c(1,1), lend=2, tcl=-0.15, mar=c(3,3,1,1)+0.1, mgp=c(1.1,0.15,0))

plot(fit, col=mycol[2:1], lwd=2, conf.int=TRUE,
     ylim=c(0,1), xlab="Survival time (days)", ylab="Probability of survival (RFS)",
     cex.axis=0.9, cex.lab=0.9)

max.x <- max(xx$time, na.rm=TRUE)
min.x <- min(xx$time, na.rm=TRUE)
legend(min.x, 0.2, mytex, lwd=2, col=mycol, bty="n", cex=0.9)
text(max.x*0.9, 0.85, myp, cex=0.9, adj=c(1,0))
text(max.x*0.9, 0.80, myh, cex=0.9, adj=c(1,0))
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.7] KM plot for Plasma.cells_C12.ES -- Never-Smokers (with CI)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list=ls())
library(survival)
#-----------------------------
# Paths
output <- "/LUNG_scRNAseq/Images/New_Plots/"
myinf1 <- "/Data/Okayama_GSE31210_Data.rda"
myFig  <- paste0(output, "Fig4_KM_Okayama_PlasmaC12_Nonsmoker_CI.pdf")

#-----------------------------
# Load data
load(myinf1)
data <- ORI.data
info <- ORI.info

comxx <- intersect(row.names(data), row.names(info))
data <- data[comxx,]
info <- info[comxx,]

se  <- which(colnames(data) == "Plasma.cells_C12.ES")
mys <- as.numeric(data[,se])
mydat <- cbind(info, mys)
colnames(mydat)[1:2] <- c("time", "event")

#-----------------------------
# Filter for never-smokers
cat("Unique smoking.status categories:\n")
print(unique(as.character(mydat$smoking.status)))

# Adjust this string if needed (check printed categories)
mydat <- mydat[which(as.character(mydat$smoking.status) == " Never-smoker" & mydat$time > 0), ]

if(nrow(mydat) < 2){
  stop("No or too few never-smoker samples after filtering. Check 'smoking.status' values.")
}

#-----------------------------
# KM + COX
mycat <- ifelse(mydat$mys > median(mydat$mys, na.rm=TRUE), 1, 0)
xx <- cbind(mycat, mydat)
mytex <- paste(c("High-Score","Low-Score")," (n=",
               c(sum(xx$mycat==1,na.rm=TRUE),sum(xx$mycat==0,na.rm=TRUE)),")",sep="")

fit.cox <- coxph(Surv(time,event) ~ mycat, data=xx)
tmp <- summary(fit.cox)$coefficients
myp <- tmp[5]
myh <- round(tmp[2],3)
myp <- ifelse(myp<0.001,formatC(myp,format="e",digits=0),signif(myp,1))
myp <- paste("P=",myp,sep="")
myh <- paste("HR=",myh,sep="")

fit <- survfit(Surv(time,event) ~ mycat, data=xx)
mycol <- c('lightcoral','skyblue')

#-----------------------------
pdf(myFig, width=4, height=3)
par(mfrow=c(1,1), lend=2, tcl=-0.15, mar=c(3,3,1,1)+0.1, mgp=c(1.1,0.15,0))

plot(fit, col=mycol[2:1], lwd=2, conf.int=TRUE,
     ylim=c(0,1), xlab="Survival time (days)", ylab="Probability of survival (RFS)",
     cex.axis=0.9, cex.lab=0.9)

max.x <- max(xx$time, na.rm=TRUE)
min.x <- min(xx$time, na.rm=TRUE)
legend(min.x, 0.2, mytex, lwd=2, col=mycol, bty="n", cex=0.9)
text(max.x*0.9, 0.85, myp, cex=0.9, adj=c(1,0))
text(max.x*0.9, 0.80, myh, cex=0.9, adj=c(1,0))
dev.off()














