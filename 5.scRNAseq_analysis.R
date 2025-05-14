# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [1] define marker gene profiles based on scRNA-seq data for infer cell ratio
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.1] define marker gene signatures with weights
# Chao data: /mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(ggplot2)
library(NMF)
library(reshape2)
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Immune_cells_data/cd45pos_nsclc_Immune_cells.rds"
myinf2 = "/mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/NSCLC_output/cd45negc_celltype_data.rds"

myoutf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/LungCancer_CellMarker_AmosscRNAseq_4base.txt"

CD45pos = readRDS(myinf1)
CD45neg = readRDS(myinf2)

info.pos = CD45pos@meta.data
info.neg = CD45neg@meta.data
# se = which(info.neg$Cell_type %in% c("AT2", "Ciliated cells", "Clara", "BASC"))
# info.neg$Cell_type[se] = "Epithelial"

se = which(info.neg$Cell_type1 %in% c("BASC-C8",  "BASC-C10"))
info.neg$Cell_type1[se] = "BASC"

info.pos = cbind(info.pos[, 1:6], Cell_type= info.pos$Cell_type1)
info.neg = cbind(info.neg[, 1:6], Cell_type= info.neg$Cell_type1) # Chao used info.neg$Cell_type I used info.neg$Cell_type1 with making all BASC in to one group
info = rbind(info.pos, info.neg)

dat.pos <- GetAssayData(CD45pos, layer = "count")
dat.pos = as.matrix(dat.pos)
dat.neg <- GetAssayData(CD45neg, layer = "count")
dat.neg = as.matrix(dat.neg)
data = cbind(dat.pos, dat.neg)
head(data)
xx = apply(data, 2, sum)
for(k in 1:ncol(data)) { # Combines datasets and normalizes each cell by total read count, scaling to counts per million (CPM)
  cat("\r", k)
  data[,k] = data[,k]*1e6/xx[k]
}
xx = apply(data>0, 1, sum) # Retains genes expressed in at least 100 cells and applies a log2 transformation
se = which(xx>=100)
data = data[se,]
dim(data)
data = log2(data+1)

xx = info$Cell_type
xx = gsub("φ", "ac", xx)
xx = gsub("γδ", "GamDel", xx)
xx = gsub("\\|", " ", xx)
xx = gsub(" ", ".", xx)
info$Cell_type = xx
table(info$Cell_type)

#-------------------------
myfac = info$Cell_type
xx = sort(table(myfac))
mycel = names(xx)[xx>=100] # Identifies cell types with at least 100 cells for downstream analysis

[1] "DC_C18"                "Mast_C17"              "Prolifirating.Mac_C16"
[4] "Fibroblasts"           "CD8_C15"               "Endothelial.cells"
[7] "GamDelT_C14"           "B.cells_C13"           "Plasma.cells_C12"
[10] "Mac.DC_C11"            "CD4_C10"               "Plasma.cells_C9"
[13] "Mac_C8"                "CD4_C7"                "CD4_C6"
[16] "Neutrophil_C5"         "NK.cells_C4"           "Mono_C3"
[19] "B.cells_C2"            "CD8_C1"                "CD8_C0"
[22] "Epithelial"


# assign each cell with associated gene expression (gene expression is mean of all cells in the 22 single cell type)
# Generating Pseudo-Cells and their expression data for 100 samples*expression for each cell type.
mydat = NULL
mylab = NULL
ncel = 100
nsiz = 20

for(k in 1:length(mycel))
{
  mylab = c(mylab, rep(mycel[k], ncel)) # mylab representing the repeat the colname for each cell for 100 times
  subset = which(myfac==mycel[k])
  subdat = data[, subset]
  for(i in 1:ncel)
  {
    cat("\r\r\r", k, "-->", i)
    se = sample(1:ncol(subdat), nsiz)
    tmp = apply(subdat[, se], 1, mean) # calculate the mean or average for each gene
    mydat = cbind(mydat, tmp)
  }
}
row.names(mydat) = row.names(data)
dim(mydat)
# For each cell type, 100 pseudo-cells are created.
# The dimensions of mydat grow to (number of genes) × (number of pseudo-cells).

#-------------------------
myfac = mylab
mycat = unique(myfac)
xx = apply(mydat>0, 1, sum)
se = which(xx>=10)
mydat = mydat[se,]

tmp = matrix(0, nrow(mydat), length(mycat))
row.names(tmp) = row.names(mydat)
colnames(tmp) = mycat
res.PV = res.TS = tmp # Differential Expression Analysis save PV p-values and TS T-scores, compare cell type 1 with all other cells.

for(k in 1:ncol(res.PV))
{
  cat("\r", k)
  se = which(myfac==mycat[k])
  dat1 = mydat[,se]
  dat2 = mydat[, -se]
  myavg1 = apply(dat1, 1, mean)
  myavg2 = apply(dat2, 1, mean)
  myvar1 = apply(dat1, 1, var)
  myvar2 = apply(dat2, 1, var)
  n1 = ncol(dat1)
  n2 = ncol(dat2)
  fc = myavg1-myavg2
  tscore = (myavg1-myavg2)/sqrt(myvar1/n1+myvar2/n2)
  df= (myvar1/n1+myvar2/n2)^2/((myvar1/n1)^2/(n1-1) + (myvar2/n2)^2/(n2-1))
  pval = pt(-abs(tscore), df)*2
  res.TS[,k] = tscore
  res.PV[,k] = pval
}

#-------------------------
# Converts p-values to a log10 scale.
# Filters insignificant results (pval > 0.001) and caps highly significant values


# res.TS is the matrix of t-scores for the same differential expression test.
# A negative t-score indicates that the gene is less expressed in the given cell type compared to others.
# Genes with negative t-scores are not of interest here, so their p-values are set to 0.

# Only genes with a -log10 p-value of 3 or higher are retained (equivalent to a p-value < 0.001).
# This ensures only highly significant results are considered.

PV = -log10(res.PV)
for(k in 1:ncol(PV))
{
  cat("\r", k)
  tmp = PV[,k]
  tmp[res.TS[,k]<0]=0
  tmp[tmp<3] = 0
  tmp[tmp>100] = 100
  PV[,k] = tmp
}

# Divides all values in the PV matrix by 100, scaling them to a range between 0 and 1:
#   Values close to 1 indicate very significant results.
# Values close to 0 indicate less significant or filtered-out results.

PV = PV/100
dim(PV)
# [1] 20223    22
 write.table(PV, myoutf1, sep="\t", quote=F)
# calculated cell type expression data genes expression of all immune and non-immune cells
heatmap(PV)

library(pheatmap)
pdf("/home/u251079/scRNA/Lung_scRNA/Heatmap_sig_Pseudobulk_PV.pdf", width = 8, height = 10)
pheatmap(PV, 
         scale = "none",         
         cluster_rows = T,    
         cluster_cols = T,  
         show_rownames = FALSE,
         show_colnames = F,
         treeheight_col= 0,treeheight_row=0, legend = F, fontsize = 8)    # Show column names
dev.off() 



# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [2] BASE analysis
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[2.1] Okayama_GSE31210
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/PRECOG/processed/GSE31210.HGU133Plus2_EntrezCDF.MAS5.pcl"
myinf2 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/LungCancer_CellMarker_AmosscRNAseq_4base.txt"

myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Okayama_GSE31210_AmoscRNAMarkers_Base.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)
reg = read.table(myinf2, sep="\t", header=T, row.names=1)

source("/home/u251079/BLCA_code/base5.R")

##
xx = base5(data, reg, perm=1000, myoutf, median.norm=F)

# data=log2(data+1)
# xx = rowSums(data) <= 0.5
# data=data[!xx,]
# 
# color_palette <- colorRampPalette(c("red","white","blue"))(100)
# pheatmap(data, 
#          scale = "none",         
#          cluster_rows = T,    
#          cluster_cols = T,  
#          show_rownames = FALSE,
#          color = color_palette,
#          show_colnames = F,
#          treeheight_col= 0,treeheight_row=0, legend = F)




# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1.2] BCM_LC_RNAseq # BCM bulk RNAseq data
# rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/BCM_prot_gene_TPM_Symbol.csv"
myinf2 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/LungCancer_CellMarker_AmosscRNAseq_4base.txt"

myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/BCM_LC_RNAseq_AmoscRNAMarkers_Base.txt"

data = read.csv(myinf1, header=T, row.names=1, check.names=F)

data = log2(data+1)
dim(data)
# Use this code in downstream to extract T samples
# se = grep("T", colnames(data)) 
# data = data[,se]

reg = read.table(myinf2, sep="\t", header=T, row.names=1)

xx = base5(data, reg, perm=1000, myoutf, median.norm=F)

 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1.3] Sato_GSE41271
# rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Sato_GSE41271/Sato_GSE41271_Lung_Symbol_expr.txt"
myinf2 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/LungCancer_CellMarker_AmosscRNAseq_4base.txt"

myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Sato_GSE41271_AmoscRNAMarkers_Base.txt"

mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

#source("~/WorSpa/system/myRprogram/base5.R")
data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=F)



# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1.4]  Rousseaux_GSE30219 (n=85)
#https://link.springer.com/article/10.1007/s13402-023-00883-w
# A total of 611 patients with LUAD were included from the GSE30219 (n=85), GSE50081 (n=128), and GSE72094 (n=398) cohorts in the Gene Expression Omnibus (GEO) database
# GSE30219 have 85 lung adenocarcinoma samples, selected based on !Sample_characteristics_ch1: histology: ADC, have patients with relapse and non-relapse data with RFS information 
# /mount/ictr1/chenglab/venu/scRNAseq_lung/cellranger_scran/Other_data
Chao # /mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Rousseaux_GSE30219
file_path1 <- "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Rousseaux_GSE30219/Rousseaux_GSE30219_GPL570_symbol.txt"
myinf2 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Rousseaux_GSE30219_AmoscRNAMarkers_Base.txt"

mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
data <- read.table(file_path1, sep="\t", header=T, row.names=1)
head(data)

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=F)

# Clinical info
file_path2 <- "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Rousseaux_GSE30219/Sample_info.txt"
info = read.table(file_path2, sep="\t", header=T, row.names=1, quote="")






# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1.5] Der_data_GSE50081 
this data has 128 patients with Lung adenocarcinoma samples with relapse and RFS information with gene expression data
!Sample_characteristics_ch1 histology: adenocarcinoma
# /mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Der_GSE50081
file_path1 <- "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Der_GSE50081/Der_GSE50081_GPL570_symbol.txt"
myinf2 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Der_GSE50081_AmoscRNAMarkers_Base.txt"

mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
data <- read.table(file_path1, sep="\t", header=T, row.names=1)
head(data)

##
reg = mywt
xx = base5(data, reg, perm=1000, myoutf, median.norm=F)

# Clinical info
file_path2 <- "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Der_GSE50081/Sample_info.txt"
info = read.table(file_path2, sep="\t", header=T, row.names=1, quote="")

# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2.1.6] Schabath_data_GSE72094 
This data has OS, time of event, days, smoking, KRAS mutation-associated gene expression, p53 and STK11 mutations, proliferation and immune surveillance in lung adenocarcinoma
kras_status	egfr_status	stk11_status	tp53_status
# /mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Schabath_GSE72094
gse72094_data <- read.table(paste0(file_path, "GSE72094_series_matrix.txt"), header = TRUE, sep = "\t")

head(gse72094_data)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.2.1] relpase-free survival or disease free survival
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Okayama_GSE31210_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Clinical/GSE31210_Clinical_info.txt"
#myoutf = "/home/u251079/scRNA/Lung_scRNA/Okayama_GSE31210_AmoscRNAMarkers_Base.txt"
data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

# color_palette <- colorRampPalette(c("#CD2626","white","#698B69"))(100)
# pheatmap(data, 
#          scale = "none",         
#          cluster_rows = T,    
#          cluster_cols = T,  
#          show_rownames = FALSE,
#          color = color_palette,
#          show_colnames = F,
#          treeheight_col= 0,treeheight_row=0)






info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$days.before.relapse.censor)
e.surv = ifelse(info$relapse ==" relapsed", 1, 0)
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]

# info=info[info$t.surv <= 1825,]
# info$t.surv
library(survival)

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

info_out<-"/home/u251079/scRNA/Lung_scRNA/Okayama.txt"
Okayama<-merge(data, info, by = "row.names")
row.names(Okayama)<-Okayama$Row.names
Okayama$Row.names<-NULL
write.table(Okayama, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)


survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1, lb1, ub1)
xx = res[order(res[,2]),]
xx[1:20,]


name  coxph.pval1 coxph.qval1       hr1
9       Plasma.cells_C12.ES 5.502311e-05 0.001210508 0.5728493
6      Endothelial.cells.ES 8.444327e-04 0.009288759 0.6203990
22            Epithelial.ES 6.795041e-03 0.049830303 0.7208717
19            B.cells_C2.ES 9.414795e-03 0.051781373 0.7715797
2               Mast_C17.ES 2.087167e-02 0.091835368 0.7926445
11               CD4_C10.ES 3.542043e-02 0.129874921 0.8533657
3  Prolifirating.Mac_C16.ES 6.612741e-02 0.207828987 1.2424546
15                CD4_C6.ES 9.092785e-02 0.250051580 0.8931429
7            GamDelT_C14.ES 1.511977e-01 0.369594384 1.1786485
20                CD8_C1.ES 2.992099e-01 0.658261784 0.9004017
1                 DC_C18.ES 3.741874e-01 0.748374854 0.8885928
10            Mac.DC_C11.ES 4.385858e-01 0.759255677 1.1098391
8            B.cells_C13.ES 5.056100e-01 0.759255677 0.9183909
4            Fibroblasts.ES 5.703378e-01 0.759255677 0.9578591
18               Mono_C3.ES 5.764718e-01 0.759255677 1.0727918
16         Neutrophil_C5.ES 5.841697e-01 0.759255677 1.0742327
5                CD8_C15.ES 5.866976e-01 0.759255677 0.9440011
14                CD4_C7.ES 7.177926e-01 0.877302055 0.9837151
13                Mac_C8.ES 8.979174e-01 0.988509757 0.9833984
21                CD8_C0.ES 9.231322e-01 0.988509757 0.9916137


## without re-normalization
name  coxph.pval1 coxph.qval1       hr1
9       Plasma.cells_C12.ES 0.0001318596 0.002900910 0.5939028
15                CD4_C6.ES 0.0006684360 0.006493808 0.6936179
6      Endothelial.cells.ES 0.0008855193 0.006493808 0.6301796
2               Mast_C17.ES 0.0012192206 0.006581071 0.8059368
8            B.cells_C13.ES 0.0014956979 0.006581071 0.6838933


#@@@@@@@@@@@@@@@@@@@@@. PLOTS
library(ggplot2)
rownames(xx)
names(xx)
xx$name
xx$logPval <- -log10(xx$coxph.qval1)
xx$Significance <- ifelse(xx$coxph.qval1 < 0.05, "Significant", "Not Significant")
Okayama1=xx #[!c(xx$name %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]
Okayama1$name <- gsub("\\.ES", "", Okayama1$name)

Plot_Okayama<-ggplot(Okayama1, aes(x = hr1, y = logPval, color = Significance)) +
  geom_point(size = 2) +  # Main scatter points
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.7) +  # Red horizontal line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.7) +  # Vertical line at HR=1
  scale_color_manual(values = c("Significant" = "steelblue", "Not Significant" = "black")) +  # Color settings
  geom_text(data = subset(Okayama1, Significance == "Significant"), 
            aes(label = name), vjust = -0.5, size = 4, check_overlap = TRUE) +  # Label only significant points
  theme_classic() + scale_x_continuous(limits = c(0.5, max(Okayama1$hr1)+0.1)) + 
  scale_y_continuous(limits = c(0, 3.5)) + 
  labs(title = NULL,
       x = "Hazard Ratio",
       y = "-log10(adj p-value)") +
  theme(legend.position = c(0.85,0.8))

Plot_Okayama
min(Okayama1$hr1)
max(Okayama1$hr1)

#ggsave("/home/u251079/scRNA/Lung_scRNA/Okayama_volcano.png", plot = Plot_Okayama, width = 6, height = 5, dpi = 1200)

pdf("/home/u251079/scRNA/Lung_scRNA/Okayama_volcano.pdf", width = 5, height =4)
print(Plot_Okayama)
dev.off() 






#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.2.1.2] compare relapsed with non-relapsed
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Okayama_GSE31210_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Clinical/GSE31210_Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

## 
library(ROCR)
se = which(info$relapse== " relapsed")
dat1 = data[se, ]
se = which(info$relapse== " not relapsed")
dat2 = data[se, ]

res = matrix(0, ncol(data), 5)
row.names(res) = colnames(data)
colnames(res) = c("Avg.R", "Avg.NR", "TS", "PV.t", "PV.w")
res = as.data.frame(res)
cat = c(rep("R", nrow(dat1)), rep("NR", nrow(dat2)))
myauc = rep(0, ncol(data))

res[,1] = apply(dat1, 2, mean)
res[,2] = apply(dat2, 2, mean)
for(k in 1:ncol(data))
{
  cat("\r", k)
  
  tmp = t.test(dat1[,k], dat2[,k])
  res[k,3] = tmp$statistic
  res[k,4] = tmp$p.value
  tmp = wilcox.test(dat1[,k], dat2[,k])
  res[k,5] = tmp$p.value
  score = c(dat1[,k], dat2[,k])
  pred <- prediction(score, cat)
  auc.perf <- performance(pred, measure = "auc")
  myauc[k] = as.numeric(auc.perf@y.values)
}
myres = cbind(myauc, res)
xx = myres[order(myres[,6]), ]
xx[1:20, ]
xx[1:20, c(1,6)]

myauc        Avg.R      Avg.NR          TS
Endothelial.cells.ES     0.3252315   1.77634193   2.1789166 -4.08406255
Plasma.cells_C12.ES      0.3341049 -12.62359133 -12.0883547 -3.99682814
Mast_C17.ES              0.3464506  -2.61137269  -1.9018124 -2.61487801
B.cells_C2.ES            0.3737461   0.56216660   1.0951550 -2.54449478
Epithelial.ES            0.3863812   4.76732825   5.1224142 -2.40306197
B.cells_C13.ES           0.3937114   0.62816130   0.7345755 -0.67978563
CD4_C6.ES                0.3947724   0.46265560   0.9692033 -1.71550306
CD8_C15.ES               0.4059606   0.22613309   0.3354540 -0.60965891
NK.cells_C4.ES           0.4111690   0.23218547   0.2455129 -0.07272991
CD4_C10.ES               0.4140625   0.83262752   1.3641112 -2.04595228
Prolifirating.Mac_C16.ES 0.5847801  16.72830775  16.4817633  1.62477300
GamDelT_C14.ES           0.5847801   5.49211051   5.3022852  1.34427209
CD8_C0.ES                0.4313272   0.56396848   0.6013689 -0.17055433
CD4_C7.ES                0.4316165  -0.01515664   0.1160159 -0.30251637
CD8_C1.ES                0.4485918   0.84893269   1.0619667 -1.05399451
Mac.DC_C11.ES            0.5445602  11.35306003  11.2155966  1.00326452
Fibroblasts.ES           0.4595872   1.79943579   1.9660010 -0.73368738
Plasma.cells_C9.ES       0.4630594   6.43095150   6.4577636 -0.15661871
Mono_C3.ES               0.5340471   9.76677224   9.6988148  0.45009505
DC_C18.ES                0.4689429  12.88977915  13.0187145 -0.86194847
PV.t         PV.w
Endothelial.cells.ES     6.694482e-05 4.305816e-05
Plasma.cells_C12.ES      1.088591e-04 1.033085e-04
Mast_C17.ES              9.932676e-03 3.260667e-04
B.cells_C2.ES            1.235056e-02 3.130385e-03
Epithelial.ES            1.786103e-02 7.841174e-03
B.cells_C13.ES           4.980369e-01 1.287447e-02
CD4_C6.ES                8.904218e-02 1.380166e-02
CD8_C15.ES               5.432134e-01 2.777479e-02
NK.cells_C4.ES           9.421466e-01 3.766159e-02
CD4_C10.ES               4.322716e-02 4.434891e-02
Prolifirating.Mac_C16.ES 1.069659e-01 4.729135e-02
GamDelT_C14.ES           1.811144e-01 4.729135e-02
CD8_C0.ES                8.648800e-01 1.081443e-01
CD4_C7.ES                7.628382e-01 1.096394e-01
CD8_C1.ES                2.942630e-01 2.292077e-01
Mac.DC_C11.ES            3.176155e-01 2.973737e-01
Fibroblasts.ES           4.645901e-01 3.446657e-01
Plasma.cells_C9.ES       8.758944e-01 3.877541e-01
Mono_C3.ES               6.535622e-01 4.260573e-01
DC_C18.ES                3.906940e-01 4.678648e-01

#@@@@@@@@@@@@@@@@@@@@ PLOT  Eventhough low AUC it has strongr prediction power 
library(ggplot2)
Okayama2=xx[!c(rownames(xx) %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]
Okayama_auc<- data.frame(Feature = rownames(Okayama2), AUC = Okayama2[, 1], PV=Okayama2$PV.t)

# Plot the AUC values
plot2_Okayama<-ggplot(Okayama_auc, aes(x = reorder(Feature, AUC), y = AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(PV,4)), vjust = 0.5, hjust = 1, size = 4, color="white") + 
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "AUC in Okayama patient samples", x = "Feature", y = "AUC") +
  theme_minimal()
print(plot2_Okayama)

# ROC_AUC curves 
xx=xx[!c(rownames(xx) %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]
xx_significant = xx#[xx[, 5] < 0.05, ]

library(RColorBrewer)

plot(NULL, xlim=c(0, 1), ylim=c(0, 1), xlab="False Positive Rate", ylab="True Positive Rate", main="Okayama ROC Curves with AUC")
dark_colors <- brewer.pal(nrow(xx_significant), "Dark2")  # "Dark2" provides a good set of dark colors
legend_labels <- list()

for(k in 1:nrow(xx_significant)) {
  score = c(dat1[,k], dat2[,k])
  pred <- prediction(score, cat)
  perf <- performance(pred, "tpr", "fpr")
  lines(perf@x.values[[1]], perf@y.values[[1]], col=dark_colors[k], lwd=2)
  auc_value = round(myauc[k], 2)
  pv_value = round(xx_significant[k, "PV.t"], 5)  # rounding PV.t value

  legend_labels[[k]] <- paste(rownames(xx_significant)[k], "AUC =", auc_value, "PV.t =", pv_value)
}
legend("bottomright", legend = unlist(legend_labels), col=dark_colors, lwd=2, cex=0.35, box.lwd=0)



library(ROCR)
se_relapsed <- which(info$relapse == " relapsed")
se_not_relapsed <- which(info$relapse == " not relapsed")
dat_relapsed <- data[se_relapsed, ]
dat_not_relapsed <- data[se_not_relapsed, ]

labels <- c(rep(1, nrow(dat_relapsed)), rep(0, nrow(dat_not_relapsed)))
combined_data <- rbind(dat_relapsed, dat_not_relapsed)
roc_curves <- list()
auc_values <- numeric(length(colnames(data)))
names(auc_values) <- colnames(data)

for (i in seq_along(colnames(data))) {
  cell_type <- colnames(data)[i]
  cat("Processing:", cell_type, "\n")
  scores <- combined_data[, cell_type]
  pred <- prediction(scores, labels)
  perf <- performance(pred, "tpr", "fpr")
  roc_curves[[cell_type]] <- perf
  auc_perf <- performance(pred, measure = "auc")
  auc_values[cell_type] <- auc_perf@y.values[[1]]
}

# Plot all ROC curves
plot(roc_curves[[1]], col = 1, lwd = 2, main = "ROC Curves for All Cell Types", xlab = "False Positive Rate (FPR)", ylab = "True Positive Rate (TPR)")
for (i in 2:length(roc_curves)) {
  plot(roc_curves[[i]], col = i, lwd = 2, add = TRUE)
}
legend_labels <- paste0(names(roc_curves), " (AUC = ", round(auc_values, 3), ")")
legend("bottomright", legend = legend_labels, col = 1:length(roc_curves), lwd = 2, cex = 0.35)

abline(a = 0, b = 1, lty = 2, col = "black")








#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.2.2] relpase-free survival
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/BCM_LC_RNAseq_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/Clinical_info_20231111.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}



info = read.table(myinf2, sep="\t", header=T, row.names=3, quote="")
# Standardize the names
info$Histologic.type..WHO.2015.edition <- tolower(info$Histologic.type..WHO.2015.edition) # Convert to lowercase
info$Histologic.type..WHO.2015.edition <- trimws(info$Histologic.type..WHO.2015.edition) # Remove extra spaces

# Replace inconsistent names
info$Histologic.type..WHO.2015.edition <- recode(info$Histologic.type..WHO.2015.edition,
                                                 "adenocarcinoma" = "Adenocarcinoma",
                                                 "mucinous adenocarcinoma" = "Mucinous Adenocarcinoma",
                                                 "muncinous adenocarcinoma" = "Mucinous Adenocarcinoma",
                                                 "squamos cell carcinoma" = "Squamous Cell Carcinoma",
                                                 "squamos cell carcinoma" = "Squamous Cell Carcinoma",
                                                 "squamous cell carcinoma" = "Squamous Cell Carcinoma",
                                                 "squamous cell carcinoma" = "Squamous Cell Carcinoma",
                                                 "squamous cell carcinoma" = "Squamous Cell Carcinoma"
)

# Verify the cleaned categories
table(info$Histologic.type..WHO.2015.edition)
se = grep("Adenocarcinoma|adenocarcinoma", info$Histologic.type..WHO.2015.edition)
info = info[se,]
se = grep("N", rownames(info))
se = grep("T", rownames(info))
se = which(info$Pstage!="0") # Normal samples 
info = info[se,]
xx = tolower(info$Pstage)
xx[grep("ii", xx)] = "II"
xx[grep("i", xx)] = "I"
info$Pstage = xx
t.surv = info$OS.month_09232023 # RFS or DFS
e.surv = info$Death # Recurrence 0 or 1
info = cbind(t.surv, e.surv, info)

#@@@@@@@@@@@
library(survival)
info$Recurrence_Status <- ifelse(info$DateRecurrence != "", "Recurrence", "Non-Recurrence")
data$patient_status<-ifelse(row.names(data) %in% row.names(info), info$Recurrence_Status, "Normal")
str(data)
library(ggpubr)
# Convert patient_status to a factor with "Normal" as reference
data$patient_status <- factor(data$patient_status, levels = c("Normal", "Non-Recurrence", "Recurrence"))
es_scores <- colnames(data)[grepl("\\.ES$", colnames(data))]
p_box_list <- lapply(es_scores, function(es) {
  ggboxplot(data, 
            x = "patient_status", 
            y = es, 
            fill = "patient_status", 
            palette = c("gray", "lightblue", "lightpink")) +
    stat_compare_means(method = "wilcox.test", ref.group = "Normal", label = "p") + 
    labs(title = es, y = "Signature Scores", x = "Patient Status") +
    theme(plot.title = element_text(size = 12), legend.position = "none")  
})
p_box_combined <- ggarrange(plotlist = p_box_list, ncol = 5, nrow = ceiling(length(es_scores) / 5))
p_box_combined
#@@@@@@@@@@@


# Tumor samples
se = grep("N", rownames(info))
data = data[se,]
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

#@@@@@@@@@@@
info_out<-"/home/u251079/scRNA/Lung_scRNA/BCMLC.txt"
BCMLC<-merge(data, info, by = "row.names")
row.names(BCMLC)<-BCMLC$Row.names
BCMLC$Row.names<-NULL
write.table(BCMLC, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)
#@@@@@@@@@@@

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1,  lb1,  ub1)
xx = res[order(res[,2]),]
xx

name coxph.pval1 coxph.qval1       hr1
2               Mast_C17.ES 0.005869773   0.1291350 0.5967699
15                CD4_C6.ES 0.038026994   0.2463960 0.6758567
19            B.cells_C2.ES 0.039240320   0.2463960 0.6791505
3  Prolifirating.Mac_C16.ES 0.057201387   0.2463960 1.6444648
11               CD4_C10.ES 0.061646617   0.2463960 0.6951698
5                CD8_C15.ES 0.067198920   0.2463960 0.6887156
8            B.cells_C13.ES 0.090701129   0.2625696 0.7009899
14                CD4_C7.ES 0.095479869   0.2625696 0.6918552
7            GamDelT_C14.ES 0.121732168   0.2878691 1.4572083
21                CD8_C0.ES 0.130849590   0.2878691 0.7507327
17           NK.cells_C4.ES 0.177053285   0.3396037 0.8068485
20                CD8_C1.ES 0.185238406   0.3396037 0.7546740
9       Plasma.cells_C12.ES 0.472079537   0.7533452 0.9097356
6      Endothelial.cells.ES 0.500120381   0.7533452 0.8899001
1                 DC_C18.ES 0.513644474   0.7533452 0.9032741
13                Mac_C8.ES 0.589694381   0.7767172 1.1277455
12       Plasma.cells_C9.ES 0.604735644   0.7767172 1.0982935
16         Neutrophil_C5.ES 0.651046507   0.7767172 0.9267111
18               Mono_C3.ES 0.670801258   0.7767172 1.0897766
22            Epithelial.ES 0.727341721   0.8000759 1.1381290
#@@@@@@@@@@@@ PLOTS
# Volcano plot for Cox p-values vs Hazard Ratios
xx$logPval <- -log10(xx$coxph.qval1)
xx$Significance <- ifelse(xx$coxph.qval1 < 0.05, "Significant", "Not Significant")

BCMLC1=xx #[!c(xx$name %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]

Plot1_BCMLC<-ggplot(BCMLC1, aes(x = reorder(name, hr1), y = hr1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lb1, ymax = ub1), width = 0.2) +
  geom_text(aes(label = format(coxph.pval1, scientific = TRUE, digits = 2)), 
            position = position_dodge(width = 0.5), vjust = -0.5) +  # Position text slightly above points
  coord_flip() +
  labs(title = "BCMLC data", 
       x = NULL, 
       y = "Hazard Ratio (HR)") +
  theme_minimal()
print(Plot1_BCMLC)


Plot_BCMLC<-ggplot(BCMLC1, aes(x = hr1, y = logPval, color = Significance)) +
  geom_point(size = 4) +  # Main scatter points
  theme_classic()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +  # Red horizontal line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 1) +  # Vertical line at HR=1
  scale_color_manual(values = c("Significant" = "darkred", "Not Significant" = "navyblue")) +  # Color settings
  geom_text(aes(label = name), vjust = -0.5, size = 4, check_overlap = TRUE) +  # Label only significant points
  theme_classic() + scale_x_continuous(limits = c(0.4, max(BCMLC1$hr1)+0.1)) + 
  scale_y_continuous(limits = c(0, 2)) + 
  labs(title = "BCMLC_data",
       x = "Hazard Ratio",
       y = "-log10(P-value)") +
  theme(legend.position = "top")

Plot_BCMLC





#------------------------------
## 
library(ROCR)
se = which(info$DateRecurrence != "")
dat1 = data[se, ]
dat2 = data[-se, ]

res = matrix(0, ncol(data), 5)
row.names(res) = colnames(data)
colnames(res) = c("Avg.R", "Avg.NR", "TS", "PV.t", "PV.w")
res = as.data.frame(res)
cat = c(rep("R", nrow(dat1)), rep("NR", nrow(dat2)))
myauc = rep(0, ncol(data))

res[,1] = apply(dat1, 2, mean)
res[,2] = apply(dat2, 2, mean)
for(k in 1:ncol(data))
{
  cat("\r", k)
  
  tmp = t.test(dat1[,k], dat2[,k])
  res[k,3] = tmp$statistic
  res[k,4] = tmp$p.value
  tmp = wilcox.test(dat1[,k], dat2[,k])
  res[k,5] = tmp$p.value
  score = c(dat1[,k], dat2[,k])
  pred <- prediction(score, cat)
  auc.perf <- performance(pred, measure = "auc")
  myauc[k] = as.numeric(auc.perf@y.values)
}
myres = cbind(myauc, res)
xx = myres[order(myres[,6]), ]



myauc      Avg.R     Avg.NR            TS
Plasma.cells_C12.ES      0.2842713 -1.9455010 -0.3885687 -2.1763781897
CD8_C15.ES               0.3318903 -1.4044547 -0.5876850 -1.8691468307
CD8_C0.ES                0.3391053 -1.4132005 -0.6545089 -1.5675138698
CD4_C6.ES                0.3665224 -1.4797466 -0.7917932 -1.3786692858

PV.t       PV.w
Plasma.cells_C12.ES      0.04715402 0.02357308
CD8_C15.ES               0.08263955 0.07796962
CD8_C0.ES                0.13912575 0.09167711
CD4_C6.ES                0.18855160 0.16213817
#@@@@@@@@@@@@@@@@@@@@@@@@@@. PLOTS
### 1. **Bar Plot of AUC values**
library(ggplot2)
xx=xx[!c(rownames(xx) %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]
# Create a data frame for the AUC plot
auc_data <- data.frame(Feature = rownames(xx), AUC = xx[, 1], PV=xx$PV.t)

# Plot the AUC values
plot2_BCMLC<-ggplot(auc_data, aes(x = reorder(Feature, AUC), y = AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(PV,4)), vjust = 0.5, hjust = 1, size = 4, color="white") + 
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "AUC in BCMLC patient samples", x = "Feature", y = "AUC") +
  theme_minimal()
print(plot2_BCMLC)

# Load necessary libraries
library(RColorBrewer)
library(ROCR)

# Filter out specific rows from the dataset
xx = xx[!rownames(xx) %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES"), ]
xx_significant = xx  # If you want to filter by significance, comment the line below
# xx_significant = xx[xx[, 5] < 0.05, ]


plot(NULL, xlim=c(0, 1), ylim=c(0, 1), xlab="False Positive Rate", ylab="True Positive Rate", main="Okayama ROC Curves with AUC")
dark_colors <- brewer.pal(nrow(xx_significant), "Dark2")  
legend_labels <- list()

# Loop through each significant cell type to plot its ROC curve
for(k in 1:nrow(xx_significant)) {
  score = c(dat1[,k], dat2[,k])  
  pred <- prediction(score, cat)  
  perf <- performance(pred, "tpr", "fpr") 
  lines(perf@x.values[[1]], perf@y.values[[1]], col=dark_colors[k], lwd=2)
  auc_value = round(myauc[k], 2)
  pv_value = round(xx_significant[k, "PV.t"], 5)  

  legend_labels[[k]] <- paste(rownames(xx_significant)[k], "AUC =", auc_value, "PV.t =", pv_value)
}

legend("bottomright", 
       legend = unlist(legend_labels),  
       col = dark_colors,               
       lwd = 2,                        
       cex = 0.35,                      
       box.lwd = 0)                     





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.2.3] relpase-free survival
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Sato_GSE41271_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Sato_GSE41271/Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}


info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = info$time.recurrence
e.surv = ifelse(info[, "recurrence"]==" Y", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[info$t.surv>0, ]
se = which(info$histology==" Adenocarcinoma")
info = info[se,]


##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data)		## 442

# info_out<-"/home/u251079/scRNA/Lung_scRNA/sato.txt"
# sato<-merge(data, info, by = "row.names")
# row.names(sato)<-sato$Row.names
# sato$Row.names<-NULL
# write.table(sato, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)



library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1,  lb1,  ub1)
xx = res[order(res[,2]),]
xx
name coxph.pval1 coxph.qval1       hr1
9       Plasma.cells_C12.ES  0.02283280   0.4221557 0.8738371
8            B.cells_C13.ES  0.05321276   0.4221557 0.8088541
19            B.cells_C2.ES  0.05756668   0.4221557 0.5250537
14                CD4_C7.ES  0.08255324   0.4540428 0.8376568
10            Mac.DC_C11.ES  0.19029826   0.8169339 1.1013396
5                CD8_C15.ES  0.24279807   0.8169339 0.8663824
4            Fibroblasts.ES  0.31049939   0.8169339 0.8701463
11               CD4_C10.ES  0.32627055   0.8169339 0.9122546
21                CD8_C0.ES  0.33902264   0.8169339 0.9247809
15                CD4_C6.ES  0.37599663   0.8169339 0.9123444
12       Plasma.cells_C9.ES  0.40846693   0.8169339 0.8980004
2               Mast_C17.ES  0.46989134   0.8244919 0.9471397
13                Mac_C8.ES  0.50635466   0.8244919 1.0806022
17           NK.cells_C4.ES  0.55500245   0.8244919 1.0584720

#@@@@@@@@@@@@@@@ PLOTS:
xx$logPval <- -log10(xx$coxph.pval1)
xx$Significance <- ifelse(xx$coxph.pval1 < 0.05, "Significant", "Not Significant")

sato1=xx[!c(xx$name %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]

Plot1_sato<-ggplot(sato1, aes(x = reorder(name, hr1), y = hr1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lb1, ymax = ub1), width = 0.2) +
  geom_text(aes(label = format(coxph.pval1, scientific = TRUE, digits = 2)), 
            position = position_dodge(width = 0.5), vjust = -0.5) +  # Position text slightly above points
  coord_flip() +
  labs(title = "Sato data", 
       x = NULL, 
       y = "Hazard Ratio (HR)") +
  theme_minimal()
print(Plot1_sato)



Plot_sato<-ggplot(sato1, aes(x = hr1, y = logPval, color = Significance)) +
  geom_point(size = 4) +  # Main scatter points
  geom_hline(yintercept = -log10(0.025), linetype = "dashed", color = "red", size = 1) +  # Red horizontal line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 1) +  # Vertical line at HR=1
  scale_color_manual(values = c("Significant" = "darkred", "Not Significant" = "navyblue")) +  # Color settings
  geom_text(data = subset(sato1, Significance == "Significant"), 
            aes(label = name), vjust = -0.5, size = 4, check_overlap = TRUE) +  # Label only significant points
  theme_classic() + scale_x_continuous(limits = c(0.5, max(sato1$hr1)+0.1)) + 
#  scale_y_continuous(limits = c(0, 5)) + 
  labs(title = "Sato_data",
       x = "Hazard Ratio",
       y = "-log10(P-value)") +
  theme(legend.position = "top")

Plot_sato

#------------------
library(ROCR)
se = which(info$recurrence == " Y")
dat1 = data[se, ]
se = which(info$recurrence == " N")
dat2 = data[se, ]

res = matrix(0, ncol(data), 5)
row.names(res) = colnames(data)
colnames(res) = c("Avg.R", "Avg.NR", "TS", "PV.t", "PV.w")
res = as.data.frame(res)
cat = c(rep("R", nrow(dat1)), rep("NR", nrow(dat2)))
myauc = rep(0, ncol(data))

res[,1] = apply(dat1, 2, mean)
res[,2] = apply(dat2, 2, mean)
for(k in 1:ncol(data))
{
  cat("\r", k)
  
  tmp = t.test(dat1[,k], dat2[,k])
  res[k,3] = tmp$statistic
  res[k,4] = tmp$p.value
  tmp = wilcox.test(dat1[,k], dat2[,k])
  res[k,5] = tmp$p.value
  score = c(dat1[,k], dat2[,k])
  pred <- prediction(score, cat)
  auc.perf <- performance(pred, measure = "auc")
  myauc[k] = as.numeric(auc.perf@y.values)
}
myres = cbind(myauc, res)
xx = myres[order(myres[,6]), ]
xx[1:20, ]

myauc       PV.w
B.cells_C13.ES           0.4043534 0.03008281
DC_C18.ES                0.5619718 0.16006980
Neutrophil_C5.ES         0.4413572 0.18375153
B.cells_C2.ES            0.4464789 0.22510303
CD8_C0.ES                0.4489117 0.24691934
GamDelT_C14.ES           0.5454545 0.30300062
CD4_C7.ES                0.4546735 0.30436659
CD4_C10.ES               0.4588988 0.35173406
CD8_C1.ES                0.4760563 0.58789677
Epithelial.ES            0.4782330 0.62237163
Mac.DC_C11.ES            0.5189501 0.66823642
Plasma.cells_C9.ES       0.4834827 0.70887749
Mast_C17.ES              0.5158771 0.71971756
Fibroblasts.ES           0.4854033 0.74156565
Mac_C8.ES                0.4862996 0.75698546
CD8_C15.ES               0.5134443 0.76140925
NK.cells_C4.ES           0.4870679 0.77028014
Prolifirating.Mac_C16.ES 0.4939821 0.89252617
Endothelial.cells.ES     0.5040973 0.92707540
CD4_C6.ES                0.4959027 0.92707540

library(ggplot2)
xx=xx[!c(rownames(xx) %in% c("Endothelial.cells.ES","Epithelial.ES", "Fibroblasts.ES")),]
# Create a data frame for the AUC plot
auc_data <- data.frame(Feature = rownames(xx), AUC = xx[, 1], PV=xx$PV.t)

# Plot the AUC values
plot2_sato<-ggplot(auc_data, aes(x = reorder(Feature, AUC), y = AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(PV,4)), vjust = 0.5, hjust = 1, size = 4, color="white") + 
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "AUC in Sato patient samples", x = "Feature", y = "AUC") +
  theme_minimal()
print(plot2_sato)


# ROC_AUC curves 
xx_significant = xx[xx[, 5] < 0.05, ]

library(RColorBrewer)

plot(NULL, xlim=c(0, 1), ylim=c(0, 1), xlab="False Positive Rate", ylab="True Positive Rate", main="Sato ROC Curves with AUC")
dark_colors <- brewer.pal(nrow(xx_significant), "Dark2")  # "Dark2" provides a good set of dark colors
legend_labels <- list()

for(k in 1:nrow(xx_significant)) {
  score = c(dat1[,k], dat2[,k])
  pred <- prediction(score, cat)
  perf <- performance(pred, "tpr", "fpr")
  lines(perf@x.values[[1]], perf@y.values[[1]], col=dark_colors[k], lwd=2)
  auc_value = round(myauc[k], 2)
  pv_value = round(xx_significant[k, "PV.t"], 5)  # rounding PV.t value
  
  legend_labels[[k]] <- paste(rownames(xx_significant)[k], "AUC =", auc_value, "PV.t =", pv_value)
}
legend("bottomright", legend = unlist(legend_labels), col=dark_colors, lwd=2, cex=1, box.lwd=0)






# 
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [3] CelRat analysis
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [3.1.1] Okayama_GSE31210
# rm(list=ls())
# myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/PRECOG/processed/GSE31210.HGU133Plus2_EntrezCDF.MAS5.pcl"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# # myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Okayama_GSE31210_AmoscRNAMarkers_BaseRat.txt"
# 
# data = read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)
# reg = read.table(myinf2, sep="\t", header=T, row.names=1)
# 
# source("~/WorSpa/system/myRprogram/BaseRat3_scRNAmarker_weighted.R")
# 
# xx = BaseRat3_pwPath(data, reg, perm=1000, myoutf, median.norm=F)
# 
# 
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [3.2.1] relpase-free survival
# rm(list=ls())
# myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Okayama_GSE31210_AmoscRNAMarkers_BaseRat.txt"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Clinical/GSE31210_Clinical_info.txt"
# 
# data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
# mystd = apply(abs(data), 2, sd)
# for(k in 1:ncol(data))
# {
#   data[,k] = data[,k]/mystd[k]
# }
# 
# info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# t.surv = as.numeric(info$days.before.relapse.censor)
# e.surv = ifelse(info$relapse ==" relapsed", 1, 0)
# info = cbind(t.surv, e.surv, info)
# tag = !is.na(info$t.surv)
# info = info[tag==1,]
# 
# library(survival)
# 
# comxx = intersect(row.names(data), row.names(info))
# data = data[comxx,]
# info = info[comxx,]
# 
# survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
# hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
# for(k in 1:ncol(data))
# {
#   cat("\r", k)
#   mytf = as.numeric(data[,k])
#   xx = cbind(mytf, info)
#   xx = xx[xx[, "t.surv"]>0,]
#   mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
#   mycox = summary(mycox)$table
#   survreg.pval1[k] = mycox["mytf", "p"]
#   mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
#   mycox = summary(mycox)
#   coxph.pval1[k] = mycox$coefficients[5]
#   tmp = mycox$conf.int
#   hr1[k] = tmp[1]
#   lb1[k] = tmp[3]
#   ub1[k] = tmp[4]
# }
# survreg.qval1 = p.adjust(survreg.pval1, "BH")
# coxph.qval1 = p.adjust(coxph.pval1, "BH")
# 
# name = colnames(data)
# res = data.frame(name, coxph.pval1, coxph.qval1, hr1 )
# xx = res[order(res[,2]),]
# xx[1:20,]
# 
# name  coxph.pval1  coxph.qval1
# 59                         Mono_C3__VS__Mast_C17 1.774288e-06 0.0004057550
# 57                   Neutrophil_C5__VS__Mast_C17 4.863609e-06 0.0004057550
# 54                          Mac_C8__VS__Mast_C17 6.350173e-06 0.0004057550
# 44           Prolifirating.Mac_C16__VS__Mast_C17 6.415099e-06 0.0004057550
# 23                          Mast_C17__VS__DC_C18 8.841915e-06 0.0004474009
# 129               Mono_C3__VS__Endothelial.cells 4.825416e-05 0.0020347170
# 66  Endothelial.cells__VS__Prolifirating.Mac_C16 5.860768e-05 0.0021000010
# 9                               Plasma.cells_C12 6.640319e-05 0.0021000010
# 124                Mac_C8__VS__Endothelial.cells 8.307028e-05 0.0023351978
# 127         Neutrophil_C5__VS__Endothelial.cells 1.109469e-04 0.0028069558
# 235                B.cells_C2__VS__Neutrophil_C5 1.563101e-04 0.0035184329
# 27                 Endothelial.cells__VS__DC_C18 1.803990e-04 0.0035184329
# 135            Plasma.cells_C12__VS__GamDelT_C14 1.962428e-04 0.0035184329
# 233               NK.cells_C4__VS__Neutrophil_C5 2.017379e-04 0.0035184329
# 69   Plasma.cells_C12__VS__Prolifirating.Mac_C16 2.086027e-04 0.0035184329
# 48                     GamDelT_C14__VS__Mast_C17 2.244307e-04 0.0035488109
# 226                    Neutrophil_C5__VS__CD4_C6 3.178434e-04 0.0047302572
# 87             Plasma.cells_C12__VS__Fibroblasts 3.423417e-04 0.0048118024
# 79         B.cells_C2__VS__Prolifirating.Mac_C16 3.899297e-04 0.0051922221
# 192                   Neutrophil_C5__VS__CD4_C10 4.338911e-04 0.0054252558
# hr1
# 59  1.8686524
# 57  1.7549115
# 54  1.8718612
# 44  1.7960520
# 23  0.5327898
# 129 1.7837283
# 66  0.6272203
# 9   0.5783210
# 124 1.7912396
# 127 1.6426165
# 235 0.6171255
# 27  0.5526097
# 135 0.6229251
# 233 0.5935679
# 69  0.6055614
# 48  1.5416848
# 226 1.5927163
# 87  0.6334139
# 79  0.6158208
# 192 1.5526189






#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4] Apply the signatures for prognostic prediction
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#------------------------------------
# [4.1.1] Shedden
# rm(list=ls())
# myinf1 = "/mount/ictr1//chenglab/cc59/PubDat/Cancer/Lung/Shedden/Shedden_GSE68465_Lung_Symbol_expr.txt"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# # myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Shedden_AmoscRNAMarkers_Base.txt"
# 
# mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# 
# source("~/WorSpa/system/myRprogram/base5.R")
# data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
# 
# ##
# reg = mywt
# xx = base5(data, reg, perm=1000, myoutf, median.norm=F)
# 
# 
# #------------------------------------
# [4.1.2] Takeuchi_GSE11969
# 
# rm(list=ls())
# myinf1 = "/mount/ictr1//chenglab/cc59/PubDat/Cancer/Lung/Takeuchi_GSE11969/Takeuchi_GSE11969_Lung_Symbol_expr.txt"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# # myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Takeuchi_GSE11969_AmoscRNAMarkers_Base.txt"
# 
# mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# 
# source("~/WorSpa/system/myRprogram/base5.R")
# data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
# 
# ##
# reg = mywt
# xx = base5(data, reg, perm=1000, myoutf, median.norm=F)
# 
# EOF
# 
# 
# #------------------------------------
# [4.1.3] Tomida_GSE13213
# 
# rm(list=ls())
# myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/PRECOG/processed/GSE13213.GPL6480.matrix.txt"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Tomida_GSE13213_AmoscRNAMarkers_Base.txt"
# 
# mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# 
# source("~/WorSpa/system/myRprogram/base5.R")
# data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
# 
# ##
# reg = mywt
# xx = base5(data, reg, perm=1000, myoutf, median.norm=F)
# 
# 
# 
# #------------------------------------
# [4.1.4] Tang_GSE42127
# 
# rm(list=ls())
# myinf1 = "/mount/ictr1//chenglab/cc59/PubDat/Cancer/Lung/Tang_GSE42127/Tang_GSE42127_Lung_expr_Symbol.txt"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# #myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Tang_GSE42127_AmoscRNAMarkers_Base.txt"
# 
# mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# 
# source("~/WorSpa/system/myRprogram/base5.R")
# data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
# 
# ##
# reg = mywt
# xx = base5(data, reg, perm=1000, myoutf, median.norm=F)
# 
# 
# 
# #------------------------------------
# [4.1.5] TCGA LUAD
# rm(list=ls())
# myinf1 = "/mount/ictr1//chenglab/cc59/PubDat/Dataset/Firehose/RNAseqv2_Hiseq/processed/LUAD_RNAseqv2_Tumor_Symbol.rda"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# # myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/TCGA_LUAD_AmoscRNAMarkers_Base.txt"
# 
# mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# 
# source("~/WorSpa/system/myRprogram/base5.R")
# load(file= myinf1)
# data = mydata	## 262 samples
# tag = apply(data>0,1, sum)
# data = data[tag>10,]
# data = log10(data+1)
# 
# ##
# reg = mywt
# xx = base5(data, reg, perm=1000, myoutf, median.norm=F)
# 
# 
# 
# #------------------------------------
# [4.1.6] Gentles_GSE67639
# rm(list=ls())
# myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Gentles_GSE67639_meta/GSE67639_Gentles_Lung_expr.txt"
# myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/SingleCell/human/LungCancer_GSE131907/LungCancer_CellMarker_AmosscRNAseq_4base.txt"
# 
# # myoutf = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Gentles_GSE67639_AmoscRNAMarkers_Base.txt"
# 
# mywt = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
# 
# source("~/WorSpa/system/myRprogram/base5.R")
# data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1, stringsAsFactors=F)
# dim(data)		## [1] 11094  1106
# myavg = apply(data, 1, mean, na.rm=T)
# mystd = apply(data, 1, sd, na.rm=T)
# sum(is.na(data))/nrow(data)/ncol(data) 	## 0.101421
# xx = apply(is.na(data), 1, sum)
# range(xx)
# 
# for(k in 1:ncol(data))
# {
#   cat("\r", k)
#   tmp = data[,k]
#   tmp = (tmp-myavg)/mystd
#   tmp[is.na(tmp)] = 0
#   data[,k] = tmp
# }
# 
# ## don't do median normalization, don't do quantile 
# 
# ##
# reg = mywt
# xx = base5(data, reg, perm=1000, myoutf, quantile.norm=F, median.norm=F)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++
[4.2.1] Shedden_GSE68465
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Shedden_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1//chenglab/cc59/PubDat/Cancer/Lung/Shedden/Clinical_info.txt"


data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}


info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
info = info[info[, "vital_status"]!="", ]
t.surv = as.numeric(as.character(info[, "months_to_last_contact_or_death"]))
xx = as.character(info[, "vital_status"])
e.surv = ifelse(xx==" Dead", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[!is.na(info$t.surv), ]
info = info[!is.na(info$t.surv), ]
#info<-info[info$t.surv <=62,]
##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data)		## 442

# info_out<-"/home/u251079/scRNA/Lung_scRNA/Shedden.txt"
# Shedden<-merge(data, info, by = "row.names")
# row.names(Shedden)<-Shedden$Row.names
# Shedden$Row.names<-NULL
# write.table(Shedden, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)

library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1, lb1, ub1 )
xx = res[order(res[,2]),]
xx

name coxph.pval1 coxph.qval1       hr1
2               Mast_C17.ES  0.01273372   0.1027717 0.8616190
19            B.cells_C2.ES  0.01736271   0.1027717 0.8538938
3  Prolifirating.Mac_C16.ES  0.02634186   0.1027717 1.1756005
5                CD8_C15.ES  0.03047960   0.1027717 0.9308235
8            B.cells_C13.ES  0.03107731   0.1027717 0.8625643
14                CD4_C7.ES  0.03136701   0.1027717 0.8608757
15                CD4_C6.ES  0.03270008   0.1027717 0.8832716
21                CD8_C0.ES  0.03972626   0.1092472 0.9038068
7            GamDelT_C14.ES  0.05377005   0.1314379 1.1650685
12       Plasma.cells_C9.ES  0.09112794   0.2004815 0.8887349
9       Plasma.cells_C12.ES  0.10813470   0.2162694 0.8963888
6      Endothelial.cells.ES  0.13012670   0.2385656 0.8923897
4            Fibroblasts.ES  0.22665585   0.3521145 0.9233902
20                CD8_C1.ES  0.23699751   0.3521145 0.9393546
11               CD4_C10.ES  0.24007808   0.3521145 0.9286825
13                Mac_C8.ES  0.29374006   0.4038926 1.0791569
22            Epithelial.ES  0.45589955   0.5899877 0.9491326
17           NK.cells_C4.ES  0.53103968   0.6275113 0.9707465
10            Mac.DC_C11.ES  0.54194162   0.6275113 0.9536274
18               Mono_C3.ES  0.58021139   0.6382325 1.0411556
1                 DC_C18.ES  0.76499111   0.8014193 0.9784672
16         Neutrophil_C5.ES  0.92497928   0.9249793 1.0070136

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS
[1]. Volcano plots
library(ggplot2)
library(ggrepel)
xx <- xx[!(xx$name %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES")), ]
xx$significant <- ifelse(xx$coxph.pval1 <= 0.05, "Sig", "NS")
pval_threshold <- 0.05
x_min <- min(xx$hr1, na.rm = TRUE)
x_max <- max(xx$hr1, na.rm = TRUE)
y_min <- min(-log10(xx$coxph.pval1), na.rm = TRUE)
y_max <- max(-log10(xx$coxph.pval1), na.rm = TRUE)
xx_significant <- xx[xx$significant == "Sig", ]

Shedden_plot<- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +  # Points, colored based on q-value significance
  scale_color_manual(values = c("black", "red")) +  # Red for significant, black for non-significant
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at HR = 1
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 1) +  # Horizontal line for p-value threshold
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits based on the data range
  scale_y_continuous(limits = c(y_min, y_max)) +  # Set y-axis limits based on the data range
  geom_text_repel(data = xx_significant, aes(label = name),  # Add labels only for significant points
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Shedden data: Survival Analysis") +  # Title and axis labels
  theme_classic() +  # Clean minimal theme
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(Shedden_plot)

Shedden_plot <- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "darkred", linewidth = 1) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  geom_text_repel(aes(label = name),
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Shedden: Survival Analysis") +
  theme_classic() +
  theme(legend.position = "none")

# Print the plot
print(Shedden_plot)



library(ggplot2)
library(ggrepel)
xx <- xx[!(xx$name %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES")), ]

# Define p-value threshold
pval_threshold <- 0.05

# Define axis limits
x_min <- min(xx$hr1, na.rm = TRUE)
x_max <- max(xx$hr1, na.rm = TRUE)
y_min <- min(-log10(xx$coxph.pval1), na.rm = TRUE)
y_max <- max(-log10(xx$coxph.pval1), na.rm = TRUE)

# Create the Volcano plot
Shedden_plot <- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 3, alpha = 0.7) +  # Color significant points
  scale_color_manual(values = c("black", "red")) +  # Black for non-significant, red for significant
  geom_vline(xintercept = 1, linetype = "dashed", color = "blue", linewidth = 1) +  # Vertical line at HR = 1
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 1) +  # Horizontal line at p = 0.05
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits
  scale_y_continuous(limits = c(y_min, y_max)) +  # Set y-axis limits
  geom_text_repel(aes(label = name), size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +  # Label points
  labs(
    x = "Hazard Ratio (HR)",
    y = "-Log10(P-value)",
    title = "Shedden: Survival Analysis (Volcano Plot)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Remove legend
  )

# Print the plot
print(Shedden_plot)


[2]. KM plot

library(survival)
library(survminer)

# Example for a significant gene (e.g., "Mast_C17.ES")
gene <- "Plasma.cells_C12.ES"
gene_expression <- data[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")

# Fit survival data
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = info)

# Plot Kaplan-Meier curve
km_plot <- ggsurvplot(fit, data = info, pval = TRUE, risk.table = F,
                      title = paste("KM_plot in Shedden_data","(",gene,")"),
                      legend.labs = c("High Expression", "Low Expression"),
                      xlab="Time in months",
                      palette = c("red", "blue"))

print(km_plot)


[3].Hazard Ratio Forest Plot
# Subset significant results
significant_res <- res[res$coxph.pval1 < 0.05, ]

# Create a forest plot
forest_plot <- ggplot(significant_res, aes(x = name, y = hr1, ymin = lb1, ymax = ub1)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Gene", y = "Hazard Ratio (95% CI)", title = "Forest Plot of Hazard Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@
info_out<-"/home/u251079/scRNA/Lung_scRNA/Shedden.txt"
Shedden <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)



#+++++++++++++++++++++++++++++++++
[4.2.2] Takeuchi_GSE11969
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Takeuchi_GSE11969_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1//chenglab/cc59/PubDat/Cancer/Lung/Takeuchi_GSE11969/GSE11969_163_Patient_info_for_GEOarchive.txt"


data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

info = read.table(myinf2, sep="\t", header=T, quote="")
info = info[!is.na(info[, "Annotation"]), ]
info = info[info[, "Status"]!="", ]
row.names(info) = paste("GSM", 302996:303144, sep="")
t.surv = as.numeric(info[, "Survival"])
xx = as.character(info[, "Status"])
e.surv = ifelse(xx=="Dead", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[info$Histology=="AD", ]

##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data)			## 90

info_out<-"/home/u251079/scRNA/Lung_scRNA/Takeuchi.txt"
Takeuchi<-merge(data, info, by = "row.names")
row.names(Takeuchi)<-Takeuchi$Row.names
Takeuchi$Row.names<-NULL
write.table(Takeuchi, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)



library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1 )
xx = res[order(res[,2]),]
xx

name  coxph.pval1 coxph.qval1       hr1
7            GamDelT_C14.ES 0.0003447879 0.007585335 2.0238466
3  Prolifirating.Mac_C16.ES 0.0054732673 0.060205940 1.3715389
2               Mast_C17.ES 0.0277249927 0.203316613 0.6159931
17           NK.cells_C4.ES 0.0504629503 0.277546227 0.7181858
5                CD8_C15.ES 0.1031536456 0.406432793 0.7273209
15                CD4_C6.ES 0.1362254962 0.406432793 0.7554207
9       Plasma.cells_C12.ES 0.1608719815 0.406432793 0.8372322
1                 DC_C18.ES 0.1677935036 0.406432793 1.2490010
8            B.cells_C13.ES 0.1677994017 0.406432793 0.8118796
21                CD8_C0.ES 0.1847421787 0.406432793 0.7811005
19            B.cells_C2.ES 0.2164924539 0.432984908 0.8119816
14                CD4_C7.ES 0.2561087255 0.469532663 0.8213846
11               CD4_C10.ES 0.3721763701 0.629836934 0.8574025
10            Mac.DC_C11.ES 0.4173921195 0.655901902 0.8477215
20                CD8_C1.ES 0.4771813186 0.699865934 0.8861757
13                Mac_C8.ES 0.6009049069 0.773420716 1.0898990
6      Endothelial.cells.ES 0.6303110538 0.773420716 0.8931173
12       Plasma.cells_C9.ES 0.6327987679 0.773420716 0.9276555
22            Epithelial.ES 0.7331140057 0.848868849 1.0532034
16         Neutrophil_C5.ES 0.8521891449 0.937408059 0.9641370
18               Mono_C3.ES 0.9089925462 0.952277906 1.0209161
4            Fibroblasts.ES 0.9926170845 0.992617085 1.0016882
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS
[1]. Volcano plots
library(ggplot2)
library(ggrepel)
xx <- xx[!(xx$name %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES")), ]
xx$significant <- ifelse(xx$coxph.pval1 <= 0.05, "Sig", "NS")
pval_threshold <- 0.05
x_min <- min(xx$hr1, na.rm = TRUE)
x_max <- max(xx$hr1, na.rm = TRUE)
y_min <- min(-log10(xx$coxph.pval1), na.rm = TRUE)
y_max <- max(-log10(xx$coxph.pval1), na.rm = TRUE)
xx_significant <- xx[xx$significant == "Sig", ]

Takeuchi_plot<- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +  # Points, colored based on q-value significance
  scale_color_manual(values = c("black", "red")) +  # Red for significant, black for non-significant
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at HR = 1
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 1) +  # Horizontal line for p-value threshold
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits based on the data range
  scale_y_continuous(limits = c(y_min, y_max)) +  # Set y-axis limits based on the data range
  geom_text_repel(data = xx_significant, aes(label = name),  # Add labels only for significant points
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Takeuchi data: Survival Analysis") +  # Title and axis labels
  theme_classic() +  # Clean minimal theme
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(Takeuchi_plot)

Takeuchi_plot <- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "darkred", linewidth = 1) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  geom_text_repel(aes(label = name),
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Takeuchi: Survival Analysis") +
  theme_classic() +
  theme(legend.position = "none")


print(Takeuchi_plot)






[2]. KM plot

library(survival)
library(survminer)

# Example for a significant gene (e.g., "Mast_C17.ES")
gene <- "Plasma.cells_C12.ES"
gene_expression <- data[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = info)

km_plot <- ggsurvplot(fit, data = info, pval = TRUE, risk.table = TRUE, 
                      title = paste("KM Takeuchi_data", "(",gene,")"),
                      legend.labs = c("High Expression", "Low Expression"),
                      palette = c("red", "blue"))

print(km_plot)


[3].Hazard Ratio Forest Plot
# Subset significant results
significant_res <- res[res$coxph.pval1 < 0.05, ]

# Create a forest plot
forest_plot <- ggplot(significant_res, aes(x = name, y = hr1, ymin = lb1, ymax = ub1)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Gene", y = "Hazard Ratio (95% CI)", title = "Forest Plot of Hazard Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)


[4]. Heatmap 
library(pheatmap)

# Subset significant genes
significant_genes <- res$name[res$coxph.pval1 < 0.05]
significant_data <- data[, significant_genes]

# Annotate with survival status
annotation <- data.frame(Survival_Status = info$e.surv)
rownames(annotation) <- rownames(info)

# Plot heatmap
pheatmap(significant_data, annotation_row = annotation, 
         show_rownames = FALSE, show_colnames = TRUE, 
         scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Heatmap of Significant Genes")


[5]. Correlation plot
library(corrplot)

# Calculate correlations
cor_matrix <- cor(data[, significant_genes], use = "complete.obs")

# Plot correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Heatmap of Significant Genes")

[6].Survival Regression Coefficients Plot
# Combine survreg results
survreg_results <- data.frame(name = colnames(data), pval = survreg.pval1, qval = survreg.qval1)

# Plot coefficients
coef_plot <- ggplot(survreg_results, aes(x = name, y = -log10(pval), label = name)) +
  geom_point(aes(color = qval < 0.05), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  geom_text(data = subset(survreg_results, qval < 0.05), aes(label = name), vjust = 1.5, hjust = 1.5, size = 3) +
  labs(x = "Gene", y = "-Log10(P-value)", title = "Survival Regression Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(coef_plot)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@
info_out<-"/home/u251079/scRNA/Lung_scRNA/Takeuchi.txt"
Takeuchi <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)


#+++++++++++++++++++++++++++++++++
[4.2.3] Tomida_GSE13213
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Okayama_GSE31210_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Clinical/GSE31210_Clinical_info.txt"

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$days.before.death.censor)
e.surv = ifelse(info$death==" dead", 1, 0)
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]
dim(info)

data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data)		## 117

info_out<-"/home/u251079/scRNA/Lung_scRNA/Tomida.txt"
Tomida<-merge(data, info, by = "row.names")
row.names(Tomida)<-Tomida$Row.names
Tomida$Row.names<-NULL
write.table(Tomida, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)


library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survreg(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1 )
xx = res[order(res[,2]),]
xx

name coxph.pval1 coxph.qval1       hr1
9       Plasma.cells_C12.ES  0.01545082   0.2233063 0.6436630
6      Endothelial.cells.ES  0.02030057   0.2233063 0.6262576
19            B.cells_C2.ES  0.03618270   0.2653398 0.7509976
22            Epithelial.ES  0.07978468   0.3635810 0.7505370
3  Prolifirating.Mac_C16.ES  0.08263206   0.3635810 1.3120274
10            Mac.DC_C11.ES  0.11626473   0.3951834 1.3584795
2               Mast_C17.ES  0.12574018   0.3951834 0.8187326
18               Mono_C3.ES  0.18808053   0.4095803 1.2368545
15                CD4_C6.ES  0.19975430   0.4095803 0.8908532
8            B.cells_C13.ES  0.20062612   0.4095803 0.7864445
11               CD4_C10.ES  0.20479015   0.4095803 0.8800721
16         Neutrophil_C5.ES  0.23049933   0.4225821 1.2284572
20                CD8_C1.ES  0.30783663   0.5209543 0.8701757
7            GamDelT_C14.ES  0.40998613   0.6288957 1.1408809
17           NK.cells_C4.ES  0.42879253   0.6288957 0.8829786
13                Mac_C8.ES  0.53219433   0.6944507 1.1113986
5                CD8_C15.ES  0.53662097   0.6944507 0.9131943
14                CD4_C7.ES  0.69006429   0.8434119 0.9757352
1                 DC_C18.ES  0.76018858   0.8802184 0.9467598
4            Fibroblasts.ES  0.92843236   0.9847017 1.0096482
12       Plasma.cells_C9.ES  0.97710948   0.9847017 1.0050629
21                CD8_C0.ES  0.98470166   0.9847017 1.0022366


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS
[1]. Volcano plots
xx <- xx[!(xx$name %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES")), ]
xx$significant <- ifelse(xx$coxph.pval1 <= 0.05, "Sig", "NS")
library(ggplot2)
library(ggrepel)
pval_threshold <- 0.05
x_min <- min(xx$hr1, na.rm = TRUE)
x_max <- max(xx$hr1, na.rm = TRUE)
y_min <- min(-log10(xx$coxph.pval1), na.rm = TRUE)
y_max <- max(-log10(xx$coxph.pval1), na.rm = TRUE)
xx_significant <- xx[xx$significant == "Sig", ]

Tomida_plot<- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +  # Points, colored based on q-value significance
  scale_color_manual(values = c("black", "red")) +  # Red for significant, black for non-significant
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at HR = 1
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 1) +  # Horizontal line for p-value threshold
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits based on the data range
  scale_y_continuous(limits = c(y_min, y_max)) +  # Set y-axis limits based on the data range
  geom_text_repel(data = xx_significant, aes(label = name),  # Add labels only for significant points
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Tomida: Survival Analysis") +  # Title and axis labels
  theme_classic() +  # Clean minimal theme
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(Tomida_plot)





[2]. KM plot

library(survival)
library(survminer)

# Example for a significant gene (e.g., "Mast_C17.ES")
gene <- "Plasma.cells_C12.ES"
gene_expression <- data[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")

# Fit survival data
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = info)

km_plot <- ggsurvplot(fit, data = info, pval = TRUE, risk.table = T, 
                      title = paste("Tomida data KM Curve for", "(",gene,")"),
                      legend.labs = c("High Expression", "Low Expression"),
                      xlab = "Time in days",
                      palette = c("red", "blue"))


print(km_plot)


[3].Hazard Ratio Forest Plot
# Subset significant results
significant_res <- res[res$coxph.pval1 < 0.05, ]

# Create a forest plot
forest_plot <- ggplot(significant_res, aes(x = name, y = hr1, ymin = lb1, ymax = ub1)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Gene", y = "Hazard Ratio (95% CI)", title = "Forest Plot of Hazard Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)


[4]. Heatmap 
library(pheatmap)

# Subset significant genes
significant_genes <- res$name[res$coxph.pval1 < 0.05]
significant_data <- data[, significant_genes]

# Annotate with survival status
annotation <- data.frame(Survival_Status = info$e.surv)
rownames(annotation) <- rownames(info)

# Plot heatmap
pheatmap(significant_data, annotation_row = annotation, 
         show_rownames = FALSE, show_colnames = TRUE, 
         scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Heatmap of Significant Genes")


[5]. Correlation plot
library(corrplot)

# Calculate correlations
cor_matrix <- cor(data[, significant_genes], use = "complete.obs")

# Plot correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Heatmap of Significant Genes")

[6].Survival Regression Coefficients Plot
# Combine survreg results
survreg_results <- data.frame(name = colnames(data), pval = survreg.pval1, qval = survreg.qval1)

# Plot coefficients
coef_plot <- ggplot(survreg_results, aes(x = name, y = -log10(pval), label = name)) +
  geom_point(aes(color = qval < 0.05), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  geom_text(data = subset(survreg_results, qval < 0.05), aes(label = name), vjust = 1.5, hjust = 1.5, size = 3) +
  labs(x = "Gene", y = "-Log10(P-value)", title = "Survival Regression Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(coef_plot)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@
info_out<-"/home/u251079/scRNA/Lung_scRNA/Tomida.txt"
Tomida <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.2.4] Gentles_GSE67639_meta
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Gentles_GSE67639_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Gentles_GSE67639_meta/Clinical_info.txt"    


data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

info = read.table(myinf2, sep="\t", header=T, row.names=1)
t.surv = as.numeric(info[, "OS_Time"])
e.surv = info[, "OS_Status"]
info = cbind(t.surv, e.surv, info)
table(info$Histology_MD)
##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data)		## 58


info_out<-"/home/u251079/scRNA/Lung_scRNA/Gentles.txt"
Gentles<-merge(data, info, by = "row.names")
row.names(Gentles)<-Gentles$Row.names
Gentles$Row.names<-NULL
write.table(Gentles, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)

library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  mytag = ifelse(mytf>median(mytf), 1, 0)
  xx = cbind(mytf, mytag, info)
  xx = xx[xx[, "t.surv"]>0,]
  mycox = survdiff(Surv(t.surv, e.surv)~mytag, xx) 
  survreg.pval1[k] = pchisq(mycox$chisq, 1, lower.tail=F)
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1 )
xx = res[order(res[,2]),]
xx[1:20,]

name  coxph.pval1  coxph.qval1       hr1
8            B.cells_C13.ES 2.215380e-07 4.873835e-06 0.8768122
19            B.cells_C2.ES 5.251799e-06 5.776979e-05 0.9001946
11               CD4_C10.ES 5.356253e-05 3.927919e-04 0.9080497
2               Mast_C17.ES 7.751190e-05 4.263154e-04 0.9178333
5                CD8_C15.ES 2.541759e-04 1.118374e-03 0.8852225
15                CD4_C6.ES 3.556084e-04 1.303897e-03 0.8883631
14                CD4_C7.ES 4.270926e-04 1.342291e-03 0.9247380
17           NK.cells_C4.ES 1.823741e-03 5.015288e-03 0.9320803
21                CD8_C0.ES 6.275488e-03 1.534008e-02 0.9210992
12       Plasma.cells_C9.ES 1.044424e-02 2.297734e-02 0.9292731
20                CD8_C1.ES 1.611454e-02 3.222907e-02 0.9314662
6      Endothelial.cells.ES 8.376330e-02 1.428166e-01 0.9543542
22            Epithelial.ES 8.439162e-02 1.428166e-01 0.9564694
1                 DC_C18.ES 1.855575e-01 2.915904e-01 0.9680648
3  Prolifirating.Mac_C16.ES 2.049394e-01 3.005777e-01 1.0334801
7            GamDelT_C14.ES 3.764040e-01 5.175555e-01 1.0221576
4            Fibroblasts.ES 6.487781e-01 8.395951e-01 0.9857787
16         Neutrophil_C5.ES 7.295585e-01 8.896078e-01 1.0125578
13                Mac_C8.ES 7.682977e-01 8.896078e-01 0.9921776
9       Plasma.cells_C12.ES 8.146948e-01 8.961642e-01 0.9897927
18               Mono_C3.ES 8.922538e-01 9.347421e-01 0.9963356
10            Mac.DC_C11.ES 9.493907e-01 9.493907e-01 1.0016470
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS
[1]. Volcano plots
xx <- xx[!(xx$name %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES")), ]
xx$significant <- ifelse(xx$coxph.pval1 <= 0.05, "Sig", "NS")
library(ggplot2)
library(ggrepel)
pval_threshold <- 0.05
x_min <- min(xx$hr1, na.rm = TRUE)
x_max <- max(xx$hr1, na.rm = TRUE)
y_min <- min(-log10(xx$coxph.pval1), na.rm = TRUE)
y_max <- max(-log10(xx$coxph.pval1), na.rm = TRUE)
xx_significant <- xx[xx$significant == "Sig", ]

Gentles_plot<- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +  # Points, colored based on q-value significance
  scale_color_manual(values = c("black", "red")) +  # Red for significant, black for non-significant
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at HR = 1
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 1) +  # Horizontal line for p-value threshold
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits based on the data range
  scale_y_continuous(limits = c(y_min, y_max)) +  # Set y-axis limits based on the data range
  geom_text_repel(data = xx_significant, aes(label = name),  # Add labels only for significant points
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Gentles data: Survival Analysis") +  # Title and axis labels
  theme_classic() +  # Clean minimal theme
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(Gentles_plot)

Gentles_plot <- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "darkred", linewidth = 1) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  geom_text_repel(aes(label = name),
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Gentles: Survival Analysis") +
  theme_classic() +
  theme(legend.position = "none")

# Print the plot
print(Gentles_plot)


[2]. KM plot

library(survival)
library(survminer)

# Example for a significant gene (e.g., "Mast_C17.ES")
gene <- "Plasma.cells_C12.ES"
gene_expression <- data[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")

# Fit survival data
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = info)

# Plot Kaplan-Meier curve
km_plot <- ggsurvplot(fit, data = info, pval = TRUE, risk.table = TRUE, 
                      title = paste("Gentles data KM Curve for", "(",gene,")"),
                      legend.labs = c("High Expression", "Low Expression"),
                      xlab="Time in months (?)",
                      palette = c("red", "blue"))

print(km_plot)


[3].Hazard Ratio Forest Plot
# Subset significant results
significant_res <- res[res$coxph.pval1 < 0.05, ]

# Create a forest plot
forest_plot <- ggplot(significant_res, aes(x = name, y = hr1, ymin = lb1, ymax = ub1)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Gene", y = "Hazard Ratio (95% CI)", title = "Forest Plot of Hazard Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)


[4]. Heatmap 
library(pheatmap)

# Subset significant genes
significant_genes <- res$name[res$coxph.pval1 < 0.05]
significant_data <- data[, significant_genes]

# Annotate with survival status
annotation <- data.frame(Survival_Status = info$e.surv)
rownames(annotation) <- rownames(info)

# Plot heatmap
pheatmap(significant_data, annotation_row = annotation, 
         show_rownames = FALSE, show_colnames = TRUE, 
         scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Heatmap of Significant Genes")


[5]. Correlation plot
library(corrplot)

# Calculate correlations
cor_matrix <- cor(data[, significant_genes], use = "complete.obs")

# Plot correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Heatmap of Significant Genes")

[6].Survival Regression Coefficients Plot
# Combine survreg results
survreg_results <- data.frame(name = colnames(data), pval = survreg.pval1, qval = survreg.qval1)

# Plot coefficients
coef_plot <- ggplot(survreg_results, aes(x = name, y = -log10(pval), label = name)) +
  geom_point(aes(color = qval < 0.05), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  geom_text(data = subset(survreg_results, qval < 0.05), aes(label = name), vjust = 1.5, hjust = 1.5, size = 3) +
  labs(x = "Gene", y = "-Log10(P-value)", title = "Survival Regression Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(coef_plot)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@
info_out<-"/home/u251079/scRNA/Lung_scRNA/Gentles.txt"
Gentles<- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.2.5] Tang
rm(list=ls())
myinf1 = "/home/u251079/scRNA/Lung_scRNA/Lung_output/Tang_GSE42127_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Tang_GSE42127/Clinical_info.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = info$OS
e.surv = ifelse(info$status=="D", 1, 0)
info = cbind(t.surv, e.surv, info)
info$age = as.numeric(info$age)

##
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
dim(data) 		## 63

info_out<-"/home/u251079/scRNA/Lung_scRNA/Tang.txt"
Tang<-merge(data, info, by = "row.names")
row.names(Tang)<-Tang$Row.names
Tang$Row.names<-NULL
write.table(Tang, file = info_out, sep = "\t", quote = FALSE, row.names = TRUE)


library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  xx = xx[xx[, "t.surv"]>0,]
  se = which(xx$adjuvant_chemo=="FALSE")
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx[se,]) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  se = which(xx$adjuvant_chemo=="TRUE")
  mycox = coxph(Surv(t.surv, e.surv)~mytf, xx[se,]) 
  mycox = summary(mycox)
  coxph.pval2[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr2[k] = tmp[1]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = colnames(data)
res = data.frame(name, coxph.pval1, coxph.qval1, hr1)
xx = res[order(res[,2]),]
xx

name coxph.pval1 coxph.qval1       hr1
8            B.cells_C13.ES  0.02807035   0.2637639 0.7258822
19            B.cells_C2.ES  0.05167501   0.2637639 0.7644313
10            Mac.DC_C11.ES  0.05864134   0.2637639 0.9062961
14                CD4_C7.ES  0.06521799   0.2637639 0.7960578
5                CD8_C15.ES  0.07478691   0.2637639 0.7826381
11               CD4_C10.ES  0.08900778   0.2637639 0.7929163
22            Epithelial.ES  0.09611803   0.2637639 0.7903376
21                CD8_C0.ES  0.10327630   0.2637639 0.8163440
15                CD4_C6.ES  0.10790340   0.2637639 0.8055265
20                CD8_C1.ES  0.13713369   0.3016941 0.8290407
17           NK.cells_C4.ES  0.15964656   0.3192931 0.8296899
2               Mast_C17.ES  0.24906273   0.4401896 0.8947263
4            Fibroblasts.ES  0.26011202   0.4401896 1.1254541
1                 DC_C18.ES  0.28923520   0.4545125 0.8796175
13                Mac_C8.ES  0.65316028   0.9517445 0.9422714
6      Endothelial.cells.ES  0.74390209   0.9517445 0.9662604
16         Neutrophil_C5.ES  0.80841299   0.9517445 0.9718486
12       Plasma.cells_C9.ES  0.84326719   0.9517445 0.9828238
7            GamDelT_C14.ES  0.84999978   0.9517445 1.0237732
18               Mono_C3.ES  0.88538044   0.9517445 0.9816537
9       Plasma.cells_C12.ES  0.95087677   0.9517445 0.9940783
3  Prolifirating.Mac_C16.ES  0.95174448   0.9517445 0.9919072


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS
[1]. Volcano plots
xx <- xx[!(xx$name %in% c("Endothelial.cells.ES", "Epithelial.ES", "Fibroblasts.ES")), ]
xx$significant <- ifelse(xx$coxph.pval1 <= 0.05, "Sig", "NS")
library(ggplot2)
library(ggrepel)
pval_threshold <- 0.05
x_min <- min(xx$hr1, na.rm = TRUE)
x_max <- max(xx$hr1, na.rm = TRUE)
y_min <- min(-log10(xx$coxph.pval1), na.rm = TRUE)
y_max <- max(-log10(xx$coxph.pval1), na.rm = TRUE)
xx_significant <- xx[xx$significant == "Sig", ]

Tang_plot<- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +  # Points, colored based on q-value significance
  scale_color_manual(values = c("black", "red")) +  # Red for significant, black for non-significant
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at HR = 1
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red", linewidth = 1) +  # Horizontal line for p-value threshold
  scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits based on the data range
  scale_y_continuous(limits = c(y_min, y_max)) +  # Set y-axis limits based on the data range
  geom_text_repel(data = xx_significant, aes(label = name),  # Add labels only for significant points
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Tang data: Survival Analysis") +  # Title and axis labels
  theme_classic() +  # Clean minimal theme
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(Tang_plot)
Tang_plot <- ggplot(xx, aes(x = hr1, y = -log10(coxph.pval1), label = name)) +
  geom_point(aes(color = coxph.qval1 < 0.05), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "darkred", linewidth = 1) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  geom_text_repel(aes(label = name),
                  size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20) +
  labs(x = "Hazard Ratio (HR)", y = "-Log10(P-value)", title = "Tang: Survival Analysis") +
  theme_minimal() +
  theme(legend.position = "none")

print(Tang_plot)

[2]. KM plot

library(survival)
library(survminer)

# Example for a significant gene (e.g., "Mast_C17.ES")
gene <- "Plasma.cells_C12.ES"
gene_expression <- data[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")

# Fit survival data
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = info)

# Plot Kaplan-Meier curve
km_plot <- ggsurvplot(fit, data = info, pval = TRUE, risk.table = TRUE, 
                      title = paste("Tang data KM Curve for","(",gene,")"),
                      legend.labs = c("High Expression", "Low Expression"),
                      xlab="Time in months (?)",
                      palette = c("red", "blue"))

print(km_plot)


[3].Hazard Ratio Forest Plot
# Subset significant results
significant_res <- res[res$coxph.pval1 < 0.05, ]

# Create a forest plot
forest_plot <- ggplot(significant_res, aes(x = name, y = hr1, ymin = lb1, ymax = ub1)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Gene", y = "Hazard Ratio (95% CI)", title = "Forest Plot of Hazard Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)


[4]. Heatmap 
library(pheatmap)

# Subset significant genes
significant_genes <- res$name[res$coxph.pval1 < 0.05]
significant_data <- data[, significant_genes]

# Annotate with survival status
annotation <- data.frame(Survival_Status = info$e.surv)
rownames(annotation) <- rownames(info)

# Plot heatmap
pheatmap(significant_data, annotation_row = annotation, 
         show_rownames = FALSE, show_colnames = TRUE, 
         scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Heatmap of Significant Genes")


[5]. Correlation plot
library(corrplot)

# Calculate correlations
cor_matrix <- cor(data[, significant_genes], use = "complete.obs")

# Plot correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Heatmap of Significant Genes")

[6].Survival Regression Coefficients Plot
# Combine survreg results
survreg_results <- data.frame(name = colnames(data), pval = survreg.pval1, qval = survreg.qval1)

# Plot coefficients
coef_plot <- ggplot(survreg_results, aes(x = name, y = -log10(pval), label = name)) +
  geom_point(aes(color = qval < 0.05), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  geom_text(data = subset(survreg_results, qval < 0.05), aes(label = name), vjust = 1.5, hjust = 1.5, size = 3) +
  labs(x = "Gene", y = "-Log10(P-value)", title = "Survival Regression Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(coef_plot)





