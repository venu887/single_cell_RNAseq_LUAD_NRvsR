# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [1] define marker gene profiles based on scRNA-seq data for infer cell ratio
# #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1.1] define marker gene signatures with weights
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(ggplot2)
library(NMF)
library(reshape2)
rm(list=ls())
myinf1 = "Immune_cells_data/cd45pos_nsclc_Immune_cells.rds"
myinf2 = "NSCLC_output/cd45negc_celltype_data.rds"
myoutf1 = "LungCancer_CellMarker_AmosscRNAseq_4base.txt"
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

# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [2] BASE analysis
# #@@@@@@@@@@@@@@@
# [2.1] BCM_LC_RNAseq # BCM bulk RNAseq data
# In-House data
myinf1 = "BCM_prot_gene_TPM_Symbol.csv"
myinf2 = "LungCancer_CellMarker_AmosscRNAseq_4base.txt"
myoutf = "BCM_LC_RNAseq_AmoscRNAMarkers_Base.txt"
data = read.csv(myinf1, header=T, row.names=1, check.names=F)
data = log2(data+1)
reg = read.table(myinf2, sep="\t", header=T, row.names=1)
data = base5(data, reg, perm=1000, myoutf, median.norm=F)

# #@@@@@@@@@@@@@@@
#[2.2] Okayama_GSE31210
rm(list=ls())
myinf1 = "GSE31210.HGU133Plus2_EntrezCDF.MAS5.pcl"
myinf2 = "LungCancer_CellMarker_AmosscRNAseq_4base.txt"
myoutf = "Okayama_GSE31210_AmoscRNAMarkers_Base.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)
reg = read.table(myinf2, sep="\t", header=T, row.names=1)
# Can freely download from online papers 
source("base5.R")
data = base5(data, reg, perm=1000, myoutf, median.norm=F)


 
#@@@@@@@@@@@@@@@ Relpase-free survival @@@@@@@@@@@@@@@
[3] relpase-free survival or disease free survival
Figure 4 
[4.1] Vocano plot to show univariable Cox results with RFS
rm(list=ls())
myinf1 = "BCM_LC_RNAseq.rda"
load(myinf1)
data = ORI.data
info = ORI.info
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

PV1 = PV2 = rep(0, ncol(data))
HR1 = HR2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
	cat("\r", k)
	mytf = as.numeric(data[,k])
	xx = cbind(mytf, info)
	xx = xx[xx[, "t.surv"]>0,]
	mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
	mycox = summary(mycox)
	PV1[k] = mycox$coefficients[5]
	tmp = mycox$conf.int
	HR1[k] = tmp[1]
	mycox = coxph(Surv(t.surv, e.surv)~mytf + Age + Sex + Smoking + Stage, xx) 
	mycox = summary(mycox)
	PV2[k] = mycox$coefficients["mytf", 5]
	tmp = mycox$conf.int
	HR2[k] = tmp["mytf", 1]
}
QV1 = p.adjust(PV1, method="BH")
QV2 = p.adjust(PV2, method="BH")

res = data.frame(PV1, HR1, PV2, HR2)
row.names(res) = colnames(data)
res = res[order(res$PV1),]
xx = row.names(res)
xx = gsub("\\.ES", "", xx)
xx = gsub("\\.cells", "", xx)
row.names(res) = xx
se = which(!row.names(res) %in% c("Endothelial", "Epithelial", "Fibroblasts"))
res = res[se,]

data = res
data$FDR = p.adjust(data$PV1, method="BH")
se = which(data$FDR<0.05)
range(log2(data$HR1))
 
#---------------------------------------
data$category= ifelse(data$PV1>0.05, "NS", ifelse(data$HR1>1, "Hazardous", "Protective"))
data$label =  rep("", nrow(data))
data$label[1:2] = row.names(data)[1:2]
lab.sub = which(c("Hazardous", "NS", "Protective") %in% data$category)
myy.intercept = -log10(0.05)
#---------------------------------------
mygg <- ggplot(data, aes(x = log2(HR1), y = -log10(PV1), col = category, label = label)) +
  geom_vline(xintercept = 0, col = "black", linetype = 'dashed') +
  geom_hline(yintercept = myy.intercept, col = "black", linetype = 'dashed') + 
  geom_point(size = 2, show.legend = FALSE) + 
  scale_color_manual(values = c("lightcoral", "grey", "skyblue")[lab.sub], 					
                     labels = c("Hazardous", "Not significant", "Protective")[lab.sub]) + 	
  guides(color="none") +	
  coord_cartesian(ylim = c(0, 2.5), xlim = c(-0.8, 0.8)) + 
  labs(color = '',
       x = "Log2(HR)", y = "-Log10(P-value)") + 							
  geom_text_repel(max.overlaps = Inf) + 									
  theme_classic()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.2] Vocano plot to show difference between R and NR 
rm(list=ls())
myinf1 = "BCM_LC_RNAseq.rda"
myFig = "Fig3_Volcano_BCM_RvsNR.pdf"
load(myinf1)
data = ORI.data
info = ORI.info
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
table(info$Recurrence)
se = which(info$Recurrence == 1)
dat1 = data[se, ]
se = which(info$Recurrence == 0)
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
res$AUC = myauc
res$Dif = res$Avg.R - res$Avg.NR
subset = c("Dif", "AUC", "PV.w")
res = res[, subset]
xx = row.names(res)
xx = gsub("\\.ES", "", xx)
xx = gsub("\\.cells", "", xx)
row.names(res) = xx
se = which(!row.names(res) %in% c("Endothelial", "Epithelial", "Fibroblasts"))
res = res[se,]
data = res
colnames(data) = c("Diff", "AUC", "PV1")
data = data[order(data$PV1),]

range(data$Diff)

data$category= ifelse(data$PV1>0.05, "NS", ifelse(data$Diff>0, "Hazardous", "Protective"))
data$label =  rep("", nrow(data))
data$label[1] = row.names(data)[1]
lab.sub = which(c("Hazardous", "NS", "Protective") %in% data$category)
myy.intercept = -log10(0.05)

mygg <- ggplot(data, aes(x = Diff, y = -log10(PV1), col = category, label = label)) +
  geom_vline(xintercept = 0, col = "black", linetype = 'dashed') +
  geom_hline(yintercept = myy.intercept, col = "black", linetype = 'dashed') + 
  geom_point(size = 2, show.legend = FALSE) + 
  scale_color_manual(values = c("lightcoral", "grey", "skyblue")[lab.sub], 					
                     labels = c("Hazardous", "Not significant", "Protective")[lab.sub]) + 	
  guides(color="none") +	
  coord_cartesian(ylim = c(0, 2.5), xlim = c(-1.5, 1.5)) + 
  labs(color = '', 
       x = "Difference (R-NR)", y = "-Log10(P-value)") + 
  geom_text_repel(max.overlaps = Inf) + 				
  theme_classic()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.3] KM plot for Plasma.cells_C12.ES
rm(list=ls())
myinf1 = "BCM_LC_RNAseq.rda"
myFig = "Fig3_KM_BCM_PlasmaC12.pdf"
load(myinf1)
data = ORI.data
info = ORI.info
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))
mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")
tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")
tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (months)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")
max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.4] ForestPlot  for Plasma.cells_C12.ES
rm(list=ls())
myinf1 = "BCM_LC_RNAseq.rda"

myFig = "Fig3_ForestPlot_BCM_PlasmaC12.pdf"
load(myinf1)
data = ORI.data
info = ORI.info

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
xx$mytag = ifelse(xx$mys>median(xx$mys), 1, 0)
se = which(xx$Stage=="")
xx$Stage[se] = NA
mycox = coxph(Surv(t.surv, e.surv)~mytag + Age + Sex + Smoking + Stage, xx) 
mycox = summary(mycox)
tmp = mycox$coefficients
PV1 = tmp[,5]
HR1 = tmp[,2]
tmp = mycox$conf.int
HR.lo = tmp[,3]
HR.hi = tmp[,4]
res = data.frame(PV1, HR1, HR.lo, HR.hi)
row.names(res) = row.names(tmp)

library(tidyverse)
library(forestplot)
ids2 = row.names(res)
ids2 = c("Plasma_C12 (High vs. Low)", "Age", "Sex (Male vs. Female)", "Smoking (Yes vs. No)", "Stage II (vs. I)", "Stage III (vs. I)")
mhr = res$HR1
mup = res$HR.hi
mlw = res$HR.lo
myp = res$PV1
myp = as.character(ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1)))
hr = paste(formatC(mhr, 2, format = 'f'), ' (', formatC(mlw, 2, format = 'f') ,' to ',formatC(mup, 2, format = 'f'),')' , sep="")

data = tibble(Feature = ids2, mean  = mhr, lower = mlw, upper = mup,  hr=hr, coxphP = myp)
header <- tibble(Feature = c("Feature"), hr='HR (95% CI)',coxphP = c('P-value'))
data = bind_rows(header, data)

limits = c(0.08, 8)  
xticks <- c(0.08,1,8)  
xtlab <- rep(c(TRUE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab

cairo_pdf(myFig, width = 5, height = 2.5) 
data %>%  
  forestplot(labeltext = c(Feature, hr,coxphP),
             graph.pos = 2,
             zero = 1,
             is.summary=c(TRUE, rep(FALSE, nrow(data)-1)), 
             xlog = FALSE,
             clip = limits,
             xlab = "HR",
             ci.vertices = TRUE,
             hrzl_lines= TRUE,
             boxsize = 0.2, 
             line.margin = 0.1,
             mar = unit(rep(0.1, times = 4), "mm"),
             xticks = xticks,
             graphwidth = "auto",
             graphhight = "auto",
             align = rep("l", 3),
             txt_gp=fpTxtGp(label = gpar(cex = 0.9),
                            title = gpar(cex = 0.9),
                            ticks = gpar(cex = 0.9),
                            xlab = gpar(cex = 0.9)),
             col = fpColors(box = "black",
                            line = "black")) |>
  fp_set_zebra_style("#EFEFEF")
dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.5] Boxplot  for Plasma.cells_C12.ES (Recurrence --> Smoking --> Stage ---> age)
rm(list=ls())
myinf1 = "BCM_LC_RNAseq.rda"
myFigA = "Fig3_BoxPlot_BCM_PlasmaC12_RvsNR.pdf"
myFigB = "Fig3_BoxPlot_BCM_PlasmaC12_Smoking.pdf"
myFigC = "Fig3_BoxPlot_BCM_PlasmaC12_Stage.pdf"
load(myinf1)
data = ORI.data
info = ORI.info

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
mydat = cbind(mys, info)

mycat = ifelse(mydat$Recurrence==1, "R", "NR")
data = cbind(mydat, mycat)
se = which(!is.na(data$mycat))
data = data[se,]

xx1 = data$mys[data$mycat=="R"]
xx2 = data$mys[data$mycat=="NR"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")

mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(mycat, mys), y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
mygg

pdf(file = myFigA, width = 2.2, height = 3) 
mygg
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$Smoking %in% c("Y", "N"))
xx = mydat[se,]
mycat = ifelse(xx$Smoking=="Y", "Smoker", "Non-smoker")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]

xx1 = data$mys[data$mycat=="Smoker"]
xx2 = data$mys[data$mycat=="Non-smoker"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(mycat, mys), y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
#mygg

pdf(file = myFigB, width = 2.2, height = 3) 
mygg
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$Stage %in% c("I", "II", "III"))
xx = mydat[se,]
mycat = ifelse(xx$Stage=="I", "Early", "Late")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]
data$mycat = as.factor(as.character(data$mycat))

xx1 = data$mys[data$mycat=="Early"]
xx2 = data$mys[data$mycat=="Late"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=mycat, y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
#mygg

pdf(file = myFigC, width = 2.2, height = 3) 
mygg
dev.off()



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
# Create Age Groups (example: median split)
mydat$Age_Group <- ifelse(mydat$Age >= median(mydat$Age, na.rm = TRUE), "Old", "Young")
mydat$Age_Group <- factor(mydat$Age_Group, levels = c("Young", "Old"))

# mycat = ifelse(mydat$Recurrence==1, "R", "NR")

data = mydat
se = which(!is.na(data$Age_Group))
data = data[se,]

xx1 = data$mys[data$Age_Group=="Young"]
xx2 = data$mys[data$Age_Group=="Old"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(Age_Group, mys), y=mys,  fill=Age_Group)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
mygg

pdf(file = myFigA_sub, width = 2.2, height = 3) 
mygg
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4] Okayama data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Figure 5 (Okayama result)
[5.1] Vocano plot to show univariable Cox results with RFS
rm(list=ls())
myinf1 = "Okayama_GSE31210_Data.rda"
myFig = "Fig4_Volcano_Okayama_Prognosis.pdf"
load(myinf1)
data = ORI.data
info = ORI.info
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
PV1 = PV2 = rep(0, ncol(data))
HR1 = HR2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
	cat("\r", k)
	mytf = as.numeric(data[,k])
	xx = cbind(mytf, info)
	xx = xx[xx[, "t.surv"]>0,]
	mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
	mycox = summary(mycox)
	PV1[k] = mycox$coefficients[5]
	tmp = mycox$conf.int
	HR1[k] = tmp[1]
	mycox = coxph(Surv(t.surv, e.surv)~mytf + age + gender + smoking.status + pstage.iorii, xx) 
	mycox = summary(mycox)
	PV2[k] = mycox$coefficients["mytf", 5]
	tmp = mycox$conf.int
	HR2[k] = tmp["mytf", 1]
}
QV1 = p.adjust(PV1, method="BH")
QV2 = p.adjust(PV2, method="BH")
res = data.frame(PV1, HR1, PV2, HR2)
row.names(res) = colnames(data)
res = res[order(res$PV1),]
xx = row.names(res)
xx = gsub("\\.ES", "", xx)
xx = gsub("\\.cells", "", xx)
row.names(res) = xx
se = which(!row.names(res) %in% c("Endothelial", "Epithelial", "Fibroblasts"))
res = res[se,]
data = res
data$FDR = p.adjust(data$PV1, method="BH")
se = which(data$FDR<0.05)
length(se)
range(log2(data$HR1))
#---------------------------------------
data$category= ifelse(data$PV1>0.05, "NS", ifelse(data$HR1>1, "Hazardous", "Protective"))
data$label =  rep("", nrow(data))
data$label[1:3] = row.names(data)[1:3]
lab.sub = which(c("Hazardous", "NS", "Protective") %in% data$category)
myy.intercept = -log10(0.05)

mygg <- ggplot(data, aes(x = log2(HR1), y = -log10(PV1), col = category, label = label)) +
  geom_vline(xintercept = 0, col = "black", linetype = 'dashed') +
  geom_hline(yintercept = myy.intercept, col = "black", linetype = 'dashed') + 
  geom_point(size = 2, show.legend = FALSE) + 
  scale_color_manual(values = c("lightcoral", "grey", "skyblue")[lab.sub], 						## c("lightcoral", "grey", "skyblue")
                     labels = c("Hazardous", "Not significant", "Protective")[lab.sub]) + 		## c("Hazardous", "Not significant", "Protective")
  guides(color="none") +	
  coord_cartesian(ylim = c(0, 5), xlim = c(-0.9, 0.8)) + 
  labs(color = '', #legend_title, 
       x = "Log2(HR)", y = "-Log10(P-value)") + 
#  scale_x_continuous(breaks = seq(-10, 10, 2)) + 
#  ggtitle('put title here') + 													# Plot title 
  geom_text_repel(max.overlaps = Inf) + 											# To show all labels 
  theme_classic()
pdf(file = myFig, width = 3.5, height = 4) 
mygg
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.2] KM plot for Plasma.cells_C12.ES
rm(list=ls())
myinf1 = "Okayama_GSE31210_Data.rda"
myFig = "Fig5_KM_Okayama_PlasmaC12.pdf"
load(myinf1)
data = ORI.data
info = ORI.info

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (months)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.3] ForestPlot  for Plasma.cells_C12.ES
rm(list=ls())
myinf1 = "Okayama_GSE31210_Data.rda"
myFig = "Fig5_ForestPlot_Okayama_PlasmaC12.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

library(survival)
se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
xx$mytag = ifelse(xx$mys>median(xx$mys), 1, 0)
xx$smoking = ifelse(xx$smoking.status ==" Ever-smoker", 1 , 0)
mycox = coxph(Surv(t.surv, e.surv)~mys + age + gender + smoking + pstage.iorii, xx) 
mycox = summary(mycox)
tmp = mycox$coefficients
PV1 = tmp[,5]
HR1 = tmp[,2]
tmp = mycox$conf.int
HR.lo = tmp[,3]
HR.hi = tmp[,4]
res = data.frame(PV1, HR1, HR.lo, HR.hi)
row.names(res) = row.names(tmp)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(forestplot)

ids2 = row.names(res)
ids2 = c("Plasma_C12 score", "Age", "Sex (Male vs. Female)", "Smoking (Ever vs. Never)", "Stage II (vs. I)")
mhr = res$HR1
mup = res$HR.hi
mlw = res$HR.lo
myp = res$PV1
myp = as.character(ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1)))
hr = paste(formatC(mhr, 2, format = 'f'), ' (', formatC(mlw, 2, format = 'f') ,' to ',formatC(mup, 2, format = 'f'),')' , sep="")

data = tibble(Feature = ids2, mean  = mhr, lower = mlw, upper = mup,  hr=hr, coxphP = myp)
header <- tibble(Feature = c("Feature"), hr='HR (95% CI)',coxphP = c('P-value'))
data = bind_rows(header, data)

limits = c(0.4, 5)  
xticks <- c(0.4,1,5)  
xtlab <- rep(c(TRUE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab

cairo_pdf(myFig, width = 5, height = 2.5) 
data %>%  
  forestplot(labeltext = c(Feature, hr,coxphP),
             graph.pos = 2,
             zero = 1,
             is.summary=c(TRUE, rep(FALSE, nrow(data)-1)), 
             xlog = FALSE,
             clip = limits,
             xlab = "HR",
             ci.vertices = TRUE,
             hrzl_lines= TRUE,
             boxsize = 0.2, 
             line.margin = 0.1,
             mar = unit(rep(0.1, times = 4), "mm"),
             xticks = xticks,
             graphwidth = "auto",
             graphhight = "auto",
             align = rep("l", 3),
             txt_gp=fpTxtGp(label = gpar(cex = 0.9),
                            title = gpar(cex = 0.9),
                            ticks = gpar(cex = 0.9),
                            xlab = gpar(cex = 0.9)),
             col = fpColors(box = "black",
                            line = "black")) |>
  fp_set_zebra_style("#EFEFEF")
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.3] Boxplot  for Plasma.cells_C12.ES (Recurrence --> Smoking --> Stage)
rm(list=ls())
myinf1 = "Okayama_GSE31210_Data.rda"

myFigA = "Fig5_BoxPlot_Okayama_PlasmaC12_RvsNR.pdf"
myFigB = "Fig5_BoxPlot_Okayama_PlasmaC12_Smoking.pdf"
myFigC = "Fig5_BoxPlot_Okayama_PlasmaC12_Stage.pdf"

load(myinf1)
data = ORI.data
info = ORI.info
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
mydat = cbind(mys, info)

se = which(mydat$relapse %in% c(" relapsed", " not relapsed"))
xx = mydat[se,]
mycat = ifelse(xx$relapse== " relapsed", "R", "NR")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]

xx1 = data$mys[data$mycat=="R"]
xx2 = data$mys[data$mycat=="NR"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")

mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(mycat, mys), y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
mygg

pdf(file = myFigA, width = 2.2, height = 3) 
mygg
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$smoking.status %in% c(" Ever-smoker", " Never-smoker"))
xx = mydat[se,]
mycat = ifelse(xx$smoking.status==" Ever-smoker", "Smoker", "Non-smoker")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]

xx1 = data$mys[data$mycat=="Smoker"]
xx2 = data$mys[data$mycat=="Non-smoker"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(mycat, mys), y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)

pdf(file = myFigB, width = 2.2, height = 3) 
mygg
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$pstage.iorii %in% c(" I", " II"))
xx = mydat[se,]
mycat = ifelse(xx$pstage.iorii==" I", "I", "II")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]
data$mycat = as.factor(as.character(data$mycat))

xx1 = data$mys[data$mycat=="I"]
xx2 = data$mys[data$mycat=="II"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")

mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=mycat, y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
#mygg

pdf(file = myFigC, width = 2.2, height = 3) 
mygg
dev.off()




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[5.1] Vocano plot to show univariable Cox results with RFS (Okayama result)
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_Volcano_Okayama_Prognosis.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------
library(survival)

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

PV1 = PV2 = rep(0, ncol(data))
HR1 = HR2 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
	cat("\r", k)
	mytf = as.numeric(data[,k])
	xx = cbind(mytf, info)
	xx = xx[xx[, "t.surv"]>0,]
	mycox = coxph(Surv(t.surv, e.surv)~mytf, xx) 
	mycox = summary(mycox)
	PV1[k] = mycox$coefficients[5]
	tmp = mycox$conf.int
	HR1[k] = tmp[1]

	mycox = coxph(Surv(t.surv, e.surv)~mytf + age + gender + smoking.status + pstage.iorii, xx) 
	mycox = summary(mycox)
	PV2[k] = mycox$coefficients["mytf", 5]
	tmp = mycox$conf.int
	HR2[k] = tmp["mytf", 1]
}
QV1 = p.adjust(PV1, method="BH")
QV2 = p.adjust(PV2, method="BH")

res = data.frame(PV1, HR1, PV2, HR2)
row.names(res) = colnames(data)
res = res[order(res$PV1),]
xx = row.names(res)
xx = gsub("\\.ES", "", xx)
xx = gsub("\\.cells", "", xx)
row.names(res) = xx
se = which(!row.names(res) %in% c("Endothelial", "Epithelial", "Fibroblasts"))
res = res[se,]


data = res
data$FDR = p.adjust(data$PV1, method="BH")
se = which(data$FDR<0.05)
length(se)

range(log2(data$HR1))
 
#---------------------------------------
data$category= ifelse(data$PV1>0.05, "NS", ifelse(data$HR1>1, "Hazardous", "Protective"))
data$label =  rep("", nrow(data))
data$label[1:3] = row.names(data)[1:3]
lab.sub = which(c("Hazardous", "NS", "Protective") %in% data$category)
myy.intercept = -log10(0.05)
#---------------------------------------
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

mygg <- ggplot(data, aes(x = log2(HR1), y = -log10(PV1), col = category, label = label)) +
  geom_vline(xintercept = 0, col = "black", linetype = 'dashed') +
  geom_hline(yintercept = myy.intercept, col = "black", linetype = 'dashed') + 
  geom_point(size = 2, show.legend = FALSE) + 
  scale_color_manual(values = c("lightcoral", "grey", "skyblue")[lab.sub], 						## c("lightcoral", "grey", "skyblue")
                     labels = c("Hazardous", "Not significant", "Protective")[lab.sub]) + 		## c("Hazardous", "Not significant", "Protective")
  guides(color="none") +	
  coord_cartesian(ylim = c(0, 5), xlim = c(-0.9, 0.8)) + 
  labs(color = '', #legend_title, 
       x = "Log2(HR)", y = "-Log10(P-value)") + 
#  scale_x_continuous(breaks = seq(-10, 10, 2)) + 
#  ggtitle('put title here') + 													# Plot title 
  geom_text_repel(max.overlaps = Inf) + 											# To show all labels 
  theme_classic()

pdf(file = myFig, width = 3.5, height = 4) 
mygg
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.2] KM plot for Plasma.cells_C12.ES
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (days)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.3] ForestPlot  for Plasma.cells_C12.ES
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_ForestPlot_Okayama_PlasmaC12.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

library(survival)
se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
xx$mytag = ifelse(xx$mys>median(xx$mys), 1, 0)
xx$smoking = ifelse(xx$smoking.status ==" Ever-smoker", 1 , 0)
mycox = coxph(Surv(t.surv, e.surv)~mys + age + gender + smoking + pstage.iorii, xx) 
mycox = summary(mycox)
tmp = mycox$coefficients
PV1 = tmp[,5]
HR1 = tmp[,2]
tmp = mycox$conf.int
HR.lo = tmp[,3]
HR.hi = tmp[,4]
res = data.frame(PV1, HR1, HR.lo, HR.hi)
row.names(res) = row.names(tmp)
res

                         PV1       HR1     HR.lo     HR.hi
mys             0.0009664662 0.6227619 0.4700760 0.8250417
age             0.0129855199 1.0488086 1.0101067 1.0889932
gender male     0.9916204870 0.9963250 0.5011871 1.9806245
smoking         0.9378498725 0.9724101 0.4813119 1.9645917
pstage.iorii II 0.0002403648 2.6737766 1.5818401 4.5194715

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(forestplot)

ids2 = row.names(res)
ids2 = c("Plasma_C12 score", "Age", "Sex (Male vs. Female)", "Smoking (Ever vs. Never)", "Stage II (vs. I)")
mhr = res$HR1
mup = res$HR.hi
mlw = res$HR.lo
myp = res$PV1
myp = as.character(ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1)))
hr = paste(formatC(mhr, 2, format = 'f'), ' (', formatC(mlw, 2, format = 'f') ,' to ',formatC(mup, 2, format = 'f'),')' , sep="")

data = tibble(Feature = ids2, mean  = mhr, lower = mlw, upper = mup,  hr=hr, coxphP = myp)
header <- tibble(Feature = c("Feature"), hr='HR (95% CI)',coxphP = c('P-value'))
data = bind_rows(header, data)
max(mhr)
max(mup)
min(mlw)
data

limits = c(0.4, 5)  
xticks <- c(0.4,1,5)  
xtlab <- rep(c(TRUE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab

cairo_pdf(myFig, width = 5, height = 2.5) 
data %>%  
  forestplot(labeltext = c(Feature, hr,coxphP),
             graph.pos = 2,
             zero = 1,
             is.summary=c(TRUE, rep(FALSE, nrow(data)-1)), 
             xlog = FALSE,
             clip = limits,
             xlab = "HR",
             ci.vertices = TRUE,
             hrzl_lines= TRUE,
             boxsize = 0.2, 
             line.margin = 0.1,
             mar = unit(rep(0.1, times = 4), "mm"),
             xticks = xticks,
             graphwidth = "auto",
             graphhight = "auto",
             align = rep("l", 3),
             txt_gp=fpTxtGp(label = gpar(cex = 0.9),
                            title = gpar(cex = 0.9),
                            ticks = gpar(cex = 0.9),
                            xlab = gpar(cex = 0.9)),
             col = fpColors(box = "black",
                            line = "black")) |>
  fp_set_zebra_style("#EFEFEF")
dev.off()

##------------------
cd ~/Dropbox/02__Project/MS/MS_BCMLC/Draft/Figures/Raw
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_ForestPlot_Okayama_PlasmaC12.pdf ./


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.3] Boxplot  for Plasma.cells_C12.ES (Recurrence --> Smoking --> Stage)
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFigA = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_BoxPlot_Okayama_PlasmaC12_RvsNR.pdf"
myFigB = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_BoxPlot_Okayama_PlasmaC12_Smoking.pdf"
myFigC = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_BoxPlot_Okayama_PlasmaC12_Stage.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

library(survival)
se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
mydat = cbind(mys, info)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$relapse %in% c(" relapsed", " not relapsed"))
xx = mydat[se,]
mycat = ifelse(xx$relapse== " relapsed", "R", "NR")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]

xx1 = data$mys[data$mycat=="R"]
xx2 = data$mys[data$mycat=="NR"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(mycat, mys), y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
mygg

pdf(file = myFigA, width = 2.2, height = 3) 
mygg
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$smoking.status %in% c(" Ever-smoker", " Never-smoker"))
xx = mydat[se,]
mycat = ifelse(xx$smoking.status==" Ever-smoker", "Smoker", "Non-smoker")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]

xx1 = data$mys[data$mycat=="Smoker"]
xx2 = data$mys[data$mycat=="Non-smoker"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="l")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=reorder(mycat, mys), y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
#mygg

pdf(file = myFigB, width = 2.2, height = 3) 
mygg
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

se = which(mydat$pstage.iorii %in% c(" I", " II"))
xx = mydat[se,]
mycat = ifelse(xx$pstage.iorii==" I", "I", "II")
data = cbind(xx, mycat)
se = which(!is.na(data$mycat))
data = data[se,]
data$mycat = as.factor(as.character(data$mycat))

xx1 = data$mys[data$mycat=="I"]
xx2 = data$mys[data$mycat=="II"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))
myp = wilcox.test(xx1, xx2, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


mycol = c('lightcoral', 'lightblue')[2:1]
mygg <- ggplot(data, aes(x=mycat, y=mys,  fill=mycat)) + 
  geom_boxplot(size=1) +
  scale_fill_manual(values=mycol) +
  labs(title='',x= '', y = 'Plasma_C12 score')+
  theme_classic() + 
  theme(legend.position="none")

posx = 0.5
posy = 1
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=10)))
mygg <- mygg + annotation_custom(mytxt)
#mygg

pdf(file = myFigC, width = 2.2, height = 3) 
mygg
dev.off()


##------------------
cd ~/Dropbox/02__Project/MS/MS_BCMLC/Draft/Figures/Raw
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_BoxPlot_Okayama_PlasmaC12_RvsNR.pdf ./
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_BoxPlot_Okayama_PlasmaC12_Smoking.pdf ./
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_BoxPlot_Okayama_PlasmaC12_Stage.pdf ./


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.4] KM plot for Plasma.cells_C12.ES  -- In Stage I
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_StageI.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

se = which(mydat$pstage.iorii ==" I")
mydat = mydat[se,]
mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (days)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

##------------------
cd ~/Dropbox/02__Project/MS/MS_BCMLC/Draft/Figures/Raw
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_StageI.pdf ./



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.5] KM plot for Plasma.cells_C12.ES  -- In Stage II
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_StageII.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

se = which(mydat$pstage.iorii ==" II")
mydat = mydat[se,]
mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (days)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

##------------------
cd ~/Dropbox/02__Project/MS/MS_BCMLC/Draft/Figures/Raw
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_StageII.pdf ./



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.6] KM plot for Plasma.cells_C12.ES  -- In Smokers
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_Smoker.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

se = which(mydat$smoking.status ==" Ever-smoker")
mydat = mydat[se,]
mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (days)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

se = which(mydat$smoking.status ==" Ever-smoker")
mydat = mydat[se,]
mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (days)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

##------------------
cd ~/Dropbox/02__Project/MS/MS_BCMLC/Draft/Figures/Raw
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_Smoker.pdf ./



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.7] KM plot for Plasma.cells_C12.ES  -- In Never-Smokers
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/data/Okayama_GSE31210_Data.rda"

myFig = "/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_Nonsmoker.pdf"

load(myinf1)
data = ORI.data
info = ORI.info

#---------------------------------------

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

se = which(colnames(data)=="Plasma.cells_C12.ES")
mys = as.numeric(data[,se])
xx = cbind(mys, info)
xx = xx[xx[, "t.surv"]>0,]
mydat = cbind(info, mys)
colnames(mydat)[1:2] = c("time", "event")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(survival)

pdf(myFig, width=4, height= 3)
par(mfrow=c(1, 1))
par(lend=2)
par(tcl= -0.15)  # 
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

se = which(mydat$smoking.status ==" Never-smoker")
mydat = mydat[se,]
mycat = ifelse(mydat$mys> median(mydat$mys), 1, 0)
xx = cbind(mycat, mydat)
mytex = paste(c("High-Score", "Low-Score"), " (n=", c(sum(xx$mycat==1, na.rm=T), sum(xx$mycat==0, na.rm=T)), ")", sep="")

tmp <- coxph(Surv(time, event)~mycat, data=xx)
tmp = summary(tmp)$coefficients
myp = tmp[5]
myh = round(tmp[2],3)
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
myp = paste("P=", myp, sep="")
myh = paste("HR=", myh, sep="")

tmp <- survfit(Surv(time, event)~mycat, data=xx)
mycol = c('lightcoral', 'skyblue')

plot(tmp, col=mycol[2:1], ylim=c(0, 1), xlab="Survival time (days)", ylab="Probability of survival (RFS)", cex.axis=0.9, cex.lab=0.9, main="")

max.x = max(xx$time, na.rm=T)
min.x = min(xx$time, na.rm=T)
posx = min.x + 0*(max.x - min.x)/10
posy = 0.2
legend(posx, posy,  mytex, lwd=2, col=mycol, bty="n", cex=0.9)
posx = max.x - 1*(max.x - min.x)/10
posy = 0.85
adjy = 0.05
text(posx, posy, pos=3, myp, cex=0.9)
text(posx, posy+adjy, pos=1, myh, cex=0.9)

dev.off()

##------------------
cd ~/Dropbox/02__Project/MS/MS_BCMLC/Draft/Figures/Raw
scp chaoc@curie.dldcc.bcm.edu:/mount/ictr1/chenglab/cc59/WorSpa/BCM/Publication/Y2025_scRNAseq_LUAD/figure/Fig4_KM_Okayama_PlasmaC12_Nonsmoker.pdf ./




