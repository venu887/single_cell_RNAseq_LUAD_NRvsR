rm(list = ls()) 
#https://pmc.ncbi.nlm.nih.gov/articles/PMC10654435/
#https://pmc.ncbi.nlm.nih.gov/articles/PMC10442056/
# Forest plot:https://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups
# NEW data GSE90623 https://www.oncotarget.com/article/19161/text/
# New data: https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.794293/full
rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(survival)
library(survminer)
library(patchwork)
library(ggpubr)
library(forestplot)
library(forestmodel)
output<-"/home/u251079/scRNA/Lung_scRNA/"
# OTHER LUNG datasets
"/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung"



#@@@@@@@@@@@@@@@ 1.Okayama
info_out<-"/home/u251079/scRNA/Lung_scRNA/Okayama.txt"
Okayama <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(Okayama)
Okayama$gender <- trimws(Okayama$gender)
Okayama$Stage<-trimws(Okayama$pstage.iorii)
Okayama$Stage <- factor(Okayama$Stage, levels = c("I", "II"))
names(Okayama)[1:28]<- gsub("\\.ES$", "", names(Okayama)[1:28])


Okayama2<-Okayama
#Okayama2=Okayama[Okayama$t.surv <= 1825,]
med_age<-median(Okayama2$age)
Okayama2$Age<-ifelse(Okayama2$age <= med_age, "age≤61", "age≥61")
Okayama2$Age <- factor(Okayama2$Age, levels = c("age≤61", "age≥61"))
str(Okayama2)


#@@@@@@  KM PLOT
gene <- "Plasma.cells_C12"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)
km_plot <- ggsurvplot(fit, data = Okayama2, pval = TRUE,  pval.method = TRUE,
                      title = gene,
                      legend.title = gene,
                      legend.labs = c("High scores", "Low scores"),
                      xlab="Time (days)",
                      palette = c("#6E8B3D", "orangered"))
print(km_plot)

km_plot$plot<-km_plot$plot+ theme(plot.title = element_text(size = 12),
                                  legend.position = c(0.95, 0.05), 
                                  legend.justification = c(1, 0))

gene <- "Clara.C5"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)

km_plot1 <- ggsurvplot(fit, data = Okayama2, pval = TRUE, pval.method = TRUE,
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       xlab = "Time (days)",
                       palette = c("#6E8B3D", "orangered"))
print(km_plot1)
km_plot1$plot <- km_plot1$plot + theme(plot.title = element_text(size = 12),
                                       legend.position = c(0.95, 0.05), 
                                       legend.justification = c(1, 0))


gene <- "Clara.C3"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)

km_plot2 <- ggsurvplot(fit, data = Okayama2, pval = TRUE, pval.method = TRUE,
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       xlab = "Time (days)",
                       palette = c("#6E8B3D", "orangered"))

print(km_plot2)
km_plot2$plot <- km_plot2$plot + theme(plot.title = element_text(size = 12),legend.position = c(0.95, 0.05), 
                                       legend.justification = c(1, 0))

gene <- "Clara.C3"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)

km_plot2 <- ggsurvplot(fit, data = Okayama2, pval = TRUE, pval.method = TRUE,
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       xlab = "Time (days)",
                       palette = c("#6E8B3D", "orangered"))
print(km_plot2)
km_plot2$plot <- km_plot2$plot + theme(plot.title = element_text(size = 12),
                                       legend.position = c(0.95, 0.05), 
                                       legend.justification = c(1, 0))


gene <- "Endothelial.cells.C4"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)

km_plot3 <- ggsurvplot(fit, data = Okayama2, pval = TRUE, pval.method = TRUE,
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       xlab = "Time (days)",
                       palette = c("#6E8B3D", "orangered"))
print(km_plot3)
km_plot3$plot <- km_plot3$plot + theme(plot.title = element_text(size = 12),
                                       legend.position = c(0.95, 0.05), 
                                       legend.justification = c(1, 0))


gene <- "Ciliated.cells.C7"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)

km_plot4 <- ggsurvplot(fit, data = Okayama2, pval = TRUE, pval.method = TRUE,
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       xlab = "Time (days)",
                       palette = c("#6E8B3D", "orangered"))
print(km_plot4)
km_plot4$plot <- km_plot4$plot + theme(plot.title = element_text(size = 12),
                                       legend.position = c(0.95, 0.05), 
                                       legend.justification = c(1, 0))


gene <- "Mac.DC_C11"
gene_expression <- Okayama2[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")

fit <- survfit(Surv(t.surv, e.surv) ~ group, data = Okayama2)

km_plot5 <- ggsurvplot(fit, data = Okayama2, pval = TRUE, pval.method = TRUE,
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       xlab = "Time (days)",
                       palette = c("#6E8B3D", "orangered"))

print(km_plot5)
km_plot5$plot <- km_plot5$plot + theme(plot.title = element_text(size = 12),
                                       legend.position = c(0.95, 0.05), 
                                       legend.justification = c(1, 0))

plots<-km_plot$plot+km_plot1$plot+km_plot2$plot+km_plot3$plot+km_plot4$plot+km_plot5$plot
plots




pdf(paste0(output, "Okayama_RFS_all.pdf"), width = 12, height =8)
print(plots)
dev.off() 


#@@@@@@@@ Boxplot 
#for Gender
p1 <- ggboxplot(Okayama2, x = "gender", y = gene, fill = "gender",
                palette = c("#66CCEE", "#EE6677"), 
                add = "jitter") +
  stat_compare_means(method = "wilcox.test", label ="p", label.x = 1.5) + 
  labs(title = NULL, y = "Signature scores", x = "Gender")+
  theme(plot.title = element_text(size = 12),  # Set title size to 12
        legend.position = "none")  # Remove the legend


# for Age
p2 <- ggboxplot(Okayama2, x = "Age", y = gene, fill = "Age",
                palette = c("#66CCEE", "#EE6677"), 
                add = "jitter") +
  stat_compare_means(method = "wilcox.test", label ="p", label.x = 1.5) + 
  labs(title = NULL, y = "Signature scores", x = "Age Group")+
  theme(plot.title = element_text(size = 12),  # Set title size to 12
        legend.position = "none")  # Remove the legend

# for Stage
p3 <- ggboxplot(Okayama2, x = "pathological.stage", y = gene, fill = "pathological.stage",
                palette = c("#66CCEE", "#EE6677", "yellow"), 
                add = "jitter") +
  stat_compare_means(method = "wilcox.test", label ="p", label.x = 1.5, ref.group = "IA", label.y = -9) + 
  labs(title = NULL, y = "Signature scores", x = "pathological.stage")+
  theme(plot.title = element_text(size = 12),  # Set title size to 12
  legend.position = "none")  # Remove the legend

p3

# For males and females
Okayama_male <- subset(Okayama2, gender == "male")
Okayama_female <- subset(Okayama2, gender == "female")
# Boxplot for Males (Stage I vs. II)
p_male <- ggboxplot(Okayama_male, x = "Stage", y = gene, fill = "Stage",
                    palette = c("purple", "orange"), 
                    add = "jitter") +
  stat_compare_means(method = "wilcox.test", label ="p", label.x = 1.5) +
  labs(title = NULL, y = "Signature scores", x = "Males_Stage")+
  theme(plot.title = element_text(size = 12),  # Set title size to 12
        legend.position = "none")  # Remove the legend


# For Females (Stage I vs. II)
p_female <- ggboxplot(Okayama_female, x = "Stage", y = gene, fill = "Stage",
                      palette = c("blue", "orange"), 
                      add = "jitter") +
  stat_compare_means(method = "wilcox.test", label ="p", label.x = 1.5) +
  labs(title = NULL, y = "Signature scores", x = "Females_Stages")+
  theme(plot.title = element_text(size = 12),  # Set title size to 12
        legend.position = "none")  # Remove the legend

# Boxplot for Females (Age ≤61 vs. ≥61)
p_female_age <- ggboxplot(Okayama_female, x = "Age", y = gene, fill = "Age",
                          palette = c("blue", "orange"), 
                          add = "jitter") +
  stat_compare_means(method = "wilcox.test", label ="p", label.x = 1.5) +
  labs(title = NULL, 
       y = "Signature scores", x = "Females_Age")+
  theme(plot.title = element_text(size = 12),  # Set title size to 12
        legend.position = "none")  # Remove the legend


combined_plot <- ggarrange(p1, p2, p3, p_male, p_female, p_female_age, 
                                   ncol = 3, nrow = 2, 
                                   common.legend = TRUE, legend = "none") +
  ggtitle("Plasma Cells Expression in Okayama Data") +  # Add title to the whole plot
  theme(plot.title = element_text(hjust = 0.5, size = 14))
combined_plot
# Save the plot
# pdf(filename = paste0(output, "Okayama_Box.pdf"), plot = combined_plot, width = 10, height = 6, dpi = 1200)
pdf(paste0(output, "Okayama_Box.pdf"), width = 10, height =8)
print(combined_plot)
dev.off() 



#@@@@@@@ FOREST Plots
library(forestmodel)
library(patchwork)
# fit.coxph <- coxph(Surv(t.surv, e.surv) ~Plasma.cells_C12 +gender+Age+Stage, data=Okayama2)
# p1<-ggforest(fit.coxph, data=Okayama2, main = "Okayama_Multivariate")
# print(p1)
# Plot 1: Okayama Multivariate
genes <- c("Plasma.cells_C12", "Clara.C5", "Clara.C3", 
           "Endothelial.cells.C4", "Ciliated.cells.C7", "Mac.DC_C11")
p1 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 + gender + Age + Stage, data = Okayama2)) 
p2 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Clara.C5 + gender + Age + Stage, data = Okayama2)) 
p3 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Clara.C3  + gender + Age + Stage, data = Okayama2)) 
p4 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Endothelial.cells.C4 + gender + Age + Stage, data = Okayama2)) 
p5 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Ciliated.cells.C7 + gender + Age + Stage, data = Okayama2)) 
p6 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Mac.DC_C11 + gender + Age + Stage, data = Okayama2)) 


Okayama_plot <- p1 + p2+ p3+p4 +p5+p6+ plot_layout(nrow = 3, ncol = 2)+
  plot_annotation(title = "Okayama Multivariate Analysis")
print(Okayama_plot)

# pdf(filename = paste0(output, "Okayama_Forest_Plot.pdf"), 
#        plot = Okayama_plot, width = 17, height = 5.5, dpi = 1200)

pdf(paste0(output, "Okayama_Forest_Plot.pdf"), width = 17, height =8)
print(Okayama_plot)
dev.off() 

library(survival)
library(ggplot2)
library(dplyr)

# Define models
models <- list(
  "Plasma cells C12" = coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 + gender + Age + Stage, data = Okayama2),
  "Clara C5" = coxph(Surv(t.surv, e.surv) ~ Clara.C5 + gender + Age + Stage, data = Okayama2),
  "Clara C3" = coxph(Surv(t.surv, e.surv) ~ Clara.C3 + gender + Age + Stage, data = Okayama2),
  "Endothelial cells C4" = coxph(Surv(t.surv, e.surv) ~ Endothelial.cells.C4 + gender + Age + Stage, data = Okayama2),
  "Ciliated cells C7" = coxph(Surv(t.surv, e.surv) ~ Ciliated.cells.C7 + gender + Age + Stage, data = Okayama2),
  "Mac.DC C11" = coxph(Surv(t.surv, e.surv) ~ Mac.DC_C11 + gender + Age + Stage, data = Okayama2)
)

# Extract HR, CI, P-value, and Sample Size
model_results <- lapply(names(models), function(name) {
  model <- models[[name]]
  summary_model <- summary(model)
  
  HR <- summary_model$coefficients[1, "exp(coef)"]
  CI_lower <- summary_model$conf.int[1, "lower .95"]
  CI_upper <- summary_model$conf.int[1, "upper .95"]
  p_value <- summary_model$coefficients[1, "Pr(>|z|)"]
  num_samples <- model$n
  
  data.frame(Model = name, HR = HR, CI_lower = CI_lower, CI_upper = CI_upper, 
             P_value = p_value, Samples = num_samples)
})

forest_data <- bind_rows(model_results)
forest_data <- forest_data %>%
  mutate(HR_label = sprintf("HR: %.2f", HR),
         P_label = sprintf("P: %.3g", P_value))

forest_data_sorted <- forest_data[order(forest_data$HR, decreasing = TRUE), ]
forest_data_sorted$Model <- factor(forest_data_sorted$Model, levels = forest_data_sorted$Model)
forest_plot <- ggplot(forest_data_sorted, aes(y = Model, x = HR, xmin = CI_lower, xmax = CI_upper)) +
  geom_point(color = "#00688B", size = 3) +
  geom_errorbarh(height = 0.2, color = "#009ACD") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10(limits = c(0.45, 1.8), breaks = c(0.45, 0.5, 1, 1.5, 1.8)) +
  theme_minimal() +
  labs(title = "Okayama Multivariate Analysis (Gender, Age, Stage)", 
       x = "Log2(Hazard Ratio) with 95% CI", y = NULL) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(aes(label = HR_label), hjust = -0.3, vjust = -0.5, size = 4, color = "black") +  # HR text
  geom_text(aes(label = P_label), hjust = 1.3, vjust = -0.5, size = 4, color = "black")  # P-value text

print(forest_plot)

ggsave(paste0(output,"Okayama_multi_Plot.pdf"), plot = forest_plot, width = 7, height = 4, dpi = 1200)







#@@@@@@@@@@@@@@@@@@@@@@ BCMLC
# BCMLC is significant only after removing patiens with out Pstage information.
rm(list = ls())
info_out<-"/home/u251079/scRNA/Lung_scRNA/BCMLC.txt"
output <- "/home/u251079/scRNA/Lung_scRNA/"
BCMLC <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(BCMLC)
names(BCMLC)[1:20] <- gsub("\\.ES", "", names(BCMLC)[1:20])
BCMLC$DateRecurrence <- as.character(BCMLC$DateRecurrence)
BCMLC$Recurrence_Status <- ifelse(BCMLC$DateRecurrence != "", "Recurrence", "Non-Recurrence")
BCMLC$Recurrence_Status <- factor(BCMLC$Recurrence_Status, levels = c("Non-Recurrence", "Recurrence"))
table(BCMLC$Recurrence_Status)


BCMLC$Race <- toupper(trimws(BCMLC$Race))
BCMLC$Race <- ifelse(BCMLC$Race == "A", "ASIAN",
                     ifelse(BCMLC$Race == "B", "BLACK",
                            ifelse(BCMLC$Race == "W", "WHITE", BCMLC$Race)))

BCMLC$Race[BCMLC$Race == ""] <- NA
BCMLC$Pstage[BCMLC$Pstage == ""] <- NA
BCMLC$Sex <- ifelse(BCMLC$Sex %in% c("M", "m"), "Male", "Female")
Age_med <- median(BCMLC$Age.at.Op)
BCMLC$Age <- ifelse(BCMLC$Age.at.Op <= Age_med, "Age ≤66.7", "Age >66.7")
BCMLC$Age <- factor(BCMLC$Age, levels = c("Age ≤66.7", "Age >66.7"))
# BCMLC1 <- BCMLC[!is.na(BCMLC$Pstage), ]
# BCMLC2 <- BCMLC[BCMLC$t.surv<60, ]


# THis information is already used by Chao, Clinical_info_20231111.txt

myinf3 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/Clinical_info_20231202.txt"
info2 = read.table(myinf3, sep="\t", header=T, quote="", fill = TRUE)
info2$LID
info2$LID <- gsub(" ", "", info2$LID)
xx<-grep("Adenocarcinoma",info2$Histologic.type)
info3<-info2[xx,]
info3$LID
info3$LID <- gsub("LC67/LC324", "LC67", info3$LID)

BCMLC$sample <- paste0("LC", gsub("T$", "", rownames(BCMLC)))
BCMLC$Rdeath<-ifelse(BCMLC$sample %in% info3$LID, info3$Recurrence.or.death, NA)
BCMLC$Rtime<-ifelse(BCMLC$sample %in% info3$LID, info3$RFS.or.DFS, NA)

myinf2 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/Clinical_info_20231212.txt"
info1 = read.table(myinf2, sep="\t", header=T, quote="")
info1$LC
names(info1)



BCMLC$smoking<-ifelse(BCMLC$sample %in% info1$LC, info1$SmokingHx, NA)
table(BCMLC$smoking)
BCMLC$smoking<-trimws(BCMLC$smoking)

table(BCMLC$Recurrence_Status)
#@@@@@@@@@@@@@@@@@@@@@@@ Box plot based on recurrence status
genes<-colnames(BCMLC)[1:28]
genes<-c("Plasma.cells_C12", "B.cells_C13") # significant that patients with non-recurrent compared with recurrent
p_box_list <- lapply(genes, function(gene) {
  ggboxplot(BCMLC, 
            x = "Recurrence_Status",  
            y = gene, fill = "Recurrence_Status",
            palette = c("#FFC125", "#CD5555")) +
    stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +  # Compare smokers vs non-smokers
    labs(title = gene, 
         y = "Signature Scores", 
         x = "Recurrence Status") +
    theme(plot.title = element_text(size = 12))
})
p_box_combined1 <- ggarrange(plotlist = p_box_list, ncol = 2, nrow = 1)
p_box_combined1


library(ggpubr)
#@@@@@@@@@@@@@@@@@@@@@@@ Recurrent patiens with Smoking and non smoking
genes<-colnames(BCMLC)[1:28]
table(BCMLC$Recurrence_Status)
genes<-c("Clara.C5", "Clara.C3") # "BASC" cells are two clusters with sample size less than 100 cells in each cluster
p_box_list <- lapply(genes, function(gene) {
  ggboxplot(BCMLC, 
            x = "Recurrence_Status",  
            y = gene, 
            fill = "smoking", 
            palette = c("lightblue", "lightpink")) +
    stat_compare_means(method = "wilcox.test", aes(group = smoking), label = "p") +  # Compare smokers vs non-smokers
    labs(title = gene, 
         y = "Signature Scores", 
         x = "Recurrence Status") +
    theme(plot.title = element_text(size = 12))
})
p_box_combined1 <- ggarrange(plotlist = p_box_list, ncol = 2, nrow = 1)
p_box_combined1

#@@@@@@@@@@@@@@@.   Signature scores differentiating patients into smokers and non-smokers
genes <- colnames(BCMLC)[1:28]
sig_genes<-c("GaDelT_C14","Clara.C3","Clara.C5", "Ciliated.cells.C7", "Plasma.cells_C12")
genes<-sig_genes
p_box_list <- lapply(genes, function(gene) {
  ggboxplot(BCMLC, 
            x = "smoking",  
            y = gene, 
            fill = "smoking", 
            palette = c("#8DEEEE", "#EEA2AD")) +
    stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +
    labs(title = gene, 
         y = "Signature Scores", 
         x = "smoking") +
    theme(plot.title = element_text(size = 12), legend.position = "none")
})

p_box_combined <- ggarrange(plotlist = p_box_list, ncol = 3, nrow = 2)
p_box_combined <- annotate_figure(p_box_combined, 
                                  top = text_grob("Signature scores differentiating patients into smokers and non-smokers", 
                                                  face = "bold", size = 14))

print(p_box_combined)

# Define cell types
genes <- c("Plasma.cells_C12", "Clara.C5", "Clara.C3", 
           "Endothelial.cells.C4", "Ciliated.cells.C7", "Mac.DC_C11")

gene <- c("Plasma.cells_C12") # "B.cells_C13"
BCMLC$Rtime <- as.numeric(BCMLC$Rtime)
BCMLC$Rtime<-BCMLC$Rtime/30
BCMLC$Recurrence_S<-ifelse(BCMLC$Recurrence_Status %in% "Recurrence", 1, 0)
BCMLC_NA <- BCMLC[!is.na(BCMLC$Rtime) & !is.na(BCMLC$Recurrence_Status), ]
gene_expression <- BCMLC_NA[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC_NA$Rtime, BCMLC_NA$Recurrence_S)
fit <- survfit(surv_object ~ group)
plot_km <- ggsurvplot(fit, 
                      data = BCMLC_NA,
                      pval = TRUE, pval.method = TRUE,         
                      title = paste0(gene, " (Recurrence=1, Non-Recurrence=0)"),
                      legend.title = gene,
                      legend.labs = c("High scores", "Low scores"),    
                      palette = c("#6E8B3D", "orangered"),
                      xlab = "Time (Months)",
                      ylab = "Recurrence-Free Survival Probability")

plot_km$plot <- plot_km$plot + theme(plot.title = element_text(size = 12),                                        
                                     legend.position = c(0.95, 0.05),                                         
                                     legend.justification = c(1, 0))
plot_km$plot

genes <- c("Plasma.cells_C12", "Clara.C5", "Clara.C3", 
           "Endothelial.cells.C4", "Ciliated.cells.C7", "Mac.DC_C11")
plots <- list()  
for (gene in genes) {
  gene_expression <- BCMLC_NA[, gene]
  median_expr <- median(gene_expression, na.rm = TRUE)
  group <- factor(ifelse(gene_expression > median_expr, "High", "Low"), levels = c("High", "Low"))
  surv_object <- Surv(BCMLC_NA$Rtime, BCMLC_NA$Recurrence_S)
  fit <- survfit(surv_object ~ group)
  plot_km <- ggsurvplot(fit, 
                        data = BCMLC_NA,
                        pval = TRUE, pval.method = TRUE,         
                        title = paste0(gene, " (Recurrence=1, Non-Recurrence=0)"),
                        legend.title = gene,
                        legend.labs = c("High scores", "Low scores"),    
                        palette = c("#6E8B3D", "orangered"),
                        xlab = "Time (Months)",
                        ylab = "Recurrent free Survival")
  plot_km$plot <- plot_km$plot + 
    theme(plot.title = element_text(size = 12),                                        
          legend.position = c(0.95, 0.05),                                         
          legend.justification = c(1, 0))
  plots[[gene]] <- plot_km$plot
}

p_combined <- ggarrange(plotlist = plots, ncol = 3, nrow = 2)
print(p_combined)







# 2. Kaplan-Meier survival plot for Plasma cells C12
genes <- c("Plasma.cells_C12", "Clara.C5", "Clara.C3", 
           "Endothelial.cells.C4", "Ciliated.cells.C7", "Mac.DC_C11")
gene <-genes[1]
gene_expression <- BCMLC[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
fit <- survfit(surv_object ~ group)
plot_km <- ggsurvplot(fit, 
                      data = BCMLC,
                      pval = TRUE, pval.method = T,         
                      title = gene,
                      legend.title = gene,
                      legend.labs = c("High scores", "Low scores"),    
                      palette = c("#6E8B3D", "orangered"),
                      xlab = "Time (Months)")
plot_km$plot<-plot_km$plot+theme(plot.title = element_text(size = 12),                                        
                                 legend.position = c(0.95, 0.05),                                         
                                 legend.justification = c(1, 0))
plot_km$plot


gene <-genes[2]
gene_expression <- BCMLC[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
fit <- survfit(surv_object ~ group)
plot_km2 <- ggsurvplot(fit, 
                      data = BCMLC,
                      pval = TRUE, pval.method = T,         
                      title = gene,
                      legend.title = gene,
                      legend.labs = c("High scores", "Low scores"),  
                      palette = c("#6E8B3D", "orangered"),
                      xlab = "Time (Months)")
plot_km2$plot<-plot_km2$plot+theme(plot.title = element_text(size = 12),
                                   legend.position = c(0.95, 0.05),
                                   legend.justification = c(1, 0))
plot_km2$plot




gene <-genes[3]
gene_expression <- BCMLC[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
fit <- survfit(surv_object ~ group)
plot_km3 <- ggsurvplot(fit, 
                      data = BCMLC,
                      pval = TRUE, pval.method = T,         
                      title = gene,
                      legend.title = gene,
                      legend.labs = c("High scores", "Low scores"),
                      palette = c("#6E8B3D", "orangered"),
                      xlab = "Time (Months)")
plot_km3$plot<-plot_km3$plot+theme(plot.title = element_text(size = 12),
                                   legend.position = c(0.95, 0.05),
                                   legend.justification = c(1, 0))
plot_km3$plot


gene <-genes[4]
gene_expression <- BCMLC[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
fit <- survfit(surv_object ~ group)
plot_km4 <- ggsurvplot(fit, 
                       data = BCMLC,
                       pval = TRUE, pval.method = T,         
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),
                       palette = c("#6E8B3D", "orangered"),
                       xlab = "Time (Months)")
plot_km4$plot<-plot_km4$plot+theme(plot.title = element_text(size = 12),
                                   legend.position = c(0.95, 0.05), 
                                   legend.justification = c(1, 0))
plot_km4$plot


gene <-genes[5]
gene_expression <- BCMLC[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
fit <- survfit(surv_object ~ group)
plot_km5 <- ggsurvplot(fit, 
                      data = BCMLC,
                      pval = TRUE, pval.method = T,         
                      title = gene,
                      legend.title = gene,
                      legend.labs = c("High scores", "Low scores"),   
                      palette = c("#6E8B3D", "orangered"),
                      xlab = "Time (Months)")
plot_km5$plot<-plot_km5$plot+theme(plot.title = element_text(size = 12),
                                   legend.position = c(0.95, 0.05),
                                   legend.justification = c(1, 0))
plot_km5$plot


gene <-genes[6]
gene_expression <- BCMLC[, gene]
median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
fit <- survfit(surv_object ~ group)
plot_km6 <- ggsurvplot(fit, 
                       data = BCMLC,
                       pval = TRUE, pval.method = T,         
                       title = gene,
                       legend.title = gene,
                       legend.labs = c("High scores", "Low scores"),   
                       palette = c("#6E8B3D", "orangered"),
                       xlab = "Time (Months)")
plot_km6$plot<-plot_km6$plot+theme(plot.title = element_text(size = 12),
                                   legend.position = c(0.95, 0.05),
                                   legend.justification = c(1, 0))
plot_km6$plot

p_combined<-plot_km$plot+plot_km2$plot+plot_km3$plot+plot_km4$plot+plot_km5$plot+plot_km6$plot
print(p_combined)

pdf(paste0(output, "BCMLC_OS_all.pdf"), width = 12, height =8)
print(p_combined)
dev.off() 

#@@@@@@@@@@@@@@@ RFS
genes <- c("Plasma.cells_C12", "B.cells_C2.ES", "CD4_C10", "GaDelT_C14", "B.cells_C13", "CD8_C1.ES",
           "NK.cells_C4.ES", "CD4_C6.ES", "CD8_C15", "Mono_C3.ES", "Neutrophil_C5.ES", "DC_C18","Mac_C8",
           "Prolifirating.Mac_C16", "Mac.DC_C11","Fibroblasts.C6", "Endothelial.cells.C4", "AT2.C2", "Clara.C5",
           "AT2.C1","Clara.C3") 

plots <- list()  # Store plots
for (i in seq_along(genes)) {
  gene <- genes[i]
  gene_expression <- BCMLC[, gene]
  median_expr <- median(gene_expression, na.rm = TRUE)
  group <- ifelse(gene_expression > median_expr, "High", "Low")
  
  surv_object <- Surv(BCMLC$t.surv, BCMLC$e.surv)
  fit <- survfit(surv_object ~ group)
  
  plot_km <- ggsurvplot(fit, 
                        data = BCMLC,
                        pval = TRUE, pval.method = TRUE,         
                        title = gene,
                        legend.title = gene,
                        legend.labs = c("High scores", "Low scores"),    
                        palette = c("#6E8B3D", "orangered"),
                        xlab = "Time (Months)")
  
  plot_km$plot <- plot_km$plot + theme(plot.title = element_text(size = 12),                                        
                                       legend.position = c(0.95, 0.05),                                         
                                       legend.justification = c(1, 0))
  
  plots[[i]] <- plot_km$plot
}
p_combined <- wrap_plots(plots, ncol = 3, nrow = 2)  
print(p_combined)



surv_by_smoking <- survfit(Surv(t.surv, e.surv) ~ smoking, data = BCMLC)
summary(surv_by_smoking)
library(survminer)
p<-ggsurvplot(surv_by_smoking, data = BCMLC, pval = TRUE, pval.method = T,         
              legend.title = "Smoking",
              legend.labs = c("Non-Smokers (19)", "Smokers (55)"),
              palette = c("#6E8B3D", "orangered"), xlab = "Time (Months)", 
              title = "BCMLC Survival by Smoking Status")
p$plot<-p$plot+theme(plot.title = element_text(size = 12),
                     legend.position = c(0.95, 0.05),
                     legend.justification = c(1, 0))
p







# 3. Forest plot 1
# Success for stage-1 samples
BCMLC2 <- BCMLC[BCMLC$Pstage %in% "I", ]
p1 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 +Age +Sex+  Recurrence_Status, BCMLC2))
p1

BCMLC1 <- BCMLC[!is.na(BCMLC$Pstage),]
table(BCMLC1$Pstage)
surv_by_smoking <- survfit(Surv(t.surv, e.surv) ~ smoking, data = BCMLC1)
summary(surv_by_smoking)
library(survminer)
p<-ggsurvplot(surv_by_smoking, data = BCMLC1, pval = TRUE, 
           legend.title = "Smoking",
           legend.labs = c("Non-Smokers (17)", "Smokers (53)"),
           palette = c("#6E8B3D", "orangered"), xlab = "Time (Months)", 
           title = "BCMLC Survival by Smoking Status")
p$plot<-p$plot+theme(plot.title = element_text(size = 12),
             legend.position = c(0.95, 0.05),
             legend.justification = c(1, 0))
p





p11 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 +Age +Sex+ Recurrence_Status+Pstage, BCMLC1))
p2 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Clara.C5 +Pstage+Recurrence_Status, BCMLC1))
p3 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Clara.C3  +Age+smoking, BCMLC1))
p4 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Endothelial.cells.C4+Age +smoking, BCMLC))
p5 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Ciliated.cells.C7 +Age +smoking, BCMLC))
p6 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Mac.DC_C11+Age +smoking, BCMLC))

BCMLC_plot <- p1 + p2+ p3+p4 +p5+p6+ plot_layout(nrow = 3, ncol = 2)+
  plot_annotation(title = "BCMLC_Multivariate Analysis")
print(BCMLC_plot)

#pdf(filename = paste0(output, "BCMLC_forest.pdf"), plot = p_forest1, width = 10, height = 4)
pdf(paste0(output, "BCMLC_forest.pdf"), width = 10, height = 4)
print(p_forest1)
dev.off() 
# 4. Forest plot 2
p_forest2 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 + Age, BCMLC_F))

#pdf(filename = paste0(output, "BCMLC_forest_F.pdf"), plot = p_forest2, width = 10, height = 3)
pdf(paste0(output, "BCMLC_forest_F.pdf"), width = 10, height = 3)
print(p_forest2)
dev.off() 




#@@@@@@@@@@@@@@@@@@ Plasma cells on smokers multivariate analysis and stages. 
BCMLC_Y <- BCMLC[BCMLC$smoking %in% "Y", ]
BCMLC_Y <- BCMLC_Y[!is.na(BCMLC_Y$Pstage),]
p_1 <- forest_model(coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 +Age+ Recurrence_Status, BCMLC_Y))
p_1

library(ggpubr)
pstage_counts <- table(BCMLC_Y$Pstage)
new_labels <- paste0(names(pstage_counts), "\n(n=", pstage_counts, ")")
ps <- ggviolin(BCMLC_Y, 
               x = "Pstage",  
               y = "Plasma.cells_C12", 
               fill = "Pstage",
               palette = c("#00688B", "#8B3A62"), 
               add = "boxplot", add.params = list(fill = "white")) +  # Add boxplot inside violin
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +
  labs(title = "Plasma_cells_C12 in Smokers", 
       y = "Signature Scores", 
       x = "Pstage") +
  theme(plot.title = element_text(size = 12), 
        legend.position = "none") +
  scale_x_discrete(labels = new_labels)  # Apply new labels with sample sizes

ps





#@@@@@@@@@@@@@@ BCMLC clinical data:
rm(list=ls())
myinf2 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/Clinical_info_20231111.txt"

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
info$Race <- tolower(info$Race)
info$Race <- trimws(info$Race)
info$Race <- ifelse(info$Race %in% c("a", "asian"), "Asian",
                    ifelse(info$Race %in% c("b", "black"), "Black",
                           ifelse(info$Race %in% c("w", "white"), "White",
                                  info$Race)))

xx = tolower(info$Pstage)
xx[grep("ii", xx)] = "II"
xx[grep("i", xx)] = "I"
info$Pstage = xx
t.surv = info$OS.month_09232023 # RFS or DFS
e.surv = info$Death # Recurrence 0 or 1
info = cbind(t.surv, e.surv, info)
table(info$Race)
table(info$DateRecurrence)
# normal samples
se = grep("N", rownames(info))
info_N<-info[se,]
table(info_N$Sex)
table(info_N$Pstage)
table(info_N$Race)
table(info_N$e.surv)
table(info_N$DateRecurrence)

# Tumor samples
se = grep("T", rownames(info))
info_T<-info[se,]
table(info_T$Sex)
table(info_T$Pstage)
table(info_T$Race)
table(info_T$e.surv)
table(info_T$DateRecurrence)









#$$$$$$$$$$$$$$$$$$$$$     DEG ANALYSIS BCMLC
myinf1 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/BCM_prot_gene_TPM_Symbol.csv"
data = read.csv(myinf1, header=T, row.names=1, check.names=F)

data = log2(data+1)
dim(data)
se = grep("T", colnames(data))
data = data[,se]

myinf2 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/Clinical_info_20231111.txt"
info = read.table(myinf2, sep="\t", header=T, row.names=3, quote="")
se = grep("Adenocarcinoma|adenocarcinoma", info$Histologic.type..WHO.2015.edition)
info = info[se,]
se = which(info$Pstage!="0")
info = info[se,]
xx = tolower(info$Pstage)
xx[grep("ii", xx)] = "II"
xx[grep("i", xx)] = "I"
info$Pstage = xx
t.surv = info$OS.month_09232023
e.surv = info$Death
info$sample<-paste0("LC",info$Patient.id)
names(info)
myinf2 = "/mount/ictr1/chenglab/cc59/PriDat/Amos_Data/BCM_LC_RNAseq/Clinical_info_20231212.txt"
info1 = read.table(myinf2, sep="\t", header=T, quote="")
info1$LC
names(info1)


info$smoking<-ifelse(info$sample %in% info1$LC, info1$SmokingHx, NA)
table(info$smoking)
info$smoking<-trimws(info$smoking)





comxx = intersect(colnames(data), row.names(info))
data = data[,comxx]
info = info[comxx,]
info$DateRecurrence <- as.character(info$DateRecurrence)
info$Recurrence_Status <- ifelse(info$DateRecurrence != "", "Rec", "NonRec")
info$Recurrence_Status <- factor(info$Recurrence_Status, levels = c("NonRec", "Rec"))
table(info$Recurrence_Status)


library(limma)
library(edgeR)
library(dplyr)

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


deg_filtered <- deg_results[deg_results$adj.P.Val < 0.05, ]




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





#@@@@@@@@@@@@@@@@@@@@@@@@@ Rousseaux_GSE30219 (n=85)
file_path1 <- "/home/u251079/scRNA/Lung_scRNA/Lung_output/Rousseaux_GSE30219_AmoscRNAMarkers_Base.txt"
data <- read.table(file_path1, sep="\t", header=T, row.names=1)
head(data)
head(data)
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}



# Clinical info
file_path2 <- "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Rousseaux_GSE30219/Sample_info.txt"
info = read.table(file_path2, sep="\t", header=T, row.names=1, quote="")
data1<-cbind(data, info)
names(data1)
str(data1)
library(dplyr)

data1 <- data1 %>%
  rename(relapse_event = relapse..event.1..no.event.0.)
names(data1)

data_filtered <- subset(data1, !is.na(relapse_event))
data_filtered <-data_filtered [data_filtered$histology %in% " ADC",]
str(data_filtered)
p <- ggplot(data_filtered, aes(x = relapse_event, y = Plasma.cells_C12.ES, fill = relapse_event)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p") + 
  scale_fill_manual(values = c("0" = "lightblue", "1" = "darkorange")) +  
  theme_minimal() +
  labs(x = "Recurrence Status", y = "Plasma Cells_C12 Signature Score", fill = "Recurrence") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")

print(p)

library(ggplot2)

# Convert pt.stage to a factor with the correct order
data_filtered$pt.stage <- factor(data_filtered$pt.stage, levels = c(" T1", " T2", " T3"))

# Create the box plot
ggplot(data_filtered, aes(x = pt.stage, y = Plasma.cells_C12.ES, fill = pt.stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot without outlier points
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add jittered points
  scale_fill_manual(values = c("#6E8B3D", "orangered", "lightblue")) +  # Custom colors
  theme_minimal() +
  labs(x = "Pathological Stage", y = "Plasma Score (C12)", title = "Plasma Score C12 by Pathological Stage") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))





#@@@@@@@@@@@@@@@@@@@@@@@@@ Der_data_GSE50081 
this data has 128 patients with Lung adenocarcinoma samples with relapse and RFS information with gene expression data
!Sample_characteristics_ch1 histology: adenocarcinoma
# /mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Der_GSE50081
file_path1 <-"/home/u251079/scRNA/Lung_scRNA/Lung_output/Der_GSE50081_AmoscRNAMarkers_Base.txt"
data <- read.table(file_path1, sep="\t", header=T, row.names=1)
head(data)
nn = ncol(data)/2
data = data[, 1:nn]
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}


# Clinical info
file_path2 <- "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Der_GSE50081/Sample_info.txt"
info = read.table(file_path2, sep="\t", header=T, row.names=1, quote="")
xx<-intersect(rownames(data), rownames(info))

data1<-cbind(data, info)
names(data1)
str(data1)

data1$recurrence <- factor(str_trim(data1$recurrence), levels = c("Y", "N"))
data_filtered <- subset(data1, !is.na(recurrence))

p <- ggplot(data_filtered, aes(x = recurrence, y = Plasma.cells_C12.ES, fill = recurrence)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p") + 
  scale_fill_manual(values = c("N" = "lightblue", "Y" = "darkorange")) +  
  theme_minimal() +
  labs(x = "Recurrence Status", y = "Plasma Cells_C12 Signature Score", fill = "Recurrence") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")

print(p)











#@@@@@@@@@@@@@@@@@@@@@@@@@@@
# OVERALL SURVIVAL
#@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
#@@@@@@@@@@@@@@@@@@@@@@@@@@@ Tomida data  
library("survival")
library("dplyr")
info_out <- "/home/u251079/scRNA/Lung_scRNA/Tomida.txt"
output<-"/home/u251079/scRNA/Lung_scRNA/"
Tomida <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(Tomida)
# Tomida <- Tomida[Tomida$days.before.death.censor <= 1830, ]
Tomida$gender <- factor(Tomida$gender)
Tomida$smoking.status <- factor(Tomida$smoking.status)
Tomida$pathological.stage <- factor(Tomida$pathological.stage)

Tomida$T_stage <- trimws(Tomida$pathological.stage)
Tomida$T_stage <- factor(ifelse(Tomida$T_stage %in% c("IA", "IB"), "I", Tomida$T_stage), levels = c("I", "II"))
med_age<-median(Tomida$age)
Tomida$Age<-ifelse(Tomida$age <=med_age, "age≤61", "age≥61")
Tomida$Plasma.cells_C12 <- gsub("\\.ES$", "", Tomida$Plasma.cells_C12.ES)
Tomida$Plasma.cells_C12<-as.numeric(Tomida$Plasma.cells_C12)

# fit.coxph_Tob <- coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12+gender +Age +T_stage, data = Tomida)
# summary(fit.coxph_Tob)
# p1_Tob <- ggforest(fit.coxph_Tob, data = Tomida, main = "Tomida Multivariate (5 Years)")
# print(p1_Tob)
p1<-forest_model(coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12+gender +Age +T_stage, data = Tomida))

# KM plot
gene <- "Plasma.cells_C12"
gene_expression <- as.numeric(Tomida[[gene]])
median_expr <- median(gene_expression, na.rm = TRUE)
Tomida$group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(as.numeric(Tomida$t.surv), as.numeric(Tomida$e.surv))
fit <- survfit(surv_object ~ group, data = Tomida)
cox_model <- coxph(surv_object ~ group, data = Tomida)
cox_pvalue <- summary(cox_model)$sctest[3] 
p2 <- ggsurvplot(fit, 
                 data = Tomida,
                 pval = TRUE,             
                 pval.method = TRUE,     
                 conf.int = FALSE,        
                 risk.table = FALSE,     
                 ggtheme = theme_classic(),
                 palette = c("darkblue", "darkred"),
                 title = NULL,
                 legend.title = "Signature scores",
                 legend.labs = c("Low", "High"),
                 xlab = "Time (Days)",
                 ylab = "Survival Probability")


p2$plot<-p2$plot+theme(axis.text.x = element_text(color = "black"),  
                       axis.text.y = element_text(color = "black"),  
                       axis.title.x = element_text(color = "black"),  
                       axis.title.y = element_text(color = "black"),
                       legend.position = c(0.8, 0.25),  
                       egend.justification = c("left", "center"),
                       axis.title = element_text(size = 14, face = "bold"))
print(p2$plot)

# Boxplot for Plasma Cells C12 Expression by Gender
p3 <- ggplot(Tomida, aes(x = gender, y = Plasma.cells_C12, fill = gender)) +
  geom_boxplot() +
  theme_classic() +  # Or use theme_classic2() if you prefer
  labs(title = NULL, x = "Gender", y = "Signature scores") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  scale_fill_manual(values = c("#FF9999",  "lightblue")) +  
  theme(plot.title = element_text(size = 12),  # Set title size to 12
    legend.position = "none",  # Hide legend
    axis.title.x = element_text(color = "black"),  # Black color for x-axis title
    axis.title.y = element_text(color = "black"),
          axis.text.x = element_text(color="black",size = 12),
          axis.text.y = element_text(color="black",size = 12),
          axis.title = element_text(size = 14, face = "bold"))

print(p3)


#@@@@@@@@@@@@@@ T stages 
Tomida_Tstage <- Tomida[!is.na(Tomida$T_stage), ]

# Create the plot
p4 <- ggplot(Tomida_Tstage, aes(x = T_stage, y = Plasma.cells_C12, fill = T_stage)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # Violin plot with transparency
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9) +  # Add boxplot on top
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +  # Center the p-value label
  scale_fill_manual(values = c("#66B3FF",  "#CC99FF")) +  # Custom colors
  theme_classic() +
  labs(title = NULL,
       x = "Tumor Stage", 
       y = "Signature scores")+
  theme(legend.position = "none",  # Remove legend if unnecessary
        axis.text.x = element_text(color="black",size = 12),
        axis.text.y = element_text(color="black",size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  

print(p4)


# Ensure T_stage is a factor with proper levels
Tomida_Tstage$T_stage <- factor(Tomida_Tstage$T_stage, levels = c("I", "II"))
Tomida_Tstage1 <- Tomida_Tstage[Tomida_Tstage$T_stage %in% "I",]
Tomida_Tstage1$Plasma_cells_Group <- ifelse(Tomida_Tstage1$Plasma.cells_C12 > median(Tomida_Tstage1$Plasma.cells_C12, na.rm = TRUE), 
                                            "High", "Low")
surv_object <- Surv(Tomida_Tstage1$t.surv, Tomida_Tstage1$e.surv)
km_fit <- survfit(surv_object ~ Plasma_cells_Group, data = Tomida_Tstage1)
p5 <- ggsurvplot(km_fit, 
                 data = Tomida_Tstage1,
                 pval = TRUE,  
                 pval.method = TRUE,  
                 conf.int = FALSE,
                 risk.table = FALSE,  
                 ggtheme = theme_classic(),  
                 palette = c("darkblue", "darkred"),  
                 title = NULL,
                 legend.title = "Signature scores",
                 legend.labs = c("Low Score", "High Score"),
                 xlab = "Time (Days)",
                 ylab = "Survival Probability")

p5$plot <- p5$plot + theme(axis.text.x = element_text(color = "black"),  
                           axis.text.y = element_text(color = "black"),  
                           axis.title.x = element_text(color = "black"),  
                           axis.title.y = element_text(color = "black"),
                           legend.position = c(0.8, 0.25),  
                           egend.justification = c("left", "center"),
                           axis.title = element_text(size = 14, face = "bold"))
print(p5)
library(gridExtra)
p2_plot<-p2$plot
p5_plot <- p5$plot 
(p2_plot| p1) / (p3 | p4 | p5_plot) 

pdf(paste0(output, "Tomida_Plots_Arranged.pdf"), width = 20, height = 8)  
print((p2_plot | p1) / (p3 | p4 | p5_plot))  
dev.off()

pdf(paste0(output, "Tomida_A.pdf"), width = 5, height = 4)
print(p2$plot)
dev.off()

pdf(paste0(output,  "Tomida_B.pdf"), width = 10, height = 3.5)
print(p1)
dev.off()

pdf(paste0(output, "Tomida_C.pdf"), width = 4, height = 3)
print(p3)
dev.off()

pdf(paste0(output, "Tomida_D.pdf"), width = 4, height = 3.5)
print(p4)
dev.off()

pdf(paste0(output,  "Tomida_E.pdf"), width = 5, height = 4)
print(p5$plot)
dev.off()









#@@@@@@@@@@@@@@@@@@@@ Shedden data
# for male did not significant anything
# for females overall survival is significant KM plot
rm(list = ls())
output<-"/home/u251079/scRNA/Lung_scRNA/"
info_out<-"/home/u251079/scRNA/Lung_scRNA/Shedden.txt"
Shedden <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
Shedden$smoking_history<-as.character(Shedden$smoking_history)
Shedden$T_stage <- sub(".*(T[1-4]).*", "\\1", Shedden$disease_stage)
Shedden$T_stage <- factor(Shedden$T_stage, levels = c("T1", "T2", "T3", "T4"))
table(Shedden$T_stage)
Shedden$smoking_history <- trimws(Shedden$smoking_history)
Shedden$smoking_history <- factor(Shedden$smoking_history)
table(Shedden$smoking_history)
table(Shedden$Sex)
Shedden$Plasma.cells_C12 <- gsub("\\.ES$", "", Shedden$Plasma.cells_C12.ES)
Shedden$Plasma.cells_C12<-as.numeric(Shedden$Plasma.cells_C12)

# Kaplan-Meier Survival Analysis (Plasma Cells C12 Expression)
#@@@@@@@@@@@@@@@@@@@@ 5 ;years 
Shedden<-Shedden[Shedden$t.surv <=62,]
Shedden$smoking_history<-as.character(Shedden$smoking_history)
Shedden$T_stage <- sub(".*(T[1-4]).*", "\\1", Shedden$disease_stage)
Shedden$T_stage <- factor(Shedden$T_stage, levels = c("T1", "T2", "T3", "T4"))
table(Shedden$T_stage)
Shedden$Sex <- factor(Shedden$Sex)
Shedden$smoking_history <- factor(Shedden$smoking_history)
table(Shedden$smoking_history)
summary(Shedden$Plasma.cells_C12)

# Kaplan-Meier Survival Analysis (Plasma Cells C12 Expression)
gene <- "Plasma.cells_C12"
gene_expression <- as.numeric(Shedden[[gene]])
median_expr <- median(gene_expression, na.rm = TRUE)
Shedden$group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(as.numeric(Shedden$t.surv), as.numeric(Shedden$e.surv))
fit <- survfit(surv_object ~ group, data = Shedden)
cox_model <- coxph(surv_object ~ group, data = Shedden)
cox_pvalue <- summary(cox_model)$sctest[3] 
p1 <- ggsurvplot(fit, 
                 data = Shedden,
                 pval = TRUE,             
                 pval.method = TRUE,     
                 conf.int = FALSE,        
                 risk.table = FALSE,     
                 ggtheme = theme_classic(),
                 palette = c("darkblue", "darkred"),
                 title = NULL,
                 legend.title = "Signature scores",
                 legend.labs = c("Low", "High"),
                 xlab = "Time (months)",
                 ylab = "Survival Probability")


p1$plot<-p1$plot+theme(axis.text.x = element_text(color = "black"),  
                        axis.text.y = element_text(color = "black"),  
                        axis.title.x = element_text(color = "black"),  
                        axis.title.y = element_text(color = "black"),
                       legend.position = c(0.7, 0.8), 
                       legend.justification = c("left", "center"))

p1
pdf(file.path(output, "Shedden_KM.pdf"), width = 4.5, height = 3)
print(p1$plot)
dev.off()

# Tumor Stage 
Shedden_Tstage<-Shedden[!Shedden$T_stage %in% NA,]
p2 <- ggplot(Shedden_Tstage, aes(x = T_stage, y = Plasma.cells_C12, fill = T_stage)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("darkblue", "darkred", "darkgreen", "darkorange")) +
  theme_classic() +
  scale_y_continuous(limits = c(NA, max(Shedden_Tstage$Plasma.cells_C12, na.rm = TRUE) + 2)) +  # Corrected scale
  stat_compare_means(method = "anova", label.x = 0.8, label.y = max(Shedden_Tstage$Plasma.cells_C12) + 1) +  
  stat_compare_means(method = "wilcox.test", ref.group = "T1", label = "p", label.y = max(Shedden_Tstage$Plasma.cells_C12) +0.5) +  # Fixed ref.group
  labs(title = NULL, x = "Stage", y = "Signature scores") +
  theme(plot.title = element_text(size = 10),  
        legend.position = "none",axis.text.x = element_text(color = "black"),  
        axis.text.y = element_text(color = "black"),  
        axis.title.x = element_text(color = "black"),  
        axis.title.y = element_text(color = "black"))


pdf(file.path(output, "Shedden_Tstage.pdf"), width = 4.5, height = 3)
print(p2)
dev.off()

Shedden_Tstage$T_stage_group <- factor(Shedden_Tstage$T_stage, levels = c("T1", "T2", "T3", "T4"))
surv_object <- Surv(as.numeric(Shedden_Tstage$t.surv), as.numeric(Shedden_Tstage$e.surv))
fit <- survfit(surv_object ~ T_stage_group, data = Shedden_Tstage)
cox_model <- coxph(surv_object ~ T_stage_group, data = Shedden_Tstage)
cox_pvalue <- summary(cox_model)$sctest[3] 
p3 <- ggsurvplot(fit, 
                 data = Shedden_Tstage,
                 pval = TRUE,             
                 pval.method = TRUE,     
                 conf.int = FALSE,        
                 risk.table = FALSE,     
                 ggtheme = theme_classic(),
                 palette = c("darkblue", "darkred", "darkgreen", "darkorange"),
                 title = NULL,
                 legend.title = "T Stage",
                 legend.labs = c("T1", "T2", "T3", "T4"))

p3$plot<-p3$plot+theme(axis.text.x = element_text(color = "black"),  
                       axis.text.y = element_text(color = "black"),  
                       axis.title.x = element_text(color = "black"),  
                       axis.title.y = element_text(color = "black"),
                       legend.position = c(0.85, 0.75), 
                       legend.justification = c("left", "center"))
p3


pdf(file.path(output, "Shedden_KM_T_stages.pdf"), width = 4.5, height = 3)
print(p3$plot)
dev.off()


#@@@@@@@@@@@@@@@  Females
Shedden$Sex <- factor(trimws(Shedden$Sex))
str(Shedden$Sex)
Shedden_female <- Shedden[Shedden$Sex == "Female", ]

# Kaplan-Meier Survival Analysis (Plasma Cells C12 Expression)
gene <- "Plasma.cells_C12"
gene_expression <- as.numeric(Shedden_female[[gene]])
median_expr <- median(gene_expression, na.rm = TRUE)
Shedden_female$group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(as.numeric(Shedden_female$t.surv), as.numeric(Shedden_female$e.surv))
fit <- survfit(surv_object ~ group, data = Shedden_female)
cox_model <- coxph(surv_object ~ group, data = Shedden_female)
cox_pvalue <- summary(cox_model)$sctest[3] 
p4 <- ggsurvplot(fit, 
                 data = Shedden_female,
                 pval = TRUE,             
                 pval.method = TRUE,     
                 conf.int = FALSE,        
                 risk.table = FALSE,     
                 ggtheme = theme_classic(),
                 palette = c("darkblue", "darkred"),
                 title = NULL,
                 legend.title = "Signature scores",
                 legend.labs = c("Low", "High"),
                 xlab = "Time (months)",
                 ylab = "Survival Probability")


p4$plot <- p4$plot + theme(axis.text.x = element_text(color = "black"),  
                           axis.text.y = element_text(color = "black"),  
                           axis.title.x = element_text(color = "black"),  
                           axis.title.y = element_text(color = "black"),
                           legend.position = c(0.75, 0.80),  
                           legend.justification = c("left", "center"))  # Corrected


p4

pdf(file.path(output, "Shedden_female_KM.pdf"), width = 4.5, height = 3)
print(p4$plot)
dev.off()

#@@@@@@@@@@@@@@@  Males
# for males KM plot not significant,
# stages not significant 
# smoking history not significant

#@@@@@@@@@@@@@@@ Smoking History
# KM is not significant
# Box Wilcox is not significant between current vs never  and current vs never
# all combination with smoking is failed. 
# Even KM plots, Box plots













#@@@@@@@@@@@@@@@@@@@@@@@@@@@ KM plot not significant.   even for 5 Years of survival 
# BUT 5 years patients EGFR wt vs Het is highly significant overall survival. 
rm(list = ls())
info_out<-"/home/u251079/scRNA/Lung_scRNA/Takeuchi.txt"
output<-"/home/u251079/scRNA/Lung_scRNA/"
Takeuchi <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(Takeuchi)
table(Takeuchi$pStage)


Takeuchi$Stage <- factor(
  ifelse(Takeuchi$pStage %in% c("IA", "IB"), "I",
         ifelse(Takeuchi$pStage %in% c("IIA", "IIB"), "II", "III")),
  levels = c("I", "II", "III")
)

table(Takeuchi$Stage)
Takeuchi$Plasma.cells_C12 <- gsub("\\.ES$", "", Takeuchi$Plasma.cells_C12.ES)
Takeuchi$Plasma.cells_C12<-as.numeric(Takeuchi$Plasma.cells_C12)
#Takeuchi<-Takeuchi[Takeuchi$t.surv <= 1830,]



gene_expression <- "Plasma.cells_C12" 
# Boxplot for Plasma.cells_C12 expression by EGFR Status
p_egfr_box <- ggplot(Takeuchi, aes(x = EGFR.status, y = Takeuchi[[gene_expression]], fill = EGFR.status)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = paste(gene_expression, "Expression by EGFR Status"),
       x = "EGFR Status",
       y = "Expression Level") +
  scale_fill_manual(values = c("Mut" = "blue", "Wt" = "red")) +
  stat_compare_means(method = "wilcox.test", label = "p") 

print(p_egfr_box)

# KM PLOT
gene <- "Plasma.cells_C12"
gene_expression <- as.numeric(Takeuchi[[gene]])
median_expr <- median(gene_expression, na.rm = TRUE)
Takeuchi$group <- ifelse(gene_expression > median_expr, "High", "Low")
surv_object <- Surv(as.numeric(Takeuchi$t.surv), as.numeric(Takeuchi$e.surv))
fit <- survfit(surv_object ~ group, data = Takeuchi)
cox_model <- coxph(surv_object ~ group, data = Takeuchi)
cox_pvalue <- summary(cox_model)$sctest[3] 
p_km <- ggsurvplot(fit, 
                   data = Takeuchi,
                   pval = TRUE,
                   pval.method = TRUE,
                   conf.int = FALSE,  
                   risk.table = FALSE, 
                   ggtheme = theme_classic(),
         #          palette = c("darkblue", "darkred", "darkgreen"),  # Adjusted to 3 colors
         #          title = "KM Survival Curve by Stage in Takeuchi",
          #         legend.title = "Stage",  # Adjusted to 'Stage'
          #        legend.labs = c("I", "II", "III"),  # Adjusted to 3 stages
                   xlab = "Time (Days)",
                   ylab = "Survival Probability")

print(p_km)




# Boxplot for Plasma.cells_C12 expression by K.ras.Status
p_kras_box <- ggplot(Takeuchi, aes(x = K.ras.Status, y = Takeuchi[[gene_expression]], fill = K.ras.Status)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = paste(gene_expression, "Expression by K.ras Status"),
       x = "K.ras Status",
       y = "Expression Level") +
  scale_fill_manual(values = c("Mut" = "green", "Wt" = "orange")) +
  stat_compare_means(method = "wilcox.test", label = "p")
print(p_kras_box)

# Boxplot for Plasma.cells_C12 expression by p53.Status
p_p53_box <- ggplot(Takeuchi, aes(x = p53.Status, y = Takeuchi[[gene_expression]], fill = p53.Status)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = paste(gene_expression, "Expression by p53 Status"),
       x = "p53 Status",
       y = "Expression Level") +
  scale_fill_manual(values = c("Mut" = "purple", "Wt" = "yellow")) +
  stat_compare_means(method = "wilcox.test", label = "p")
print(p_p53_box)



# Boxplot for Plasma.cells_C12 expression by Sex
p_sex_box <- ggplot(Takeuchi, aes(x = Sex, y = Takeuchi[[gene_expression]], fill = Sex)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = paste(gene_expression, "Expression by Sex"),
       x = "Sex",
       y = "Expression Level") +
  scale_fill_manual(values = c("M" = "blue", "F" = "pink")) +
  stat_compare_means(method = "wilcox.test", label = "p")

# Print Sex boxplot
print(p_sex_box)











#@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
info_out<-"/home/u251079/scRNA/Lung_scRNA/Tang.txt"
Tang <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(Tang)

table(Tang$stage)
Tang$Stage <- factor(
  ifelse(Tang$stage %in% c("IA", "IB"), "I",
         ifelse(Tang$stage %in% c("IIA", "IIB"), "II", "III")),
  levels = c("I", "II", "III")
)
table(Tang$histology)
Tang<-Tang[Tang$histology %in% "Adenocarcinoma",]
Tang$Plasma.cells_C12 <- gsub("\\.ES$", "", Tang$Plasma.cells_C12.ES)
Tang$Plasma.cells_C12<-as.numeric(Tang$Plasma.cells_C12)


# Categorize Plasma.cells_C12 into low and high expression groups
threshold <- median(Tang$Plasma.cells_C12, na.rm = TRUE)  # Use median as threshold
Tang$Plasma_cells_C12_group <- ifelse(Tang$Plasma.cells_C12 > threshold, "High", "Low")
Tang$Plasma_cells_C12_group <- factor(Tang$Plasma_cells_C12_group, levels = c("Low", "High"))
table(Tang$Plasma_cells_C12_group)
surv_object_tang <- Surv(Tang$t.surv, Tang$e.surv)
fit_plasma_c12 <- survfit(surv_object_tang ~ Plasma_cells_C12_group, data = Tang)
p_km_plasma_c12 <- ggsurvplot(fit_plasma_c12, 
                              data = Tang,
                              pval = TRUE,                 # Display p-value
                              pval.method = TRUE,          # Method for p-value
                              conf.int = FALSE,            # No confidence interval
                              risk.table = FALSE,          # Hide risk table
                              ggtheme = theme_classic(),   # Use classic theme
                              palette = c("blue", "red"),  # Color for Low and High expression
                              title = "Kaplan-Meier Survival Curve for Plasma.cells_C12 Expression",
                              legend.title = "Plasma.cells_C12",
                              legend.labs = c("Low Expression", "High Expression"),
                              xlab = "Time (Days)",
                              ylab = "Survival Probability")

print(p_km_plasma_c12)


# Boxplot for Plasma.cells_C12 expression by Gender
p_gender_box <- ggplot(Tang, aes(x = gender, y = Plasma.cells_C12, fill = gender)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Plasma.cells_C12 Expression by Gender",
       x = "Gender",
       y = "Expression Level") +
  scale_fill_manual(values = c("M" = "blue", "F" = "pink")) +
  stat_compare_means(method = "wilcox.test", label = "p")

# Print the Gender boxplot
print(p_gender_box)



# Boxplot for Plasma.cells_C12 expression by Stage
p_stage_box <- ggplot(Tang, aes(x = Stage, y = Plasma.cells_C12, fill = Stage)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Plasma.cells_C12 Expression by Stage",
       x = "Stage",
       y = "Expression Level") +
  scale_fill_manual(values = c("I" = "blue", "II" = "red", "III" = "green")) +
  stat_compare_means(method = "wilcox.test", label = "p")
print(p_stage_box)




#@@@@@@@@@@@@@@@@@@@@@@@@@@@
info_out<-"/home/u251079/scRNA/Lung_scRNA/Gentles.txt"
Gentles<- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(Gentles)
table(Gentles$Gender)
table(Gentles$Histology_progno)
Gentles$Plasma.cells_C12 <- gsub("\\.ES$", "",Gentles$Plasma.cells_C12.ES)
Gentles$Plasma.cells_C12<-as.numeric(Gentles$Plasma.cells_C12)
Gentles<-Gentles[Gentles$t.surv <= 60,]

# Create High/Low groups based on the median expression level
gene_expression <- "Plasma.cells_C12"
gene_expression_values <- Gentles[[gene_expression]]
median_expr <- median(gene_expression_values, na.rm = TRUE)
Gentles$gene_group <- ifelse(gene_expression_values > median_expr, "High", "Low")
surv_object <- Surv(Gentles$t.surv, Gentles$e.surv)
fit_gene <- survfit(surv_object ~ gene_group, data = Gentles)
p_gene_km <- ggsurvplot(fit_gene, 
                        data = Gentles,
                        pval = TRUE,
                        pval.method = TRUE,
                        conf.int = FALSE,  
                        risk.table = FALSE, 
                        ggtheme = theme_classic(),
                        palette = c("blue", "red"),  # Adjust colors for High and Low
                        title = paste("Kaplan-Meier Survival Curve for", gene_expression),
                        legend.title = "Gene Expression",
                        legend.labs = c("Low", "High"),
                        xlab = "Time (Days)",
                        ylab = "Survival Probability")

print(p_gene_km)






















# RUFF Non-significant 
#@@@@@@@@@@@@@@@@ Sato data
rm(list=ls())
info_out<-"/home/u251079/scRNA/Lung_scRNA/sato.txt"
sato <- read.table(info_out, header = TRUE, sep = "\t", row.names = 1)
str(sato)
names(sato)[1:28]<-gsub("\\.ES$", "", names(sato)[1:28])
sato$recurrence <- trimws(sato$recurrence)  # Remove leading/trailing spaces
sato$recurrence <- factor(sato$recurrence, levels = c("N", "Y"),
                                 labels = c("Non-Recurrence", "Recurrence"))

library(ggpubr)

p_sato <- ggboxplot(sato, 
                    x = "recurrence", 
                    y = "Plasma.cells_C12", 
                    fill = "recurrence", 
                    palette = c("blue", "orange"), 
                    add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +
  labs(title = "Plasma Cells in Sato: Non-Rec vs Rec", 
       y = "Signature Scores", 
       x = "Recurrence Status")

print(p_sato)


# sato=sato[sato$t.surv <= 1825,]
sato$t.surv
table(sato$final.patient.stage)
sato$final.patient.stage <- trimws(sato$final.patient.stage)
sato <- sato %>%
  mutate(Tstage = case_when(
    final.patient.stage %in% c("IA", "IB") ~ "I",
    final.patient.stage %in% c("IIA", "IIB") ~ "II",
    final.patient.stage %in% c("IIIA", "IIIB", "IA Vs IIIB") ~ "III",
    final.patient.stage == "IV" ~ "IV",
    TRUE ~ final.patient.stage  # Ensure any unrecognized values remain unchanged
  ))

table(sato$Tstage)

sato$Sex <- trimws(sato$gender)
sato <- sato[sato$Sex != "", ]
table(sato$Sex)
sato$Sex<- factor(sato$Sex, levels = c("M", "F"))
table(sato$Sex)



sato$race <- trimws(sato$race)
sato$race <- ifelse(sato$race %in% c("A", "ASIAN"), "Asian",
                    ifelse(sato$race %in% c("B", "BLACK"), "Black",
                           ifelse(sato$race %in% c("W", "WHITE"), "White", sato$race)))

table(sato$race)


str(sato)
Age_med<-median(sato$age)
sato$Age<-ifelse(sato$age <= Age_med, "age≤63.76", "age≥63.76")
sato$Age <- factor(sato$Age, levels = c("age≤63.76", "age≥63.76"))




gene <- "Plasma.cells_C12"
gene_expression <- sato[, gene]

# Calculate Z-scores
sato[, gene] <- (gene_expression - mean(gene_expression, na.rm=TRUE)) / sd(gene_expression, na.rm=TRUE)

median_expr <- median(gene_expression, na.rm = TRUE)
group <- ifelse(gene_expression > median_expr, "High", "Low")
fit <- survfit(Surv(t.surv, e.surv) ~ group, data = sato)
km_plot <- ggsurvplot(fit, data = sato, pval = TRUE, 
                      title = paste("RFS KM plot of Plasma cells in Sato data"),
                      legend.labs = c("High Expression", "Low Expression"),
                      xlab="Time in months",
                      palette = c("darkblue", "darkred"))
print(km_plot)



# Boxplot for Plasma Cells C12 by Gender
p1 <- ggboxplot(sato, x = "Sex", y = "Plasma.cells_C12", fill = "Sex",
                palette = c("#66CCEE", "#EE6677"), 
                add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +
  labs(title = "Plasma Cells C12 by Gender", y = "Signature Scores", x = "Gender")
p1

# Boxplot for Plasma Cells C12 by Age
p2 <- ggboxplot(sato, x = "Age", y = "Plasma.cells_C12", fill = "Age",
                palette = c("#66CCEE", "#EE6677"), 
                add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +
  labs(title = "Plasma Cells C12 by Age", y = "Signature Scores", x = "Age Group")
p2

# Boxplot for Plasma Cells C12 by Tstage
table(sato$Tstage)
sato$Tstage<- factor(sato$Tstage, levels = c("I","II","III" ,"IV"))
my_comparisons <- list( c("I","II"), c("II","III"), c("III" ,"IV") , c("I" ,"IV"))
p4 <- ggboxplot(sato, x = "Tstage", y = "Plasma.cells_C12", fill = "Tstage",
                palette = c("#66CCEE", "#EE6677", "#FFB94C", "#8B3A62"), 
                add = "jitter") +
  stat_compare_means(comparisons = my_comparisons)+
  labs(title = "Plasma Cells C12 by Tstage", y = "Signature Scores", x = "Tstage")
p4

# Boxplot for Plasma Cells C12 by Tobacco History
sato1<-sato[!is.na(sato$tobacco.history),]
p5 <- ggboxplot(sato1, x = "tobacco.history", y = "Plasma.cells_C12", fill = "tobacco.history",
                palette = c("#66CCEE", "#EE6677"), 
                add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5) +
  labs(title = "Plasma Cells C12 by Tobacco History", y = "Signature Scores", x = "Tobacco History")
p5


fit.coxph <- coxph(Surv(t.surv, e.surv) ~Plasma.cells_C12 +Age+ Sex+tobacco.history +Tstage+race, data=sato)
p1 <- ggforest(fit.coxph, data=sato, main = "sato Multivariate (5 Years)")
print(p1)


sato_M <- sato[sato$Sex == "M", ]
fit.coxph_M <- coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 + Age+tobacco.history +Tstage, data = sato_M)
summary(fit.coxph_M)
p1_M <- ggforest(fit.coxph_M, data = sato_M, main = "Male sato Multivariate (5 Years)")
print(p1_M)


sato_F <- sato[sato$Sex == "F", ]
fit.coxph_F <- coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12 +Age+tobacco.history +Tstage, data = sato_F)
summary(fit.coxph_F)
p1_F <- ggforest(fit.coxph_F, data = sato_F, main = "Female sato Multivariate (5 Years)")
print(p1_F)


sato$tobacco.history <- trimws(sato$tobacco.history)
sato_Tob <- sato[sato$tobacco.history == "Y", ]
sato_Tob <- sato_Tob[!sato_Tob$Tstage == c("IV"), ]

fit.coxph_Tob <- coxph(Surv(t.surv, e.surv) ~ Plasma.cells_C12  +Tstage, data = sato_Tob)
summary(fit.coxph_Tob)
p1_Tob <- ggforest(fit.coxph_Tob, data = sato_Tob, main = "Tobako sato Multivariate (5 Years)")
print(p1_Tob)











GSE72094 (Schabath_GSE72094): This dataset involves the study of genomic instability-related genes in lung adenocarcinoma, which may be pertinent to recurrence analysis. 
GSE4115 (Spira_GSE4115): This dataset focuses on gene expression profiles in lung cancer, potentially including samples with varying recurrence statuses. 
GSE30219 (Rousseaux_GSE30219): This dataset includes gene expression profiles from lung cancer samples, potentially relevant to recurrence studies. 
GSE41271 (Sato_GSE41271): This dataset contains gene expression data from lung adenocarcinoma samples, which may include both recurrent and non-recurrent cases. 
GSE37745 (Botling_GSE37745): This dataset provides gene expression data from lung cancer samples, which may encompass both recurrent and non-recurrent cases. 
GSE4115 (Spira_GSE4115): This dataset focuses on gene expression profiles in lung cancer, potentially including samples with varying recurrence statuses. 


#GSE157011 (Bueno_GSE157011): This dataset contains gene expression data from early-stage squamous lung cancer and may also include samples with recurrence and non-recurrence status. However, it’s primarily focused on early-stage cancers, so recurrence details should be verified in the specific study context.
GSE72094 (Schabath_GSE72094): In this study of lung adenocarcinoma, recurrence status might be inferred based on genomic instability and associated factors, but checking the publication details would be useful for confirming recurrence information.
GSE41271 (Sato_GSE41271): This dataset contains expression data from lung adenocarcinoma, including samples with different clinical outcomes, such as recurrence.


https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.794293/full
https://ascopubs.org/doi/full/10.1200/CCI.22.00062
https://pubmed.ncbi.nlm.nih.gov/38025809/
https://academic.oup.com/jnci/article/107/6/djv059/870811
https://www.oncotarget.com/article/19161/text/
https://tlcr.amegroups.org/article/view/80004/pdf
https://academic.oup.com/bib/article/25/3/bbae080/7638265#447211070
https://www.sciencedirect.com/science/article/pii/S0169500223009510
https://pmc.ncbi.nlm.nih.gov/articles/PMC10654435/
https://www.nature.com/articles/s41586-022-05672-3
https://www.researchgate.net/publication/379024175_A_new_model_using_deep_learning_to_predict_recurrence_after_surgical_resection_of_lung_adenocarcinoma
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0300442
https://cardiothoracicsurgery.biomedcentral.com/articles/10.1186/s13019-021-01476-0
https://www.sciencedirect.com/science/article/pii/S2001037022001106
https://www.cell.com/cell-reports/fulltext/S2211-1247(22)00841-5
https://www.tandfonline.com/doi/full/10.2217/fon-2023-0024
https://www.jtcvs.org/article/S0022-5223(24)00440-9/abstract
https://www.cancerbiomed.org/content/18/3/734
https://www.nature.com/articles/s41598-023-42851-2
https://www.nature.com/articles/s41598-021-02946-0
https://www.nature.com/articles/s41598-024-56867-9
https://www.jtcvs.org/article/S0022-5223(24)00440-9/abstract




# Cancer statistics
https://acsjournals.onlinelibrary.wiley.com/doi/10.3322/caac.21654






https://www.pnas.org/doi/full/10.1073/pnas.1902651116






rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/Lung/ImmGen/data/Botling_GSE37745_AmoscRNAMarkers_BaseRat.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Cancer/Lung/Botling_GSE37745/Sample_info.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = which(info$histology ==" adeno")
info = info[se,]
t.surv = as.numeric(info$days.to.recurrence...to.last.visit)
e.surv = ifelse(info$recurrence ==" yes", 1, 0)
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]

library(survival)

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]

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
xx[1:20,]
se = grep("_VS_", xx[,1])
xx[-se,]

volcano_df <- xx[-se, ]

# Add log2HR and -log10pval
volcano_df <- volcano_df %>%
  mutate(log2HR = log2(hr1),
         negLog10P = -log10(coxph.pval1),
         Sig = case_when(
           coxph.pval1 < 0.01 & hr1 < 1 ~ "Protective",
           coxph.pval1 < 0.01 & hr1 > 1 ~ "Risk",
           TRUE ~ "NS"
         ))

# Volcano plot
library(ggrepel)
plot_1 <- ggplot(volcano_df, aes(x = log2HR, y = negLog10P, label = name)) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = c("Protective" = "#1E90FF", "Risk" = "red", "NS" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = subset(volcano_df, coxph.pval1 < 0.01),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  xlim(-1.2, 0.5) +
  ylim(0, 4) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.9)) +
  labs(x = "log2(Hazard Ratio)", y = "-log10(p-value)", color = "Type") +
  ggtitle("Botling_GSE37745 data RFS")


plot_1

se = grep("_VS_", colnames(data))
data_surv<-data[,-se]
str(info)





library(survival)
library(survminer)
info$Plasma_C9_group <- ifelse(data_surv$`Plasma.cells_C9` > median(data_surv$`Plasma.cells_C9`, na.rm = TRUE),
                               "High Plasma_C9", "Low Plasma_C9")

surv_obj <- Surv(time = info$t.surv, event = info$e.surv)
fit <- survfit(surv_obj ~ Plasma_C9_group, data = info)

km_plot <- ggsurvplot(fit,
                      data = info,
                      pval = TRUE,pval.method = T,
                      risk.table = F,
                      palette = c("darkblue", "darkred"),
                      legend.title = "Plasma.cells_C9",
                      legend.labs = c("High scores", "Low scores"))


km_plot



library(ggplot2)
info$recurrence <- trimws(info$recurrence)  
info$recurrence <- factor(info$recurrence, levels = c("no", "yes"),
                          labels = c("Non-recurrent", "Recurrent"))

plot_2<-ggplot(info, aes(x = recurrence, y = data_surv$`Plasma.cells_C9`, fill = recurrence)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Non-recurrent" = "#56B4E9", "Recurrent" = "#D55E00")) +
  labs(title = "Plasma.cells_C9 Expression by Recurrence Status",
       x = "Recurrence Status",
       y = "Plasma.cells_C9 signature scores") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format", 
                     label.y = 10)+
  theme_classic() +
  theme(legend.position = "none")

plot_2









rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/Lung/ImmGen/data/TCGA_LUAD_AmoscRNAMarkers_Base.txt"
myinf2 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/Firehose/done/Thorsson_2018_TCGA_immunelandscape.csv"


data = read.table(myinf1, sep="\t", header=T, row.names=1,  quote="")
mystd = apply(abs(data), 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/mystd[k]
}

#----------------------------
imm= read.table(myinf2, sep=",", header=T, row.names=1, stringsAsFactors=F)
se = which(imm$TCGA.Study=="LUAD")
imm = imm[se,]
se = c(4:31, 36:63)
imm = imm[,se]


comxx = intersect(row.names(data), row.names(imm))
data = data[comxx,]
imm = imm[comxx,]

myx = as.numeric(data[, "Plasma.cells_C9.ES"])
xx = cor(imm,  myx, method="s", use="pair")

xx = xx[order(xx[,1], decreasing=T), ]
xx







