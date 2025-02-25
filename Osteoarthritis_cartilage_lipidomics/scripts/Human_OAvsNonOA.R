# 1. Packages ##################################################

library(plyr)
library(tidyverse)
library(data.table)
library(pheatmap)
library(ggplot2)
library(ggbreak)
library(ggfortify)
library(broom)
library(ggrepel)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(sva)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)


# 2.  Sample annotation ##################################################
# 2.1 Group annotation
anno_table <- read.csv(file = "~/Sample information.csv",
                       header = TRUE,na.strings = ".",check.names = FALSE)

anno_table$gender2 <- ifelse(anno_table$Gender == "F", "Female", "Male")


# 2.2 Remove outliers
# Remove OA-7 (too young age)
anno_table <- anno_table %>% filter(Number != "OA-7")

# Remove OA-4 (outlier in PCA)
anno_table <- anno_table %>% filter(Number != "OA-4")


# 2.3 Adjust format
row.names(anno_table) <- anno_table$Number

anno_table <- anno_table %>% select(c("Age","gender2","Batch","Group"))

anno_table$Batch <- factor(anno_table$Batch)

colnames(anno_table)[2] <- "Gender"


# 2.4 Sample descriptive statistics
OA <- subset(anno_table, Group == "OA")
length(which(OA$Gender == "Female"))
length(which(OA$Gender == "Male"))
min(OA$Age)
max(OA$Age)
mean(OA$Age)
sd(OA$Age)

NonOA <- subset(anno_table, Group == "Non-OA")
length(which(NonOA$Gender == "Female"))
length(which(NonOA$Gender == "Male"))
min(NonOA$Age)
max(NonOA$Age)
mean(NonOA$Age)
sd(NonOA$Age)



# 3.  Read datasets and filter QC RSD within batch ##############################################

# @3.1 Batch 1 ###########################################
conc1 <- read.csv(file = "~/Batch1_concentration.csv",
                  header = TRUE,na.strings = ".",check.names = F)

# Adjust format
rownames(conc1) <- conc1[,1]
conc1 <- conc1[,-1]
rownames(conc1)[6] <- "NonOA-1"
conc1 <- conc1[-1,]

# Remove OA-4
conc1 <- conc1[-which(row.names(conc1) %like% "OA-4"),]

# Remove RSD QC > 25% within batch 1
qc1 <- conc1[which(row.names(conc1) %like% "QC"),]

RSD_1 <- apply(qc1, 2, function(x) sd(x)/mean(x))

length(which(RSD_1 > 0.25))

conc1_qc <- conc1[,-which(RSD_1 > 0.25)]


# @3.2 Batch 2 ##############################################
conc2 <- read.csv(file = "~/Batch2_concentration.csv",
                  header = TRUE,na.strings = ".",check.names = FALSE)

# Adjust format
conc2 <- conc2[c(5,4,6,7,3,2,10),]
rownames(conc2) <- c("NonOA-2","NonOA-3","OA-5","OA-6","OA-7","QC1-00011","QC2-00011")
conc2 <- conc2[,-(1:2)]

# Remove OA-7
conc2 <- conc2[-which(row.names(conc2) %like% "OA-7"),]

# Remove RSD QC > 25% within batch 2
qc2 <- conc2[which(row.names(conc2) %like% "QC"),]

RSD_2 <- apply(qc2, 2, function(x) sd(x)/mean(x))

length(which(RSD_2 > 0.25))

conc2_qc <- conc2[,-which(RSD_2 > 0.25)]


# @3.3 Batch 3 ##############################################
conc3 <- read.csv(file = "~/Batch3_concentration.csv",
                  header = TRUE,na.strings = ".",check.names = FALSE)

# Remove P211, P212 (surgery one week before)
# Adjust format
conc3 <- conc3[c(10,12,7,3,5,11,8,6,13,2,9,15),]
rownames(conc3) <- c("OA-8","OA-9","OA-10","OA-11",
                     "NonOA-4","NonOA-5","NonOA-6","NonOA-7","NonOA-8",
                     "QC1-00014","QC2-00014","QC3-00005")
conc3 <- conc3[,-1]

# Remove RSD QC > 25% within batch 3
qc3 <- conc3[which(row.names(conc3) %like% "QC"),]

RSD_3 <- apply(qc3, 2, function(x) sd(x)/mean(x))

length(which(RSD_3 > 0.25))

conc3_qc <- conc3[,-which(RSD_3 > 0.25)]


# @3.4 Batch 4 ##############################################
conc4 <- read.csv(file = "~/Batch4_concentration.csv",
                  header = TRUE,na.strings = ".",check.names = FALSE)

# Adjust format
conc4 <- conc4[c(5,4,7,9,12,3,13,10,11,6,2,8,14),]
rownames(conc4) <- c("NonOA-9","NonOA-10","NonOA-11","NonOA-12","NonOA-13",
                     "OA-12","OA-13","OA-14","OA-15","OA-16",
                     "QC1-00016","QC2-00016","QC3-00007")
conc4 <- conc4[,-1]

# Remove RSD QC > 25% within batch 4
qc4 <- conc4[which(row.names(conc4) %like% "QC"),]

RSD_4 <- apply(qc4, 2, function(x) sd(x)/mean(x))

length(which(RSD_4 > 0.25))

conc4_qc <- conc4[,-which(RSD_4 > 0.25)]



# 4.  Merge tables ################################################
a <- intersect(intersect(intersect(colnames(conc1_qc),colnames(conc2_qc)),colnames(conc3_qc)),colnames(conc4_qc))

conc_bind_all <- rbind.fill(conc1_qc, conc2_qc, conc3_qc, conc4_qc)
row.names(conc_bind_all) <- c(row.names(conc1_qc), row.names(conc2_qc), row.names(conc3_qc), row.names(conc4_qc))

conc_bind <- conc_bind_all[,a]

colnames(conc_bind)[colnames(conc_bind) %like% "DAG"] <- str_replace(colnames(conc_bind)[colnames(conc_bind) %like% "DAG"], "/", "_")
colnames(conc_bind)[colnames(conc_bind) %like% "PC"] <- str_replace(colnames(conc_bind)[colnames(conc_bind) %like% "PC"], "/", "_")
colnames(conc_bind)[colnames(conc_bind) %like% "PE"] <- str_replace(colnames(conc_bind)[colnames(conc_bind) %like% "PE"], "/", "_")
colnames(conc_bind)[colnames(conc_bind) %like% "O-"] <- str_replace(colnames(conc_bind)[colnames(conc_bind) %like% "O-"], "_", "/")
colnames(conc_bind)[colnames(conc_bind) %like% "P-"] <- str_replace(colnames(conc_bind)[colnames(conc_bind) %like% "P-"], "_", "/")



# 5.  Inter-batch QC RSD ############################################
# Remove RSD QC > 25% 
conc_bind_qc <- conc_bind[row.names(conc_bind) %like% "QC",]

RSD_all <- apply(conc_bind_qc, 2, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))

length(which(RSD_all > 0.25))

conc_qc <- conc_bind[,-which(RSD_all > 0.25)]



# 6.  Missing values ################################################
conc_qc_sample <- conc_qc[-which(row.names(conc_qc) %like% "QC"),]

# Check variables with missing values > 45%
length(which(colMeans(is.na(conc_qc_sample[,])) > 0.45))

# Remove variables with missing value > 45%
conc_rmNA <- conc_qc[,-which(colMeans(is.na(conc_qc_sample[,])) > 0.45)]

conc_rmNA_sample <- conc_rmNA[-which(row.names(conc_rmNA) %like% "QC"),]

conc_rmNA_qc <- conc_rmNA[which(row.names(conc_rmNA) %like% "QC"),]

# Replace NA by 1/2 of the minimum positive value of each variable
f_na = function(x) {
  x[is.na(x)] = min(x, na.rm=TRUE)/2
  x
}

conc_rpNA <- apply(conc_rmNA_sample,2,f_na)



# 7.  Reorder and normalize to cell number ##########################
# Reorder
conc_order <- conc_rpNA[match(row.names(anno_table),row.names(conc_rpNA)),]

# Normalize to cell number
conc_cell <- conc_order/(40*2.5)



# Composition sample wise
comp <- as.data.frame(conc_cell)
FFA <- rowSums(comp %>% select(starts_with("FFA")))
CE <- comp %>% select(starts_with("CE"))
CE <- rowSums(CE[,-which(colnames(CE) %like% "CER")])
DAG <- rowSums(comp %>% select(starts_with("DAG")))
TAG <- rowSums(comp %>% select(starts_with("TAG")))
PE <- rowSums(comp %>% select(starts_with("PE")))
PC <- rowSums(comp %>% select(starts_with("PC")))
SM <- rowSums(comp %>% select(starts_with("SM")))
CER <- rowSums(comp %>% select(starts_with("CER")))
DCER <- rowSums(comp %>% select(starts_with("DCER")))
HCER <- rowSums(comp %>% select(starts_with("HCER")))
LCER <- rowSums(comp %>% select(starts_with("LCER")))

comp_class <- rbind(FFA,CE,DAG,TAG,PE,PC,SM,CER,DCER,HCER,LCER)

comp_class_long <- melt(comp_class)

ggplot(comp_class_long, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="fill", stat="identity",width=0.9) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.title = element_text(size = 9, color = "black", face = "bold"),
        axis.text = element_text(size = 9, color = "black"),
        axis.line = element_line(size = 0.5), # line size = font size/22
        axis.ticks = element_line(size = 0.5)) +
  labs(x="",y="Proportions of lipid classes", fill="Lipid class")



# 8. Log transformation ##########################################
conc_cell_log <- log2(conc_cell)

conc_cell_log_t <- t(conc_cell_log)



# 9. Check quality ##############################################
# @9.1 Technical variability ############################
techVar <- apply(conc_rmNA_qc, 2, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100)

techVar <- data.frame(var = techVar, class = str_sub(colnames(conc_rmNA_qc), 1,3),
                      row.names = colnames(conc_rmNA_qc))

ggplot(techVar, aes(x=class, y=var)) +
  geom_boxplot()  +
  geom_point() +
  ylim(c(0,100)) # check x axis tick labels

ggplot(techVar, aes(x=class, y=var)) +
  geom_boxplot(outlier.size = 1, color="#848080", fill="#E6E4E4") +
  scale_x_discrete(labels = c("CE","CER","DAG","DCER","FFA","HCER","LCER","PC","PE","SM","TAG"),
                   guide = guide_axis(angle = 45)) +
  ylim(c(0,100)) +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(size = 0.5), # line size = font size/22
        axis.ticks = element_line(size = 0.5),
        panel.background = element_rect(fill = "white")) +
  labs(x="Lipid class", y="Coefficient of variation (%)")

ggsave("technical variation.png", width = 11, height = 8, units = "cm", dpi = 300)


# @9.3 Heatmap ##############################################
pheatmap(conc_cell_log_t,
         scale = "row",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#00aedb","white","#d11141"))(50),
         cellwidth = 10,
         cellheight = 0.5,
         display_numbers = FALSE,
         border = FALSE,
         show_rownames = FALSE,
         fontsize = 10, fontsize_row = 8, fontsize_col = 8,
         annotation_col = anno_table) 


# @9.4 PCA #################################################
# log-transformed scaled data
# Calculate PC
PCA.PC_qc <- prcomp(conc_cell_log, scale. = TRUE)

PCA_qc <- as.data.frame(conc_cell_log)
PCA_qc$Batch <- anno_table$Batch
PCA_qc$Age <- anno_table$Age
PCA_qc$Gender <- anno_table$Gender

autoplot(PCA.PC_qc, data = PCA_qc, colour = "Batch",
         frame = T, frame.type = "norm", frame.level=0.95, size = 2) +
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  theme_test()+
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_beforeCombat_batch.png", width = 1096, height = 871, units = "px", dpi = 300)

autoplot(PCA.PC_qc, data = PCA_qc, colour = "Age") +
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  geom_text_repel(aes(label = PCA_qc$Batch)) +
  theme_test()+
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_beforeCombat_age.png", width = 1100, height = 800, units = "px", dpi = 300)

autoplot(PCA.PC_qc, data = PCA_qc, colour = "Gender",
         frame = T, frame.type = "norm", frame.level=0.95, size = 2) +
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  scale_fill_manual(values = c("#56B4E9","#D55E00")) + 
  scale_color_manual(values = c("#56B4E9","#D55E00")) +
  geom_text_repel(aes(label = PCA_qc$Batch)) +
  theme_test()+
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_beforeCombat_gender.png", width = 1300, height = 950, units = "px", dpi = 300)



# 10. Batch correction #######################################
modcombat = model.matrix(~as.factor(Group), data = anno_table)

combat_data = ComBat(dat = conc_cell_log_t, batch = anno_table$Batch, mod = modcombat, par.prior = TRUE, prior.plots = T)

combat_data_t <- t(combat_data)

combat_data_t <- as.data.frame(combat_data_t)

conc <- 2^combat_data_t

# write.csv(combat_data_t, file = "human_combat_data_t.csv")
# write.csv(conc, file = "Table S2_hACs lipidomics.csv")


  
# 11. Lipid annotation ###################################################
# Create annotation table
anno_table_lpd <- data.frame(Class = word(colnames(combat_data_t),1,sep = "\\("),
                             row.names = colnames(combat_data_t))

anno_table_lpd$Class[which(colnames(combat_data_t) %like% "TAG")] <- "TAG"

unique(anno_table_lpd$Class)

anno_table_lpd$Class <- factor(anno_table_lpd$Class,
                               levels = c("FFA","CE","DAG","TAG","PE","PC","LPE","LPC","SM","CER","DCER","HCER","LCER"),
                               ordered = TRUE)

# Create lipid palette
pal_custom <- c("#31992A","#93D073", 
                "#88BBD9","#1F78B4",
                "#E5423C","#F39b7F","#987884","#DAE7E6",
                "#EEE48D","#83EADD","#B295C7","#EB47A3","#F4CCE0")
names(pal_custom) <- c("FFA","CE","TAG","DAG","PE","PC","LPE","LPC","SM","CER","DCER","HCER","LCER")



# 12. Lipid number proportion  #######################################
countClass <- anno_table_lpd %>% count(Class)
countClass <- countClass[(order(countClass$n, decreasing = TRUE)),]
countClass$perc <- countClass$n/sum(countClass$n)*100



# 13.  Multivariate analysis###################################

# @13.1 PCA #################################################
# ****13.1.1 Sample PCA ###################################
PCA.PC_sample <- prcomp(combat_data_t, scale. = TRUE)

?prcomp

fviz_pca_ind(PCA.PC_sample,
             habillage = anno_table$Group,
             legend.title = "Group",
             repel = TRUE,
             addEllipses= TRUE)

# Contributions of variables to PC1
fviz_contrib(PCA.PC_sample, choice = "var", axes = 1, top = 30)
# Contributions of variables to PC2
fviz_contrib(PCA.PC_sample, choice = "var", axes = 2, top = 30)

# ggplot2
PCAdata_sample <- combat_data_t

PCAdata_sample <- as.data.frame(PCAdata_sample)

PCAdata_sample$Group <- anno_table$Group

PCAdata_sample$Batch <- anno_table$Batch

PCAdata_sample$Age <- anno_table$Age

PCAdata_sample$Gender <- anno_table$Gender

# PCA group
autoplot(PCA.PC_sample, data = PCAdata_sample, colour = "Group",
         frame = TRUE, frame.type = "norm", frame.level=0.95, size = 2) +
  scale_fill_manual(values = c("#848080","#E69F00")) + 
  scale_color_manual(values = c("#848080","#E69F00")) +
  geom_hline(yintercept=0,linetype=3, color = "grey") +
  geom_vline(xintercept=0,linetype=3, color = "grey") +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca1.png", dpi = 300)

# PCA age
autoplot(PCA.PC_sample, data = PCAdata_sample, colour = "Age") +
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  geom_text_repel(aes(label = PCAdata_sample$Batch)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_afterCombat_age.png", width = 1085, height = 800, unit = "px", dpi = 300)

# PCA batch
autoplot(PCA.PC_sample, data = PCAdata_sample, colour = "Batch") +
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  theme_test()+
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_batch.png", dpi = 300)

# PCA gender
autoplot(PCA.PC_sample, data = PCAdata_sample, colour = "Gender",
         frame = TRUE, frame.type = "norm", frame.level=0.95, size = 2) +
  scale_fill_manual(values = c("#56B4E9","#D55E00")) + 
  scale_color_manual(values = c("#56B4E9","#D55E00")) +
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  geom_text_repel(aes(label = PCAdata_sample$Batch)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_afterCombat_gender.png", width = 1300, height = 950, unit = "px", dpi = 300)


# ****13.1.2 Lipid PCA ####################################################
PCA.PC_lipid <- prcomp(combat_data)

PCAdata_lipid <- as.data.frame(combat_data)
PCAdata_lipid$Class <- anno_table_lpd$Class

autoplot(PCA.PC_lipid, data = PCAdata_lipid, colour = "Class",
         frame = FALSE, frame.type = "norm", frame.level=0.95, size = 2)+
  geom_hline(yintercept=0,linetype=3, color = "grey") +
  geom_vline(xintercept=0,linetype=3, color = "grey") +
  scale_color_manual(values= pal_custom) +
  theme_test()



# @13.2 Clustering #######################################################
# For checking
pheatmap(combat_data,
         scale = "row",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#00aedb","white","#d11141"))(100),
         cellwidth = 11,
         cellheight = 0.8,
         display_numbers = FALSE,
         border = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize = 10, fontsize_row = 8, fontsize_col = 8,
         angle_col = 45,
         annotation_col = anno_table,
         annotation_row = anno_table_lpd,
         annotation_names_row = FALSE)


# For presenting
# pheatmap
pheatmap(combat_data,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#00aedb","white","#d11141"))(100),
         cellwidth = 11,
         cellheight = 0.8,
         display_numbers = FALSE,
         border = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         angle_col = 45,
         annotation_row = anno_table_lpd,
         annotation_names_row = FALSE,
         gaps_col = c(13))


# ComplexHeatmap
?scale # column scaling
hmdata <- t(scale(t(combat_data), center = TRUE, scale = TRUE))
min(hmdata)
max(hmdata)

col_fun <- colorRamp2(c(-4,-2,0,2,4), c("#51446E","#937BB1","white","#EE9591","#F55C45"))

row_anno <- rowAnnotation(
  `Lipid class` = anno_table_lpd$Class,
  show_annotation_name = F,
  annotation_legend_param = list(labels_gp = gpar(fontsize = 11),
                                 title_gp = gpar(fontsize = 11, fontface = "bold")),
  col = list(
    `Lipid class`= c("FFA" = "#31992A",
                     "CE" = "#93D073",
                     "TAG" = "#88BBD9",
                     "DAG" = "#1F78B4",
                     "PE" = "#E5423C",
                     "PC" = "#F39b7F",
                     "SM" = "#EEE48D",
                     "CER" = "#83EADD",
                     "DCER" = "#B295C7",
                     "HCER" = "#EB47A3",
                     "LCER" = "#F4CCE0")
    )
  )

png(file="heatmap1.png", width = 1100, height = 1300, res = 300)
heatmap1 <- Heatmap(hmdata,
                    name = "Z-score",
                    col = col_fun,
                    cluster_columns = F,
                    column_split = anno_table$Group,
                    column_gap = unit(2,"mm"),
                    border = T,
                    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_row_names = F,
                    show_column_names = F,
                    left_annotation = row_anno,
                    heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                    title_gp = gpar(fontsize = 11, fontface = "bold"),
                                    direction = "vertical")
                    #width = unit(8,"cm"),
                    #height = unit(8,"cm")
                    )
draw(heatmap1)
dev.off()



# 14.  Differential lipids #######################################################

# @14.1 Significant test #######################################################
anno_table$Group[15:27] <- rep("NonOA",13)

design <- model.matrix(~0+factor(anno_table$Group))
colnames(design) <- levels(factor(anno_table$Group))
row.names(design) <- row.names(anno_table)

contrast.matrix <- makeContrasts(paste0(unique(anno_table$Group), collapse = "-"), levels = design)

fit <- lmFit(combat_data, design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, adjust.method ="BH", number=1000)

write.csv(output, "human_output_original.csv")



# @14.2 Volcano plot #######################################################
logic1 <- (output$logFC > 0 & output$adj.P.Val< 0.05)
logic2 <- (output$logFC < 0 & output$adj.P.Val < 0.05)

output$regulation <- ifelse(logic1,"Up-regulated",ifelse(logic2,"Down-regulated","Non-significant"))
output$regulation <- factor(output$regulation, levels = c("Non-significant", "Down-regulated", "Up-regulated", ordered = TRUE))

# Mark down-regulated lipids
down <- subset(output,regulation == "Down-regulated")
down_top40 <- down[order(down$logFC),][1:40,]

write.csv(down,file = "C:/Users/u0134161/Documents/Qiongfei Zhou/07 Raw data/1 AC lipidomics/Human/20220214_batch1-4/combat_down_p0.05.csv")

# Mark up-regulated lipids
up <- subset(output,regulation == "Up-regulated")
up <- up[order(up$logFC),]

diff <- rbind(down,up)
diff <- diff[order(diff$logFC),]

write.csv(diff,file = "diff.csv")

oddFA1 <- diff[(row.names(diff) %like% "17:0"),]
oddFA2 <- diff[(row.names(diff) %like% "15:0"),]
oddFA <- rbind(oddFA1,oddFA2)
write.csv(oddFA,file = "combat_diff_oddFA.csv")

# Plot
FC_expression <- expression(bold(paste(Log[2], "FC (OA/non-OA)")))
p_expression <- expression(paste(bold(-Log[10]), bold("(BH "), bolditalic("P"), bold("-value)")))

ggplot(output,aes(x = logFC,y = -log10(adj.P.Val))) + 
  geom_hline(yintercept = -log10(0.05), color = "#8A8989", linetype = 3) +
  geom_point(cex = 3, aes(fill = regulation), shape= 21, color = "black", stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = c("#B3B3B3","#BCADCF","#F2796B")) + 
  xlim(c(-1.8, 1.8)) +
  ylim(c(0,7)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        #legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x = FC_expression, y = p_expression, fill = "")

ggsave("volcano.png", dpi = 300)

# Plot with class
output_volcano <- merge(anno_table_lpd, output, by = "row.names")
row.names(output_volcano) <- output_volcano$Row.names

output_volcano$Color <- as.character(output_volcano$Class)

output_volcano$Color[which(output_volcano$regulation == "Non-significant")] <- "Non-significant"

output_volcano$Color <- factor(output_volcano$Color, levels = c("Non-significant","LCER","DCER","SM","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

pal_volcano <- c(FFA="#31992A",CE="#93D073", DAG="#1F78B4",TAG="#88BBD9",PE="#E5423C",PC="#F39b7F",SM="#EEE48D",DCER="#B295C7",LCER="#F4CCE0")

volcano_text <- output_volcano[c("PE(O-18:0/18:2)","TAG53:2-FA18:2", "TAG53:3-FA17:0", "PC(18:0_22:4)", "PE(P-18:1/22:4)", "PC(18:0_20:4)"),]

ggplot(output_volcano,aes(x = logFC,y = -log10(adj.P.Val))) + 
  geom_hline(yintercept = -log10(0.05), color = "#8A8989", linetype = 3) +
  geom_point(size = 1, aes(color = Color), alpha = 0.7) +
  scale_color_manual(values = c(pal_volcano,c("Non-significant"="#B3B3B3")),
                    limits = c("FFA","CE", "DAG", "TAG", "PE", "PC", "SM","DCER","LCER", "Non-significant")) +
  xlim(c(-1.8, 1.8)) +
  ylim(c(0,7)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")) +
  labs(x = FC_expression, y = p_expression, color = "") +
  geom_text_repel(data = volcano_text,
                  aes(x = logFC,y = -log10(adj.P.Val), label = row.names(volcano_text)), 
                  size = 3, box.padding = unit(0.5, "lines"))

ggsave("volcano_class_transparent.png", dpi = 300, width = 11, height = 7, units = "cm", bg = "transparent")



# @14.3 Bubble plot #######################################################
bubble <- output["logFC"]
bubble$adjp <- -log10(output$adj.P.Val)

bubble <- bubble[match(row.names(anno_table_lpd),row.names(bubble)),]

bubble$Class <- anno_table_lpd$Class

ggplot(bubble, aes(x = logFC, y = Class)) + 
  geom_vline(xintercept = c(-1, 0, 1), color = "#8A8989", linetype = 3) +
  geom_point(aes(size = adjp, fill = Class), shape= 21, color = "black", stroke = 0.8, alpha = 0.7) +
  scale_size(range = c(1.5, 7), name = p_expression) +
  scale_fill_manual(values = pal_custom) +
  scale_y_discrete(limits = rev(levels(anno_table_lpd$Class))[-c(6, 7)]) + #remove LPC LPE label
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white")) +
  labs(x = FC_expression, y = "Lipid class") +
  guides(fill = "none") 

ggsave("bubble.png", dpi = 300)

ggplot(bubble, aes(x = logFC, y = Class)) + 
  geom_vline(xintercept = c(-1, 0, 1), color = "#8A8989", linetype = 3) +
  geom_point(aes(size = adjp, fill = Class), shape= 21, color = "black", stroke = 0.8, alpha = 0.7) +
  scale_size(range = c(1.5, 7), name = p_expression) +
  scale_fill_manual(values = pal_custom) +
  scale_y_discrete(limits = rev(levels(anno_table_lpd$Class))[-c(6, 7)]) + #remove LPC LPE label
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position = "bottom"
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent"),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(fill = "transparent")
        ) +
  labs(x = FC_expression, y = "Lipid class") +
  guides(fill = "none") 
ggsave("bubble_transperent by r.png", width = 1206, height = 1181, units = "px", bg = "transparent", dpi = 300)

ggplot(bubble, aes(x = logFC, y = Class)) +
  geom_vline(xintercept = c(-1, 0, 1), color = "#8A8989", linetype = 3) +
  geom_boxplot(aes(fill = Class), width = 0.5, alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = pal_custom) + 
  scale_y_discrete(limits = rev(levels(anno_table_lpd$Class))[-c(6, 7)]) + #remove LPC LPE label
  xlim(c(-1.8,1.2)) +
  theme_test() +
  theme(axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        #legend.title = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x ="log2(FC)", y = "Lipid class", fill = "LipidClass")

ggsave("all lipid boxplot.png", dpi = 300)


# @14.4 Proportion of differential lipid class##############################
down_class <- merge(down,anno_table_lpd, by = "row.names")
row.names(down_class) <- down_class$Row.names
unique(down_class$Class)
down_class <- down_class %>% count(Class)
down_class <- down_class[order(down_class$n, decreasing = TRUE),]
down_class$prec <- down_class$n/sum(down_class$n)*100

up_class <- merge(up,anno_table_lpd, by = "row.names")
row.names(up_class) <- up_class$Row.names
up_class %>% count(Class)


# @14.5 Heatmap top regulated lipids ########################################
# Down top 20
hmdata_down <- hmdata[(row.names(hmdata) %in% row.names(down_top20)),]

hmdata_down <- hmdata_down[match(row.names(down_top20),row.names(hmdata_down)),]

min(-log10(down_top20$adj.P.Val))
max(-log10(down_top20$adj.P.Val))

col_down_p <- colorRamp2(row.names(hmdata_down), colorRampPalette(c("#00aedb","white","#d11141"))(20))

row_anno_down <- rowAnnotation(
  `log2(FC)` = anno_barplot(
    down_top20$logFC, 
    baseline = -1,
    #gp = gpar(
      #fill = col_down_p
      #)
    )
  )

png(file="heatmap2.png", width = 1500, height = 1200, res = 300)
heatmap2 <- Heatmap(hmdata_down,
                    name = "Z-score",
                    col = col_fun,
                    cluster_columns = F,
                    cluster_rows = F,
                    column_split = anno_table$Group,
                    column_gap = unit(1,"mm"),
                    border = T,
                    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                    show_row_names = T,
                    show_column_names = F,
                    row_names_gp = gpar(fontsize = 11),
                    right_annotation = row_anno_down,
                    heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                                title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                direction = "horizontal")
                    #width = unit(8,"cm"),
                    #height = unit(8,"cm")
)

draw(heatmap2)
dev.off()

# Down top 40
hmdata_down <- hmdata[(row.names(hmdata) %in% row.names(down_top40)),]

hmdata_down <- hmdata_down[match(row.names(down_top40),row.names(hmdata_down)),]

min(-log10(down_top40$adj.P.Val))
max(-log10(down_top40$adj.P.Val))

png(file="heatmap_down40.png", width = 1200, height = 1400, res = 300)
heatmap_down40 <- Heatmap(hmdata_down,
                    name = "Z-score",
                    col = col_fun,
                    cluster_columns = F,
                    cluster_rows = F,
                    column_split = anno_table$Group,
                    column_gap = unit(1,"mm"),
                    border = T,
                    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_row_names = T,
                    show_column_names = F,
                    row_names_gp = gpar(fontsize = 11),
                    #right_annotation = row_anno_down,
                    heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                                title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                direction = "horizontal")
                    #width = unit(8,"cm"),
                    #height = unit(8,"cm")
)

draw(heatmap_down40)
dev.off()

down_top40$FA <- row.names(down_top40)

ggplot(down_top40,aes(x = logFC, y = FA, fill = -log10(adj.P.Val))) +
  geom_bar(stat="identity",width=0.7) +
  scale_y_discrete(limits = rev(row.names(down_top40)), position = "right") +
  scale_x_continuous(breaks = c(-1,0)) +
  scale_fill_gradient(low = "#C9F3D4", high = "#35D35D", limits = c(1,6.5)) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  theme_test() +
  theme(axis.ticks = element_line(size = 0.5),
        legend.position="bottom",
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        #legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x= FC_expression, y="", fill = p_expression)

ggsave("bar_downtop40.png", width = 709, height = 1721, units = "px", dpi = 300)
ggsave("bar_downtop40_2.png", width = 1331, height = 1721, units = "px", dpi = 300)


# Top
hmdata_up <- hmdata[(row.names(hmdata) %in% row.names(up)),]

hmdata_up <- hmdata_up[match(row.names(up),row.names(hmdata_up)),]

hmdata_up <- hmdata_up[rev(1:nrow(hmdata_up)),]

min(-log10(up$adj.P.Val))
max(-log10(up$adj.P.Val))

png(file="heatmap_up.png", width = 1200, height = 500, res = 300)
heatmap_up <- Heatmap(hmdata_up,
                          name = "Z-score",
                          col = col_fun,
                          cluster_columns = F,
                          cluster_rows = F,
                          column_split = anno_table$Group,
                          column_gap = unit(1,"mm"),
                          border = T,
                          column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                          show_row_names = T,
                          show_column_names = F,
                          row_names_gp = gpar(fontsize = 11),
                          heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                                      title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                      direction = "horizontal"))

draw(heatmap_up)
dev.off()

up$FA <- row.names(up)

ggplot(up,aes(x = logFC, y = FA, fill = -log10(adj.P.Val))) +
  geom_bar(stat="identity",width=0.7) +
  scale_y_discrete(limits = row.names(up), position = "right") +
  scale_x_continuous(limits = c(0,1), breaks = c(0,1)) +
  scale_fill_gradient(low = "#C9F3D4", high = "#35D35D", limits = c(1,6.5)) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  theme_test() +
  theme(axis.ticks = element_line(size = 0.5),
        legend.position="bottom",
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        #legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x= FC_expression, y="", fill = p_expression)

ggsave("bar_up.png", width = 700, height = 1000, units = "px", dpi = 300)



# 16.   Deconstruct lipids ########################################################
# @16.1 Batch-corrected data ############################
DAG <- combat_data_t %>% select(starts_with("DAG"))
colnames(DAG) <- paste(colnames(DAG),"_1", sep = "")

PC <- combat_data_t %>% select(starts_with("PC"))
colnames(PC) <- paste(colnames(PC),"_1", sep = "")

PE <- combat_data_t %>% select(starts_with("PE"))
colnames(PE) <- paste(colnames(PE),"_1", sep = "")

combat_data_FA <- cbind(combat_data_t,DAG,PC,PE)


# @16.2 Differential data ################################
DAG_FC <- output[which(row.names(output) %like% "DAG"),]

row.names(DAG_FC) <- paste(row.names(DAG_FC),"_1", sep = "")

PC_FC <- output[which(row.names(output) %like% "PC"),]

row.names(PC_FC) <- paste(row.names(PC_FC),"_1", sep = "")

PE_FC <- output[which(row.names(output) %like% "PE"),]

row.names(PE_FC) <- paste(row.names(PE_FC),"_1", sep = "")

output_FA <- rbind(output,DAG_FC,PC_FC,PE_FC)


# @16.3 Fatty acyl annotation ######################################################
anno_table_lpd$FA <- str_sub(row.names(anno_table_lpd),-5,-2)

anno_table_lpd$FA[which(row.names(anno_table_lpd) %like% "TAG")] <- str_sub(row.names(anno_table_lpd)[which(row.names(anno_table_lpd) %like% "TAG")],-4,-1)

DAG_chain1st <- data.frame(row.names = colnames(DAG), Class = rep("DAG", ncol(DAG)), FA = str_sub(colnames(DAG), -12,-9))
PE_chain1st <- data.frame(row.names = colnames(PE), Class = rep("PE", ncol(PE)), FA = str_sub(colnames(PE), -12,-9))
PC_chain1st <- data.frame(row.names = colnames(PC), Class = rep("PC", ncol(PC)), FA = str_sub(colnames(PC), -12,-9))

anno_table_lpd_FA <- rbind(anno_table_lpd, DAG_chain1st, PE_chain1st, PC_chain1st)


anno_table_lpd_FA$Length <- str_sub(anno_table_lpd_FA$FA,1,2)
anno_table_lpd_FA$Unsaturation <- str_sub(anno_table_lpd_FA$FA,4,4)


anno_table_lpd_FA$LengthClass <- ifelse(anno_table_lpd_FA$Length >21, "VLCFA", 
                                        ifelse(anno_table_lpd_FA$Length < 13, "MCFA", "LCFA"))
anno_table_lpd_FA$UnsatClass <- ifelse(anno_table_lpd_FA$Unsaturation == 0, "SFA", 
                                       ifelse(anno_table_lpd_FA$Unsaturation == 1, "MUFA", "PUFA"))


output_FA <- merge(anno_table_lpd_FA, output_FA, by = "row.names")
row.names(output_FA) <- output_FA$Row.names

output_FA$FA <- paste("C", output_FA$FA, sep = "")

min(output_FA$logFC)
max(output_FA$logFC)

# 17. Analysis on fatty acyls ######################################################

# @ 17.1 Certain fatty acyls #####################################################
# **** 17.1.1 FA in all the lipids ###########################################
ggplot(output_FA, aes(x = logFC, y = FA)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.7, outlier.shape = NA) +
  geom_jitter(size = 0.7, aes(color = -log10(adj.P.Val)), height = 0.3) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  # scale_color_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$FA)))) +
  # scale_fill_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$FA)))) +
  scale_y_discrete(limits = rev) +
  xlim(c(-1.7,1.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position="none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Fatty acid chain", x = FC_expression, color = p_expression)

ggsave("allFA_FC+p_2.png", dpi = 300, width = 821, height = 1481, units = "px")

# Legend
ggplot(output_FA, aes(x = logFC, y = FA)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.7, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.3) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = rev) +
  xlim(c(-1.7,1.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position="bottom",
        panel.background = element_rect(fill = "white")) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  labs(color = p_expression)

ggsave("legend.png", dpi = 300)


# **** 17.1.2 FA in significantly regulated lipids ##########################
output_FA_sig <- subset(output_FA, regulation == "Up-regulated" | regulation == "Down-regulated")

# ATTENTION!!!!!!!!
unique(output_FA_sig$Class)

pal_custom_sig <- c(FFA="#31992A",CE="#93D073", DAG="#1F78B4",TAG="#88BBD9",PE="#E5423C",PC="#F39b7F",SM="#EEE48D",DCER="#B295C7",LCER="#F4CCE0")

ggplot(output_FA_sig, aes(x = logFC, y = FA, color = Class)) + 
  geom_hline(yintercept = unique(output_FA_sig$FA), linetype = 3, color = "grey") +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_jitter(aes(color = Class), size = 1.1, alpha = 0.7, height = 0.3) +
  scale_color_manual(values = pal_custom_sig) + 
  scale_y_discrete(limits = rev) +
  xlim(c(-1.7,1.2)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Fatty acid chain", color = "Lipid class")

ggsave("allFA_sig_2.png", dpi = 300, width = 1081, height = 1481, units = "px")



# **** 17.1.3 The number of individual FA #####################################################
count_FA <- output_FA %>% group_by(regulation,Class) %>% count(FA)

count_FA$Color <- as.character(count_FA$Class)

count_FA$Color[which(count_FA$regulation == "Non-significant")] <- "Non-significant"

count_FA$Color <- factor(count_FA$Color, levels = c("Non-significant","LCER","DCER","SM","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA, aes(x = n, y = FA, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")),
                    limits = c("FFA","CE", "DAG", "TAG", "PE", "PC", "SM","DCER","LCER", "Non-significant")) +
  scale_y_discrete(limits = rev) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Fatty acid chain", fill = "Lipid class")

count_FA[which(count_FA$FA == "18:1"),]

ggsave("allFA_n_2.png", dpi = 300, width = 1193, height = 1481, units = "px")


# **** 17.1.4 FA PCA  #####################################################
PCA.PC_indiFA <- prcomp(t(combat_data_FA))

PCAdata_indiFA <- as.data.frame(t(combat_data_FA))
PCAdata_indiFA$FA <- anno_table_lpd_FA$FA

autoplot(PCA.PC_indiFA, data = PCAdata_indiFA, colour = "FA",
         frame = FALSE, frame.type = "norm", frame.level=0.95, size = 2)+
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  scale_color_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$FA))))



# @ 17.2 The length of fatty acyls #####################################################

# **** 17.2.1 FA in all the lipids #####################################################
# Group by carbon number
ggplot(output_FA, aes(x = logFC, y = Length)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.7, aes(color = -log10(adj.P.Val)), height = 0.4) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = rev) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Carbon number of fatty acid chain", x = FC_expression)

ggsave("length_FC+p_2.png", dpi = 300, width = 706, height = 1265, units = "px")

# Group by length classification
ggplot(output_FA, aes(x = logFC, y = LengthClass)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.7, aes(color = -log10(adj.P.Val)), height = 0.4) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = c("VLCFA", "LCFA", "MCFA")) +
  xlim(c(-1.8,1.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Length of fatty acid chain", x = FC_expression, color = p_expression)

ggsave("lengthClass_FC+p_2.png", dpi = 300, width = 706, height = 1040, units = "px")



# **** 17.2.2 FA in significantly regulated lipids #####################################################
# Group by carbon number
ggplot(output_FA_sig, aes(x = logFC, y = Length, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1.1, alpha = 0.7, height = 0.4) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = rev) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Carbon number of fatty acid chain", color = "Lipid class")

ggsave("length_sig_2.png", dpi = 300, width = 993, height = 1265, units = "px")

# Group by length classification
ggplot(output_FA_sig, aes(x = logFC, y = LengthClass, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1.1, alpha = 0.7, height = 0.5) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = c("VLCFA", "LCFA", "MCFA")) +
  xlim(c(-1.8,1.2)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Length of fatty acid chain", color = "Lipid class")

ggsave("lengthClass_sig_3.png", dpi = 300, width = 1081, height = 1043, units = "px")


 # **** 17.2.3 The amount of individual FA #####################################################
# Group by carbon number
count_FA_length <- output_FA %>% group_by(regulation,Class) %>% count(Length)

count_FA_length$Color <- as.character(count_FA_length$Class)

count_FA_length$Color[which(count_FA_length$regulation == "Non-significant")] <- "Non-significant"

count_FA_length$Color <- factor(count_FA_length$Color, levels = c("Non-significant","LCER","DCER","SM","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_length, aes(x = n, y = Length, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")),
                    limits = c("FFA","CE", "DAG", "TAG", "PE", "PC", "SM","DCER","LCER", "Non-significant")) +
  scale_y_discrete(limits = rev) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Carbon number of fatty acid chain", fill = "Lipid class")

ggsave("length_n_2.png", dpi = 300, width = 1143, height = 1265, units = "px")


# Group by length classification
count_FA_lengthClass <- output_FA %>% group_by(regulation,Class) %>% count(LengthClass)

count_FA_lengthClass$Color <- as.character(count_FA_lengthClass$Class)

count_FA_lengthClass$Color[which(count_FA_lengthClass$regulation == "Non-significant")] <- "Non-significant"

count_FA_lengthClass$Color <- factor(count_FA_lengthClass$Color, levels = c("Non-significant","LCER","DCER","SM","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_lengthClass, aes(x = n, y = LengthClass, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")),
                    limits = c("FFA","CE", "DAG", "TAG", "PE", "PC", "SM","DCER","LCER", "Non-significant")) +
  scale_y_discrete(limits = c("VLCFA", "LCFA", "MCFA")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Length of fatty acid chain", fill = "Lipid class")

count_FA_lengthClass[which(count_FA_lengthClass$LengthClass == "LCFA"),]

ggsave("lengthClass_n_2.png", dpi = 300, width = 1150, height = 1040, units = "px")



# **** 17.2.4 FA length PCA#####################################################
# Group by carbon number
PCAdata_indiFA$Length <- anno_table_lpd_FA$Length

autoplot(PCA.PC_indiFA, data = PCAdata_indiFA, colour = "Length",
         frame = FALSE, frame.type = "norm", frame.level=0.95, size = 2)+
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  scale_color_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$Length))))

# Group by length classification
PCAdata_indiFA$LengthClass <- anno_table_lpd_FA$LengthClass

autoplot(PCA.PC_indiFA, data = PCAdata_indiFA, colour = "LengthClass",
         frame = FALSE, frame.type = "norm", frame.level=0.95, size = 2)+
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  scale_color_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$LengthClass))))




# @ 17.3 The Unsaturation of fatty acyls #####################################################

# **** 17.3.1 FA in all the lipids #####################################################
# Group by bouble bond number
ggplot(output_FA, aes(x = logFC, y = Unsaturation)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.7, aes(color = -log10(adj.P.Val)), height = 0.4) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = rev) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Double bond number of fatty acid chain", x = FC_expression)

ggsave("sat_FC+p_2.png", dpi = 300, width = 656, height = 1265, units = "px")


# Group by saturation classification
ggplot(output_FA, aes(x = logFC, y = UnsatClass)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.7, aes(color = -log10(adj.P.Val)), height = 0.4) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = c("PUFA", "MUFA", "SFA")) +
  xlim(c(-1.8,1.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Unsaturation of fatty acid chain", x = FC_expression)

ggsave("satClass_FC+p_2.png", dpi = 300, width = 706, height = 1040, units = "px")


# **** 17.3.2 FA in significantly regulated lipids #####################################################
# Group by double bond number
ggplot(output_FA_sig, aes(x = logFC, y = Unsaturation, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1.1, alpha = 0.7, height = 0.4) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = rev) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Double bond number of fatty acid chain", color = "Lipid class")

ggsave("sat_sig_2.png", dpi = 300, width = 968, height = 1265, units = "px")


# Group by saturation classification
ggplot(output_FA_sig, aes(x = logFC, y = UnsatClass, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1.1, alpha = 0.7, height = 0.4) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = c("PUFA", "MUFA", "SFA")) +
  xlim(c(-1.8,1.2)) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Unsaturation of fatty acid chain", color = "Lipid class")

ggsave("satClass_sig_2.png", dpi = 300, width = 1081, height = 1043, units = "px")



# **** 17.3.3 The amount of individual FA #####################################################
# Group by double bond number
count_FA_Unsat <- output_FA %>% group_by(regulation,Class) %>% count(Unsaturation)

count_FA_Unsat$Color <- as.character(count_FA_Unsat$Class)

count_FA_Unsat$Color[which(count_FA_Unsat$regulation == "Non-significant")] <- "Non-significant"

count_FA_Unsat$Color <- factor(count_FA_Unsat$Color, levels = c("Non-significant","LCER","DCER","SM","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_Unsat, aes(x = n, y = Unsaturation, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")),
                    limits = c("FFA","CE", "DAG", "TAG", "PE", "PC", "SM","DCER","LCER", "Non-significant")) +
  scale_y_discrete(limits = rev) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Double bond number of fatty acid chain", fill = "Lipid class")

ggsave("Unsat_n_2.png", dpi = 300, width = 1143, height = 1265, units = "px")


# Group by saturation classification
count_FA_UnsatClass <- output_FA %>% group_by(regulation,Class) %>% count(UnsatClass)

count_FA_UnsatClass$Color <- as.character(count_FA_UnsatClass$Class)

count_FA_UnsatClass$Color[which(count_FA_UnsatClass$regulation == "Non-significant")] <- "Non-significant"

count_FA_UnsatClass$Color <- factor(count_FA_UnsatClass$Color, levels = c("Non-significant","LCER","DCER","SM","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_UnsatClass, aes(x = n, y = UnsatClass, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")),
                    limits = c("FFA","CE", "DAG", "TAG", "PE", "PC", "SM","DCER","LCER", "Non-significant")) +
  scale_y_discrete(limits = c("PUFA", "MUFA", "SFA")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Unsaturation of fatty acid chain", fill = "Lipid class")

count_FA_UnsatClass %>% group_by(UnsatClass) %>% summarise(sum=sum(n))
count_FA_UnsatClass %>% group_by(UnsatClass, regulation) %>% summarise(sum=sum(n))

ggsave("UnsatClass_n_2.png", dpi = 300, width = 1143, height = 1040, units = "px")


# **** 17.3.4 FA Unsaturation PCA#####################################################
# Group by double bond number
PCAdata_indiFA$Unsaturation <- anno_table_lpd_FA$Unsaturation

autoplot(PCA.PC_indiFA, data = PCAdata_indiFA, colour = "Unsaturation",
         frame = FALSE, frame.type = "norm", frame.level=0.95, size = 2)+
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  scale_color_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$Unsaturation))))

# Group by saturation classification
PCAdata_indiFA$UnsatClass <- anno_table_lpd_FA$UnsatClass

autoplot(PCA.PC_indiFA, data = PCAdata_indiFA, colour = "UnsatClass",
         frame = FALSE, frame.type = "norm", frame.level=0.95, size = 2)+
  geom_hline(yintercept=0,linetype=3) +
  geom_vline(xintercept=0,linetype=3) +
  scale_color_manual(values = colorRampPalette(c("#aee9fe","blue"))(length(unique(anno_table_lpd_FA$UnsatClass))))

ggsave("test.png",dpi = 300)


# 18. Analysis on fatty acyls group by classes ######################################
ggplot(output_FA, aes(x = Length, y = logFC)) +
  geom_hline(yintercept = c(-1, 0, 1), color = "#8A8989", linetype = 3) +
  geom_point(aes(size = -log10(adj.P.Val), fill = Class), shape= 21, color = "black", stroke = 0.8, alpha = 0.6) +
  scale_size(range = c(1, 5)) +
  scale_fill_manual(values= pal_custom) +
  facet_wrap(~Class) +
  ylim(c(-1.8, 1.1)) +
  theme_test() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x ="Carbon number of fatty acyls", y = FC_expression, size = p_expression) +
  guides(fill = "none")

ggsave("length_class.png", dpi = 300)    


ggplot(output_FA, aes(x = Unsaturation, y = logFC)) +
  geom_hline(yintercept = c(-1, 0, 1), color = "#8A8989", linetype = 3) +
  geom_point(aes(size = -log10(adj.P.Val), fill = Class), shape= 21, color = "black", stroke = 0.8, alpha = 0.6) +
  scale_size(range = c(1, 5)) +
  scale_fill_manual(values= pal_custom) +
  facet_wrap(~Class) +
  ylim(c(-1.8, 1.1)) +
  theme_test() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x ="Double bond number of fatty acyls", y = FC_expression, size = p_expression) +
  guides(fill = "none")

ggsave("sat_class.png", dpi = 300)
