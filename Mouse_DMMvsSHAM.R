# 1.  Packages ##################################################

library(plyr)
library(tidyverse)
library(data.table)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(broom)
library(ggrepel)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(limma)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(ggpubr)
library(car)


# 2.  Sample annotation #########################################

# Read annotation table 
anno_table <- read.csv(file = "~/Sample_annotation.csv",
                       header = TRUE,check.names = FALSE)


row.names(anno_table) <- anno_table[,1]

anno_table$Time <- c(rep("W2",20),rep("W4",20),rep("W12",20))

anno_table$Group <- rep(c(rep("DMM",10),rep("Sham",10)),3)

anno_table$Side <- rep(c("R","L"),30)

# Technical issue with 4_DMM_R_3413, remove it
anno_table <- anno_table[-which(row.names(anno_table) %like% "4_DMM_R_3413"),]

# For right knees only
anno_table <- subset(anno_table, Side == "R")



# 3. Read datasets and QC RSD within batch ######################
# @3.1 Batch 1 ##################################################
# Already remove 4_DMM_R_3413 in csv
conc1 <- read.csv(file = "~/MGLeuvenAC_Batch1.csv",
                  header = TRUE, check.names = FALSE, na.strings = 0)

# Adjust format
row.names(conc1) <- conc1[,1]
conc1 <- conc1[,-1]
conc1 <- conc1[,-which(colnames(conc1) %like% "PA")]
conc1 <- conc1[,-which(colnames(conc1) %like% "PG")]
conc1 <- conc1[,-which(colnames(conc1) %like% "PI")]
conc1 <- conc1[,-which(colnames(conc1) %like% "PS")]

# Remove RSD QC > 25% within batch 1
qc1 <- conc1[which(row.names(conc1) %like% "QC"),]

RSD_1 <- apply(qc1, 2, function(x) sd(x)/mean(x))

length(which(RSD_1 > 0.25))

conc1_qc <- conc1[,-which(RSD_1 > 0.25)]


# @3.2 Batch 2 ###########################################
conc2 <- read.csv(file = "~/MGLeuvenAC_Batch2.csv",
                  header = TRUE, check.names = FALSE, na.strings = 0)

# Adjust format
row.names(conc2) <- conc2[,1]
conc2 <- conc2[,-1]
conc2 <- conc2[,-which(colnames(conc2) %like% "PA")]
conc2 <- conc2[,-which(colnames(conc2) %like% "PG")]
conc2 <- conc2[,-which(colnames(conc2) %like% "PI")]
conc2 <- conc2[,-which(colnames(conc2) %like% "PS")]

# Remove RSD QC > 25% within batch 2
qc2 <- conc2[which(row.names(conc2) %like% "QC"),]

RSD_2 <- apply(qc2, 2, function(x) sd(x)/mean(x))

length(which(RSD_2 > 0.25))

conc2_qc <- conc2[,-which(RSD_2 > 0.25)]



# 4.  Merge tables ################################################
a <- intersect(colnames(conc1_qc),colnames(conc2_qc))

conc_bind <- rbind.fill(conc1_qc, conc2_qc)
row.names(conc_bind) <- c(row.names(conc1_qc), row.names(conc2_qc))

conc_bind <- conc_bind[,a]

conc_bind <- conc_bind[-which(row.names(conc_bind) %like% "Blank"),]


# Lipid annotation update 2024.07.11 #################################
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


# For right knees only ####
conc_qc <- conc_qc[-which(row.names(conc_qc) %like% "L"),]



# 6.   Missing values ################################################
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



# 7.   Weight normalization and reorder #############################
conc_rpNA <- as.data.frame(conc_rpNA)

conc_rpNA$Batch <- c(rep(1,13),rep(2,16))

conc_norm1 <- conc_rpNA[match(row.names(anno_table),row.names(conc_rpNA)),]

anno_table$Batch <- conc_norm1$Batch
anno_table$Batch <- factor(anno_table$Batch)
anno_table_W12 <- subset(anno_table, Time == "W12")

conc_norm2 <- conc_norm1[,-ncol(conc_norm1)]

# nmol/mg=concLipidyzer/(40*MassTissue)
conc_norm2$weight <- anno_table$weight

conc_norm3 <- sweep(conc_norm2[,1:(ncol(conc_norm2)-1)], 1, unlist(conc_norm2$weight),  "/")

conc_norm4 <- conc_norm3/40

conc_norm <- conc_norm4[(row.names(conc_norm4) %like% "12_"),]

write.csv(conc_norm4, file = "Table S4_mouse articular cartilage lipidomics_time course.csv")



# 8. Lipid annotation ###################################################
# Create annotation table
anno_table_lpd <- data.frame(Class = word(colnames(conc_norm), 1, sep = "\\("),
                             row.names = colnames(conc_norm))

anno_table_lpd$Class[which(colnames(conc_norm) %like% "TAG")] <- "TAG"

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



# 9. Lipid number proportion #################################################
lpd_class_prop <- anno_table_lpd %>% count(Class)
lpd_class_prop$prop <- lpd_class_prop$n/ncol(conc_norm)*100
lpd_class_prop[order(lpd_class_prop$n, decreasing = TRUE),]



# 10. Technical variation ##################################################
techVar <- apply(conc_rmNA_qc, 2, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100)

techVar <- data.frame(var = techVar, class = str_sub(colnames(conc_rmNA_qc), 1,3),
                      row.names = colnames(conc_rmNA_qc))

ggplot(techVar, aes(x=class, y=var)) +
  geom_boxplot()  +
  geom_point() +
  ylim(c(0,100)) # check x axis tick labels

ggplot(techVar, aes(x = class, y = var)) +
  geom_boxplot(outlier.size = 1, color = "#848080", fill = "#E6E4E4") +
  scale_x_discrete(labels = c("CE","CER","DAG","DCER","FFA","HCER","LCER","LPC","LPE","PC","PE","SM","TAG"),
                   guide = guide_axis(angle = 45)) +
  ylim(c(0,100)) +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(size = 0.5), # line size = font size/22
        axis.ticks = element_line(size = 0.5),
        panel.background = element_rect(fill = "white")) +
  labs(x = "Lipid class", y = "Coefficient of variation (%)")

ggsave("technical variation.png", width = 11, height = 8, units = "cm", dpi = 300)



# 11. Multivariate analysis ######################################

# @11.1 PCA ######################################################
PCAdata <- log2(conc_norm)

# Add subset information
PCAdata$Time <- as.factor(anno_table_W12$Time)
PCAdata$Group <- as.factor(anno_table_W12$Group)
PCAdata$Batch <- as.factor(anno_table_W12$Batch)

# Calculate PC
PCA.PC <- prcomp(as.matrix(PCAdata[,1:ncol(conc_norm)]), scale. = TRUE)


autoplot(PCA.PC, data = PCAdata, colour = "Group",
         frame = TRUE, frame.type = "norm", frame.level = 0.95, size = 2) +
  geom_hline(yintercept=0,linetype=3, color = "grey") +
  geom_vline(xintercept=0,linetype=3, color = "grey") +
  scale_color_manual(values = c("#00BDC4", "#F8766D"), breaks = c("Sham","DMM")) +
  scale_fill_manual(values = c("#00BDC4", "#F8766D"), breaks = c("Sham","DMM")) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white"))

ggsave("pca_W12.png", dpi = 300, width = 1109, height = 818, units = "px")


# @11.2 Heatmap #################################################
HMdata <- PCAdata[,1:ncol(conc_norm)]

HMdata <- t(scale(HMdata, center = TRUE, scale = TRUE))

min(HMdata)
max(HMdata)

col_fun <- colorRamp2(c(-4,-2,0,2,4), c("#51446E","#937BB1","white","#EE9591","#F55C45"))

anno_table_W12$Group <- factor(anno_table_W12$Group, levels = c("Sham", "DMM"), ordered = TRUE)

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
                     "LPE" = "#987884",
                     "LPC" = "#DAE7E6",
                     "SM" = "#EEE48D",
                     "CER" = "#83EADD",
                     "DCER" = "#B295C7",
                     "HCER" = "#EB47A3",
                     "LCER" = "#F4CCE0")
  )
)

png(file="heatmap_W12.png", width = 1000, height = 1300, res = 300)
heatmap1 <- Heatmap(HMdata,
                    name = "Z-score",
                    col = col_fun,
                    cluster_columns = F,
                    column_split = anno_table_W12$Group,
                    column_gap = unit(2,"mm"),
                    border = T,
                    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_row_names = F,
                    show_column_names = F,
                    left_annotation = row_anno,
                    heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                                title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                direction = "vertical"))
draw(heatmap1)
dev.off()



# 12.  Differential analysis DMM vs.Sham #######################################
# For right knees only
FC_expression <- expression(bold(paste(Log[2], "FC (DMM/sham)")))
p_expression <- expression(paste(bold(-Log[10]), bold("(BH "), bolditalic("P"), bold("-value)")))

# Rebuttal: to show the kinetics of lipid profiling 2024.12.20 #################
# @ 12.1 Time point W2 #########################################################
df_W2 <- conc_norm4[grep("^2_", rownames(conc_norm4)),]

df_W2_log <- t(log2(df_W2))

design_W2 <- model.matrix(~0+factor(c(rep("DMM",5),rep("Sham",5))))
colnames(design_W2) <- levels(factor(c(rep("DMM",5),rep("Sham",5))))
row.names(design_W2) <- row.names(df_W2)

contrast.matrix_W2 <- makeContrasts(paste0(unique(c("DMM","Sham")), collapse = "-"), levels = design_W2)

fit_W2 <- lmFit(df_W2_log, design_W2)

fit2_W2 <- contrasts.fit(fit_W2, contrast.matrix_W2)
fit2_W2 <- eBayes(fit2_W2)

output_W2 <- topTable(fit2_W2, adjust.method="BH", number=1000)

logic1 <- (output_W2$logFC > 0 & output_W2$adj.P.Val< 0.05)
logic2 <- (output_W2$logFC < 0 & output_W2$adj.P.Val < 0.05)

output_W2$regulation <- ifelse(logic1,"Up-regulated",ifelse(logic2,"Down-regulated","Non-significant"))
output_W2$regulation <- factor(output_W2$regulation, levels = c("Non-significant", "Down-regulated", "Up-regulated", ordered = TRUE))


# Volcano plot
ggplot(output_W2,aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) + 
  geom_hline(yintercept = -log10(0.05), color = "#8A8989", linetype = 3) +
  geom_point(cex = 1, alpha = 0.6) +
  scale_color_manual(values = c("#C6C6C5","#6AB6E3","#E5423C")) + 
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        plot.title = element_text(size = 11, color = "black", face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white")) +
  labs(title = "W2", x = FC_expression, y = p_expression, color = "")

ggsave("volcano_W2.png", width = 1300, height = 900, units = "px", dpi = 300)


# @ 12.2 Time point W4 ###################################################
df_W4 <- conc_norm4[grep("^4_", rownames(conc_norm4)),]

df_W4_log <- t(log2(df_W4))

design_W4 <- model.matrix(~0+factor(c(rep("DMM",4),rep("Sham",5))))
colnames(design_W4) <- levels(factor(c(rep("DMM",4),rep("Sham",5))))
row.names(design_W4) <- row.names(df_W4)

contrast.matrix_W4 <- makeContrasts(paste0(unique(c("DMM","Sham")), collapse = "-"), levels = design_W4)

fit_W4 <- lmFit(df_W4_log, design_W4)

fit2_W4 <- contrasts.fit(fit_W4, contrast.matrix_W4)
fit2_W4 <- eBayes(fit2_W4)

output_W4 <- topTable(fit2_W4, adjust.method="BH", number=1000)

logic3 <- (output_W4$logFC > 0 & output_W4$adj.P.Val< 0.05)
logic4 <- (output_W4$logFC < 0 & output_W4$adj.P.Val < 0.05)

output_W4$regulation <- ifelse(logic3,"Up-regulated",ifelse(logic4,"Down-regulated","Non-significant"))
output_W4$regulation <- factor(output_W4$regulation, levels = c("Non-significant", "Down-regulated", "Up-regulated", ordered = TRUE))

# Volcano plot
ggplot(output_W4,aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) + 
  geom_hline(yintercept = -log10(0.05), color = "#8A8989", linetype = 3) +
  geom_point(cex = 1, alpha = 0.6) +
  scale_color_manual(values = c("#C6C6C5","#6AB6E3","#E5423C")) + 
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        plot.title = element_text(size = 11, color = "black", face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white")) +
  labs(title = "W4", x = FC_expression, y = p_expression, color = "")

ggsave("volcano_W4.png", width = 1300, height = 900, units = "px", dpi = 300)


# @ 12.3 Time point W12 ###################################################
df_W12 <- conc_norm

df_W12_log <- t(log2(df_W12))

design_W12 <- model.matrix(~0+factor(c(rep("DMM",5),rep("Sham",5))))
colnames(design_W12) <- levels(factor(c(rep("DMM",5),rep("Sham",5))))
row.names(design_W12) <- row.names(df_W12)

contrast.matrix_W12 <- makeContrasts(paste0(unique(c("DMM","Sham")), collapse = "-"), levels = design_W12)

fit_W12 <- lmFit(df_W12_log, design_W12)

fit2_W12 <- contrasts.fit(fit_W12, contrast.matrix_W12)
fit2_W12 <- eBayes(fit2_W12)

output_W12 <- topTable(fit2_W12, adjust.method="BH", number=1000)

logic5 <- (output_W12$logFC > 0 & output_W12$adj.P.Val< 0.05)
logic6 <- (output_W12$logFC < 0 & output_W12$adj.P.Val < 0.05)

output_W12$regulation <- ifelse(logic5,"Up-regulated", ifelse(logic6, "Down-regulated", "Non-significant"))
output_W12$regulation <- factor(output_W12$regulation, levels = c("Non-significant", "Down-regulated", "Up-regulated", ordered = TRUE))

down_W12 <- subset(output_W12, regulation == "Down-regulated")
down_W12 <- down_W12[order(down_W12$logFC),]

down_W12_top40 <- down_W12[order(down_W12$logFC),][1:40,]
write.csv(down_W12_top40,file = "~/down_W12_top40.csv")

up_W12 <- subset(output_W12,regulation == "Up-regulated")
up_W12 <- up_W12[rev(order(up_W12$logFC)),]

diff <- rbind(down_W12,up_W12)
diff <- diff[order(diff$logFC),]

write.csv(diff,file = "diff.csv")

# Volcano plot
ggplot(output_W12,aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) + 
  geom_hline(yintercept = -log10(0.05), color = "#8A8989", linetype = 3) +
  geom_point(cex = 1, alpha = 0.6) +
  scale_color_manual(values = c("#C6C6C5","#6AB6E3","#E5423C")) + 
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        plot.title = element_text(size = 11, color = "black", face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white")) +
  labs(title = "W12", x = FC_expression, y = p_expression, color = "")

ggsave("volcano_W12.png", width = 1300, height = 900, units = "px", dpi = 300)


# Bubble plot
bubble <- output_W12["logFC"]
bubble$adjp <- -log10(output_W12$adj.P.Val)

bubble <- bubble[match(row.names(anno_table_lpd),row.names(bubble)),]

bubble$Class <- anno_table_lpd$Class

ggplot(bubble, aes(x = logFC, y = Class)) + 
  geom_vline(xintercept = c(-1, 0, 1), color = "#8A8989", linetype = 3) +
  geom_point(aes(size = adjp, fill = Class), shape= 21, color = "black", stroke = 0.8, alpha = 0.7) +
  scale_size(range = c(1.5, 7), name = p_expression) +
  scale_fill_manual(values = pal_custom) +
  scale_y_discrete(limits = rev(levels(anno_table_lpd$Class))) +
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


# @ 12.4 Differentially regulated lipid proportion ###########################
down_class <- subset(output_W12,regulation == "Down-regulated")
up_class <- subset(output_W12,regulation == "Up-regulated")

down_class$Class <- str_sub(row.names(down_class), 1,3)
down_class %>% count(Class)

up_class$Class <- str_sub(row.names(up_class), 1,3)
up_class %>% count(Class)


# @ 12.5 Heatmap top regulated lipids ########################################
# Down top 40
HMdata_W12_down <- HMdata[(row.names(HMdata) %in% row.names(down_W12_top40)),]

HMdata_W12_down <- HMdata_W12_down[match(row.names(down_W12_top40),row.names(HMdata_W12_down)),]

min(-log10(down_W12_top40$adj.P.Val))
max(-log10(down_W12_top40$adj.P.Val))

# Heatmap
png(file="heatmap_W12_top40.png", width = 1000, height = 1400, res = 300)
heatmap2 <- Heatmap(HMdata_W12_down,
                    name = "Z-score",
                    col = col_fun,
                    cluster_columns = F,
                    cluster_rows = F,
                    column_split = anno_table_W12$Group,
                    column_gap = unit(1,"mm"),
                    border = T,
                    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_row_names = T,
                    show_column_names = F,
                    row_names_gp = gpar(fontsize = 11),
                    heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                                title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                direction = "horizontal"))

draw(heatmap2)
dev.off()

# Bar plot
down_W12_top40$FA <- row.names(down_W12_top40)

ggplot(down_W12_top40,aes(x = logFC, y = FA, fill = -log10(adj.P.Val))) +
  geom_bar(stat="identity",width=0.7) +
  scale_y_discrete(limits = rev(row.names(down_W12_top40)), position = "right") +
  scale_x_continuous(breaks = c(-1,0, 1)) +
  scale_fill_gradient(low = "#C9F3D4", high = "#35D35D", limits = c(1,4)) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  theme_test() +
  theme(axis.ticks = element_line(size = 0.5),
        legend.position="bottom",
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x= FC_expression, y="", fill = p_expression)

ggsave("bar_top.png", width = 659, height = 1615, units = "px", dpi = 300)


# Up
HMdata_W12_up <- HMdata[(row.names(HMdata) %in% row.names(up_W12)),]

HMdata_W12_up <- HMdata_W12_up[match(row.names(up_W12),row.names(HMdata_W12_up)),]

min(-log10(up_W12$adj.P.Val))
max(-log10(up_W12$adj.P.Val))

png(file="heatmap_W12_up.png", width = 1200, height = 300, res = 300)
heatmap2 <- Heatmap(HMdata_W12_up,
                    name = "Z-score",
                    col = col_fun,
                    cluster_columns = F,
                    cluster_rows = F,
                    column_split = anno_table_W12$Group,
                    column_gap = unit(1,"mm"),
                    border = T,
                    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_row_names = T,
                    show_column_names = F,
                    row_names_gp = gpar(fontsize = 11),
                    heatmap_legend_param = list(labels_gp = gpar(fontsize = 11),
                                                title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                direction = "horizontal"))

draw(heatmap2)
dev.off()

# Bar plot
up_W12$FA <- row.names(up_W12)

ggplot(up_W12,aes(x = logFC, y = FA, fill = -log10(adj.P.Val))) +
  geom_bar(stat="identity",width=0.7) +
  scale_y_discrete(limits = rev(row.names(up_W12)), position = "right") +
  scale_x_continuous(limits = c(0, 1.5), breaks = c(0,1)) +
  scale_fill_gradient(low = "#C9F3D4", high = "#35D35D", limits = c(1,4)) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  theme_test() +
  theme(axis.ticks = element_line(size = 0.5),
        legend.position="bottom",
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x= FC_expression, y="", fill = p_expression)

ggsave("bar_up.png", width = 700, height = 600, units = "px", dpi = 300)



# ATTENTION!!! Following analyses mirror human data #########################
# 13.   Deconstruct lipids ########################################################
# @ 13.1 Matrix data ############################
DAG <- df_W12 %>% select(starts_with("DAG"))
colnames(DAG) <- paste(colnames(DAG),"_1", sep = "")

PC <- df_W12 %>% select(starts_with("PC"))
colnames(PC) <- paste(colnames(PC),"_1", sep = "")

PE <- df_W12 %>% select(starts_with("PE"))
colnames(PE) <- paste(colnames(PE),"_1", sep = "")

df_W12_FA <- cbind(df_W12,DAG,PC,PE)


# @ 13.2 Differential data ################################
DAG_FC <- output_W12[which(row.names(output_W12) %like% "DAG"),]
row.names(DAG_FC) <- paste(row.names(DAG_FC),"_1", sep = "")

PC_FC <- output_W12[which(row.names(output_W12) %like% "PC"),]
PC_FC <- PC_FC[-which(row.names(PC_FC) %like% "LPC"),]
row.names(PC_FC) <- paste(row.names(PC_FC),"_1", sep = "")

PE_FC <- output_W12[which(row.names(output_W12) %like% "PE"),]
PE_FC <- PE_FC[-which(row.names(PE_FC) %like% "LPE"),]
row.names(PE_FC) <- paste(row.names(PE_FC),"_1", sep = "")

output_W12_FA <- rbind(output_W12,DAG_FC,PC_FC,PE_FC)


# @ 13.3 Fatty acyl annotation ######################################################
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

output_W12_FA <- merge(anno_table_lpd_FA, output_W12_FA, by = "row.names")
row.names(output_W12_FA) <- output_W12_FA$Row.names

output_W12_FA$FA <- paste("C", output_W12_FA$FA, sep = "")

min(output_W12_FA$logFC)
max(output_W12_FA$logFC)



# 14. Analysis on fatty acyls ######################################################

# @ 14.1 Certain fatty acyls #####################################################
# **** 14.1.1 FA in all the lipids ###########################################
ggplot(output_W12_FA, aes(x = logFC, y = FA)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.7, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.3) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = rev) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        legend.position="none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Fatty acid chain", x = FC_expression, color = p_expression)

ggsave("allFA_FC+p.png", width = 820, height = 1480, units = "px", dpi = 300)

# Legend
ggplot(output_W12_FA, aes(x = logFC, y = FA)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.7, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.3) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = rev) +
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


# **** 14.1.2 FA in significantly regulated lipids ##########################
output_W12_FA_sig <- subset(output_W12_FA, regulation == "Up-regulated" | regulation == "Down-regulated")

# ATTENTION!!!!!!!!
unique(output_W12_FA_sig$Class)

pal_custom_sig <- c(FFA = "#31992A",
                    CE = "#93D073", 
                    DAG = "#1F78B4",
                    TAG = "#88BBD9",
                    PE = "#E5423C",
                    PC = "#F39b7F",
                    LPE = "#987884",
                    LPC = "#DAE7E6",
                    SM = "#EEE48D",
                    DCER = "#B295C7",
                    LCER = "#F4CCE0")

ggplot(output_W12_FA_sig, aes(x = logFC, y = FA, color = Class)) + 
  geom_hline(yintercept = unique(output_W12_FA_sig$FA), linetype = 3, color = "grey") +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_jitter(aes(color = Class), size = 1.3, alpha = 0.7, height = 0.3) +
  scale_color_manual(values = pal_custom_sig) + 
  scale_y_discrete(limits = rev) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Fatty acid chain", color = "Lipid class")

ggsave("allFA_sig.png", width = 1080, height = 1480, units = "px", dpi = 300)


# **** 14.1.3 The number of individual FA #####################################################
count_FA <- output_W12_FA %>% group_by(regulation,Class) %>% count(FA)

sum(count_FA$n)

count_FA$Color <- as.character(count_FA$Class)

count_FA$Color[which(count_FA$regulation == "Non-significant")] <- "Non-significant"

count_FA$Color <- factor(count_FA$Color, levels = c("Non-significant","LCER","DCER","SM","LPC","LPE","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA, aes(x = n, y = FA, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E1DCDE")), limits = rev) +
  scale_y_discrete(limits = rev) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Fatty acid chain", fill = "Lipid class")

count_FA[which(count_FA$FA == "C18:1"),]

ggsave("allFA_n.png", width = 1200, height = 1480, units = "px", dpi = 300)



# @ 14.2 The length of fatty acyls #####################################################

# **** 14.2.1 FA in all the lipids #####################################################
# Group by carbon number
ggplot(output_W12_FA, aes(x = logFC, y = Length)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.4) +
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

ggsave("length_FC+p.png", width = 720, height = 1280, units = "px", dpi = 300)

# Group by length classification
ggplot(output_W12_FA, aes(x = logFC, y = LengthClass)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.4) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = c("VLCFA", "LCFA", "MCFA")) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Length of fatty acid chain", x = FC_expression, color = p_expression)

ggsave("lengthClass_FC+p.png", width = 820, height = 1040, units = "px", dpi = 300)



# **** 14.2.2 FA in significantly regulated lipids #####################################################
# Group by carbon number
ggplot(output_W12_FA_sig, aes(x = logFC, y = Length, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1, alpha = 0.7, height = 0.4) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = rev) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Carbon number of fatty acid chain", color = "Lipid class")

ggsave("length_sig.png", width = 1080, height = 1265, units = "px", dpi = 300)

# Group by length classification
ggplot(output_W12_FA_sig, aes(x = logFC, y = LengthClass, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1, alpha = 0.7, height = 0.5) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = c("VLCFA", "LCFA", "MCFA")) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Length of fatty acid chain", color = "Lipid class")

ggsave("lengthClass_sig.png", width = 1080, height = 1040, units = "px", dpi = 300)


# **** 14.2.3 The amount of individual FA #####################################################
# Group by carbon number
count_FA_length <- output_W12_FA %>% group_by(regulation,Class) %>% count(Length)

count_FA_length$Color <- as.character(count_FA_length$Class)

count_FA_length$Color[which(count_FA_length$regulation == "Non-significant")] <- "Non-significant"

count_FA_length$Color <- factor(count_FA_length$Color, levels = c("Non-significant","LCER","DCER","SM","LPC","LPE","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_length, aes(x = n, y = Length, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")), limits = rev) +
  scale_y_discrete(limits = rev) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Carbon number of fatty acid chain", fill = "Lipid class")

ggsave("length_n.png", width = 1130, height = 1280, units = "px", dpi = 300)


# Group by length classification
count_FA_lengthClass <- output_W12_FA %>% group_by(regulation,Class) %>% count(LengthClass)

count_FA_lengthClass$Color <- as.character(count_FA_lengthClass$Class)

count_FA_lengthClass$Color[which(count_FA_lengthClass$regulation == "Non-significant")] <- "Non-significant"

count_FA_lengthClass$Color <- factor(count_FA_lengthClass$Color, levels = c("Non-significant","LCER","DCER","SM","LPC","LPE","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_lengthClass, aes(x = n, y = LengthClass, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")), limits = rev) +
  scale_y_discrete(limits = c("VLCFA", "LCFA", "MCFA")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Length of fatty acid chain", fill = "Lipid class")

count_FA_lengthClass[which(count_FA_lengthClass$LengthClass == "LCFA"),]

ggsave("lengthClass_n.png", width = 1150, height = 1040, units = "px", dpi = 300)



# @ 14.3 The Unsaturation of fatty acyls #####################################################

# **** 14.3.1 FA in all the lipids #####################################################
# Group by bouble bond number
ggplot(output_W12_FA, aes(x = logFC, y = Unsaturation)) +
  geom_vline(xintercept = 0, linetype = 3, color = "black") +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.4) +
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

ggsave("sat_FC+p.png", width = 720, height = 1280, units = "px", dpi = 300)


# Group by saturation classification
ggplot(output_W12_FA, aes(x = logFC, y = UnsatClass)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_boxplot(color = "#8A8989", width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 1, aes(color = -log10(adj.P.Val)), height = 0.4) +
  scale_colour_gradient(low = "#EAE5E7", high = "red") +
  scale_y_discrete(limits = c("PUFA", "MUFA", "SFA")) +
  theme_classic() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(y = "Unsaturation of fatty acid chain", x = FC_expression)

ggsave("satClass_FC+p.png", width = 820, height = 1040, units = "px", dpi = 300)


# **** 14.3.2 FA in significantly regulated lipids #####################################################
# Group by double bond number
ggplot(output_W12_FA_sig, aes(x = logFC, y = Unsaturation, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1, alpha = 0.7, height = 0.4) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = rev) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Double bond number of fatty acid chain", color = "Lipid class")

ggsave("sat_sig.png", width = 1080, height = 1265, units = "px", dpi = 300)


# Group by saturation classification
ggplot(output_W12_FA_sig, aes(x = logFC, y = UnsatClass, color = Class)) + 
  geom_vline(xintercept = 0, linetype =3) +
  geom_jitter(aes(color = Class), size = 1, alpha = 0.7, height = 0.4) +
  scale_color_manual(values = pal_custom_sig) +
  scale_y_discrete(limits = c("PUFA", "MUFA", "SFA")) +
  theme_test() +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x = FC_expression, y = "Unsaturation of fatty acid chain", color = "Lipid class")

ggsave("satClass_sig.png", width = 1080, height = 1040, units = "px", dpi = 300)



# **** 14.3.3 The amount of individual FA #####################################################
# Group by double bond number
count_FA_Unsat <- output_W12_FA %>% group_by(regulation,Class) %>% count(Unsaturation)

count_FA_Unsat$Color <- as.character(count_FA_Unsat$Class)

count_FA_Unsat$Color[which(count_FA_Unsat$regulation == "Non-significant")] <- "Non-significant"

count_FA_Unsat$Color <- factor(count_FA_Unsat$Color, levels = c("Non-significant","LCER","DCER","SM","LPC","LPE","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_Unsat, aes(x = n, y = Unsaturation, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")), limits = rev) +
  scale_y_discrete(limits = rev) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Double bond number of fatty acid chain", fill = "Lipid class")

ggsave("Unsat_n.png", width = 1130, height = 1280, units = "px", dpi = 300)


# Group by saturation classification
count_FA_UnsatClass <- output_W12_FA %>% group_by(regulation,Class) %>% count(UnsatClass)

count_FA_UnsatClass$Color <- as.character(count_FA_UnsatClass$Class)

count_FA_UnsatClass$Color[which(count_FA_UnsatClass$regulation == "Non-significant")] <- "Non-significant"

count_FA_UnsatClass$Color <- factor(count_FA_UnsatClass$Color, levels = c("Non-significant","LCER","DCER","SM","LPC","LPE","PC","PE","TAG","DAG","CE","FFA"), ordered = TRUE)

ggplot(count_FA_UnsatClass, aes(x = n, y = UnsatClass, fill = Color)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c(pal_custom,c("Non-significant"="#E8E3E5")), limits = rev) +
  scale_y_discrete(limits = c("PUFA", "MUFA", "SFA")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black")) +
  labs(x="Count", y = "Unsaturation of fatty acid chain", fill = "Lipid class")

count_FA_UnsatClass %>% group_by(UnsatClass) %>% summarise(sum=sum(n))
count_FA_UnsatClass %>% group_by(UnsatClass, regulation) %>% summarise(sum=sum(n))

ggsave("UnsatClass_n.png", width = 1150, height = 1040, units = "px", dpi = 300)



# 15. Analysis on fatty acyls group by classes ######################################
ggplot(output_W12_FA, aes(x = Length, y = logFC)) +
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


ggplot(output_W12_FA, aes(x = Unsaturation, y = logFC)) +
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


# 16. Lipid class dynamics ############################################
# Rebuttal: to show the kinetics of lipid profiling 2024.12.20 ########

# @16.1 Time-course lipid classes #####################################
CE <- conc_norm4 %>% select(starts_with("CE("))
CER <- conc_norm4 %>% select(starts_with("CER"))
LCER <- conc_norm4 %>% select(starts_with("LCER"))
HCER <- conc_norm4 %>% select(starts_with("HCER"))
DCER <- conc_norm4 %>% select(starts_with("DCER"))
SM <- conc_norm4 %>% select(starts_with("SM"))
FFA <- conc_norm4 %>% select(starts_with("FFA"))
PC <- conc_norm4 %>% select(starts_with("PC"))
PE <- conc_norm4 %>% select(starts_with("PE"))
LPC <- conc_norm4 %>% select(starts_with("LPC"))
LPE <- conc_norm4 %>% select(starts_with("LPE"))
DAG <- conc_norm4 %>% select(starts_with("DAG"))
TAG <- conc_norm4 %>% select(starts_with("TAG"))

conc_class <- data.frame(FFA = rowSums(FFA),
                         CE = rowSums(CE),
                         DAG = rowSums(DAG),
                         TAG = rowSums(TAG),
                         PE = rowSums(PE),
                         PC = rowSums(PC),
                         LPE = rowSums(LPE),
                         LPC = rowSums(LPC),
                         SM = rowSums(SM),
                         CER = rowSums(CER),
                         DCER = rowSums(DCER),
                         HCER = rowSums(HCER),
                         LCER = rowSums(LCER))


conc_class$Time <- factor(anno_table$Time, levels = c("W2","W4","W12"), ordered = T)
conc_class$Group <- anno_table$Group

# Perform 2-way ANOVA in R
conc_class_t <- melt(conc_class)

aov_list1 <- list()
postoc_list1 <- list()
for (i in colnames(conc_class)[1:(ncol(conc_class)-2)]) {
  
  data <- subset(conc_class_t, variable == i)
  
  mod <- lm(value ~ Group * Time, data = data, contrasts = list(Group = contr.sum, Time = contr.sum))
  
  qq <- plot(mod, which = 2, main = i) #Normality
  
  var <- plot(mod, which = 3, main = i) #Homogeneity of variances
  
  mod.aov <- car::Anova(mod, type = 3)
  aov_list1[[i]] <- mod.aov
  
  posthoc <- emmeans(mod, list(pairwise ~ Group | Time), adjust = "none")
  postoc_summary <- as.data.frame(pairs(posthoc))
  postoc_list1[[i]] <- postoc_summary

}

aov_list1

postoc_list1

# Line chart
linechart <- conc_class_t
linechart$Time <- gsub("W", "", linechart$Time)
linechart$Time <- as.numeric(linechart$Time)
linechart$Group <- factor(linechart$Group, levels = c("Sham","DMM"), ordered = T)

linechart %>%
  group_by(variable, Time, Group) %>%
  summarise(mean = mean(value), sd = sd(value), n = n(), se = sd/sqrt(n)) %>%
  ggplot(aes(x = Time, y = mean, color = Group)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.6, linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = Group), size = 2) +
  facet_wrap(~variable, scales = "free", ncol = 7) +
  scale_color_manual(values = c("#00BDC4", "#F8766D")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = c(2,4,12), limits = c(1,13)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0,0.15))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x ="Time (w)", y = "Concentration (nmol/mg)")

ggsave("time_class.png", width = 3500, height = 850, units = "px", dpi = 300)

# Perform 2-way ANOVA in GraphPad
names <- rep(1:5,6)
names <- names[-15]

for (i in colnames(conc_class)[1:(ncol(conc_class)-2)]) {
  
  data1 <- subset(conc_class_t, variable == i)
  data1$Group <- paste(data1$Group, names, sep = "_")
  data2 <- dcast(data1, Time ~ Group)
  assign(paste0("conc_",i), data2)
  write.csv(data2, paste0(paste0("conc_",i), ".csv"), row.names = FALSE)

}


# @16.2 Time-course FA #####################################
conc_FA1 <- data.frame(`C12:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "12:0"]),
                       `C14:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "14:0"]),
                       `C14:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "14:1"]),
                       `C15:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "15:0"]),
                       `C16:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "16:0"]),
                       `C16:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "16:1"]),
                       `C17:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "17:0"]),
                       `C18:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:0"]),
                       `C18:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:1"]),
                       `C18:2` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:2"]),
                       `C18:3` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:3"]),
                       `C18:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:4"]),
                       `C20:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:0"]),
                       `C20:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:1"]),
                       `C20:2` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:2"]),
                       `C20:3` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:3"]),
                       `C20:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:4"]),
                       `C20:5` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:5"]),
                       `C22:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:0"]),
                       `C22:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:1"]),
                       `C22:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:4"]),
                       `C22:5` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:5"]),
                       `C22:6` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:6"]),
                       `C24:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "24:0"]),
                       `C24:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "24:1"]),
                       `C26:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "26:0"]),
                       `C26:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "26:1"]),
                       check.names = FALSE)

conc_FA2 <- data.frame(`C12:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "12:0/12:0"]),
                       `C14:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "14:0/14:0"]),
                       `C14:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "14:1/14:1"]),
                       `C15:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "15:0/15:0"]),
                       `C16:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "16:0/16:0"]),
                       `C16:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "16:1/16:1"]),
                       `C17:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "17:0/17:0"]),
                       `C18:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:0/18:0"]),
                       `C18:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:1/18:1"]),
                       `C18:2` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:2/18:2"]),
                       `C18:3` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:3/18:3"]),
                       `C18:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:4/18:4"]),
                       `C20:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:0/20:0"]),
                       `C20:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:1/20:1"]),
                       `C20:2` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:2/20:1"]),
                       `C20:3` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:3/20:3"]),
                       `C20:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:4/20:4"]),
                       `C20:5` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:5/20:5"]),
                       `C22:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:0/22:0"]),
                       `C22:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:1/22:1"]),
                       `C22:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:4/22:4"]),
                       `C22:5` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:5/22:5"]),
                       `C22:6` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:6/22:6"]),
                       `C24:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "24:0/24:0"]),
                       `C24:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "24:1/24:1"]),
                       `C26:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "26:0/26:0"]),
                       `C26:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "26:1/26:1"]),
                       check.names = FALSE)

conc_FA3 <- data.frame(`C12:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "12:0_12:0"]),
                       `C14:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "14:0_14:0"]),
                       `C14:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "14:1_14:1"]),
                       `C15:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "15:0_15:0"]),
                       `C16:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "16:0_16:0"]),
                       `C16:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "16:1_16:1"]),
                       `C17:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "17:0_17:0"]),
                       `C18:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:0_18:0"]),
                       `C18:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:1_18:1"]),
                       `C18:2` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:2_18:2"]),
                       `C18:3` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:3_18:3"]),
                       `C18:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "18:4_18:4"]),
                       `C20:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:0_20:0"]),
                       `C20:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:1_20:1"]),
                       `C20:2` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:2_20:1"]),
                       `C20:3` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:3_20:3"]),
                       `C20:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:4_20:4"]),
                       `C20:5` = rowSums(conc_norm4[colnames(conc_norm4) %like% "20:5_20:5"]),
                       `C22:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:0_22:0"]),
                       `C22:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:1_22:1"]),
                       `C22:4` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:4_22:4"]),
                       `C22:5` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:5_22:5"]),
                       `C22:6` = rowSums(conc_norm4[colnames(conc_norm4) %like% "22:6_22:6"]),
                       `C24:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "24:0_24:0"]),
                       `C24:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "24:1_24:1"]),
                       `C26:0` = rowSums(conc_norm4[colnames(conc_norm4) %like% "26:0_26:0"]),
                       `C26:1` = rowSums(conc_norm4[colnames(conc_norm4) %like% "26:1_26:1"]),
                       check.names = FALSE)

conc_FA <- conc_FA1+conc_FA2+conc_FA3

conc_FA$Time <- factor(anno_table$Time, levels = c("W2","W4","W12"), ordered = T)
conc_FA$Group <- anno_table$Group

# Perform 2-way ANOVA in R
conc_FA_t <- melt(conc_FA)

aov_list2 <- list()
postoc_list2 <- list()
for (i in colnames(conc_FA)[1:(ncol(conc_FA)-2)]) {
  
  data <- subset(conc_FA_t, variable == i)
  
  mod <- lm(value ~ Group * Time, data = data, contrasts = list(Group = contr.sum, Time = contr.sum))
  
  qq <- plot(mod, which = 2, main = i) #Normality
  
  var <- plot(mod, which = 3, main = i) #Homogeneity of variances
  
  mod.aov <- car::Anova(mod, type = 3)
  aov_list2[[i]] <- mod.aov
  
  posthoc <- emmeans(mod, list(pairwise ~ Group | Time), adjust = "none")
  postoc_summary <- as.data.frame(pairs(posthoc))
  postoc_summary$sig <- ifelse(postoc_summary$p.value < 0.0001, "****",
                               ifelse(postoc_summary$p.value < 0.001, "***",
                                      ifelse(postoc_summary$p.value < 0.01, "**",
                                             ifelse(postoc_summary$p.value < 0.05, "*", "ns"))))
  postoc_list2[[i]] <- postoc_summary
  
}

aov_list2

postoc_list2

# Line chart
linechart_FA <- conc_FA_t
linechart_FA$Time <- gsub("W", "", linechart_FA$Time)
linechart_FA$Time <- as.numeric(linechart_FA$Time)
linechart_FA$Group <- factor(linechart_FA$Group, levels = c("Sham","DMM"), ordered = T)

linechart_FA %>%
  group_by(variable, Time, Group) %>%
  summarise(mean = mean(value), sd = sd(value), n = n(), se = sd/sqrt(n)) %>%
  ggplot(aes(x = Time, y = mean, color = Group)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.6, linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = Group), size = 2) +
  facet_wrap(~variable, scales = "free", ncol = 7) +
  scale_color_manual(values = c("#00BDC4", "#F8766D")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = c(2,4,12), limits = c(1,13)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0,0.15))) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        panel.background = element_rect(fill = "white")) +
  labs(x ="Time (w)", y = "Concentration (nmol/mg)")

ggsave("time_FA.png", width = 3500, height = 1600, units = "px", dpi = 300)


# Perform 2-way ANOVA in GraphPad
names <- rep(1:5,6)
names <- names[-15]

conc_FA_t2 <- conc_FA_t
conc_FA_t2$variable <- gsub(":", ".", conc_FA_t2$variable)

for (i in unique(conc_FA_t2$variable)) {
  
  data1 <- subset(conc_FA_t2, variable == i)
  data1$Group <- paste(data1$Group, names, sep = "_")
  data2 <- dcast(data1, Time ~ Group)
  assign(paste0("conc_",i), data2)
  write.csv(data2, paste0(paste0("conc_",i), ".csv"), row.names = FALSE)
  
}
