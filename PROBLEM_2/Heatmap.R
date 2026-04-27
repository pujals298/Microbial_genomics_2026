if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(ggrepel)
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)
if (!require("writexl")) install.packages("writexl")
library(writexl)
library(RColorBrewer)

Condiciones <- factor(c("Control", "Control", "Mutante", "Mutante"))

count_matrix <- PAO2 %>%
  inner_join(PAO3, by="Geneid") %>%
  inner_join(purM2, by="Geneid") %>%
  inner_join(purM3, by="Geneid")

View(count_matrix)

rownames(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL  

View(count_matrix)

Condiciones <- factor(c("Control", "Control", "Mutante", "Mutante"))

meta_data <- data.frame(Condition = Condiciones)
rownames(meta_data) <- colnames(count_matrix)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta_data,
                              design = ~ Condition)

dds <- DESeq(dds)

counts_norm <- counts(dds, normalized=TRUE)

gene_list <- genes_COG_classificats$`Locus Tag`
counts_degs <- counts_norm[rownames(counts_norm) %in% gene_list, ]
nrow(counts_degs)

counts_log <- log2(counts_degs + 1)
counts_scaled <- t(scale(t(counts_log)))

png("heatmap_DEGs.png", width = 4000, height = 3000, res = 150)
pheatmap(counts_scaled,
         fontsize = 3,
         scale = "none",           # ya escalamos manualmente arriba
         cluster_rows = TRUE,      # agrupa genes similares
         cluster_cols = TRUE,      # agrupa muestras similares
         show_rownames = TRUE,     # muestra nombres de genes
         show_colnames = TRUE,     # muestra nombres de muestras
         fontsize_row = 8,         # tamaño de letra de genes
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap de genes DEG")

dev.off()
