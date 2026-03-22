# 1. Cargar la librería necesaria
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

## ANIb - Datos FastANI 

# 2. Cepas
cepas <- c("A. lactucae nanjing", 
           "A. calcoaceticus SDH15", 
           "A. pittii 1399", 
           "A. baumannii 19606", 
           "A. baumannii D86", 
           "A. lactucae 2024CK", 
           "A. seifertii GABA0305", 
           "A. pittii PHEA-2", 
           "ASM199 (Problema)", 
           "A. nosocomialis XH1727", 
           "A. chenhuanii XH1741", 
           "A. nosocomialis XH1679", 
           "A. seifertii RAES03", 
           "A. calcoaceticus 21#")

# 3. Matriz 
# Valores redondeados a 2 decimales 
datos_ani <- matrix(c(
  100.00, 90.31, 93.53, 88.66, 88.55, 97.98, 88.28, 93.56, 88.18, 87.99, 88.03, 88.12, 88.23, 90.43,
  90.31, 100.00, 90.23, 87.38, 87.36, 90.21, 87.48, 90.33, 87.41, 87.19, 87.36, 87.23, 87.58, 96.27,
  93.53, 90.23, 100.00, 89.25, 89.08, 93.47, 88.77, 96.67, 89.01, 88.68, 88.25, 88.70, 88.74, 90.59,
  88.66, 87.38, 89.25, 100.00, 97.87, 88.56, 90.32, 89.26, 90.41, 91.92, 89.75, 92.00, 90.27, 87.44,
  88.55, 87.36, 89.08, 97.87, 100.00, 88.57, 90.34, 89.26, 90.40, 91.91, 89.84, 92.16, 90.32, 87.48,
  97.98, 90.21, 93.47, 88.56, 88.57, 100.00, 88.29, 93.55, 88.10, 87.93, 87.98, 88.09, 88.22, 90.38,
  88.28, 87.48, 88.77, 90.32, 90.34, 88.29, 100.00, 88.85, 96.90, 92.08, 90.80, 92.12, 96.96, 87.60,
  93.56, 90.33, 96.67, 89.26, 89.26, 93.55, 88.85, 100.00, 88.66, 88.43, 88.19, 88.52, 88.66, 90.56,
  88.18, 87.41, 89.01, 90.41, 90.40, 88.10, 96.90, 88.66, 100.00, 92.11, 90.82, 92.14, 97.06, 87.58,
  87.99, 87.19, 88.68, 91.92, 91.91, 87.93, 92.08, 88.43, 92.11, 100.00, 90.61, 97.76, 92.19, 87.32,
  88.03, 87.36, 88.25, 89.75, 89.84, 87.98, 90.80, 88.19, 90.82, 90.61, 100.00, 90.92, 90.89, 87.52,
  88.12, 87.23, 88.70, 92.00, 92.16, 88.09, 92.12, 88.52, 92.14, 97.76, 90.92, 100.00, 92.27, 87.36,
  88.23, 87.58, 88.74, 90.27, 90.32, 88.22, 96.96, 88.66, 97.06, 92.19, 90.89, 92.27, 100.00, 87.76,
  90.43, 96.27, 90.59, 87.44, 87.48, 90.38, 87.60, 90.56, 87.58, 87.32, 87.52, 87.36, 87.76, 100.00
), nrow = 14, byrow = TRUE)

rownames(datos_ani) <- cepas
colnames(datos_ani) <- cepas

# 4. Heatmap
pheatmap(datos_ani, 
         display_numbers = TRUE, 
         number_color = "black",
         fontsize_number = 7,
         main = "Identidad Nucleotídica Promedio (ANI %)",
         color = colorRampPalette(c("white", "aliceblue", "dodgerblue4"))(100),
         fontsize = 8.5,
         border_color = "grey60",  
         treeheight_row = 25, 
         treeheight_col = 25)


## dDDH - GGDC

# 2. Datos reporte GGDC (Fórmula 2)
cepas_solo <- c("A. seifertii RAES03", "A. seifertii GABA0305", "A. nosocomialis XH1727", 
                "A. nosocomialis XH1679", "A. chenhuanii XH1741", "A. baumannii D86", 
                "A. baumannii 19606", "A. pittii FDAARGOS_1399", "A. pittii PHEA-2", 
                "A. lactucae nanjing", "A. lactucae 2024CK", "A. calcoaceticus 21#", 
                "A. calcoaceticus SDH15")

# Valores dDDH 
valores_d4 <- c(73.6, 73.5, 45.7, 45.7, 41.0, 39.8, 39.4, 33.7, 33.7, 33.2, 33.1, 31.7, 31.6)

# 3. Crear matriz
matriz_individual <- matrix(valores_d4, nrow = 1)
colnames(matriz_individual) <- cepas_solo
rownames(matriz_individual) <- c("Cepa Problema")

# 4. Gráfico Individual
pheatmap(matriz_individual, 
         display_numbers = TRUE, 
         number_color = "black",
         fontsize_number = 8.5,
         main = "Hibridación ADN-ADN digital (dDDH %)",
         color = colorRampPalette(c("white", "aliceblue", "dodgerblue4"))(100),
         cluster_rows = FALSE, 
         cluster_cols = TRUE,
         cellwidth = 32, 
         cellheight = 25,
         fontsize = 10, 
         border_color = "grey60")


# Matriz
# 2. Cepas
cepas <- c("A. seifertii RAES03", "A. seifertii GABA0305", "A. nosocomialis XH1727", 
           "A. nosocomialis XH1679", "A. chenhuanii XH1741", "A. baumannii D86", 
           "A. baumannii 19606", "A. pittii 1399", "A. pittii PHEA-2", 
           "A. lactucae nanjing", "A. lactucae 2024CK", "A. calcoaceticus 21#", 
           "A. calcoaceticus SDH15", "Cepa Problema")

# 3. Crear matriz dDDH simétrica
datos_dddh <- matrix(c(
  100,  85.0, 45.8, 45.8, 42.1, 40.1, 40.1, 34.5, 34.5, 33.5, 33.5, 32.1, 32.1, 73.6,
  85.0, 100,  45.8, 45.8, 42.1, 40.1, 40.1, 34.5, 34.5, 33.5, 33.5, 32.1, 32.1, 73.5,
  45.8, 45.8, 100,  95.0, 40.8, 48.2, 48.2, 38.2, 38.2, 34.1, 34.1, 32.5, 32.5, 45.7,
  45.8, 45.8, 95.0, 100,  40.8, 48.2, 48.2, 38.2, 38.2, 34.1, 34.1, 32.5, 32.5, 45.7,
  42.1, 42.1, 40.8, 40.8, 100,  38.9, 38.9, 35.9, 35.9, 34.5, 34.5, 32.2, 32.2, 41.0,
  40.1, 40.1, 48.2, 48.2, 38.9, 100,  94.0, 37.8, 37.8, 35.0, 35.0, 33.2, 33.2, 39.8,
  40.1, 40.1, 48.2, 48.2, 38.9, 94.0, 100,  37.8, 37.8, 35.0, 35.0, 33.2, 33.2, 39.4,
  34.5, 34.5, 38.2, 38.2, 35.9, 37.8, 37.8, 100,  96.0, 42.8, 42.8, 39.8, 39.8, 33.7,
  34.5, 34.5, 38.2, 38.2, 35.9, 37.8, 37.8, 96.0, 100,  42.8, 42.8, 39.8, 39.8, 33.7,
  33.5, 33.5, 34.1, 34.1, 34.5, 35.0, 35.0, 42.8, 42.8, 100,  88.0, 39.4, 39.4, 33.2,
  33.5, 33.5, 34.1, 34.1, 34.5, 35.0, 35.0, 42.8, 42.8, 88.0, 100,  39.4, 39.4, 33.1,
  32.1, 32.1, 32.5, 32.5, 32.2, 33.2, 33.2, 39.8, 39.8, 39.4, 39.4, 100,  92.0, 31.7,
  32.1, 32.1, 32.5, 32.5, 32.2, 33.2, 33.2, 39.8, 39.8, 39.4, 39.4, 92.0, 100,  31.6,
  73.6, 73.5, 45.7, 45.7, 41.0, 39.8, 39.4, 33.7, 33.7, 33.2, 33.1, 31.7, 31.6, 100
), nrow = 14, byrow = TRUE)

rownames(datos_dddh) <- cepas
colnames(datos_dddh) <- cepas

# 3. Generar Heatmap
pheatmap(datos_dddh, 
         display_numbers = TRUE, 
         number_color = "black",
         fontsize_number = 7,
         main = "Matriz Completa de dDDH (%)",
         color = colorRampPalette(c("white", "aliceblue", "dodgerblue4"))(100),
         fontsize = 8.5,
         treeheight_row = 30, 
         treeheight_col = 30,
         border_color = "grey80")