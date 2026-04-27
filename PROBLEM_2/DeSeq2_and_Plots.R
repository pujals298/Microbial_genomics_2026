########################################
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

########################################
# 2. Directorio
setwd("C:/Users/Results featureCount")

#########################################

### Análisis de Expresión Diferencial (RNA-seq).

# 3. Leer los archivos individuales

# Explicación:

## Nota: featureCounts genera archivos con encabezado: Geneid + cepa 
## read.table: herramienta para importar datos de archivos de texto plano.
## header = TRUE: coloca a "Geneid" y "Cepas" como titulo y no parte de los datos.
## sep = "\t": separador, para decirle a R dónde termina una columna y empieza otra.
## stringsAsFactors=F: para que R no los lea como factores (y nos los agrupe), sino
## que lo vea como etiquetas únicas --> así evitamos problemas de renombrar o filtrar.

####

pao2 <- read.table("PAO2.tsv", header=T, sep="\t", stringsAsFactors=F)
View(pao2)

pao3 <- read.table("PAO3.tsv", header=T, sep="\t", stringsAsFactors=F)
View(pao3)

purm2 <- read.table("purM2.tsv", header=T, sep="\t", stringsAsFactors=F)
View(purm2)

purm3 <- read.table("purM3.tsv", header=T, sep="\t", stringsAsFactors=F)
View(purm3)

#########################################
# 4. Unir todas las tablas en una sola matriz maestra usando el Geneid como ancla

# Explicación:

## 1. Cada archivo tiene más de 5.500 filas de genes de P.aeruginosa, por lo que 
## inner_join actúa como un escáner; busca el Geneid (por ejemplo, PA0001) en
## los cuatro archivos y se asegura de que los datos de esa fila correspondan
## exactamente al mismo gen en todas las muestras.

## 2. El problema dice "RNA samples were extracted from duplicate cultures just 
## at exponential phase to create stranded RNAseq libraries", por lo que se 
## extrajeron de cultivos por duplicados (PAO2/PAO3 y purM2/purM3).
## Para que un software estadístico como DESeq2 pueda calcular si es un cambio es
## real o solo un ruido experimental, necesita ver los datos de las réplicas uno
## al lado del otro en una sola tabla (matriz).

## 3.Al usar inner_join, R solo conserva los genes que están presentes en los 
## cuatro archivos. Si un gen falta en alguna muestra, se elimina automáticamente 
## de la matriz final. Esto evita que el análisis falle más adelante por culpa de 
## datos incompletos o valores "NA" (vacíos). --> no es nuestro caso:).

####

count_matrix <- pao2 %>%
  inner_join(pao3, by="Geneid") %>%
  inner_join(purm2, by="Geneid") %>%
  inner_join(purm3, by="Geneid")

View(count_matrix)

#########################################
# 5. Convertir el Geneid en los nombres de las filas (rownames)

# Explicación:
## Esto es obligatorio para que DESeq2 funcione correctamente.

####

rownames(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL  # Eliminamos la columna repetida

View(count_matrix)


#########################################
# 6. Tabla de metadatos

# Explicación:

## 1. Asignación de Condiciones Experimentales
## La matriz de conteos tiene cuatro columnas: PAO2, PAO3, purM2 y purM3. 
## La tabla de metadatos indicará formalmente a DESeq2 que: PAO2 y PAO3 son 
## réplicas del grupo Control (Wild-Type), y que purM2 y purM3 son réplicas 
## del grupo Mutante.

## 2. Permitir la Comparación Estadística 
## Para responder si un proceso celular está "controlado" por el gen mutado 
## (Pregunta 3), R necesita calcular una división matemática: 
## (Expresión en Mutante/Expresión en Control).
## La tabla de metadatos define cuál es el numerador y cuál es el denominador 
## de esa comparación. Sin esta tabla, el software no sabría contra qué comparar 
## los datos.

## 3. Manejo de las Réplicas Biológicas
## El diseño del experimento incluye duplicados para cada condición. 
## La tabla de metadatos le dice a R: "Estos dos archivos son repeticiones del 
## mismo fenómeno".Esto permite que el algoritmo diferencie entre el ruido natural 
## (pequeñas diferencias entre PAO2 y PAO3) y el efecto real de la mutación 
## (diferencias grandes entre el bloque Control y el bloque Mutante).

## 4. Organización para Gráficos (Pregunta 2)
## Para el Volcano Plot o el PCA, el programa usará la tabla de metadatos para 
## saber de qué color pintar cada punto. Por ejemplo, pintará de azul los controles 
## y de rojo las mutantes, permitiéndote ver visualmente si los grupos se separan 
## correctamente.

####

# Definimos los grupos en el mismo orden que las columnas de la matriz
condiciones <- factor(c("Control", "Control", "Mutante", "Mutante"))

# Creamos la tabla de metadatos
meta_data <- data.frame(condition = condiciones)
rownames(meta_data) <- colnames(count_matrix)

print(meta_data)
View(meta_data)


#########################################
# 7. DESeq2

# Explicación:

## Función DESeqDataSetFromMatrix es un constructor: toma piezas sueltas y 
## las ensambla en un objeto complejo llamado dds (DESeqDataSet).

## 1. countData = count_matrix: tabla con los números puros de featureCounts.
## DESeq2 necesita los conteos brutos (raw counts) para aplicar su modelo 
## estadístico de distribución de Poisson. No usa porcentajes ni promedios, 
## sino los datos reales de cuántas lecturas (reads) se mapearon a cada gen.

## 2. colData = meta_data: tabla de metadatos donde dice que PAO2 es "Control" 
## y purM2 es "Mutante".
## Sin esto, R tiene números pero no sabe qué significan. 
## Esta parte le dice al programa: "Agrupa estas columnas porque son réplicas 
## del mismo fenómeno biológico".

## 3. design = ~ condition: fórmula del diseño experimental.
## Aquí se indica a R qué es lo que quieres comparar. Al poner ~ condition, le 
## estás diciendo: "Quiero que calcules la diferencia de expresión basándote 
## exclusivamente en la condición (Control vs. Mutante)".

## Hacemos esto para crear un contenedor único (el objeto dds) que guarda:
# Los datos (cuántos reads hay).
# La estructura (qué muestra es qué).
# La intención (qué vamos a comparar).

##### Análisis estadístico**

# Cuando se corre DESeq(dds), R realiza tres procesos principales:

## 1. Estimación de factores de tamaño (Normalización): 
# No todas las librerías de RNA-seq tienen el mismo número total de lecturas. 
# R calcula un factor para "equilibrar" las muestras y que puedas comparar la 
# réplica 2 con la 3 de forma justa.

## 2. Estimación de la dispersión: 
# Mide qué tanto varían los genes entre las réplicas del mismo grupo 
# (por ejemplo, qué tan diferentes son PAO2 y PAO3). Esto es vital para saber 
# si un cambio es biológicamente real o solo ruido del experimento.

## 3. Ajuste de un Modelo Lineal Generalizado: 
# Utiliza una prueba estadística (llamada Test de Wald) para calcular el p-valor 
# y determinar si la diferencia entre la mutante purM y la cepa silvestre es 
# significativa.

##!! Importante destacar que DESeq2 utiliza una Distribución Binominal Negativa.
# (Evolución de la distribución de Poisson).

# Naturaleza de los Conteos:Los datos de RNA-seq son conteos de eventos discretos 
# (cuántas veces una lectura "cayó" en un gen). La distribución de Poisson es la 
# forma matemática estándar para modelar cosas que se cuentan.

# El problema de la varianza: En una distribución de Poisson pura, se asume que 
# la media es igual a la varianza (Mean = Variance). Sin embargo, en biología, 
# los genes suelen variar mucho más de lo esperado --> sobredispersión.

# La solución (Binomial Negativa): DESeq2 usa este modelo porque añade un parámetro 
# extra para manejar esa variabilidad biológica típica de (micro)organismos.

#!! Todo esto para decir que:
# Sin este modelo estadístico, no se podría asegurar que la caída en los genes 
# de pigmento o biofilm se debe a la mutación en purM y no simplemente a que una 
# de las muestras salió "diferente" por azar. Este análisis da el p-valor ajustado 
# (padj), que es la validación científica para afirmar que los resultados son válidos.


#####


# Crear el objeto de DESeq2 uniendo los conteos y los metadatos
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta_data,
                              design = ~ condition)

# Correr el análisis estadístico**
dds <- DESeq(dds)

# Extraer los resultados finales
res <- results(dds)
head(res)

# Convertir el objeto de resultados (res) en un data.frame (tabla normal)
res_df <- as.data.frame(res)
View(res_df)

######
# Extraer los conteos ya normalizados (corregidos por tamaño de librería)
conteos_normalizados <- counts(dds, normalized=TRUE)

# Ver la tabla (Aquí SÍ verás PAO2, PAO3, purM2 y purM3)
View(conteos_normalizados)

#########################################
# 8. Métricas de DESeq2 y significancia

# baseMean: Es el promedio de los conteos normalizados de todas las muestras 
# (WT y Mutante). Dice qué tan "expresado" está el gen en general; si es muy 
# bajo, el gen casi no se detecta.

# log2FoldChange (LFC): Es la medida del cambio. Dice cuánto cambió la expresión 
# en la mutante purM respecto al Control (WT). Se usa una escala logarítmica de 
# base 2 porque es fácil de interpretar: 
## un LFC de 1 significa que la expresión se duplicó (2^1 = 2)
## un LFC de 2 significa que es 4 veces mayor (2^ = 4).

# lfcSE: Es el error estándar del log2FoldChange. Indica qué tanta incertidumbre 
# hay en ese cálculo.

# stat: Es el valor del estadístico (Test de Wald). Es la base para calcular 
# el p-valor.

# pvalue: La probabilidad de que el cambio observado sea por puro azar.

# padj (P-valor ajustado): Es el más importante. Al analizar miles de genes a 
# la vez, aumenta la probabilidad de encontrar "falsos positivos". 
# El padj corrige este error (usando el método Benjamini-Hochberg). 
# Se debe usar este para decidir si un gen es significativo.

#####


# Añadimos una columna que marque si el gen es significativo (padj < 0.05)
# y si sube (UP) o baja (DOWN)
res_df$diff_status <- "Not significant"
res_df$diff_status[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Overexpressed"
res_df$diff_status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Underexpressed"

View(res_df)

# ¿"Overexpressed" y "Underexpressed"?
# Esta clasificación es una forma de simplificar la respuesta biológica de la 
# bacteria ante la mutación en purM. 

# Se agregaron condiciones para "marcar" un gen como importante: que sea 
# estadísticamente confiable (padj < 0.05) y que el cambio sea relevante 
# (más del doble o menos de la mitad, |log2FC| > 1).

# Sobreexpresado: log2FoldChange > 1.
# Significa que la bacteria está produciendo más ARNm de ese gen en la mutante 
# que en la WT.
# Interpretación Biológica: Son genes que la bacteria "enciende" para 
# intentar compensar la falta de purinas o para responder al estrés de la mutación.

# Subexpresado: log2FoldChange < -1.
# Significa que la bacteria está produciendo menos ARNm de ese gen en la mutante.
# Interpretación Biológica: Si los genes de virulencia, pigmentos o biofilm están 
# subexpresados --> explicación de por qué la mutante purM perdió esas capacidades.



#########################################
# 9. PCA

# Explicación:

## PCA (Análisis de Componentes Principales): ver si las réplicas se parecen 
## entre sí y si la mutante es realmente distinta al control.

####

# Transformación para estabilizar la varianza (vst)
# Esto hace que los datos sean comparables y "suaviza" las diferencias extremas
vsd <- vst(dds, blind=FALSE)

# Generar el gráfico PCA

##All (5678 geneid)
plotPCA(vsd, intgroup="condition", ntop=nrow(vsd)) + 
  theme_minimal() + 
  geom_text_repel(aes(label=colnames(count_matrix)),  
                  show.legend = FALSE, 
                  max.overlaps = Inf, 
                  vjust=1.5,
                  size=3) +
  ggtitle("PCA: Wild-Type vs Mutante purM all")

##500*
plotPCA(vsd, intgroup="condition") + 
  theme_minimal() + 
  geom_text_repel(aes(label=colnames(count_matrix)),  
                  show.legend = FALSE, 
                  max.overlaps = Inf, 
                  vjust=1.5,
                  size=3) +
  ggtitle("PCA: Wild-Type vs Mutante purM")

####*

# Los 500* son los 500 genes con mayor varianza.
# Aunque el PCA está diseñado para encontrar los ejes de máxima variación, 
# incluir todos los genes suele ser contraproducente debido a la alta dimensionalidad 
# y el ruido técnico de los datos. 
# https://medium.com/@gunkurnia/is-pca-the-best-feature-selection-method-when-you-dont-know-your-data-60c827548983
#! Por lo que creo que es mejor escoger los 500.

#### Resultados de PCA:

## 1. ¿Qué es PC1 y PC2? Son los Componentes Principales. 
# Resúmenes de los datos:
# PC1 (Eje Horizontal): Es el cambio más grande que existe en el experimento. 
# En el gráfico, separa perfectamente al Control de la Mutante.

# PC2 (Eje Vertical): Es el segundo cambio más grande. 
# Separa un poco las réplicas entre sí (por ejemplo, PAO2 de PAO3).

## 2. ¿Qué significa la Varianza (88% y 11%)?
# La varianza es una medida de cuánta información o diferencia hay en cada eje.
# PC1: 88% variance:Significa que el 88% de todas las diferencias que R encontró 
# en los 5678 genes se explica por cosa: ser Control o ser Mutante.
#!! La mutación en purM es tan potente que domina casi por completo el 
#!! comportamiento de la bacteria.

# PC2: 11% variance: Solo el 11% de las diferencias se deben a otros factores 
# (pequeños errores de laboratorio, diferencias en el crecimiento de ese día, etc.).


######################################### POR MODIFICAR
# 9. Volcanoplot

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diff_status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("Underexpressed" = "blue", "Overexpressed" = "red", "Not significant" = "grey")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + # Líneas de Fold Change
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Línea de p-valor
  labs(title = "Volcano Plot: purM Mutant vs Wild-Type",
       x = "Log2 Fold Change",
       y = "-Log10 P-adj",
       color = "Estado") +
  geom_text_repel(data = head(res_df[order(res_df$padj), ], 10), 
                  aes(label = rownames(head(res_df[order(res_df$padj), ], 10))),
                  color = "black", size = 3)

top10_genes <- res_df %>%
  filter(!is.na(padj)) %>% 
  arrange(padj) %>%
  head(10)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diff_status)) +
  geom_point(alpha = 0.3, size = 1.5) + 
  scale_color_manual(values = c("Underexpressed" = "blue", 
                                "Overexpressed" = "red", 
                                "Not significant" = "grey85"),
                     labels = c("No significativo", "Sobreexpresado", "Subexpresado")) + 
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey60") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey60") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Perfil Transcriptómico: purM vs Wild-Type",
       x = "Log2 Fold Change",
       y = "-Log10 P-adj",
       color = "Estado del Gen") +
  geom_text_repel(data = top10_genes, 
                  aes(label = rownames(top10_genes)),
                  color = "black", 
                  size = 3.2, 
                  max.overlaps = Inf,
                  force = 3, 
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = NA) + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank()) 
