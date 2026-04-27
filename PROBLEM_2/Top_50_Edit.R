suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(grid)
})

## -----------------------------
## 0) Sanity checks
## -----------------------------
stopifnot(exists("dds"))
stopifnot(inherits(dds, "DESeqDataSet"))

## -----------------------------
## 1) Get results (DESeq2)
## -----------------------------
res <- DESeq2::results(dds)
# Si necesitas fijar la dirección explícitamente:
# res <- DESeq2::results(dds, contrast = c("condition", "Mutante", "Control"))

## Quedarse con filas válidas (y opcionalmente DEGs)
res <- res[!is.na(res$log2FoldChange) & is.finite(res$log2FoldChange), ]
res <- res[!is.na(res$padj), ]

## (Opcional pero típico) quedarte solo con DEGs por padj:
## Si NO quieres filtrar por significancia, comenta esta línea.
res <- res[res$padj < 0.05, , drop = FALSE]

if (nrow(res) < 2) stop("No hay suficientes genes en 'res' tras el filtrado (padj/log2FC).")

## -----------------------------
## 2) Pick TOP 25 up + TOP 25 down by signed log2FC
## -----------------------------
n_each <- 25

res_pos <- res[order(res$log2FoldChange, decreasing = TRUE), , drop = FALSE]
res_neg <- res[order(res$log2FoldChange, decreasing = FALSE), , drop = FALSE]

res_pos <- head(res_pos, n_each)
res_neg <- head(res_neg, n_each)

res_top <- rbind(res_pos, res_neg)

## Ordenar para que queden arriba positivos y abajo negativos (opcional)
res_top <- res_top[order(res_top$log2FoldChange, decreasing = TRUE), , drop = FALSE]

## -----------------------------
## 3) Build expression matrix (vst) and subset ONLY those genes
## -----------------------------
vst_obj <- DESeq2::vst(dds, blind = TRUE)
mat <- assay(vst_obj)

## Intersección para asegurar match de IDs
genes_requested <- rownames(res_top)
genes_present <- intersect(genes_requested, rownames(mat))

if (length(genes_present) < length(genes_requested)) {
  missing <- setdiff(genes_requested, rownames(mat))
  message("Aviso: ", length(missing), " genes del top no están en assay(vst(dds)). Ejemplos: ",
          paste(utils::head(missing, 5), collapse = ", "))
}

## Reordenar exactamente en el orden de res_top
genes_present <- genes_requested[genes_requested %in% rownames(mat)]

mat_top <- mat[genes_present, , drop = FALSE]
res_top2 <- res_top[genes_present, , drop = FALSE]

## Si quieres EXACTAMENTE 50, obliga a tenerlos (si no, para y revisa IDs)
if (nrow(mat_top) != 50) {
  stop("Tras cruzar con la matriz de expresión, quedaron ", nrow(mat_top),
       " genes (no 50). Revisa que rownames(res) y rownames(assay(vst(dds))) sean el mismo tipo de ID.")
}

## Eliminar genes con sd=0 para evitar NaNs en Z-score
row_sd <- matrixStats::rowSds(mat_top)
if (any(row_sd == 0)) {
  zero_var <- sum(row_sd == 0)
  stop("Hay ", zero_var, " genes con varianza 0 entre muestras dentro del top50. ",
       "No se pueden escalar a Z-score; elimina esos genes o usa otra transformación.")
}

## Z-score por gen (por fila)
zmat <- t(scale(t(mat_top)))
zmat[is.na(zmat)] <- 0

## -----------------------------
## 4) Sample annotation
## -----------------------------
stopifnot("Condition" %in% colnames(colData(dds)))
condition <- factor(colData(dds)$Condition, levels = c("Control", "Mutante"))

## (Opcional) Control primero
ord <- order(condition)
zmat <- zmat[, ord, drop = FALSE]
condition <- condition[ord]

cond_cols <- c(Control = "red", Mutante = "#1F51FF")
top_ha <- HeatmapAnnotation(
  groups = condition,
  col = list(groups = cond_cols),
  annotation_name_side = "right"
)

## -----------------------------
## 5) Right-side Log2FC bar
## -----------------------------
lfc <- res_top2$log2FoldChange
names(lfc) <- rownames(res_top2)

lfc_lim <- max(abs(lfc), na.rm = TRUE)
lfc_col_fun <- circlize::colorRamp2(
  c(-lfc_lim, 0, lfc_lim),
  c("#1F51FF", "white", "red")
)
lfc_mat <- matrix(lfc, ncol = 1, dimnames = list(names(lfc), "Log2FC"))

## -----------------------------
## 6) Plot (TODO a la derecha: Log2FC + nombres)
## -----------------------------
z_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#1F51FF", "white", "red"))

# Heatmap principal SIN nombres (los pondremos como anotación a la derecha del todo)
ht_main <- Heatmap(
  zmat,
  name = "Z-score",
  col = z_col_fun,
  top_annotation = top_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,      # <- clave
  show_column_names = TRUE,
  column_names_rot = 90
)

# Columna de Log2FC (a la derecha del heatmap principal)
ht_lfc <- Heatmap(
  lfc_mat,
  name = "Log2 Fold Change",
  col = lfc_col_fun,
  width = unit(10, "mm"),
  cluster_rows = FALSE,        # se alinea con el orden de ht_main en el dibujo combinado
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  rect_gp = gpar(col = NA),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", lfc_mat[i, j]), x, y, gp = gpar(fontsize = 7))
  }
)

# Nombres de genes como anotación de texto AL FINAL (extremo derecho)
gene_labels <- rownames(zmat)

ht_genes <- rowAnnotation(
  Gene = anno_text(
    gene_labels,
    which = "row",
    just = "left",
    location = 0,
    gp = gpar(fontsize = 9)
  ),
  width = max_text_width(gene_labels, gp = gpar(fontsize = 9)) + unit(2, "mm")
)

# Orden: heatmap -> Log2FC -> nombres
draw(
  ht_main + ht_lfc + ht_genes,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)

