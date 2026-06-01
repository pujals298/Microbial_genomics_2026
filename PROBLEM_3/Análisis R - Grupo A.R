# ==============================================================================
# PROBLEMA 3
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CARGA DE LIBRERÍAS Y CONFIGURACIÓN ESTÉTICA
# ------------------------------------------------------------------------------
library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyverse)  
library(pheatmap)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "white", color = "black"),
                  plot.title = element_text(face = "bold", hjust = 0.5)))

# ------------------------------------------------------------------------------
# 2. DATOS Y CONSOLIDACIÓN DE LA MATRIZ (9 MUESTRAS)
# ------------------------------------------------------------------------------
setwd("C:/Users/Constanza Pia/Downloads/Problema 3")

archivos <- c(
  "STIN"   = "A - Influx/Influx_2_SRR7497167.tabular", 
  "SWHIN"  = "A - Influx/Influx_3_SRR7500591.tabular", 
  "STLIN"  = "A - Influx/Influx_1_SRR7618093.tabular", 
  
  "STAS"   = "B - Activated Sludge (AS)/AS_3_SRR7497945.tabular",  
  "SWHAS"  = "B - Activated Sludge (AS)/AS_2_SRR7503007.tabular",  
  "STLAS"  = "B - Activated Sludge (AS)/AS_1_SRR7627523.tabular",  
  
  "STEFF"  = "C - Efflux/Efflux_3_SRR7499264.tabular", 
  "SWHEFF" = "C - Efflux/Efflux_2_SRR7503035.tabular", 
  "STLEFF" = "C - Efflux/Efflux_1_SRR7630418.tabular"
)

leer_metagenomica <- function(ruta_archivo, nombre_muestra) {
  primera_linea <- colnames(read_delim(ruta_archivo, delim = "\t", n_max = 1, show_col_types = FALSE))
  
  if ("name" %in% primera_linea) {
    read_delim(ruta_archivo, delim = "\t", show_col_types = FALSE) %>%
      select(name, new_est_reads) %>%
      rename(!!nombre_muestra := new_est_reads, Especie = name)
  } else {
    read_delim(ruta_archivo, delim = "\t", col_names = FALSE, show_col_types = FALSE) %>%
      select(X6, X2) %>% 
      mutate(X6 = trimws(X6)) %>%
      rename(!!nombre_muestra := X2, Especie = X6)
  }
}

lista_tablas <- imap(archivos, ~leer_metagenomica(.x, .y))

tabla_abundancia <- lista_tablas %>% 
  reduce(full_join, by = "Especie") %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  group_by(Especie) %>%
  summarise(across(everything(), sum), .groups = 'drop')

palabras_a_eliminar <- c("Homo", "Primates", "Mammalia", "Metazoa", "root", 
                         "cellular organisms", "Hominidae", "Homininae", 
                         "Hominoidea", "Catarrhini", "Simiiformes", "Eutheria", 
                         "Theria", "Vertebrata", "Chordata", "Craniata", 
                         "Bilateria", "Eukaryota", "Amniota", "Boreoeutheria",
                         "Deuterostomia", "Dipnotetrapodomorpha", "Euarchontoglires",
                         "Eumetazoa", "Euteleostomi", "Gnathostomata", "Haplorrhini",
                         "Opisthokonta", "Sarcopterygii", "Teleostomi", "Tetrapoda")

tabla_abundancia_filtrada <- tabla_abundancia %>%
  filter(!str_detect(Especie, paste(palabras_a_eliminar, collapse = "|")))

matriz_comunidad <- tabla_abundancia_filtrada %>%
  column_to_rownames(var = "Especie") %>%
  as.matrix()

# ------------------------------------------------------------------------------
# 3. CREACIÓN DE METADATOS INTEGRALES Y OBJETO PHYLOSEQ
# ------------------------------------------------------------------------------
metadatos_reales <- data.frame(
  Nicho = rep(c("Influx", "AS", "Efflux"), each = 3),
  Planta = rep(c("ST", "SWH", "STL"), times = 3),
  row.names = colnames(matriz_comunidad)
)

OTU <- otu_table(matriz_comunidad, taxa_are_rows = TRUE)
SAMPLES <- sample_data(metadatos_reales)
pseq <- phyloseq(OTU, SAMPLES)

# ------------------------------------------------------------------------------
# 4. ALFA DIVERSIDAD TOTAL CON SIGNIFICANCIA 
# ------------------------------------------------------------------------------
df_alpha <- estimate_richness(pseq, measures = c("Observed", "Shannon", "Simpson"))
df_alpha <- cbind(df_alpha, sample_data(pseq)) 

df_alpha_long <- df_alpha %>%
  gather(key = "Metrica", value = "Valor", Observed, Shannon, Simpson)

df_alpha_long$Nicho <- factor(as.character(df_alpha_long$Nicho), levels = c("Influx", "AS", "Efflux"))

print(kruskal.test(Shannon ~ Nicho, data = df_alpha))
print(kruskal.test(Shannon ~ Planta, data = df_alpha))

# --- GRÁFICO ALFA 1: ORGANIZADO POR NICHO ---
comparaciones_nicho_reales <- list(c("Influx", "AS"), c("AS", "Efflux"), c("Influx", "Efflux"))

plot_significancia_nicho_ordenado <- ggplot(df_alpha_long, aes(x = Nicho, y = Valor, fill = Nicho)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.05, size = 3.5, aes(shape = Planta), color = "black", stroke = 0.8) +
  facet_wrap(~Metrica, scales = "free_y") +
  scale_fill_manual(values = c("Influx" = "#E76F51", "AS" = "#2A9D8F", "Efflux" = "#E9C46A")) + 
  stat_compare_means(comparisons = comparaciones_nicho_reales, 
                     method = "wilcox.test", 
                     label = "p.signif") + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold")) +
  labs(title = "Alpha Diversity Across Ecological Niches", 
       x = "Ecological Niche", 
       y = "Index Value")

print(plot_significancia_nicho)
ggsave("alpha_diversity_niche.png", plot_significancia_nicho_ordenado, width = 11, height = 6, dpi = 300)

# --- GRÁFICO ALFA 2: ORGANIZADO POR PLANTA DE ORIGEN ---
comparaciones_plantas <- list(c("ST", "STL"), c("ST", "SWH"), c("STL", "SWH"))

plot_significancia_planta <- ggplot(df_alpha_long, aes(x = Planta, y = Valor, fill = Planta)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.08, size = 3.5, aes(shape = Nicho), color = "black", stroke = 0.8) +
  facet_wrap(~Metrica, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  stat_compare_means(comparisons = comparaciones_plantas, 
                     method = "wilcox.test", 
                     label = "p.signif") + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(title = "Alpha Diversity Across Wastewater Plants", 
       x = "Wastewater Treatment Plant", 
       y = "Index Value")

print(plot_significancia_planta)
ggsave("alpha_diversity_plant.png", plot_significancia_planta, width = 11, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 5. BETA DIVERSIDAD GLOBAL: ORDENACIÓN Y PERMANOVA COMPLETA (*)
# ------------------------------------------------------------------------------
pseq_rel <- transform_sample_counts(pseq, function(x) x / sum(x))
pcoa_bray_gala <- ordinate(pseq_rel, method = "PCoA", distance = "bray")

plot_beta_gala <- plot_ordination(pseq_rel, pcoa_bray_gala, color = "Planta", shape = "Nicho") +
  
  geom_point(size = 5, stroke = 1.2, alpha = 0.9) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Principal Coordinate Analysis (PCoA) of Beta Diversity",
       x = paste0("PCoA 1 [", round(pcoa_bray_gala$values$Relative_eig[1]*100, 1), "%]"),
       y = paste0("PCoA 2 [", round(pcoa_bray_gala$values$Relative_eig[2]*100, 1), "%]"),
       color = "Wastewater Plant",
       shape = "Ecological Niche") +
  
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "grey40"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

print(plot_beta)

ggsave("beta_diversity_EN.png", plot_beta_gala, width = 8, height = 5.5, dpi = 300)

# --- ANÁLISIS ESTADÍSTICO DE ADONIS ---
dist_bray <- distance(pseq_rel, method = "bray")
df_metadatos_9 <- data.frame(sample_data(pseq_rel))

permanova_9_muestras <- adonis2(dist_bray ~ Nicho * Planta, data = df_metadatos_9, permutations = 999)
print(permanova_9_muestras)

permanova_9_muestras_ok <- adonis2(dist_bray ~ Nicho + Planta, data = df_metadatos_9, permutations = 999)
print(permanova_9_muestras_ok)

# ==============================================================================
# MATCH DISTANCE 
# ==============================================================================

dist_matrix <- as.matrix(dist_bray)
metadatos <- data.frame(sample_data(pseq_rel)) %>% 
  rownames_to_column(var = "Sample_ID")

muestras <- colnames(dist_matrix)
pares <- t(combn(muestras, 2))

df_match <- data.frame(
  Muestra1 = pares[,1],
  Muestra2 = pares[,2],
  value = sapply(1:nrow(pares), function(i) dist_matrix[pares[i,1], pares[i,2]])
)

df_match <- df_match %>%
  left_join(metadatos %>% select(Sample_ID, Nicho), by = c("Muestra1" = "Sample_ID")) %>%
  rename(Nicho1 = Nicho) %>%
  left_join(metadatos %>% select(Sample_ID, Nicho), by = c("Muestra2" = "Sample_ID")) %>%
  rename(Nicho2 = Nicho) %>%
  mutate(Tipo_Comparacion = if_else(Nicho1 == Nicho2, 
                                    "Same Niche\n(Intra-group)", 
                                    "Different Niches\n(Inter-group)")) %>%
  mutate(Tipo_Comparacion = factor(Tipo_Comparacion, 
                                   levels = c("Same Niche\n(Intra-group)", 
                                              "Different Niches\n(Inter-group)")))

comparacion_match <- list(c("Same Niche\n(Intra-group)", "Different Niches\n(Inter-group)"))

plot_match_palo <- ggplot(df_match, aes(x = Tipo_Comparacion, y = value, fill = Tipo_Comparacion)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "grey30", width = 0.35, size = 0.7) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', 
               dotsize = 0.7,       
               color = "grey20",    
               alpha = 0.8,         
               stackratio = 1.1) +        
  scale_fill_manual(values = c("Same Niche\n(Intra-group)" = "#99B898", 
                               "Different Niches\n(Inter-group)" = "#FECEAB")) + 
  stat_compare_means(comparisons = comparacion_match,
                     method = "wilcox.test", 
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")),
                     size = 5, 
                     color = "grey20",
                     tip_length = 0.02) + 
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2), expand = c(0,0)) +
  theme_classic(base_size = 14) + 
  theme(axis.text.x = element_text(face = "bold", color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 12)),
        axis.line = element_line(color = "grey40"),
        legend.position = "none") +
  labs(x = "Ecological Comparison Type",
       y = "Bray-Curtis Dissimilarity (0 to 1)")

print(plot_match)
ggsave("match_distance_EN.png", plot_match_palo, width = 6.5, height = 5.5, dpi = 300)

# ==============================================================================
# PIPELINE GLOBAL-- Diseño poster
# ==============================================================================
df_global_crudo <- psmelt(pseq_rel)

top_30_bacterias_reales <- df_global_crudo %>%
  group_by(OTU) %>%
  summarise(Abundancia_Media = mean(Abundance)) %>%
  slice_max(order_by = Abundancia_Media, n = 30) %>% 
  pull(OTU)

df_global_top30_puro <- df_global_crudo %>%
  filter(OTU %in% top_30_bacterias_reales) %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  mutate(Especie_Final = factor(OTU))

df_global_top30_puro$Nicho <- factor(df_global_top30_puro$Nicho, levels = c("Influx", "AS", "Efflux"))

paleta_limpia_poster <- c(
  "#1F77B4", "#3498DB", "#5DADE2", "#AED6F1", 
  "#2ECC71", "#27AE60", "#117A65", "#16A085", 
  "#E74C3C", "#C0392B", "#962D22", "#D98880", 
  "#9B59B6", "#8E44AD", "#BB8FCE", "#D2B4DE", 
  "#F1C40F", "#F39C12", "#F5B041", "#F9E79F", 
  "#34495E", "#ff4d6d", "#84a98c", "#bb8588", 
  "#E67E22", "#b56576", "#EDBB99", "#a7c957", 
  "#1ABC9C", "#48C9B0"                        
)

plot_global_bacterias_puras <- ggplot(df_global_top30_puro, aes(x = Sample, y = Abundance, fill = Especie_Final)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.15, alpha = 0.95) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_fill_manual(values = paleta_limpia_poster) + 
  facet_wrap(~Nicho, scales = "free_x") + 
  labs(x = "Wastewater Treatment Plant Samples", 
       y = "Relative Abundance (%)",
       fill = "Top 30 Microbial Taxa") +
  theme_classic(base_size = 14) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
        axis.text.x = element_text(face = "bold", color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "grey50"),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        legend.text = element_text(size = 8, face = "italic"),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 2, byrow = FALSE))

print(plot_global_bacterias_puras)
ggsave("global_microbial.png", plot_global_bacterias_puras, width = 16, height = 10, dpi = 300)