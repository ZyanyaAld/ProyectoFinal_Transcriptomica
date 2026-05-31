# =========================================================
# RNA-seq Differential Expression Analysis with DESeq2
# STAR + featureCounts
# Batch correction included
# =========================================================

# =========================
# Cargar librerías
# =========================

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(matrixStats)
library(limma)
library(ggforce)

# =========================
# Crear directorios
# =========================

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/deseq2/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/deseq2/tables", recursive = TRUE, showWarnings = FALSE)

# =========================
# Leer metadata
# =========================

metadata <- read.csv(
  "data/metadata/sample_info.csv",
  row.names = 1
)

# =========================
# Renombrar condiciones
# =========================

metadata$condition <- factor(
  metadata$condition,
  levels = c("Control", "TreatmentA", "TreatmentB"),
  labels = c(
    "GC",
    "A1G",
    "MG"
  )
)

# Batch / replicate
metadata$replicate <- factor(metadata$replicate)

# Etiquetas completas
condition_labels <- c(
  GC = "Ground Control",
  A1G = "Artificial Earth-gravity",
  MG = "Microgravity"
)

# Revisar metadata
print(metadata)

# =========================
# Leer matriz de counts
# =========================

counts <- read.delim(
  "results/featureCounts/star_counts.txt",
  comment.char = "#"
)

# Extraer columnas de conteos
countdata <- counts[, 7:ncol(counts)]

# Asignar nombres
colnames(countdata) <- rownames(metadata)

# Asignar gene IDs
rownames(countdata) <- counts$Geneid

# =========================
# Verificar consistencia
# =========================

all(colnames(countdata) == rownames(metadata))

# =========================
# Crear objeto DESeq2
# =========================

dds <- DESeqDataSetFromMatrix(
  countData = round(countdata),
  colData = metadata,
  design = ~ replicate + condition
)

# =========================
# Filtrar genes poco expresados
# =========================

dds <- dds[rowSums(counts(dds)) > 10, ]

# =========================
# Correr DESeq2
# =========================

dds <- DESeq(dds)

# =========================
# Umbrales de significancia
# =========================

FDR <- 0.05
LFC <- 0.5

# =========================
# Función de filtrado
# =========================

get_sig_genes <- function(res, FDR, LFC) {
  res <- as.data.frame(res)

  res$significant <- !is.na(res$padj) & res$padj < FDR
  res$up <- res$significant & res$log2FoldChange > LFC
  res$down <- res$significant & res$log2FoldChange < -LFC

  return(res)
}

# =========================
# Aplicar a tus resultados
# =========================

res_A1G_vs_GC <- get_sig_genes(res_A1G_vs_GC, FDR, LFC)
res_MG_vs_GC <- get_sig_genes(res_MG_vs_GC, FDR, LFC)
res_MG_vs_A1G <- get_sig_genes(res_MG_vs_A1G, FDR, LFC)

# =========================
# Leer anotación genes
# =========================

gene_annot <- read.delim(
  "data/reference/mm39_gene_names.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(gene_annot) <- c(
  "gene_id",
  "gene_name"
)

# quitar versiones
gene_annot$gene_id <- sub(
  "\\..*",
  "",
  gene_annot$gene_id
)

# =========================
# Función anotación
# =========================

annotate_results <- function(
  res_df,
  gene_annot
) {
  res_df$gene_id <- rownames(res_df)

  # quitar versiones
  res_df$gene_id <- sub(
    "\\..*",
    "",
    res_df$gene_id
  )

  # unir anotación
  res_df <- left_join(
    res_df,
    gene_annot,
    by = "gene_id"
  )

  # usar símbolo si existe
  res_df$label <- ifelse(
    is.na(res_df$gene_name),
    res_df$gene_id,
    res_df$gene_name
  )

  return(res_df)
}

# =========================
# Anotar resultados
# =========================

res_A1G_annot <- annotate_results(
  as.data.frame(res_A1G_vs_GC),
  gene_annot
)

res_MG_annot <- annotate_results(
  as.data.frame(res_MG_vs_GC),
  gene_annot
)

res_MGA1G_annot <- annotate_results(
  as.data.frame(res_MG_vs_A1G),
  gene_annot
)

# =========================
# Guardar tablas
# =========================

write.csv(
  res_A1G_annot,
  "results/deseq2/tables/A1G_vs_GC.csv",
  row.names = FALSE
)

write.csv(
  res_MG_annot,
  "results/deseq2/tables/MG_vs_GC.csv",
  row.names = FALSE
)

write.csv(
  res_MGA1G_annot,
  "results/deseq2/tables/MG_vs_A1G.csv",
  row.names = FALSE
)

# =========================
# Transformación VST
# =========================

vsd <- vst(
  dds,
  blind = FALSE
)

# =========================
# Batch correction SOLO plots
# =========================

assay(vsd) <- removeBatchEffect(
  assay(vsd),
  batch = metadata$replicate
)

# =========================
# Heatmap distancia muestras
# =========================

sampleDists <- dist(
  t(assay(vsd))
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

rownames(sampleDistMatrix) <- rownames(metadata)
colnames(sampleDistMatrix) <- rownames(metadata)

# PDF
pdf(
  "results/deseq2/plots/sample_distance_heatmap.pdf",
  width = 8,
  height = 8
)

pheatmap(sampleDistMatrix)

dev.off()

# PNG
png(
  "results/deseq2/plots/sample_distance_heatmap.png",
  width = 1200,
  height = 1200,
  res = 200
)

pheatmap(sampleDistMatrix)

dev.off()

# =========================
# PCA normal
# =========================

pcaData <- plotPCA(
  vsd,
  intgroup = c(
    "condition",
    "replicate"
  ),
  returnData = TRUE
)

percentVar <- round(
  100 *
    attr(
      pcaData,
      "percentVar"
    )
)

p1 <- ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    color = condition,
    shape = replicate
  )
) +

  geom_point(size = 5) +

  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +

  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +

  labs(
    color = "Condition",
    shape = "Replicate"
  ) +

  theme_bw(base_size = 14)

ggsave(
  "results/deseq2/plots/PCA_plot.pdf",
  plot = p1,
  width = 7,
  height = 5
)

ggsave(
  "results/deseq2/plots/PCA_plot.png",
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

# =========================
# PCA grupos
# =========================

p2 <- ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    color = condition
  )
) +

  geom_mark_ellipse(
    aes(
      fill = condition,
      label = condition
    ),
    alpha = 0.15,
    show.legend = FALSE
  ) +

  geom_point(size = 5) +

  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +

  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +

  labs(
    color = "Condition"
  ) +

  theme_bw(base_size = 14)

ggsave(
  "results/deseq2/plots/PCA_plot_groups.pdf",
  plot = p2,
  width = 7,
  height = 5
)

ggsave(
  "results/deseq2/plots/PCA_plot_groups.png",
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

# =========================
# Heatmap genes variables
# =========================

topVarGenes <- head(
  order(
    rowVars(assay(vsd)),
    decreasing = TRUE
  ),
  30
)

mat <- assay(vsd)[
  topVarGenes,
]

# Obtener gene IDs
gene_ids <- rownames(mat)

# Quitar versiones
gene_ids_clean <- sub(
  "\\..*",
  "",
  gene_ids
)

# Match símbolos
gene_symbols <- gene_annot$gene_name[
  match(
    gene_ids_clean,
    gene_annot$gene_id
  )
]

# Reemplazar NAs
gene_symbols[
  is.na(gene_symbols)
] <- gene_ids_clean[
  is.na(gene_symbols)
]

# Asignar símbolos
rownames(mat) <- gene_symbols

# Escalar filas
mat <- t(
  scale(
    t(mat)
  )
)

annotation_col <- data.frame(
  condition = metadata$condition
)

rownames(annotation_col) <- rownames(metadata)

# PDF
pdf(
  "results/deseq2/plots/heatmap_top30.pdf",
  width = 8,
  height = 10
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 8
)

dev.off()

# PNG
png(
  "results/deseq2/plots/heatmap_top30.png",
  width = 1200,
  height = 1500,
  res = 200
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 8
)

dev.off()

# =========================================================
# EXPORTAR LISTAS PARA GO / STRING / GSEA
# =========================================================

dir.create(
  "results/deseq2/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================
# Función para exportar
# =========================

export_enrichment_files <- function(
  res_df,
  comparison_name
) {
  # quitar NAs
  res_df <- res_df[!is.na(res_df$padj), ]

  # =========================
  # genes significativos
  # =========================

  sig <- res_df %>%
    filter(padj < FDR)

  # =========================
  # upregulated
  # =========================

  up <- sig %>%
    filter(log2FoldChange > LFC)

  # =========================
  # downregulated
  # =========================

  down <- sig %>%
    filter(log2FoldChange < LFC)

  # =========================
  # background
  # =========================

  background <- res_df %>%
    filter(!is.na(label))

  # =====================================================
  # EXPORTAR PARA DAVID / STRING
  # =====================================================

  write.table(
    up$label,
    file = paste0(
      "results/deseq2/enrichment/",
      comparison_name,
      "_UP_genes.txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  write.table(
    down$label,
    file = paste0(
      "results/deseq2/enrichment/",
      comparison_name,
      "_DOWN_genes.txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  write.table(
    background$label,
    file = paste0(
      "results/deseq2/enrichment/",
      comparison_name,
      "_BACKGROUND.txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # =====================================================
  # EXPORTAR RANKING PARA GSEA
  # =====================================================

  gsea_rank <- res_df %>%
    filter(!is.na(stat)) %>%
    select(label, stat) %>%
    distinct() %>%
    arrange(desc(stat))

  write.table(
    gsea_rank,
    file = paste0(
      "results/deseq2/enrichment/",
      comparison_name,
      "_GSEA_rank.rnk"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  cat(
    "\nArchivos exportados para:",
    comparison_name,
    "\n"
  )
}

# =========================================================
# EXPORTAR TODOS LOS COMPARATIVOS
# =========================================================

export_enrichment_files(res_A1G_annot, "A1G_vs_GC")
export_enrichment_files(res_MG_annot, "MG_vs_GC")
export_enrichment_files(res_MGA1G_annot, "MG_vs_A1G")

cat("\n=========== ENRICHMENT FILES DONE ===========\n")

# =========================
# Resumen
# =========================

cat(
  "\n=========== RESUMEN ===========\n"
)

cat(
  "\nGenes significativos A1G vs GC:",
  sum(
    res_A1G_annot$padj < FDR,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Genes significativos MG vs GC:",
  sum(
    res_MG_annot$padj < FDR,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Genes significativos MG vs A1G:",
  sum(
    res_MGA1G_annot$padj < FDR,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\n=========== DONE ===========\n"
)
