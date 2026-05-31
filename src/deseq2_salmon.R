# =========================================================
# RNA-seq Differential Expression Analysis
# Salmon + tximport + DESeq2
# =========================================================

# =========================
# Librerías
# =========================

library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(matrixStats)
library(limma)
library(ggforce)

# =========================
# Directorios
# =========================

dir.create(
  "results/deseq2_salmon",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/deseq2_salmon/plots",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/deseq2_salmon/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================
# Metadata
# =========================

metadata <- read.csv(
  "data/metadata/sample_info.csv",
  row.names = 1
)

# Renombrar grupos
metadata$condition <- recode(
  metadata$condition,
  "Control" = "GC",
  "TreatmentA" = "A1G",
  "TreatmentB" = "MG"
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("GC", "A1G", "MG")
)

# IMPORTANTE
metadata$replicate <- factor(
  metadata$replicate
)

print(metadata)

# =========================
# Archivos Salmon
# =========================

samples <- rownames(metadata)

files <- file.path(
  "results/salmon",
  samples,
  "quant.sf"
)

names(files) <- samples

print(files)

# =========================
# tx2gene
# =========================

tx2gene <- read.delim(
  "data/reference/tx2gene_clean.tsv"
)

# =========================
# tximport
# =========================

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  txOut = FALSE
)

# =========================
# DESeq2 object
# =========================

dds <- DESeqDataSetFromTximport(
  txi,
  colData = metadata,
  design = ~ replicate + condition
)

# =========================
# Filtrar genes
# =========================

dds <- dds[
  rowSums(counts(dds)) > 10,
]

# =========================
# DESeq2
# =========================

dds <- DESeq(dds)

# =========================
# Umbrales (DEG criteria)
# =========================

FDR <- 0.05
LFC <- 0.5

# =========================
# Contrastes
# =========================

res_A_vs_GC <- results(
  dds,
  contrast = c(
    "condition",
    "A1G",
    "GC"
  )
)

res_MG_vs_GC <- results(
  dds,
  contrast = c(
    "condition",
    "MG",
    "GC"
  )
)

res_MG_vs_A1G <- results(
  dds,
  contrast = c(
    "condition",
    "MG",
    "A1G"
  )
)

# =========================
# Anotación genes
# =========================

gene_annot <- read.delim(
  "data/reference/mm39_gene_names.tsv",
  header = FALSE
)

colnames(gene_annot) <- c(
  "gene_id",
  "gene_name"
)

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

  res_df$gene_id <- sub(
    "\\..*",
    "",
    res_df$gene_id
  )

  res_df <- left_join(
    res_df,
    gene_annot,
    by = "gene_id"
  )

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

res_A_annot <- annotate_results(
  as.data.frame(res_A_vs_GC),
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
  res_A_annot,
  "results/deseq2_salmon/tables/A1G_vs_GC.csv",
  row.names = FALSE
)

write.csv(
  res_MG_annot,
  "results/deseq2_salmon/tables/MG_vs_GC.csv",
  row.names = FALSE
)

write.csv(
  res_MGA1G_annot,
  "results/deseq2_salmon/tables/MG_vs_A1G.csv",
  row.names = FALSE
)

# =========================================================
# VST
# =========================================================

vsd <- vst(
  dds,
  blind = FALSE
)

# =========================================================
# Batch correction SOLO plots
# =========================================================

assay(vsd) <- removeBatchEffect(
  assay(vsd),
  batch = metadata$replicate
)

# =========================================================
# PCA
# =========================================================

pcaData <- plotPCA(
  vsd,
  intgroup = c(
    "condition",
    "replicate"
  ),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

# =========================================================
# PCA normal
# =========================================================

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

  theme_bw(base_size = 14)

ggsave(
  "results/deseq2_salmon/plots/PCA_plot.png",
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

# =========================================================
# PCA elipses
# =========================================================

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

  theme_bw(base_size = 14)

ggsave(
  "results/deseq2_salmon/plots/PCA_plot_groups.png",
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

# =========================================================
# Heatmap
# =========================================================

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

gene_ids <- rownames(mat)

gene_ids_clean <- sub(
  "\\..*",
  "",
  gene_ids
)

gene_symbols <- gene_annot$gene_name[
  match(
    gene_ids_clean,
    gene_annot$gene_id
  )
]

gene_symbols[
  is.na(gene_symbols)
] <- gene_ids_clean[
  is.na(gene_symbols)
]

rownames(mat) <- gene_symbols

mat <- t(scale(t(mat)))

annotation_col <- data.frame(
  condition = metadata$condition
)

rownames(annotation_col) <- rownames(metadata)

png(
  "results/deseq2_salmon/plots/heatmap_top30.png",
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
# Enrichment files
# =========================================================

dir.create(
  "results/deseq2_salmon/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

export_enrichment_files <- function(
  res_df,
  comparison_name,
  FDR = 0.05,
  LFC = 0.5
) {
  res_df <- res_df[!is.na(res_df$padj), ]

  sig <- res_df %>%
    filter(padj < FDR)

  up <- sig %>%
    filter(log2FoldChange > LFC)

  down <- sig %>%
    filter(log2FoldChange < -LFC)

  background <- res_df %>%
    filter(!is.na(label))

  write.table(
    up$label,
    file = paste0(
      "results/deseq2_salmon/enrichment/",
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
      "results/deseq2_salmon/enrichment/",
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
      "results/deseq2_salmon/enrichment/",
      comparison_name,
      "_BACKGROUND.txt"
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  gsea_rank <- res_df %>%
    filter(!is.na(stat)) %>%
    select(label, stat) %>%
    distinct() %>%
    arrange(desc(stat))

  write.table(
    gsea_rank,
    file = paste0(
      "results/deseq2_salmon/enrichment/",
      comparison_name,
      "_GSEA_rank.rnk"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}

# =========================================================
# Exportar enrichment
# =========================================================

export_enrichment_files(res_A_annot, "A1G_vs_GC", FDR, LFC)
export_enrichment_files(res_MG_annot, "MG_vs_GC", FDR, LFC)
export_enrichment_files(res_MGA1G_annot, "MG_vs_A1G", FDR, LFC)

cat("\n=========== DONE ===========\n")
