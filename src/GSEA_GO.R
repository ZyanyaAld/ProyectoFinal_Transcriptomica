# ============================================================
# GSEA y GO enrichment - Proyecto Microgravedad
# Completamente reproducible
# ============================================================

library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)

# 1. Leer resultados de expresión diferencial
res_mg <- read.csv("results/deseq2/tables/MG_vs_GC.csv")

# 2. Crear ranking para GSEA
ranks <- res_mg$log2FoldChange
names(ranks) <- res_mg$gene_name
ranks <- sort(ranks, decreasing = TRUE)

# 3. Cargar pathways (Hallmarks y GO)
m_hallmarks <- msigdbr(species = "Mus musculus", category = "H")
hallmarks <- split(m_hallmarks$gene_symbol, m_hallmarks$gs_name)

m_go_bp <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP"
)
go_pathways <- split(m_go_bp$gene_symbol, m_go_bp$gs_name)

# 4. Ejecutar GSEA
set.seed(123)
fgsea_hallmarks <- fgsea(
  hallmarks,
  ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 1000
)

# 5. Filtrar significativos (FDR < 0.05)
sig_hallmarks <- fgsea_hallmarks[padj < 0.05][order(-NES)]

# 6. Guardar resultados
write.csv(
  sig_hallmarks,
  "results/deseq2/enrichment/GSEA_hallmarks.csv",
  row.names = FALSE
)

# 7. Gráfica de GSEA (top 10)
top_hallmarks <- sig_hallmarks[1:min(10, nrow(sig_hallmarks))]

p_gsea <- ggplot(
  top_hallmarks,
  aes(x = reorder(gsub("HALLMARK_", "", pathway), NES), y = NES, fill = NES > 0)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "GSEA - Hallmarks enriquecidos en MG vs GC",
    x = "",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("results/deseq2/plots/gsea_hallmarks.png", p_gsea, width = 8, height = 6)

# 8. GO enrichment para genes upregulated
up_genes <- readLines("results/deseq2/enrichment/MG_vs_GC_UP_genes.txt")
up_genes <- up_genes[!grepl("^#|^$", up_genes)]

# Convertir a Entrez
up_entrez <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# GO enrichment
go_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Guardar resultados GO
write.csv(
  as.data.frame(go_up),
  "results/deseq2/enrichment/GO_upregulated.csv",
  row.names = FALSE
)

# Gráfica GO (top 10)
p_go <- dotplot(
  go_up,
  showCategory = 10,
  title = "GO Biological Process - Upregulated genes"
)
ggsave("results/deseq2/plots/go_upregulated.png", p_go, width = 8, height = 6)

# 9. GO enrichment para genes downregulated
down_genes <- readLines("results/deseq2/enrichment/MG_vs_GC_DOWN_genes.txt")
down_genes <- down_genes[!grepl("^#|^$", down_genes)]

down_entrez <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

go_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

write.csv(
  as.data.frame(go_down),
  "results/deseq2/enrichment/GO_downregulated.csv",
  row.names = FALSE
)

# 10. Mensaje de éxito
cat("\n✅ Análisis completado. Archivos generados:\n")
cat("  - results/deseq2/plots/gsea_hallmarks.png\n")
cat("  - results/deseq2/plots/go_upregulated.png\n")
cat("  - results/deseq2/enrichment/GSEA_hallmarks.csv\n")
cat("  - results/deseq2/enrichment/GO_upregulated.csv\n")
cat("  - results/deseq2/enrichment/GO_downregulated.csv\n")
