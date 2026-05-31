#!/usr/bin/env Rscript
# Script para generar diagramas de Venn por separado

library(VennDiagram)

# Rutas absolutas
base_dir <- "/export/space3/users/zyanyava/2026-2/Transcriptomica/RNA_seq_Project"
star_dir <- file.path(base_dir, "results/deseq2/enrichment")
salmon_dir <- file.path(base_dir, "results/deseq2_salmon/enrichment")

# Leer archivos - UPREGULATED
star_up <- readLines(file.path(star_dir, "MG_vs_GC_UP_genes.txt"))
salmon_up <- readLines(file.path(salmon_dir, "MG_vs_GC_UP_genes.txt"))

# Leer archivos - DOWNREGULATED
star_down <- readLines(file.path(star_dir, "MG_vs_GC_DOWN_genes.txt"))
salmon_down <- readLines(file.path(salmon_dir, "MG_vs_GC_DOWN_genes.txt"))

# Limpiar (eliminar líneas vacías y comentarios)
star_up <- star_up[star_up != "" & !grepl("^#", star_up)]
salmon_up <- salmon_up[salmon_up != "" & !grepl("^#", salmon_up)]
star_down <- star_down[star_down != "" & !grepl("^#", star_down)]
salmon_down <- salmon_down[salmon_down != "" & !grepl("^#", salmon_down)]

# Calcular intersecciones
intersect_up <- length(intersect(star_up, salmon_up))
intersect_down <- length(intersect(star_down, salmon_down))

# Mostrar números en consola
cat("\n========== COMPARACIÓN STAR vs Salmon (MG vs GC) ==========\n\n")
cat("GENES UPREGULATED:\n")
cat("  STAR:", length(star_up), "\n")
cat("  Salmon:", length(salmon_up), "\n")
cat("  Intersección:", intersect_up, "\n")
cat("  Solo STAR:", length(star_up) - intersect_up, "\n")
cat("  Solo Salmon:", length(salmon_up) - intersect_up, "\n\n")
cat("GENES DOWNREGULATED:\n")
cat("  STAR:", length(star_down), "\n")
cat("  Salmon:", length(salmon_down), "\n")
cat("  Intersección:", intersect_down, "\n")
cat("  Solo STAR:", length(star_down) - intersect_down, "\n")
cat("  Solo Salmon:", length(salmon_down) - intersect_down, "\n\n")

# ============================================================
# IMAGEN 1: Venn upregulated
# ============================================================
png(
  file.path(base_dir, "docs/figures/venn_upregulated.png"),
  width = 600,
  height = 600,
  res = 120
)

draw.pairwise.venn(
  area1 = length(star_up),
  area2 = length(salmon_up),
  cross.area = intersect_up,
  category = c("STAR", "Salmon"),
  fill = c("#1B9E77", "#D95F02"),
  alpha = 0.6,
  cat.cex = 1.8,
  cex = 1.5,
  cat.fontface = "bold",
  main = "Genes upregulated (MG vs GC)",
  main.cex = 1.3,
  main.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

dev.off()
cat(
  "✓ Imagen 1 guardada:",
  file.path(base_dir, "docs/figures/venn_upregulated.png"),
  "\n"
)

# ============================================================
# IMAGEN 2: Venn downregulated
# ============================================================
png(
  file.path(base_dir, "docs/figures/venn_downregulated.png"),
  width = 600,
  height = 600,
  res = 120
)

draw.pairwise.venn(
  area1 = length(star_down),
  area2 = length(salmon_down),
  cross.area = intersect_down,
  category = c("STAR", "Salmon"),
  fill = c("#1B9E77", "#D95F02"),
  alpha = 0.6,
  cat.cex = 1.8,
  cex = 1.5,
  cat.fontface = "bold",
  main = "Genes downregulated (MG vs GC)",
  main.cex = 1.3,
  main.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

dev.off()
cat(
  "✓ Imagen 2 guardada:",
  file.path(base_dir, "docs/figures/venn_downregulated.png"),
  "\n"
)

cat("\n✅ ¡Listo! Las dos imágenes se han generado correctamente.\n")
