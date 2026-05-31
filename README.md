# Proyecto Final - Transcriptómica 2026-2

**Análisis transcriptómico de hígado de ratón en condiciones de microgravedad**

**Autora:** Zyanya Valentina Velazquez Aldrete  
**Fecha:** 31 de mayo de 2026  
**Materia:** Transcriptómica (2026-2)  
**Institución:** Universidad Nacional Autónoma de México - Licenciatura en Ciencias Genómicas

---

## Resumen

Este proyecto analiza datos de RNA-seq de hígado de ratón (*Mus musculus*) expuesto a tres condiciones experimentales: control terrestre (GC), gravedad artificial de 1g (A1G) y microgravedad real (MG). Se identificaron genes diferencialmente expresados y se realizaron análisis de enriquecimiento funcional (GSEA, GO, STRING) para caracterizar las rutas biológicas afectadas por la microgravedad.

---

## Estructura del repositorio

```
ProyectoFinal_Transcriptomica/
├── src/                   # Scripts completos del análisis
├── docs/                  # Documento final (Rmd, PDF, referencias)
├── results/               # Tablas de expresión diferencial y listas de genes
│   ├── deseq2/            # Resultados del pipeline STAR
│   ├── deseq2_salmon/     # Resultados del pipeline Salmon
│   ├── DAVID/             # Resultados de GO enrichment
│   ├── STRING/            # Imágenes de redes de interacción
│   └── GSEA/              # Resultados de GSEA (WebGestalt)
└── README.md              # Este archivo
```

---

## Datos crudos

Los datos de RNA-seq están disponibles en NCBI SRA bajo el BioProject:  
[PRJNA1005192](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1005192)

Estudio original: Kurosawa et al. (2021) - *Scientific Reports*

---

## Requisitos computacionales

- **Sistema operativo:** Linux (Rocky Linux 8.10)
- **R:** 4.3.2
- **Herramientas externas:** FastQC, MultiQC, fastp, Nextflow, STAR, featureCounts, Salmon
- **Paquetes de R:** DESeq2, tximport, ggplot2, pheatmap, clusterProfiler, STRINGdb, knitr, kableExtra

---

## Reproducibilidad

Para reproducir el análisis completo:

```bash
# Clonar el repositorio
git clone https://github.com/ZyanyaAld/ProyectoFinal_Transcriptomica.git
cd ProyectoFinal_Transcriptomica

# Ejecutar pipeline de control de calidad
nextflow run src/rna_qc.nf -profile docker

# Alineamiento con STAR y cuantificación con featureCounts
bash src/run_star.sh
bash src/run_featureCounts.sh

# Cuantificación con Salmon
bash src/run_salmon.sh

# Análisis de expresión diferencial
Rscript src/deseq2_star.R
Rscript src/deseq2_salmon.R

# Renderizar documento final
R -e "rmarkdown::render('docs/VelazquezAldrete_Zyanya_ProyectoFinal.Rmd')"
```

---

## Resultados principales

- **325 genes diferencialmente expresados** en microgravedad (228 upregulated, 97 downregulated)
- **Activación de rutas de estrés oxidativo, autofagia y metabolismo de xenobióticos**
- **Disminución de tráfico vesicular del RE y metabolismo lipídico**
- **Alta concordancia entre pipelines STAR y Salmon** (Pearson r > 0.9)

---

## Contacto

Para dudas o solicitud de archivos adicionales:  
📧 zyanyava@ciencias.unam.mx

---

## Licencia

MIT License - Copyright (c) 2026 Zyanya Valentina Velazquez Aldrete
