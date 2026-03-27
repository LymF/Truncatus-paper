# Truncatus-paper

R scripts used for the transcriptomic analyses published in:

> **Potential effect of *Wolbachia* on virus restriction in the spider mite *T. truncatus***
> *Frontiers in Microbiology*, Sec. Microbial Symbioses, Volume 16 — 28 May 2025
> https://doi.org/10.3389/fmicb.2025.1570606

Covers differential expression analysis with DESeq2, gene set enrichment analysis (GSEA) with clusterProfiler, and alpha diversity assessment.

## Overview

This repository contains the complete R analysis workflow applied to RNA-seq data from *Tetranychus truncatus*, a phytophagous spider mite. Analyses include identification of differentially expressed genes (DEGs), pathway enrichment, and ecological diversity metrics.

## Scripts

| Script | Description |
|---|---|
| `DESEQ2final.R` | Full DESeq2 differential expression pipeline — count matrix import, DESeq2 object construction, normalization, DEG calling, MA plots, heatmaps, and result export |
| `gsea.Rmd` | Gene Set Enrichment Analysis (GSEA) using clusterProfiler — ranked gene list preparation, enrichment analysis, and visualization (ridgeplots, dot plots, GSEA score plots) |
| `Alpha diversity.R` | Alpha diversity analysis — diversity index calculation and statistical comparison between sample groups |

## Dependencies

### DESeq2 analysis

```r
BiocManager::install("DESeq2")
install.packages(c("ggplot2", "pheatmap", "RColorBrewer"))
```

### GSEA analysis

```r
BiocManager::install(c("clusterProfiler", "pathview", "enrichplot"))
install.packages("ggplot2")
```

## Usage

### Differential expression (DESeq2)

The script expects a count matrix (genes × samples) and a sample metadata CSV with a `Condition` column.

```r
# Update file paths inside the script before running
source("DESEQ2final.R")
```

**Count matrix format (`counts.tsv`):**

```
gene_id    sample1    sample2    sample3
gene001    120        98         145
gene002    0          3          1
```

**Sample metadata format (`sample_info.csv`):**

```
,Condition
sample1,treated
sample2,un
sample3,treated
```

> The reference condition is set with `relevel(..., ref = "un")`. Adjust to match your experimental design.

### GSEA

Open `gsea.Rmd` in RStudio and render with:

```r
rmarkdown::render("gsea.Rmd")
```

The notebook uses organism-specific annotation databases for GO and KEGG term mapping. Replace the `organism` variable with the appropriate `OrgDb` for your target species.

## Analysis details

### DESeq2 workflow

1. Count matrix imported from gene-level quantification (e.g., STAR + featureCounts or Salmon + tximport)
2. DESeq2 object constructed with `~ Condition` design formula
3. Size factor estimation and dispersion modeling
4. Wald test for pairwise condition comparison
5. Results filtered by adjusted p-value (padj < 0.05) and log2 fold change threshold
6. Visualization: MA plot, normalized count boxplots, sample distance heatmap, PCA

### GSEA workflow

1. Genes ranked by log2 fold change from DESeq2 results
2. GSEA performed with `gseGO` or `gseKEGG` using the ranked gene list
3. Results visualized with ridgeplots, dot plots, and enrichment score plots

## Citation

If you use or adapt these scripts, please cite:

```
Potential effect of Wolbachia on virus restriction in the spider mite T. truncatus.
Frontiers in Microbiology, Sec. Microbial Symbioses, Volume 16, 28 May 2025.
https://doi.org/10.3389/fmicb.2025.1570606
```
