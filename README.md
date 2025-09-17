# Enrichment-Pipeline
Automatic R pipeline for FGSEA and classic enrichment (GO, KEGG, Reactome, Hallmark)

# FGSEA & Classic Enrichment Pipeline

This R script provides a comprehensive workflow to perform gene set enrichment analysis from RNA-seq or similar gene expression results. 

It combines:

- **FGSEA** (fast gene set enrichment analysis) for Hallmark, GO_BP, and Reactome pathways.
- **Classic enrichment analysis** using GO (BP, MF, CC), KEGG, Reactome, and Hallmark gene sets.
- **Automatic plotting**: dotplots, treeplots, cnetplots, and enrichment plots for top pathways.
- **Exporting results**: Excel files with enrichment results and PNG plots.

## Features
- Separate analysis of up- and down-regulated genes.
- Configurable parameters: adjusted p-value cutoff, log2 fold change cutoff, p/q-value cutoffs for enrichment.
- Robust handling of Excel sheets, missing data, and directory creation.
- Logs messages automatically for reproducibility and debugging.

## Requirements
- R >= 4.5.0
- Packages: `readxl`, `ggplot2`, `ggrepel`, `clusterProfiler`, `org.Mm.eg.db`, `org.Hs.eg.db`, `openxlsx`, `ReactomePA`, `msigdbr`, `enrichplot`, `fgsea`, `data.table`.

## Data frame requirements
The DEGs dataframe must contain the following column names: SYMBOL (for gene name column), log2FoldChange (for fold change column) and padj (for adjusted p-value column).

## Example usage
```r
excel_path <- "/path/to/your/data.xlsx"

run_enrichment_pipeline(
  excel_path = excel_path,
  sheet_name = "Sheet1",
  padj_cutoff = 0.05,
  log2fc_cutoff = 0,
  pval_cutoff = 0.05,
  qval_cutoff = 0.1,
  output_dir = "/path/to/output_folder"
)
