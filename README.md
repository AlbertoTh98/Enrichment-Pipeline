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
- Packages: `readxl`, `ggplot2`, `ggrepel`, `clusterProfiler`, `org.Mm.eg.db`, `org.Hs.eg.db`, `openxlsx`, `ReactomePA`, `msigdbr`, `enrichplot`, `fgsea`, `data.table`, `ggtext`, `dplyr`.

## Data frame requirements
The input DEG (Differentially Expressed Genes) file or data frame must contain the following columns:

| Column                  | Description                                                       | Accepted Names                                                                                                                               |
| :---------------------- | :---------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------- |
| **SYMBOL**              | Gene name or gene symbol.                                         | `symbol`, `gene`, `genesymbol`, `genename`                                                                                                   |
| **log2FoldChange**      | Log₂ fold change value for each gene.                             | `log2foldchange`, `logfc`, `foldchange`, `logfold`                                                                                           |
| **padj**                | Adjusted p-value (or FDR). Used to select significant DEGs.       | `padj`, `adjpval`, `adj`, `pvaladj`, `pvalueadj`, `pvalueadjusted`, `pvaladjusted`, `padjval`, `padjustedval`, `padjvalue`, `padjustedvalue` |
| *(Optional)* **pvalue** | Raw p-value column (used only if no adjusted p-values are found). | `pval`, `pvalue`, `p`, `pv`                                                                                                                  |

## Example usage
```r
excel_path <- "/path/to/your/data.xlsx"
deg_dataframe <- FindMarkers(seurat_obj, ident.1 = "Cluster3", ident.2 = NULL)

run_enrichment_pipeline(
   input_data = excel_path,                                   # or just the deg_dataframe
   sheet_name = "Sheet1",                                     # only needed if input_data is an Excel file
   padj_cutoff = 0.05,                                        # adjusted p-value cutoff for DEGs
   log2fc_cutoff = 0,                                         # log2 fold-change threshold
   pval_cutoff = 0.05,                                        # p-value cutoff for enrichment
   qval_cutoff = 0.1,                                         # q-value cutoff for enrichment
   output_dir = "/path/to/output_folder",                     # directory where results will be saved
   n_terms = 30,                                              # number of top terms to display in plots
   flexible_terms = c("mito", "resp"),                        # highlight any term containing these substrings "mito" and "resp" are examples for mitochondria or respiratory terms
   specific_terms = c("Transmembrane Transporter Complex"),   # highlight exact term(s)
   specific_genes = c("Mki67", "Arap2", "Prok2", "Syne2"),    # highlight terms containing these genes
   color_terms = "#0052DB",                                   # color used to highlight terms in plots, default value is darkred
   size_terms = 12                                            # text size for highlighted terms
)
