# -----------------------------------------------------------------------------------------------
# FGSEA & CLASSIC ENRICHMENT PIPELINE
# Alberto Atencia Rodriguez
# 30/10/2025

# -----------------------------
# Libraries & environment
rm(list=ls())
packages <- c("readxl","ggplot2","ggrepel","clusterProfiler","org.Mm.eg.db",
              "org.Hs.eg.db","openxlsx","ReactomePA","msigdbr","enrichplot",
              "fgsea","data.table","ggtext","dplyr")
invisible(lapply(packages, function(pkg){suppressMessages(library(pkg, character.only=TRUE))}))

# -----------------------------
# Helper functions
make_dir <- function(path){ if(!dir.exists(path)) dir.create(path, recursive=TRUE) }

write_log <- function(output_dir, msg){
  tryCatch({
    make_dir(output_dir)
    log_file <- file.path(output_dir,"log.txt")
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste0("[",timestamp,"] ",msg,"\n"), file=log_file, append=TRUE)
  }, error=function(e){message("Log could not be written: ", e$message)})
}

save_plot <- function(plot, filename, width=10, height=8, output_dir=NULL){
  tryCatch({
    ggsave(filename, plot=plot, width=width, height=height)
    if(!is.null(output_dir)) write_log(output_dir, paste("Saved plot:", filename))
  }, error=function(e){
    msg <- paste("Error saving plot", filename, ":", e$message)
    message(msg)
    if(!is.null(output_dir)) write_log(output_dir, msg)
  })
}

save_enrich_excel <- function(enrich_obj, filepath, output_dir=NULL){
  tryCatch({
    wb <- createWorkbook()
    addWorksheet(wb,"Results")
    df <- tryCatch(as.data.frame(enrich_obj), error=function(e) NULL)
    if(is.null(df) || nrow(df)==0){
      writeData(wb,"Results",data.frame(Message="No results"))
    } else {
      writeData(wb,"Results",df)
    }
    saveWorkbook(wb,filepath,overwrite=TRUE)
    if(!is.null(output_dir)) write_log(output_dir,paste("Saved enrichment Excel:", filepath))
  }, error=function(e){
    msg <- paste("Error saving enrichment Excel", filepath, ":", e$message)
    message(msg)
    if(!is.null(output_dir)) write_log(output_dir,msg)
  })
}

# -----------------------------
# Plotting function with specific_genes
plot_all_enrichment <- function(enrich_obj, outdir, prefix, max_show=30,
                                flexible_terms=NULL, specific_terms=NULL,
                                specific_genes=NULL,  # NEW
                                color_terms="darkred", size_terms=12){
  if(is.null(enrich_obj) || nrow(enrich_obj)==0) return(NULL)
  
  # Guardar tabla de enriquecimiento
  save_enrich_excel(enrich_obj@result, file.path(outdir, paste0(prefix,".xlsx")), outdir)
  
  # ------------------
  # Dotplot
  tryCatch({
    p_dot <- dotplot(enrich_obj, showCategory=max_show) + 
      ggtitle(paste(prefix, "Dotplot")) +
      theme(text = element_text(family="Helvetica"))
    
    labels_original <- enrich_obj@result$Description
    
    flexible_terms_lower <- if(!is.null(flexible_terms)) tolower(flexible_terms) else character(0)
    specific_terms_lower <- if(!is.null(specific_terms)) tolower(specific_terms) else character(0)
    specific_genes_upper <- if(!is.null(specific_genes)) toupper(specific_genes) else character(0)
    
    if(is.null(color_terms) || color_terms=="") color_terms <- "darkred"
    if(is.null(size_terms)) size_terms <- 12
    
    labels_custom <- sapply(seq_along(labels_original), function(i) {
      lbl <- labels_original[i]
      lbl_lower <- tolower(lbl)
      
      # Genes asociados al término
      genes_in_term <- unlist(strsplit(enrich_obj@result$geneID[i], "/"))
      
      genes_symbols <- if(all(grepl("^\\d+$", genes_in_term))) {
        # Posiblemente ENTREZID → convertir
        tryCatch({
          bitr(genes_in_term, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Mm.eg.db)$SYMBOL
        }, error=function(e) character(0))
      } else {
        # Ya son SYMBOLs
        genes_in_term
      }
      genes_symbols <- toupper(genes_symbols)
      
      highlight <- any(sapply(flexible_terms_lower, function(k) grepl(k, lbl_lower))) ||
        any(sapply(specific_terms_lower, function(k) grepl(k, lbl_lower))) ||
        any(genes_symbols %in% specific_genes_upper)
      
      if(highlight){
        paste0("<span style='font-size:", size_terms, "pt; color:", color_terms, ";'><b>", 
               tools::toTitleCase(lbl), "</b></span>")
      } else {
        paste0("<span style='font-size:", max(size_terms-2,10), "pt;'>", 
               tools::toTitleCase(lbl), "</span>")
      }
    })
    
    p_dot <- p_dot +
      scale_y_discrete(labels = labels_custom) +
      theme(axis.text.y = ggtext::element_markdown())
    
    height_dynamic <- max(8, max_show * 0.3)
    
    save_plot(p_dot,
              file.path(outdir, paste0(prefix,"_dotplot.png")),
              width = 10,
              height = height_dynamic,
              output_dir = outdir)
  }, error=function(e){write_log(outdir,paste("Error dotplot",prefix,e$message))})
  
  # ------------------
  # Treeplot
  tryCatch({
    terms <- pairwise_termsim(enrich_obj)
    p_tree <- treeplot(terms, showCategory=max_show, label_format=30) + 
      ggtitle(paste(prefix,"Treeplot"))
    
    max_label_length <- max(nchar(terms@result$Description[1:max_show]), na.rm=TRUE)
    right_margin <- max(0.1, min(0.5, max_label_length / 100))
    
    p_tree <- p_tree + 
      theme(
        plot.margin = margin(t=10, r=right_margin*100, b=10, l=10),
        plot.title = element_text(hjust=0.5)
      )
    
    save_plot(p_tree, file.path(outdir, paste0(prefix,"_treeplot.png")), 
              width = 10 + right_margin*10, height = 8, output_dir = outdir)
  }, error=function(e){write_log(outdir,paste("Error treeplot",prefix,e$message))})
  
  # ------------------
  # Cnetplot
  tryCatch({
    p_cnet <- cnetplot(enrich_obj, showCategory=max_show) + 
      ggtitle(paste(prefix,"Cnetplot"))
    save_plot(p_cnet, file.path(outdir,paste0(prefix,"_cnetplot.png")), output_dir=outdir)
  }, error=function(e){write_log(outdir,paste("Error cnetplot",prefix,e$message))})
  
  invisible(TRUE)
}


# -----------------------------
# Classic enrichment functions (GO, Reactome, KEGG, Hallmark)
perform_go_enrichment <- function(genes, universe_genes, orgdb, outdir, label, pval_cutoff, qval_cutoff,
                                  n_terms=30, flexible_terms=NULL, specific_terms=NULL,
                                  specific_genes=NULL, color_terms=NULL, size_terms=NULL){
  for(ont in c("BP","MF","CC")){
    ego <- enrichGO(gene=genes, OrgDb=orgdb, universe=universe_genes, keyType="SYMBOL",
                    ont=ont, pAdjustMethod="BH", pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff,
                    readable=TRUE)
    if(!is.null(ego) && nrow(ego)>0){
      plot_all_enrichment(ego, outdir, paste0("GO_",ont,"_",label), max_show=n_terms,
                          flexible_terms=flexible_terms, specific_terms=specific_terms,
                          specific_genes=specific_genes, color_terms=color_terms, size_terms=size_terms)
    }
  }
}

perform_reactome_enrichment <- function(genes_entrez, outdir, label, pval_cutoff, qval_cutoff,
                                        n_terms=30, flexible_terms=NULL, specific_terms=NULL,
                                        specific_genes=NULL, color_terms=NULL, size_terms=NULL){
  reactome <- enrichPathway(gene=genes_entrez, organism="mouse", pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff, readable=TRUE)
  if(!is.null(reactome) && nrow(reactome)>0){
    plot_all_enrichment(reactome, outdir, paste0("Reactome_",label), max_show=n_terms,
                        flexible_terms=flexible_terms, specific_terms=specific_terms,
                        specific_genes=specific_genes, color_terms=color_terms, size_terms=size_terms)
  }
}

perform_kegg_enrichment <- function(genes_entrez, outdir, label, pval_cutoff, qval_cutoff,
                                    n_terms=30, flexible_terms=NULL, specific_terms=NULL,
                                    specific_genes=NULL, color_terms=NULL, size_terms=NULL){
  kegg <- enrichKEGG(gene=genes_entrez, organism='mmu', pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff)
  if(!is.null(kegg) && nrow(kegg)>0){
    plot_all_enrichment(kegg, outdir, paste0("KEGG_",label), max_show=n_terms,
                        flexible_terms=flexible_terms, specific_terms=specific_terms,
                        specific_genes=specific_genes, color_terms=color_terms, size_terms=size_terms)
  }
}

perform_hallmark_enrichment <- function(genes, species, outdir, label, pval_cutoff, qval_cutoff,
                                        n_terms=30, flexible_terms=NULL, specific_terms=NULL,
                                        specific_genes=NULL, color_terms=NULL, size_terms=NULL){
  hallmark_sets <- msigdbr(species=species, category="H")
  term2gene <- hallmark_sets[, c("gs_name","gene_symbol")]
  hallmark <- enricher(gene=genes, TERM2GENE=term2gene, pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff)
  if(!is.null(hallmark) && nrow(hallmark)>0){
    plot_all_enrichment(hallmark, outdir, paste0("Hallmark_",label), max_show=n_terms,
                        flexible_terms=flexible_terms, specific_terms=specific_terms,
                        specific_genes=specific_genes, color_terms=color_terms, size_terms=size_terms)
  }
}

# -----------------------------
# Main pipeline function
run_enrichment_pipeline <- function(input_data,
                                    sheet_name = "Sheet1",
                                    padj_cutoff = 0.05, 
                                    log2fc_cutoff = 0,
                                    pval_cutoff = 0.05, 
                                    qval_cutoff = 0.1,
                                    output_dir = "output_enrichment",
                                    n_terms = 30, 
                                    flexible_terms = NULL, 
                                    specific_terms = NULL,
                                    specific_genes = NULL,  # NEW
                                    color_terms = NULL, 
                                    size_terms = NULL) {
  
  make_dir(output_dir)
  
  # ---- Load data ----
  if (is.character(input_data) && file.exists(input_data)) {
    message("📘 Reading Excel file: ", input_data)
    df <- readxl::read_excel(input_data, sheet = sheet_name) %>% as.data.frame()
  } else if (is.data.frame(input_data)) {
    message("🧬 Using provided data frame directly")
    df <- input_data
  } else {
    stop("❌ 'input_data' must be either a valid Excel file path or a data frame.")
  }
  
  # --- Normalize column names and remove duplicates ---
  orig_names <- colnames(df)
  
  make_base <- function(x) {
    x2 <- tolower(x)
    x2 <- trimws(x2)
    x2 <- sub("\\.[0-9]+$", "", x2)
    x2 <- sub("[[:punct:]]+$", "", x2)
    x2 <- gsub("[[:space:][:punct:]]+", "", x2)
    return(x2)
  }
  
  base_names <- vapply(orig_names, make_base, character(1))
  
  new_names <- orig_names
  seen <- list()
  for (i in seq_along(base_names)) {
    b <- base_names[i]
    if (!(b %in% names(seen))) {
      seen[[b]] <- 1
      new_names[i] <- orig_names[i]
    } else {
      seen[[b]] <- seen[[b]] + 1
      suffix <- seen[[b]] - 1
      new_names[i] <- paste0(orig_names[i], "_", suffix)
    }
  }
  
  colnames(df) <- new_names
  changed <- which(orig_names != new_names)
  if (length(changed) > 0) {
    msg <- paste0("Renamed duplicated-base columns:\n",
                  paste0("  '", orig_names[changed], "' -> '", new_names[changed], "'", collapse = "\n"))
    message(msg)
    write_log(output_dir, msg)
  }
  
  # ---- Flexible column detection ----
  colnames(df) <- gsub("[[:space:][:punct:]]+", "", tolower(colnames(df)))
  
  find_col <- function(df, aliases) {
    for (a in aliases) {
      if (a %in% colnames(df)) return(a)
    }
    return(NULL)
  }
  
  symbol_col <- find_col(df, c("symbol", "gene", "genesymbol", "genename"))
  padj_col   <- find_col(df, c("padj", "adjpval", "adj", "pvaladj", "pvalueadj", 
                               "pvalueadjusted", "pvaladjusted", "padjval", 
                               "padjustedval", "padjvalue", "padjustedvalue"))
  pval_col   <- find_col(df, c("pval", "pvalue", "p", "pv"))
  logfc_col  <- find_col(df, c("log2foldchange", "logfc", "foldchange", "logfold"))
  
  if (is.null(symbol_col))
    stop("❌ Could not find a 'symbol' or 'gene' column.")
  if (is.null(logfc_col))
    stop("❌ Could not find a log2FoldChange column (try naming it 'log2FoldChange' or 'logFC').")
  if (is.null(padj_col) && is.null(pval_col))
    stop("❌ No 'padj' or 'pvalue' column found. At least one is required.")
  if (is.null(padj_col) && !is.null(pval_col)) {
    warning("⚠️ No adjusted p-values found — using raw p-values instead.")
    padj_col <- pval_col
  }
  
  # ---- Internal renaming ----
  df <- df %>% rename(
    SYMBOL = !!symbol_col,
    padj = !!padj_col,
    log2FoldChange = !!logfc_col
  )
  
  # ---- Capitalize gene symbols ----
  df$SYMBOL <- sapply(df$SYMBOL, function(x)
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  )
  
  # ---- Split dataset ----
  universe_genes <- df$SYMBOL[!is.na(df$padj)]
  up_df   <- df %>% filter(padj < padj_cutoff & log2FoldChange > log2fc_cutoff)
  down_df <- df %>% filter(padj < padj_cutoff & log2FoldChange < -log2fc_cutoff)
  
  # ---- Convert to ENTREZ IDs ----
  entrez_up   <- bitr(up_df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  entrez_down <- bitr(down_df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # ---- Output directories ----
  dirs <- list(
    GO = file.path(output_dir, "GO"),
    Reactome = file.path(output_dir, "Reactome"),
    KEGG = file.path(output_dir, "KEGG"),
    Hallmark = file.path(output_dir, "Hallmark"),
    FGSEA = file.path(output_dir, "FGSEA")
  )
  lapply(dirs, make_dir)
  
  # ---- Classic enrichment ----
  if (nrow(up_df) > 0) {
    perform_go_enrichment(up_df$SYMBOL, universe_genes, org.Mm.eg.db, dirs$GO, "UP",
                          pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                          specific_genes, color_terms, size_terms)
    perform_reactome_enrichment(entrez_up$ENTREZID, dirs$Reactome, "UP",
                                pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                                specific_genes, color_terms, size_terms)
    perform_kegg_enrichment(entrez_up$ENTREZID, dirs$KEGG, "UP",
                            pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                            specific_genes, color_terms, size_terms)
    perform_hallmark_enrichment(up_df$SYMBOL, "Mus musculus", dirs$Hallmark, "UP",
                                pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                                specific_genes, color_terms, size_terms)
  }
  
  if (nrow(down_df) > 0) {
    perform_go_enrichment(down_df$SYMBOL, universe_genes, org.Mm.eg.db, dirs$GO, "DOWN",
                          pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                          specific_genes, color_terms, size_terms)
    perform_reactome_enrichment(entrez_down$ENTREZID, dirs$Reactome, "DOWN",
                                pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                                specific_genes, color_terms, size_terms)
    perform_kegg_enrichment(entrez_down$ENTREZID, dirs$KEGG, "DOWN",
                            pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                            specific_genes, color_terms, size_terms)
    perform_hallmark_enrichment(down_df$SYMBOL, "Mus musculus", dirs$Hallmark, "DOWN",
                                pval_cutoff, qval_cutoff, n_terms, flexible_terms, specific_terms,
                                specific_genes, color_terms, size_terms)
  }
  
  write_log(output_dir, "Pipeline completed successfully.")
  message("✅ Pipeline completed — results saved in: ", output_dir)
  
  return(invisible(TRUE))
}

# -----------------------------------------------------------------------------------------------
# Example usage
# -----------------------------

# It is needed to specify: 
#   - Excel_path & sheet_name (if it is changed in your excel) or data frame coming from a DE analysis
#   - Output_dir
#   - Rest of parametres are predefined
                              
# You can easily run the main function after all this code "run_enrichment_pipeline" or call it from another R script for clearness. Here are some examples:


# On a separate R script you must first call this script from your computer, then you will be able to use it from another R sript:
# source(/path/to/run_enrichment.R)

# 1. Using an Excel file
# run_enrichment_pipeline(input_data = "path/to/results.xlsx", sheet_name = "Sheet1")

# 2. Using a data frame (DEG results from a single-cell cluster)
# deg_cluster3 <- FindMarkers(seurat_obj, ident.1 = "Cluster3", ident.2 = NULL)
# run_enrichment_pipeline(input_data = deg_cluster3, output_dir = "results_cluster3")

# -----------------------------------------------------------------------------------------------
# 3. Full parameter example
# -----------------------------------------------------------------------------------------------
# Example of a complete run using all parameters.
# You can provide either:
#   - A path to an Excel/CSV file (input_data = "path/to/file.xlsx")
#   - A data frame object already loaded in R (input_data = deg_df)
#
# Note:
#   - 'sheet_name' is only required when the input is an Excel file.
#   - 'flexible_terms', 'specific_terms', and 'specific_genes' are optional
#     parameters used to highlight specific terms or genes in the plots.
# -----------------------------------------------------------------------------------------------

# run_enrichment_pipeline(
#   input_data = "/path/to/your/excel",                        # or just the dataframe
#   sheet_name = "Sheet1",                                     # only needed if input_data is an Excel file
#   padj_cutoff = 0.05,                                        # adjusted p-value cutoff for DEGs
#   log2fc_cutoff = 0,                                         # log2 fold-change threshold
#   pval_cutoff = 0.05,                                        # p-value cutoff for enrichment
#   qval_cutoff = 0.1,                                         # q-value cutoff for enrichment
#   output_dir = "/path/to/output_folder",                     # directory where results will be saved
#   n_terms = 30,                                              # number of top terms to display in plots
#   flexible_terms = c("mito", "resp"),                        # highlight any term containing these substrings
#   specific_terms = c("Transmembrane Transporter Complex"),   # highlight exact term(s)
#   specific_genes = c("Mki67", "Arap2", "Prok2", "Syne2"),    # highlight terms containing these genes
#   color_terms = "#0052DB",                                   # color used to highlight terms in plots, default value is darkred
#   size_terms = 12                                            # text size for highlighted terms
# )




