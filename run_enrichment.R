# FGSEA & CLASSIC ENRICHMENT PIPELINE
# Alberto Atencia Rodriguez
# 15/09/2025

# -----------------------------
# Libraries & environment
rm(list=ls())
packages <- c("readxl", "ggplot2", "ggrepel", "clusterProfiler", "org.Mm.eg.db", 
              "org.Hs.eg.db", "openxlsx", "ReactomePA", "msigdbr", "enrichplot", "fgsea", "data.table")
invisible(lapply(packages, function(pkg) {suppressMessages(library(pkg, character.only = TRUE))}))


# -----------------------------
# Auxiliar Functions
make_dir <- function(path){ if(!dir.exists(path)) dir.create(path, recursive=TRUE) }

write_log <- function(output_dir, msg) {
  tryCatch({
    make_dir(output_dir)
    log_file <- file.path(output_dir, "log.txt")
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste0("[", timestamp, "] ", msg, "\n"), file = log_file, append = TRUE)
  }, error=function(e){ message("Log could not be written: ", e$message) })
}

save_plot <- function(plot, filename, width=10, height=8, output_dir=NULL) {
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
    addWorksheet(wb, "Results")
    df <- tryCatch(as.data.frame(enrich_obj), error=function(e) NULL)
    if(is.null(df) || nrow(df)==0){
      writeData(wb, "Results", data.frame(Message="No results"))
    } else {
      writeData(wb, "Results", df)
    }
    saveWorkbook(wb, filepath, overwrite=TRUE)
    if(!is.null(output_dir)) write_log(output_dir, paste("Saved enrichment excel:", filepath))
  }, error=function(e){
    msg <- paste("Error saving enrichment Excel", filepath, ":", e$message)
    message(msg)
    if(!is.null(output_dir)) write_log(output_dir, msg)
  })
}

plot_all_enrichment <- function(enrich_obj, outdir, prefix, max_show=30){
  if(is.null(enrich_obj) || nrow(enrich_obj)==0) return(NULL)
  
  save_enrich_excel(enrich_obj@result, file.path(outdir, paste0(prefix,".xlsx")), outdir)
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(enrich_obj, showCategory=max_show) + ggtitle(paste(prefix,"Dotplot"))
    save_plot(p_dot, file.path(outdir,paste0(prefix,"_dotplot.png")), output_dir=outdir)
  }, error=function(e){write_log(outdir,paste("Error dotplot",prefix,e$message))})
  
  # Treeplot
  tryCatch({
    terms <- pairwise_termsim(enrich_obj)
    p_tree <- treeplot(terms, showCategory=max_show,label_format=30) + ggtitle(paste(prefix,"Treeplot"))
    save_plot(p_tree, file.path(outdir,paste0(prefix,"_treeplot.png")), output_dir=outdir)
  }, error=function(e){write_log(outdir,paste("Error treeplot",prefix,e$message))})
  
  # Cnetplot
  tryCatch({
    p_map <- cnetplot(enrich_obj, showCategory=10) + ggtitle(paste(prefix,"Cnetplot"))
    save_plot(p_map, file.path(outdir,paste0(prefix,"_cnetplot.png")), output_dir=outdir)
  }, error=function(e){write_log(outdir,paste("Error cnetplot",prefix,e$message))})
  
  invisible(TRUE)
}

# -----------------------------
# Enrichment Functions
perform_go_enrichment <- function(genes, universe_genes, orgdb, outdir, label, pval_cutoff, qval_cutoff){
  ontologies <- c("BP","MF","CC")
  for(ont in ontologies){
    ego <- enrichGO(
      gene=genes,
      OrgDb=orgdb,
      universe=universe_genes,
      keyType="SYMBOL",
      ont=ont,
      pAdjustMethod="BH",
      pvalueCutoff=pval_cutoff,
      qvalueCutoff=qval_cutoff,
      readable=TRUE
    )
    if(!is.null(ego) && nrow(ego)>0){
      plot_all_enrichment(ego,outdir,paste0("GO_",ont,"_",label))
    }
  }
}

perform_reactome_enrichment <- function(genes_entrez, outdir, label, pval_cutoff, qval_cutoff){
  reactome <- enrichPathway(gene=genes_entrez, organism="mouse", pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff, readable=TRUE)
  if(!is.null(reactome) && nrow(reactome)>0){
    plot_all_enrichment(reactome,outdir,paste0("Reactome_",label))
  }
}

perform_kegg_enrichment <- function(genes_entrez, outdir, label, pval_cutoff, qval_cutoff){
  kegg <- enrichKEGG(gene=genes_entrez, organism='mmu', pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff)
  if(!is.null(kegg) && nrow(kegg)>0){
    plot_all_enrichment(kegg,outdir,paste0("KEGG_",label))
  }
}

perform_hallmark_enrichment <- function(genes, species, outdir, label, pval_cutoff, qval_cutoff){
  hallmark_sets <- msigdbr(species=species, category="H")
  term2gene <- hallmark_sets[,c("gs_name","gene_symbol")]
  hallmark <- enricher(gene=genes, TERM2GENE=term2gene, pvalueCutoff=pval_cutoff, qvalueCutoff=qval_cutoff)
  if(!is.null(hallmark) && nrow(hallmark)>0){
    plot_all_enrichment(hallmark,outdir,paste0("Hallmark_",label))
  }
}

# -----------------------------
# FGSEA 
perform_fgsea_all_separated_excel <- function(excel_path, sheet_name="ALL GENES", output_excel, padj_cutoff=0.05, log2fc_cutoff=0){
  make_dir(dirname(output_excel))
  df <- readxl::read_excel(excel_path,sheet=sheet_name)
  df <- as.data.frame(df)
  df <- df[!is.na(df$padj) & df$padj<padj_cutoff,]
  if(nrow(df)==0) stop("No hay genes después de filtrar por padj < ", padj_cutoff)
  
  up_df <- df[df$log2FoldChange>log2fc_cutoff,]
  down_df <- df[df$log2FoldChange< -log2fc_cutoff,]
  
  sets_list <- list(
    Hallmark  = split(msigdbr(species="Mus musculus", category="H")$gene_symbol,
                      msigdbr(species="Mus musculus", category="H")$gs_name),
    GO_BP     = split(msigdbr(species="Mus musculus", category="C5", subcategory="BP")$gene_symbol,
                      msigdbr(species="Mus musculus", category="C5", subcategory="BP")$gs_name),
    Reactome  = split(msigdbr(species="Mus musculus", category="C2", subcategory="REACTOME")$gene_symbol,
                      msigdbr(species="Mus musculus", category="C2", subcategory="REACTOME")$gs_name)
  )
  
  wb <- createWorkbook()
  results_list <- list(UP=list(), DOWN=list())
  
  for(set_name in names(sets_list)){
    pathways <- sets_list[[set_name]]
    
    # UP
    if(nrow(up_df)>0){
      ranks_vector_up <- setNames(up_df$log2FoldChange, up_df$SYMBOL)
      fgseaRes_up <- fgsea(pathways=pathways, stats=ranks_vector_up, minSize=15, maxSize=500)
      fgseaRes_up <- as.data.table(fgseaRes_up)
      results_list$UP[[set_name]] <- fgseaRes_up
      addWorksheet(wb,paste0(set_name,"_UP"))
      writeData(wb,paste0(set_name,"_UP"),fgseaRes_up)
      if(nrow(fgseaRes_up)>0){
        topPath <- fgseaRes_up[order(-NES)][1,pathway]
        png(file.path(dirname(output_excel),paste0("FGSEA_",set_name,"_UP_top.png")), width=800, height=600)
        print(plotEnrichment(pathways[[topPath]], ranks_vector_up) + labs(title=paste(set_name,"UP:",topPath)))
        dev.off()
      }
    }
    
    # DOWN
    if(nrow(down_df)>0){
      ranks_vector_down <- setNames(down_df$log2FoldChange, down_df$SYMBOL)
      fgseaRes_down <- fgsea(pathways=pathways, stats=ranks_vector_down, minSize=15, maxSize=500)
      fgseaRes_down <- as.data.table(fgseaRes_down)
      results_list$DOWN[[set_name]] <- fgseaRes_down
      addWorksheet(wb,paste0(set_name,"_DOWN"))
      writeData(wb,paste0(set_name,"_DOWN"),fgseaRes_down)
      if(nrow(fgseaRes_down)>0){
        topPath <- fgseaRes_down[order(-NES)][1,pathway]
        png(file.path(dirname(output_excel),paste0("FGSEA_",set_name,"_DOWN_top.png")), width=800, height=600)
        print(plotEnrichment(pathways[[topPath]], ranks_vector_down) + labs(title=paste(set_name,"DOWN:",topPath)))
        dev.off()
      }
    }
  }
  
  saveWorkbook(wb, output_excel, overwrite=TRUE)
  message("Excel saved on: ", output_excel)
  return(results_list)
}

# -----------------------------
# Main
run_enrichment_pipeline <- function(excel_path, sheet_name="Sheet1", 
                                    padj_cutoff=0.05, log2fc_cutoff=0,
                                    pval_cutoff=0.05, qval_cutoff=0.1,
                                    output_dir="output_enrichment"){
  
  make_dir(output_dir)
  fgsea_excel <- file.path(output_dir,"FGSEA_results.xlsx")
  
  # FGSEA
  fgsea_results <- perform_fgsea_all_separated_excel(excel_path,sheet_name,fgsea_excel,padj_cutoff,log2fc_cutoff)
  
  # Classic Enrichment with Go, Reactome & Hallmark
  df <- readxl::read_excel(excel_path,sheet=sheet_name)
  df <- as.data.frame(df)
  universe_genes <- df$SYMBOL[!is.na(df$padj)]
  up_df <- df[df$padj<padj_cutoff & df$log2FoldChange>log2fc_cutoff,]
  down_df <- df[df$padj<padj_cutoff & df$log2FoldChange< -log2fc_cutoff,]
  
  entrez_up <- bitr(up_df$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  entrez_down <- bitr(down_df$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  go_dir <- file.path(output_dir,"GO")
  reactome_dir <- file.path(output_dir,"Reactome")
  kegg_dir <- file.path(output_dir,"KEGG")
  hallmark_dir <- file.path(output_dir,"Hallmark")
  lapply(list(go_dir,reactome_dir,kegg_dir,hallmark_dir),make_dir)
  
  if(nrow(up_df)>0){
    perform_go_enrichment(up_df$SYMBOL, universe_genes, org.Mm.eg.db, go_dir, "UP", pval_cutoff, qval_cutoff)
    perform_reactome_enrichment(entrez_up$ENTREZID, reactome_dir, "UP", pval_cutoff, qval_cutoff)
    perform_kegg_enrichment(entrez_up$ENTREZID, kegg_dir, "UP", pval_cutoff, qval_cutoff)
    perform_hallmark_enrichment(up_df$SYMBOL, "Mus musculus", hallmark_dir, "UP", pval_cutoff, qval_cutoff)
  }
  
  if(nrow(down_df)>0){
    perform_go_enrichment(down_df$SYMBOL, universe_genes, org.Mm.eg.db, go_dir, "DOWN", pval_cutoff, qval_cutoff)
    perform_reactome_enrichment(entrez_down$ENTREZID, reactome_dir, "DOWN", pval_cutoff, qval_cutoff)
    perform_kegg_enrichment(entrez_down$ENTREZID, kegg_dir, "DOWN", pval_cutoff, qval_cutoff)
    perform_hallmark_enrichment(down_df$SYMBOL, "Mus musculus", hallmark_dir, "DOWN", pval_cutoff, qval_cutoff)
  }
  
  message("=== PIPELINE IS FINISHED ===")
  return(list(fgsea=fgsea_results))
}

# -----------------------------
# Usage Example
# It is needed to specify excel_path, sheet_name (if it is changed in your excel) & output_dir, the rest of parametres are predefined
# You can easily run the main function after all this code "run_enrichment_pipeline" or call it from another R script for clearness.

# On a separate R:

# source(/path/to/run_enrichment.R)
# excel_path <- "/path/to/deg_excel"
# run_enrichment_pipeline(excel_path=excel_path,
#                         sheet_name="Sheet1",
#                         padj_cutoff=0.05,
#                         log2fc_cutoff=0,
#                         pval_cutoff=0.05,
#                         qval_cutoff=0.1,
#                         output_dir="/output_path")

