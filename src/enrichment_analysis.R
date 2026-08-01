#!/usr/bin/env Rscript
# ============================================================
# Gene set enrichment analysis on phase-shifted studies.
# For each of the 4 studies with significant circadian phase shifts:
#   1. Load normalized RNA-seq counts
#   2. Map sample IDs to flight/ground condition
#   3. Run limma-voom differential expression (flight vs ground)
#   4. Rank genes by t-statistic
#   5. Build Arabidopsis GO BP + KEGG gene sets from org.At.tair.db
#   6. Run fgsea (fast GSEA)
#   7. Save results + generate enrichment figure
#
# Output:
#   - /mnt/results/data/deg_results/<OSD>_deg.csv  (per-study DEG tables)
#   - /mnt/results/tables/tableS3_enrichment_results.csv
#   - /mnt/results/figures/figS8_enrichment.png/.svg
# ============================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(fgsea)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(GO.db)
  library(ggplot2)
  library(enrichplot)
  library(patchwork)
})
# Load dplyr AFTER AnnotationDbi so AnnotationDbi::select takes precedence
suppressPackageStartupMessages(library(dplyr))
# Ensure select uses AnnotationDbi version
select <- AnnotationDbi::select

# --- Configuration ---
get_path <- function(relative_path) {
  if (dir.exists(file.path(getwd(), relative_path))) {
    return(file.path(getwd(), relative_path))
  } else if (dir.exists(file.path(getwd(), "..", relative_path))) {
    return(file.path(getwd(), "..", relative_path))
  } else if (dir.exists(file.path("/workspace", relative_path))) {
    return(file.path("/workspace", relative_path))
  } else {
    return(file.path(getwd(), relative_path))
  }
}

DATA_DIR <- get_path("genelab_data")

if (dir.exists("/mnt/results")) {
  RESULTS_DIR <- "/mnt/results"
} else if (dir.exists("/results")) {
  RESULTS_DIR <- "/results"
} else {
  if (file.exists("manuscript.tex")) {
    RESULTS_DIR <- getwd()
  } else if (file.exists("../manuscript.tex")) {
    RESULTS_DIR <- file.path(getwd(), "..")
  } else {
    RESULTS_DIR <- getwd()
  }
}

DEG_DIR <- file.path(RESULTS_DIR, "data", "deg_results")
dir.create(DEG_DIR, showWarnings = FALSE, recursive = TRUE)

# Studies with significant phase shifts (RNA-seq with normalized counts)
STUDIES <- list(
  list(osd = "OSD-321", file = "GLDS-321_rna_seq_Normalized_Counts.csv",
       phase_shift = -0.43, pval = 1.8e-7, tissue = "whole_seedling", light = "dark"),
  list(osd = "OSD-193", file = "GLDS-193_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv",
       phase_shift = -0.24, pval = 0.002, tissue = "root", light = "unknown"),
  list(osd = "OSD-38", file = "GLDS-38_rna_seq_Normalized_Counts.csv",
       phase_shift = -0.52, pval = 0.026, tissue = "whole_seedling", light = "dark"),
  list(osd = "OSD-281", file = "GLDS-281_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
       phase_shift = -0.17, pval = 0.038, tissue = "root", light = "unknown")
)

# Circadian-related GO/KEGG terms to highlight
CIRCADIAN_GO <- c("GO:0007623", "GO:0042752", "GO:0042754", "GO:0042753",
                  "GO:0007622", "GO:0043153", "GO:0051592", "GO:0009266")
CIRCADIAN_KEGG <- "ath04710"

cat("=" , rep("=", 58), "\n", sep = "")
cat("GENE SET ENRICHMENT ANALYSIS\n")
cat("=" , rep("=", 58), "\n", sep = "")

# ============================================================
# Build Arabidopsis gene sets from org.At.tair.db
# ============================================================
cat("\nBuilding Arabidopsis gene sets from org.At.tair.db...\n")

# GO Biological Process gene sets
go_bp <- select(org.At.tair.db, keys = keys(org.At.tair.db),
                columns = c("TAIR", "GO", "ONTOLOGY"), keytype = "TAIR")
go_bp <- go_bp[go_bp$ONTOLOGY == "BP" & !is.na(go_bp$GO), ]

# Build gene sets: GO term -> list of TAIR IDs
go_sets <- split(go_bp$TAIR, go_bp$GO)
# Add GO term names using AnnotationDbi::Term
go_term_names <- AnnotationDbi::Term(GO.db::GOTERM)
go_term_map <- unlist(go_term_names)
names(go_sets) <- paste0(names(go_sets), " (", go_term_map[names(go_sets)], ")")

# KEGG pathway gene sets
kegg_anno <- select(org.At.tair.db, keys = keys(org.At.tair.db),
                    columns = c("TAIR", "PATH"), keytype = "TAIR")
kegg_anno <- kegg_anno[!is.na(kegg_anno$PATH), ]
kegg_sets <- split(kegg_anno$TAIR, kegg_anno$PATH)
# Add KEGG pathway names (manual for Arabidopsis)
kegg_names <- c(
  "ath04710" = "Circadian rhythm - plant",
  "ath03010" = "Ribosome",
  "ath04010" = "MAPK signaling pathway",
  "ath04075" = "Plant hormone signal transduction",
  "ath00195" = "Photosynthesis",
  "ath00196" = "Photosynthesis - antenna proteins"
)
named_kegg <- list()
for (pid in names(kegg_sets)) {
  pname <- ifelse(pid %in% names(kegg_names), kegg_names[pid], paste0("KEGG:", pid))
  named_kegg[[paste0(pid, " (", pname, ")")]] <- kegg_sets[[pid]]
}

# Combine all gene sets
all_gene_sets <- c(go_sets, named_kegg)
# Filter: min 15, max 500 genes
all_gene_sets <- all_gene_sets[sapply(all_gene_sets, length) >= 15 &
                                sapply(all_gene_sets, length) <= 500]
cat(sprintf("  Total gene sets: %d (GO BP: %d, KEGG: %d)\n",
            length(all_gene_sets), length(go_sets), length(named_kegg)))

# ============================================================
# Load metadata for sample-to-condition mapping
# ============================================================
metadata <- read.csv(file.path(DATA_DIR, "harmonized_metadata.csv"), stringsAsFactors = FALSE)

# ============================================================
# Process each study
# ============================================================
all_enrichment <- list()
all_deg <- list()

for (study in STUDIES) {
  osd_id <- study$osd
  cat(sprintf("\n--- Processing %s ---\n", osd_id))

  # Load expression data
  expr_path <- file.path(DATA_DIR, osd_id, study$file)
  if (!file.exists(expr_path)) {
    cat(sprintf("  ERROR: File not found: %s\n", expr_path))
    next
  }

  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))

  # Clean gene IDs (remove version suffixes)
  rownames(expr) <- gsub("\\.\\d+$", "", rownames(expr))
  # Filter to AGI codes
  agi_mask <- grepl("^AT[1-5]G\\d{5}$", rownames(expr))
  expr <- expr[agi_mask, ]
  # Remove duplicates
  expr <- expr[!duplicated(rownames(expr)), ]
  cat(sprintf("  After AGI filtering: %d genes\n", nrow(expr)))

  # Map samples to condition
  study_meta <- metadata[metadata$osd_id == osd_id, ]
  cat(sprintf("  Study meta rows: %d\n", nrow(study_meta)))
  sample_to_cond <- setNames(study_meta$condition, study_meta$sample_name)
  cat(sprintf("  sample_to_cond length: %d\n", length(sample_to_cond)))

  # Match columns to metadata
  common_samples <- intersect(colnames(expr), names(sample_to_cond))
  cat(sprintf("  Common samples: %d\n", length(common_samples)))
  if (length(common_samples) > 0) {
    cat(sprintf("  First 5 common: %s\n", paste(head(common_samples, 5), collapse=", ")))
    cat(sprintf("  Their conditions: %s\n", paste(as.character(sample_to_cond[common_samples[1:5]]), collapse=", ")))
  }
  if (length(common_samples) < 6) {
    cat(sprintf("  ERROR: Only %d samples matched to metadata\n", length(common_samples)))
    next
  }

  expr <- expr[, common_samples]
  conditions <- factor(sample_to_cond[common_samples], levels = c("ground_control", "flight"))
  cat(sprintf("  Samples: %d flight, %d ground\n",
              sum(conditions == "flight"), sum(conditions == "ground_control")))

  if (sum(conditions == "flight") < 2 || sum(conditions == "ground_control") < 2) {
    cat("  ERROR: Need at least 2 samples per group\n")
    next
  }

  # --- limma-voom DEG analysis ---
  cat("  Running limma-voom...\n")

  # Create DGEList and filter
  dge <- DGEList(counts = expr)
  # Filter low expression: CPM > 1 in at least min(group) samples
  cpm_matrix <- cpm(dge)
  keep <- rowSums(cpm_matrix > 1) >= min(table(conditions))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("  After filtering: %d genes\n", nrow(dge)))

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")

  # Design matrix
  design <- model.matrix(~conditions)
  colnames(design) <- c("ground_control", "flight_vs_ground")

  # Voom transform
  v <- voom(dge, design, plot = FALSE)

  # Fit linear model
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  # Extract results
  deg_results <- topTable(fit, coef = "flight_vs_ground", number = Inf, sort.by = "t")
  deg_results$TAIR <- rownames(deg_results)

  # Save DEG results
  deg_out <- file.path(DEG_DIR, paste0(osd_id, "_deg.csv"))
  write.csv(deg_results, deg_out, row.names = FALSE)
  cat(sprintf("  Saved DEG results: %d genes (%d significant at FDR<0.05)\n",
              nrow(deg_results), sum(deg_results$adj.P.Val < 0.05)))

  all_deg[[osd_id]] <- deg_results

  # --- fgsea ---
  cat("  Running fgsea...\n")

  # Rank genes by t-statistic
  ranks <- deg_results$t
  names(ranks) <- deg_results$TAIR
  ranks <- sort(ranks, decreasing = TRUE)

  # Run fgsea
  set.seed(42)
  fgsea_res <- fgsea(pathways = all_gene_sets, stats = ranks,
                     minSize = 15, maxSize = 500, nperm = 10000)

  # Add study info
  fgsea_res$osd_id <- osd_id
  fgsea_res$phase_shift <- study$phase_shift
  fgsea_res$tissue <- study$tissue
  fgsea_res$light_regime <- study$light

  # Mark circadian-related pathways
  fgsea_res$is_circadian <- grepl("GO:0007623|GO:0042752|GO:0042754|GO:0042753|circadian",
                                  fgsea_res$pathway, ignore.case = TRUE)

  all_enrichment[[osd_id]] <- fgsea_res

  # Report top results
  sig <- fgsea_res[fgsea_res$padj < 0.05, ]
  cat(sprintf("  Significant pathways (FDR<0.05): %d\n", nrow(sig)))
  if (nrow(sig) > 0) {
    circ <- sig[sig$is_circadian, ]
    if (nrow(circ) > 0) {
      cat("  Circadian-related significant pathways:\n")
      print(circ[, c("pathway", "NES", "padj")])
    }
  }
}

# ============================================================
# Combine and save enrichment results
# ============================================================
cat("\n" , rep("=", 59), "\n", sep = "")
cat("Combining results...\n")

combined_enrichment <- do.call(rbind, all_enrichment)
combined_enrichment <- as.data.frame(combined_enrichment)
# Convert list columns
combined_enrichment$leadingEdge <- sapply(combined_enrichment$leadingEdge,
                                          function(x) paste(unlist(x), collapse = ";"))

# Save full enrichment table
enrichment_out <- file.path(RESULTS_DIR, "tables", "tableS3_enrichment_results.csv")
write.csv(combined_enrichment, enrichment_out, row.names = FALSE)
cat(sprintf("Saved tableS3_enrichment_results.csv (%d rows)\n", nrow(combined_enrichment)))

# ============================================================
# Generate enrichment figure (figS8)
# ============================================================
cat("\nGenerating figS8...\n")

# Prepare data for plotting
plot_data <- combined_enrichment[combined_enrichment$padj < 0.1, ]
if (nrow(plot_data) == 0) {
  cat("  No pathways with padj < 0.1, using top 20 by padj per study\n")
  plot_data <- combined_enrichment %>%
    group_by(osd_id) %>%
    arrange(padj) %>%
    slice_head(n = 20) %>%
    ungroup()
}

# Panel A: Circadian pathway NES across studies
circadian_data <- combined_enrichment[combined_enrichment$is_circadian |
  grepl("circadian", combined_enrichment$pathway, ignore.case = TRUE), ]

# Panel B: Top enriched GO BP terms (by NES) across studies
top_pathways <- combined_enrichment %>%
  filter(padj < 0.1) %>%
  group_by(pathway) %>%
  summarise(mean_NES = mean(NES), min_padj = min(padj), n_studies = n()) %>%
  arrange(desc(abs(mean_NES))) %>%
  slice_head(n = 15)

# Create multi-panel figure
# (patchwork loaded at top)

# Panel A: NES for circadian pathways across studies
pA <- if (nrow(circadian_data) > 0) {
  ggplot(circadian_data, aes(x = osd_id, y = NES, fill = padj < 0.05)) +
    geom_col(position = "dodge") +
    facet_wrap(~pathway, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("grey70", "#FF9400"), name = "FDR < 0.05") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "Study", y = "Normalized Enrichment Score (NES)",
         title = "A  Circadian pathway enrichment") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(size = 7),
          legend.position = "right")
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No circadian pathways enriched") +
    theme_void() + labs(title = "A  Circadian pathway enrichment")
}

# Panel B: Top 15 pathways by mean |NES| (dot plot)
pB <- if (nrow(top_pathways) > 0) {
  # Get per-study NES for these pathways
  top_pathway_names <- top_pathways$pathway
  top_data <- combined_enrichment[combined_enrichment$pathway %in% top_pathway_names, ]
  # Shorten pathway names
  top_data$pathway_short <- substr(top_data$pathway, 1, 50)

  ggplot(top_data, aes(x = NES, y = reorder(pathway_short, NES),
                       color = padj < 0.05, size = -log10(padj))) +
    geom_point() +
    facet_wrap(~osd_id, nrow = 1) +
    scale_color_manual(values = c("grey60", "#0279EE"), name = "FDR < 0.05") +
    labs(x = "NES", y = "", title = "B  Top enriched pathways") +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 6),
          strip.text = element_text(size = 8),
          legend.position = "right")
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No pathways enriched") +
    theme_void() + labs(title = "B  Top enriched pathways")
}

# Panel C: GSEA running enrichment score for circadian rhythm in OSD-321
pC <- tryCatch({
  # Get ranks for OSD-321
  if ("OSD-321" %in% names(all_deg)) {
    deg321 <- all_deg[["OSD-321"]]
    ranks321 <- sort(setNames(deg321$t, deg321$TAIR), decreasing = TRUE)

    # Find circadian pathway
    circ_pathway <- grep("GO:0007623|circadian", names(all_gene_sets),
                         ignore.case = TRUE, value = TRUE)
    if (length(circ_pathway) > 0) {
      circ_pathway <- circ_pathway[1]
      # Run fgsea with plot data
      plotRes <- fgsea(all_gene_sets[circ_pathway], ranks321, nperm = 10000)

      # Create running enrichment score plot manually
      gene_order <- names(ranks321)
      pathway_genes <- all_gene_sets[[circ_pathway]]
      hits <- gene_order %in% pathway_genes

      # Compute running ES
      n <- length(gene_order)
      nh <- sum(hits)
      if (nh > 0) {
        running <- cumsum(ifelse(hits, 1/nh, -1/(n-nh)))
        running_df <- data.frame(rank = 1:n, ES = running)

        ggplot(running_df, aes(x = rank, y = ES)) +
          geom_line(color = "#0279EE", linewidth = 0.8) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          labs(x = "Gene rank", y = "Running enrichment score",
               title = paste0("C  GSEA: ", substr(circ_pathway, 1, 45), "\n(OSD-321, NES=",
                              round(plotRes$NES, 2), ", FDR=", signif(plotRes$padj, 3), ")")) +
          theme_bw(base_size = 10)
      } else {
        ggplot() + annotate("text", x=0.5, y=0.5, label="No circadian genes in ranked list") +
          theme_void() + labs(title = "C  GSEA circadian pathway (OSD-321)")
      }
    } else {
      ggplot() + annotate("text", x=0.5, y=0.5, label="Circadian pathway not found") +
        theme_void() + labs(title = "C  GSEA circadian pathway (OSD-321)")
    }
  } else {
    ggplot() + annotate("text", x=0.5, y=0.5, label="OSD-321 not available") +
      theme_void() + labs(title = "C  GSEA circadian pathway (OSD-321)")
  }
}, error = function(e) {
  ggplot() + annotate("text", x=0.5, y=0.5, label = paste("Error:", e$message)) +
    theme_void() + labs(title = "C  GSEA circadian pathway (OSD-321)")
})

# Combine panels
final_plot <- (pA / pB / pC) + plot_layout(heights = c(1.2, 1, 0.8))
final_plot <- final_plot + plot_annotation(title = "Gene set enrichment: spaceflight vs ground in phase-shifted studies",
                                            theme = theme(plot.title = element_text(size = 14, face = "bold")))

# Save
svg_path <- file.path(RESULTS_DIR, "figures", "figS8_enrichment.svg")
png_path <- file.path(RESULTS_DIR, "figures", "figS8_enrichment.png")
ggsave(svg_path, final_plot, width = 14, height = 18, dpi = 300)
ggsave(png_path, final_plot, width = 14, height = 18, dpi = 300)
cat("  Saved figS8_enrichment.svg and figS8_enrichment.png\n")

cat("\n" , rep("=", 59), "\n", sep = "")
cat("ENRICHMENT ANALYSIS COMPLETE\n")
cat("=" , rep("=", 58), "\n", sep = "")
