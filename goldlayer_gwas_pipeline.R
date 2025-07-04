# ==========================================================
# GoldLayer-305: Full Epigenomic CRISPR/gRNA Validation Pipeline
# Author: Othman Sayfuldeen Ismaeil Mohammed
# Last Update: 2025-06-19
#
# Purpose: Retrieve significant GWAS SNPs for prostate cancer,
# design CRISPR gRNAs, and validate them across multi-layer
# epigenomic data (cCREs, eQTLs, H3K27ac peaks, ChIP-seq signal).
#
# Requirements:
#   - R >= 4.1
#   - See required_pkgs at top for all dependencies
#   - Download input files as specified in comments
#
# ================================================================
# Copyright (c) 2025 Othman S. I. Mohammed
#
# Permission is hereby granted for academic, non-commercial use only.
# Any commercial use, including incorporation into proprietary software
# or services, requires explicit permission and a separate licensing agreement.
#
# For inquiries about commercial licensing, please contact: [admin@leopard.ly]
#
# All other rights reserved.
# ==========================================================
# ============================================
# 1. Installation & Library Setup (All Systems)
# ============================================

local_hg38 <- "D:/manhaj/R/rpackages/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz"  # !! Edit file paths as needed for your system !!

required_pkgs <- c(
  "gwasrapidd", "dplyr", "tidyr", "purrr", "BiocManager",
  "BSgenome.Hsapiens.UCSC.hg38", "crisprDesign", "crisprBowtie",
  "crisprScore", "crisprBase", "Rbowtie", "rmarkdown", "knitr",
  "httr", "jsonlite", "tibble", "AnnotationHub", "GenomicRanges",
  "rtracklayer", "stringr", "readr", "ggplot2"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BSgenome.Hsapiens.UCSC.hg38" && file.exists(local_hg38)) {
      install.packages(local_hg38, repos = NULL, type = "source")
    } else if (pkg %in% c("crisprDesign","crisprBowtie","crisprScore",
                          "crisprBase","BSgenome.Hsapiens.UCSC.hg38","Rbowtie",
                          "AnnotationHub","GenomicRanges","rtracklayer")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_pkgs, install_if_missing))
options(timeout = 1200)

# ============================================
# 2. Retrieve Significant Prostate-Cancer SNPs
# ============================================

fetch_prostate_assocs <- function(p_cut = 5e-8) {
  try({
    m <- get_associations(efo_id = "EFO_0001663", verbose = FALSE)
    if (nrow(m@associations) > 0) return(m)
  }, silent = TRUE)
  try({
    m <- get_associations(efo_trait = "prostate carcinoma", verbose = FALSE)
    if (nrow(m@associations) > 0) return(m)
  }, silent = TRUE)
  api <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"
  url <- sprintf("%s/traits/EFO_0001663/associations?p_upper=%g&size=5000", api, p_cut)
  pull_all <- function(u) {
    out <- list()
    repeat {
      r  <- GET(u, add_headers(Accept = "application/hal+json, application/json"))
      stop_for_status(r)
      js <- content(r, as = "parsed", type = "application/json")
      if (!"associations" %in% names(js$`_embedded`)) break
      page <- flatten(js$`_embedded`$associations)
      out[[length(out) + 1]] <- page
      nxt <- js$`_links`$`next`$href
      if (is.null(nxt) || nxt == u) break
      u <- nxt
    }
    bind_rows(out)
  }
  ss <- pull_all(url)
  if (nrow(ss) == 0)
    stop("Summary-Statistics API returned no hits – check connectivity.")
  fake <- methods::new(
    "associations",
    associations = tibble::tibble(
      association_id = paste0(ss$variant_id, ":", ss$study_accession),
      pvalue = as.numeric(ss$p_value)),
    risk_alleles = tibble::tibble(
      association_id = paste0(ss$variant_id, ":", ss$study_accession),
      variant_id = ss$variant_id)
  )
  fake
}
catalog_assocs <- fetch_prostate_assocs()
assoc_df <- catalog_assocs@associations %>% drop_na(pvalue) %>% filter(pvalue < 5e-8)
sig_assoc_ids <- assoc_df$association_id
sig_variants <- catalog_assocs@risk_alleles %>%
  filter(association_id %in% sig_assoc_ids) %>%
  pull(variant_id) %>% unique()
batch_size <- 100
batches <- split(sig_variants, ceiling(seq_along(sig_variants)/batch_size))
variant_info_list <- lapply(batches, function(vids) {
  get_variants(variant_id = vids, verbose = FALSE)
})
variant_tbl <- bind_rows(lapply(variant_info_list, function(x) x@variants))
genemap_tbl <- bind_rows(lapply(variant_info_list, function(x) x@genomic_contexts))
mapped_genes <- genemap_tbl %>%
  filter(is_mapped_gene) %>%
  select(variant_id, gene_name) %>%
  distinct() %>%
  group_by(variant_id) %>%
  summarise(mapped_gene = paste(unique(gene_name), collapse = ", "), .groups = "drop")
variant_tbl <- variant_tbl %>%
  left_join(mapped_genes, by = "variant_id") %>%
  select(variant_id, chromosome_name, chromosome_position, mapped_gene) %>%
  arrange(chromosome_name, chromosome_position)

# ============================================
# 3. Extract ±50 bp Genomic Sequences
# ============================================

hg38 <- BSgenome.Hsapiens.UCSC.hg38
variant_tbl <- variant_tbl %>%
  mutate(
    chrom = case_when(
      chromosome_name %in% c("23", "X")        ~ "chrX",
      chromosome_name %in% c("24", "Y")        ~ "chrY",
      chromosome_name %in% c("25", "MT", "M")  ~ "chrM",
      TRUE                                      ~ paste0("chr", chromosome_name)
    )
  ) %>%
  filter(chrom %in% seqnames(hg38)) %>%
  rowwise() %>%
  mutate(
    start_seq = max(1, chromosome_position - 50),
    end_seq   = min(seqlengths(hg38)[[chrom]], chromosome_position + 50)
  ) %>%
  ungroup()
variant_seqs <- Biostrings::getSeq(hg38,
                                   names  = variant_tbl$chrom,
                                   start  = variant_tbl$start_seq,
                                   end    = variant_tbl$end_seq)
names(variant_seqs) <- variant_tbl$variant_id

# ============================================
# 4. Design SpCas9 gRNAs for Each SNP Window
# ============================================

data("SpCas9", package = "crisprBase")
target_gr <- GenomicRanges::GRanges(
  seqnames   = variant_tbl$chrom,
  ranges     = IRanges::IRanges(
    start = variant_tbl$start_seq,
    end   = variant_tbl$end_seq),
  variant_id = variant_tbl$variant_id
)
guideSet <- findSpacers(
  target_gr,
  bsgenome       = hg38,
  crisprNuclease = SpCas9
)

# ============================================
# 5. Annotation and Genomic Ranges Workflow
# ============================================

ah <- AnnotationHub()
cCRE_query <- query(ah, c("encodeCcreCombined", "hg38"))
if (length(cCRE_query)) {
  cCREs <- cCRE_query[[1]]
} else {
  bb_file <- "encodeCcreCombined.bb"
  if (!file.exists(bb_file)) {
    bb_url <- "http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"
    download.file(bb_url, destfile = bb_file, mode = "wb")
  }
  cCREs <- rtracklayer::import(bb_file)
}
eqtl_file <- "D:/manhaj/R/rpackages/GTEx_Analysis_v8_eQTL/Prostate.signifpairs.txt" # !! Edit file paths as needed for your system !!
if (!file.exists(eqtl_file)) {
  stop("❌ File not found: ", eqtl_file, "\nPlace the correct file in this location and retry.")
}
gtex_eqtl <- read.table(eqtl_file, header = TRUE, sep = "\t")
if ("variant_id" %in% colnames(gtex_eqtl)) {
  parts          <- stringr::str_split_fixed(gtex_eqtl$variant_id, "_", 5)
  gtex_eqtl$chr  <- parts[, 1]
  gtex_eqtl$pos  <- as.integer(parts[, 2])
}
variant_tbl$chr <- paste0("chr", variant_tbl$chromosome_name)
variant_tbl$pos <- variant_tbl$chromosome_position
variant_gr <- GRanges(seqnames = variant_tbl$chr,
                      ranges   = IRanges(start = variant_tbl$pos, width = 1))
if (!("chr" %in% names(gtex_eqtl)) || !("pos" %in% names(gtex_eqtl))) {
  if ("variant_id" %in% names(gtex_eqtl)) {
    parts <- stringr::str_split_fixed(gtex_eqtl$variant_id, "_", 5)
    gtex_eqtl$chr <- parts[, 1]
    gtex_eqtl$pos <- as.integer(parts[, 2])
  } else if (all(c("chromosome", "position") %in% names(gtex_eqtl))) {
    gtex_eqtl$chr <- gtex_eqtl$chromosome
    gtex_eqtl$pos <- gtex_eqtl$position
  } else {
    stop("❌ Cannot extract chr and pos columns from gtex_eqtl: no recognizable columns found.")
  }
}
gtex_eqtl_gr <- GRanges(
  seqnames = gtex_eqtl$chr,
  ranges   = IRanges(start = gtex_eqtl$pos, width = 1)
)
variant_tbl$cCRE_overlap  <- overlapsAny(variant_gr, cCREs)
variant_tbl$Prostate_eQTL <- overlapsAny(variant_gr, gtex_eqtl_gr)

# ============================================
# 6. Add Azimuth On-Target Scores and Filtering
# ============================================

guideSet <- addOnTargetScores(guideSet, methods="azimuth")
hi_guides <- guideSet[guideSet$score_azimuth > 0.5]
library(dplyr)
flag_tbl <- variant_tbl %>% 
  mutate(region = paste0("region_", row_number())) %>% 
  select(region, cCRE_overlap, Prostate_eQTL)
hi_df <- as.data.frame(hi_guides) %>%  
  left_join(flag_tbl, by = "region")
n_cCRE   <- sum(hi_df$cCRE_overlap, na.rm = TRUE)
n_eQTL   <- sum(hi_df$Prostate_eQTL, na.rm = TRUE)
n_both   <- sum(hi_df$cCRE_overlap & hi_df$Prostate_eQTL, na.rm = TRUE)
saveRDS(hi_guides, file = "high_efficacy_gRNAs.rds")

# ============================================
# 7. Epigenomic Multi-Layer Validation (Gold 305)
# ============================================

goldlayer_305 <- hi_df %>% filter(cCRE_overlap & Prostate_eQTL)
cat("GoldLayer guide count:", nrow(goldlayer_305), "\n")
guides_gr <- GRanges(
  seqnames = goldlayer_305$seqnames,
  ranges   = IRanges(start = goldlayer_305$start, end = goldlayer_305$end)
)
h3_bed_path <- "D:/manhaj/GWAS/proposal 2025 may/sources/GSE105760_ENCFF341RJV_peaks_GRCh38.bed" # !! Edit file paths as needed for your system !!
h3_df <- read_tsv(
  h3_bed_path,
  col_names = c("chr", "start", "end", "name", "score", "strand",
                "signalValue", "pValue", "qValue", "peak"),
  col_types = "ciiidddddi"
)
h3_gr <- GRanges(
  seqnames = h3_df$chr,
  ranges   = IRanges(start = h3_df$start + 1, end = h3_df$end)
)
goldlayer_305$H3K27ac_peak <- overlapsAny(guides_gr, h3_gr)
h3_rep_bed <- "D:/manhaj/GWAS/proposal 2025 may/sources/GSE105760_ENCFF776KYJ_replicated_peaks_GRCh38.bed.gz" # !! Edit file paths as needed for your system !!
h3_rep_df <- read_tsv(
  h3_rep_bed, comment = "#",
  col_names = c("chr", "start", "end", "name", "score", "strand",
                "signalValue", "pValue", "qValue", "peak"),
  col_types = "ciiidddddi"
)
h3_rep_gr <- GRanges(
  seqnames = h3_rep_df$chr,
  ranges   = IRanges(start = h3_rep_df$start + 1, end = h3_rep_df$end)
)
goldlayer_305$H3K27ac_peak_rep <- overlapsAny(guides_gr, h3_rep_gr)
goldlayer_305$Epigenomic_4layer <- goldlayer_305$H3K27ac_peak_rep
n_hit <- sum(goldlayer_305$Epigenomic_4layer)
cat("GoldLayer guides in high-confidence IDR peaks =", n_hit, "of", nrow(goldlayer_305), "\n")
pval_bw <- "D:/manhaj/GWAS/proposal 2025 may/sources/GSE105760_ENCFF808YKT_signal_p-value_GRCh38.bigWig" # !! Edit file paths as needed for your system !!
p_vec <- import(pval_bw, which = guides_gr, as = "NumericList")
goldlayer_305$negLog10P <- vapply(p_vec, function(x)
  if (length(x)) max(x) else NA_real_, 0)
summary_tbl <- goldlayer_305 %>%
  summarise(
    Total_GoldLayer   = n(),
    H3K27ac_raw       = sum(H3K27ac_peak),
    H3K27ac_rep       = sum(H3K27ac_peak_rep),
    Percent_Rep       = round(100 * H3K27ac_rep / Total_GoldLayer, 1),
    Median_negLog10P  = round(median(negLog10P, na.rm = TRUE), 2),
    Guides_strong     = sum(negLog10P > 10, na.rm = TRUE)
  )
print(summary_tbl)
write_csv(goldlayer_305 %>% filter(Epigenomic_4layer),
          "Epigenomic_validated_guides_IDR.csv")
write_csv(goldlayer_305, "GoldLayer_305_fullEpi.csv")
epibar <- tibble(
  Category = c("All_GoldLayer_305", "H3K27ac_raw", "IDR_Replicated_Peak", "Signal>10"),
  Count = c(nrow(goldlayer_305),
            sum(goldlayer_305$H3K27ac_peak),
            sum(goldlayer_305$H3K27ac_peak_rep),
            sum(goldlayer_305$negLog10P > 10, na.rm = TRUE))
)
p <- ggplot(epibar, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, show.legend = FALSE) +
  labs(title = "Full Epigenomic Validation (GoldLayer 305)",
       y = "Number of Guides", x = "") +
  theme_minimal(base_size = 13)
ggsave("epigenomic_barplot_full.png", plot = p, width = 5, height = 4, dpi = 300)
top_guides <- goldlayer_305 %>%
  filter(Epigenomic_4layer) %>%
  arrange(desc(negLog10P)) %>%
  select(seqnames, start, protospacer, score_azimuth, negLog10P) %>%
  head(10)
print(top_guides)
cat("\n🔱 full_epi_validation: All computational and visual validation layers complete, tables saved.\n")
# ============================================
# 8. Stepwise Summary & Output for Validation
# ============================================
# This section provides a concise, step-by-step display of key results and diagnostics for inspection and documentation.

# 1. Show the first 10 variants and their mapped genes
cat("\n1. Top 10 Lead Variants and Mapped Genes:\n")
print(head(variant_tbl, 10))

# 2. Display the total number of lead variants after filtering
cat("\n2. Total Number of Lead Variants After Filtering:\n")
cat("Total lead variants kept:", nrow(variant_tbl), "\n")

# 3. Display the number of extracted genomic sequences (+/- 50bp)
cat("\n3. Number of Extracted Sequences (+/- 50 bp):\n")
cat("Sequences extracted:", length(variant_seqs), "\n")

# 4. Show the total number of designed gRNAs across all loci
cat("\n4. Total Number of Designed gRNAs Across All Loci:\n")
cat("gRNAs identified across all loci:", length(guideSet), "\n")

# 5. Show the first 10 gRNAs with their Azimuth scores
cat("\n5. Top 10 gRNAs with Azimuth Scores:\n")
print(as.data.frame(guideSet)[1:10, c("protospacer", "score_azimuth")])

# 6. Plot the distribution of all Azimuth scores before filtering
cat("\n6. Distribution of All Azimuth Scores (Before Filtering):\n")
hist(guideSet$score_azimuth, breaks = 30,
     main = "Distribution of All Azimuth Scores",
     xlab = "Azimuth Score", col = "lightblue", border = "white")

# 7. Display the count of high-efficacy gRNAs (Azimuth > 0.5)
cat("\n7. Number of High-Efficacy gRNAs (Azimuth > 0.5):\n")
cat("High-efficacy gRNAs (Azimuth > 0.5):", length(hi_guides), "\n")

# 8. Plot the distribution of Azimuth scores for high-efficacy gRNAs only
cat("\n8. Distribution of Azimuth Scores for Filtered High-Efficacy gRNAs (Azimuth > 0.5):\n")
hist(hi_guides$score_azimuth, breaks = 20,
     main = "Distribution of Filtered High-Efficacy gRNAs (Azimuth > 0.5)",
     xlab = "Azimuth Score", col = "seagreen", border = "white")

# 9. Show the first 6 rows of the filtered gRNAs table (with genomic flags)
cat("\n9. First 6 Rows of Filtered gRNAs Table (with Genomic Flags):\n")
print(head(hi_df))

# 10. Show overlap statistics between high-efficacy gRNAs and cCRE/eQTL features
cat("\n10. Overlap Statistics for High-Efficacy gRNAs with cCRE and eQTL:\n")
cat("High-efficacy gRNAs =", nrow(hi_df), "\n")
cat("∩ cCRE =", n_cCRE, "\n")
cat("∩ eQTL =", n_eQTL, "\n")
cat("∩ cCRE ∩ eQTL =", n_both, "\n")
# ============================================
# 9. Epigenomic Multi-Layer Validation Results
# ============================================
# This section summarizes the full results from the multi-layer epigenomic validation, as applied to the GoldLayer 305 gRNAs.

# 9.1. Display total number of GoldLayer guides (high-efficacy, cCRE+eQTL overlap)
cat("\n9.1. Number of GoldLayer Guides (cCRE & eQTL overlap):\n")
cat("GoldLayer guide count:", nrow(goldlayer_305), "\n")

# 9.2. Show number of guides overlapping raw H3K27ac peaks (ENCODE)
cat("\n9.2. Guides Overlapping H3K27ac Peaks (ENCODE raw):\n")
cat("Overlapping raw H3K27ac peaks:", sum(goldlayer_305$H3K27ac_peak), "\n")

# 9.3. Show number of guides overlapping replicated (IDR) H3K27ac peaks
cat("\n9.3. Guides Overlapping Replicated (IDR) H3K27ac Peaks:\n")
cat("Overlapping replicated H3K27ac peaks:", sum(goldlayer_305$H3K27ac_peak_rep), "\n")

# 9.4. Show percentage of GoldLayer guides confirmed by replicated IDR peaks
cat("\n9.4. Percentage of GoldLayer Guides Confirmed by Replicated Peaks:\n")
cat("Percent confirmed (IDR):",
    round(100 * sum(goldlayer_305$H3K27ac_peak_rep) / nrow(goldlayer_305), 1), "%\n")

# 9.5. Display median and range of -log10(p-value) ChIP-seq signal for all GoldLayer guides
cat("\n9.5. Median and Max -log10(p-value) Signal (H3K27ac ChIP-seq):\n")
cat("Median -log10(p):", round(median(goldlayer_305$negLog10P, na.rm = TRUE), 2), "\n")
cat("Max -log10(p):", round(max(goldlayer_305$negLog10P, na.rm = TRUE), 2), "\n")

# 9.6. Show number of GoldLayer guides with strong epigenomic evidence (signal > 10)
cat("\n9.6. GoldLayer Guides with Strong ChIP-seq Signal (-log10(p) > 10):\n")
cat("Guides with signal > 10:", sum(goldlayer_305$negLog10P > 10, na.rm = TRUE), "\n")

# 9.7. Display summary statistics table for all layers
cat("\n9.7. Summary Table for All Validation Layers:\n")
print(summary_tbl)

# 9.8. Display barplot of all validation layers (optional: ensure plot file is saved)
cat("\n9.8. Barplot of Validation Layers Saved as 'epigenomic_barplot_full.png'\n")

# 9.9. Show the top 10 epigenomically validated guides, sorted by ChIP-seq signal
cat("\n9.9. Top 10 Epigenomic-Validated GoldLayer Guides (by ChIP-seq signal):\n")
print(top_guides)

# 9.10. Export final tables (outputs confirmed in code, mention to PI)
cat("\n9.10. Final validated guides tables exported:\n")
cat(" - Epigenomic_validated_guides_IDR.csv\n")
cat(" - GoldLayer_305_fullEpi.csv\n")

cat("\n🔱 Epigenomic multi-layer validation complete: Results available in the summary table, top guide list, and exported files.\n")
