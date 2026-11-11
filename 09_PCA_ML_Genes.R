# PCA on ML-selected NULISA proteins (from 08_Clustering_Feature_Selection.R)
#
# Works with 08_Clustering_Feature_Selection.R (--panel cns | immune | both).
# With the same -i and --panel as used for 08, the selected-feature file and
# output folder are found automatically:
#   --panel both   -> <input>/Combined_CNS_Immune_panels/Without_tears/...
#   --panel cns    -> <input>/CNS_panel/Without_tears/...
#   --panel immune -> <input>/Immune_panel/Without_tears/...
#   reading  .../08_Clustering_Feature_Selection[_<Matrix>][_logLoss]/selected_features.txt
#   writing  .../09_PCA_ML_Genes[_<Matrix>][_logLoss]/
# --matrix CSF | Serum | Plasma (same as 08; CSF folders have no suffix).
# --metric logLoss (default) | ROC: same as used for 08; logLoss folders get "_logLoss".
#
# Rscript Script/09_PCA_ML_Genes.R -i CNS_immune/Results --panel cns
# Rscript Script/09_PCA_ML_Genes.R -i CNS_immune/Results --panel cns --matrix Serum
# Rscript Script/09_PCA_ML_Genes.R -i CNS_immune/Results --panel cns --metric ROC
# Rscript Script/09_PCA_ML_Genes.R -i CNS_immune/Results --panel both
# (-f / -o override the automatic paths)

suppressMessages(library(optparse))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Note: install 'ggrepel' for non-overlapping biplot labels.")
}


label_repel_geom <- function(data, mapping, ...) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = data,
      mapping = mapping,
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "italic",
      color = "gray10",
      segment.color = "gray65",
      segment.size = 0.3,
      min.segment.length = 0,
      box.padding = 0.45,
      point.padding = 0.25,
      force = 2.5,
      max.overlaps = Inf,
      ...
    )
  } else {
    geom_text(data = data, mapping = mapping, inherit.aes = FALSE, size = 3.2, fontface = "italic", color = "gray10", ...)
  }
}

option_list <- list(
  make_option(
    c("-i", "--input"), type = "character", default = "CNS_immune/Results",
    help = paste(
      "Panel base directory (as for 08), OR comma-separated paths.",
      "--panel both: cns_npq.xlsx,immune_npq.xlsx,metadata.xlsx;",
      "--panel cns/immune: npq.xlsx,metadata.xlsx"
    )
  ),
  make_option(
    c("--panel"), type = "character", default = "both",
    help = "Panel(s) used in 08: 'cns', 'immune' or 'both' [default: %default]"
  ),
  make_option(
    c("--matrix"), type = "character", default = "CSF",
    help = "Sample matrix used in 08: 'CSF', 'Serum' or 'Plasma' (case-insensitive) [default: %default]"
  ),
  make_option(
    c("--metric"), type = "character", default = "logLoss",
    help = paste(
      "Tuning metric used in 08: 'logLoss' (default; reads/writes folders with suffix '_logLoss')",
      "or 'ROC' (folders without suffix) [default: %default]"
    )
  ),
  make_option(
    c("-f", "--features"), type = "character", default = NULL,
    help = "selected_features.txt from 08 [default: derived from -i and --panel]"
  ),
  make_option(
    c("-o", "--output"), type = "character", default = NULL,
    help = "Output directory [default: derived from -i and --panel]"
  ),
  make_option(c("--require_k2"), type = "logical", default = TRUE),
  make_option(c("--pr_breaks"), type = "character", default = "0,0.5,1,2,4,Inf")
)

opt <- parse_args(OptionParser(option_list = option_list))

panel_choice <- tolower(opt$panel)
if (!panel_choice %in% c("cns", "immune", "both")) {
  stop("--panel must be 'cns', 'immune' or 'both'.")
}
panel_label <- switch(panel_choice, cns = "CNS", immune = "Immune", both = "CNS + Immune")

# Sample matrix: NPQ files use CSF/SERUM/PLASMA, metadata uses CSF/Serum/Plasma
matrix_key <- toupper(trimws(opt$matrix))
if (!matrix_key %in% c("CSF", "SERUM", "PLASMA")) {
  stop("--matrix must be 'CSF', 'Serum' or 'Plasma'.")
}
matrix_meta <- switch(matrix_key, CSF = "CSF", SERUM = "Serum", PLASMA = "Plasma")
matrix_suffix <- if (matrix_meta == "CSF") "" else paste0("_", matrix_meta)   # folder suffix, as in 08

# Tuning metric of the 08 run: logLoss folders carry "_logLoss" (after the matrix suffix), as in 08
tune_metric <- switch(tolower(trimws(opt$metric)),
  "roc" = "ROC",
  "logloss" = "logLoss",
  stop("--metric must be 'logLoss' or 'ROC'.")
)
metric_suffix <- if (tune_metric == "logLoss") "_logLoss" else ""

# Same logic as 08_Clustering_Feature_Selection.R
parse_panel_paths <- function(input_arg, panel_choice) {
  subdir <- "Without_tears"
  panel_dir <- c(CNS = "CNS_panel", Immune = "Immune_panel")
  prefixes <- switch(panel_choice, cns = "CNS", immune = "Immune", both = c("CNS", "Immune"))

  if (grepl(",", input_arg, fixed = TRUE)) {
    parts <- trimws(strsplit(input_arg, ",", fixed = TRUE)[[1]])
    if (panel_choice == "both") {
      if (length(parts) != 3L) {
        stop("--panel both with file paths needs: cns_npq.xlsx,immune_npq.xlsx,metadata.xlsx")
      }
      return(list(npq = list(CNS = parts[1], Immune = parts[2]), metadata = parts[3]))
    }
    if (length(parts) == 2L) {
      return(list(npq = setNames(list(parts[1]), prefixes), metadata = parts[2]))
    }
    if (length(parts) == 3L) {
      npq_path <- if (panel_choice == "cns") parts[1] else parts[2]
      return(list(npq = setNames(list(npq_path), prefixes), metadata = parts[3]))
    }
    stop("--panel ", panel_choice, " with file paths needs: npq.xlsx,metadata.xlsx")
  }

  npq <- lapply(prefixes, function(p) {
    file.path(input_arg, panel_dir[[p]], subdir, "00_Initialization", "npq_counts.xlsx")
  })
  names(npq) <- prefixes
  metadata <- file.path(input_arg, panel_dir[[prefixes[1]]], subdir, "00_Initialization",
                        "all_participants_IDs.xlsx")
  list(npq = npq, metadata = metadata)
}

# Panel result folder used by 08 (and now 09)
panel_result_dir <- function(input_arg, panel_choice) {
  sub_dir <- switch(panel_choice,
    both = "Combined_CNS_Immune_panels",
    cns = "CNS_panel",
    immune = "Immune_panel"
  )
  file.path(input_arg, sub_dir, "Without_tears")
}

npq_to_wide <- function(npq_path, panel_prefix) {
  read_excel(npq_path) %>%
    filter(toupper(SampleMatrixType) == matrix_key, SampleType == "Sample") %>%
    mutate(Feature = paste0(panel_prefix, "_", Target)) %>%
    select(SampleName, Feature, NPQ) %>%
    pivot_wider(names_from = Feature, values_from = NPQ) %>%
    as.data.frame()
}

panel_paths <- parse_panel_paths(opt$input, panel_choice)
is_file_input <- grepl(",", opt$input, fixed = TRUE)

features_path <- opt$features
if (is.null(features_path) || !nzchar(features_path)) {
  if (is_file_input) stop("With comma-separated --input, please give -f (selected_features.txt from 08).")
  features_path <- file.path(panel_result_dir(opt$input, panel_choice),
                             paste0("08_Clustering_Feature_Selection", matrix_suffix, metric_suffix), "selected_features.txt")
}
if (!file.exists(features_path)) {
  stop("Selected-feature file not found: ", features_path,
       "\nRun 08_Clustering_Feature_Selection.R with the same -i/--panel/--matrix/--metric first, or give -f.")
}

output_dir <- opt$output
if (is.null(output_dir) || !nzchar(output_dir)) {
  output_dir <- if (is_file_input) {
    paste0("09_PCA_ML_Genes", matrix_suffix, metric_suffix, "_", panel_choice)
  } else {
    file.path(panel_result_dir(opt$input, panel_choice), paste0("09_PCA_ML_Genes", matrix_suffix, metric_suffix))
  }
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Panel(s): ", panel_label, "; matrix: ", matrix_meta, "; 08 tuning metric: ", tune_metric)
for (p in names(panel_paths$npq)) message("Input (", p, " NPQ): ", panel_paths$npq[[p]])
message("Metadata: ", panel_paths$metadata)
message("Selected features: ", features_path)
message("Output: ", output_dir)

selected_features <- trimws(readLines(features_path, warn = FALSE))
selected_features <- selected_features[nzchar(selected_features)]
if (length(selected_features) < 2L) {
  stop("Need at least 2 features for PCA.")
}

message("Loading ", panel_label, " NPQ (", matrix_meta, " samples)...")
wide_list <- lapply(names(panel_paths$npq), function(p) npq_to_wide(panel_paths$npq[[p]], p))
common_samples <- Reduce(intersect, lapply(wide_list, function(w) w$SampleName))
npq_wide <- Reduce(
  function(a, b) inner_join(a, b, by = "SampleName"),
  lapply(wide_list, function(w) w[w$SampleName %in% common_samples, , drop = FALSE])
)
rownames(npq_wide) <- npq_wide$SampleName
npq_wide$SampleName <- NULL

missing_feat <- setdiff(selected_features, colnames(npq_wide))
if (length(missing_feat) > 0L) {
  stop("Selected features not in the ", panel_label, " NPQ matrix: ", paste(missing_feat, collapse = ", "),
       "\nCheck that --panel matches the 08 run that produced ", features_path)
}

meta <- read_excel(panel_paths$metadata) %>%
  filter(toupper(Material) == matrix_key) %>%
  mutate(
    progression_rate = suppressWarnings(as.numeric(progression_rate)),
    k2 = as.character(k2),
    across(any_of(c("slow_vital_capacity", "Nfl")), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  select(Tube_ID, k2, progression_rate, any_of(c("slow_vital_capacity", "Nfl"))) %>%
  filter(Tube_ID %in% rownames(npq_wide))

if (isTRUE(opt$require_k2)) {
  meta <- meta %>% filter(!is.na(k2), nzchar(k2))
}

expr <- npq_wide[meta$Tube_ID, selected_features, drop = FALSE]
rownames(expr) <- meta$Tube_ID
expr <- expr[stats::complete.cases(expr), , drop = FALSE]
meta <- meta[match(rownames(expr), meta$Tube_ID), , drop = FALSE]
var_ok <- apply(expr, 2, function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  !is.na(s) && s > 0
})
expr <- expr[, var_ok, drop = FALSE]
if (ncol(expr) < 2L || nrow(expr) < 3L) {
  stop("Insufficient samples or features after filtering.")
}

# ---- PCA ----
message("PCA: n = ", nrow(expr), ", p = ", ncol(expr))

pca <- stats::prcomp(expr, center = TRUE, scale. = TRUE)
var_expl <- summary(pca)$importance[2, 1:2] * 100

pr_breaks <- as.numeric(trimws(strsplit(opt$pr_breaks, ",", fixed = TRUE)[[1]]))
if (length(pr_breaks) < 3L || !is.infinite(pr_breaks[length(pr_breaks)])) {
  stop("--pr_breaks must be increasing numeric values ending with Inf.")
}
bin_labels <- vapply(seq_len(length(pr_breaks) - 1L), function(i) {
  hi <- pr_breaks[i + 1L]
  if (is.infinite(hi)) paste0("\u2265", pr_breaks[i]) else paste0("[", pr_breaks[i], ", ", hi, ")")
}, character(1L))

scores <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  k2 = factor(meta$k2, levels = c("alpha", "beta")),
  progression_rate = meta$progression_rate,
  progression_rate_bin = cut(
    meta$progression_rate,
    breaks = pr_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = bin_labels
  ),
  stringsAsFactors = FALSE
)
if ("slow_vital_capacity" %in% names(meta)) scores$slow_vital_capacity <- meta$slow_vital_capacity
if ("Nfl" %in% names(meta)) scores$Nfl <- meta$Nfl

n_bins <- nlevels(scores$progression_rate_bin)
pr_palette <- grDevices::colorRampPalette(c("#313695", "#74add1", "#fed976", "#fd8d3c", "#a50026"))(n_bins)
names(pr_palette) <- levels(scores$progression_rate_bin)

scale_k2_shape <- scale_shape_manual(values = c(alpha = 16, beta = 17), name = "k2")

plot_df <- scores %>% filter(!is.na(progression_rate_bin))

p_prog <- ggplot(plot_df, aes(x = .data$PC1, y = .data$PC2, color = .data$progression_rate_bin, shape = .data$k2)) +
  geom_point(size = 3.8, alpha = 0.92, stroke = 0.9) +
  scale_color_manual(values = pr_palette, name = "progression_rate\n(binned)", drop = FALSE) +
  scale_k2_shape +
  labs(
    title = paste0("PCA: ML-selected proteins (", panel_label, ", ", matrix_meta, " NPQ, ", ncol(expr), " proteins)"),
    x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2], 1), "%)")
  ) +
  theme_classic(base_size = 13)
ggsave(file.path(output_dir, "pca_ML_proteins_progression_k2.pdf"), p_prog, width = 8, height = 6)

if ("slow_vital_capacity" %in% names(scores)) {
  plot_df_svc <- scores %>% filter(!is.na(slow_vital_capacity))
  if (nrow(plot_df_svc) >= 3L) {
    p_svc <- ggplot(plot_df_svc, aes(x = .data$PC1, y = .data$PC2, color = .data$slow_vital_capacity, shape = .data$k2)) +
      geom_point(size = 3.6, alpha = 0.9, stroke = 0.9) +
      viridis::scale_color_viridis(option = "C", name = "slow_vital_capacity") +
      scale_k2_shape +
      labs(
        title = paste0("PCA: ML-selected proteins (", panel_label, ", ", matrix_meta, ", ", ncol(expr), " proteins)"),
        x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PC2 (", round(var_expl[2], 1), "%)")
      ) +
      theme_classic(base_size = 13)
    ggsave(file.path(output_dir, "pca_ML_proteins_slow_vital_capacity_k2.pdf"), p_svc, width = 8, height = 6)
  }
}

if ("Nfl" %in% names(scores)) {
  plot_df_nfl <- scores %>% filter(!is.na(Nfl))
  if (nrow(plot_df_nfl) >= 3L) {
    p_nfl <- ggplot(plot_df_nfl, aes(x = .data$PC1, y = .data$PC2, color = .data$Nfl, shape = .data$k2)) +
      geom_point(size = 3.6, alpha = 0.9, stroke = 0.9) +
      viridis::scale_color_viridis(option = "C", name = "Nfl") +
      scale_k2_shape +
      labs(
        title = paste0("PCA: ML-selected proteins (", panel_label, ", ", matrix_meta, ", ", ncol(expr), " proteins)"),
        x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
        y = paste0("PC2 (", round(var_expl[2], 1), "%)")
      ) +
      theme_classic(base_size = 13)
    ggsave(file.path(output_dir, "pca_ML_proteins_Nfl_k2.pdf"), p_nfl, width = 8, height = 6)
  }
}

var_cor <- stats::cor(expr, pca$x[, 1:2, drop = FALSE])
var_r <- sqrt(rowSums(var_cor^2))
var_dir <- var_cor
var_dir[var_r > 0, ] <- var_cor[var_r > 0, , drop = FALSE] / var_r[var_r > 0]
var_dir[var_r == 0, ] <- 0

# Biplot is drawn on the Nfl plot (Nfl colour; samples with missing Nfl dropped).
# Falls back to progression_rate colouring if Nfl is not in the metadata.
use_nfl_biplot <- "Nfl" %in% names(scores) && sum(!is.na(scores$Nfl)) >= 3L
bi_df <- if (use_nfl_biplot) scores %>% filter(!is.na(Nfl)) else plot_df
message("Biplot: coloured by ", if (use_nfl_biplot) "Nfl" else "progression_rate (Nfl not available)",
        ", n = ", nrow(bi_df), " of ", nrow(scores), " samples")

score_span <- max(abs(c(bi_df$PC1, bi_df$PC2)), na.rm = TRUE)
arrow_len <- 0.38 * score_span
gene_arrows <- data.frame(
  feature_label = ifelse(
    grepl("^CNS_", rownames(var_dir)),
    paste0("CNS:", sub("^CNS_", "", rownames(var_dir))),
    paste0("Imm:", sub("^Immune_", "", rownames(var_dir)))
  ),
  PC1 = var_dir[, 1] * arrow_len,
  PC2 = var_dir[, 2] * arrow_len,
  label_x = var_dir[, 1] * arrow_len * 1.06,
  label_y = var_dir[, 2] * arrow_len * 1.06,
  stringsAsFactors = FALSE
)

bi_color <- if (use_nfl_biplot) {
  list(aes(color = .data$Nfl), viridis::scale_color_viridis(option = "C", name = "Nfl"))
} else {
  list(aes(color = .data$progression_rate_bin),
       scale_color_manual(values = pr_palette, name = "progression_rate\n(binned)", drop = FALSE))
}

p_biplot <- ggplot(bi_df, aes(x = .data$PC1, y = .data$PC2, shape = .data$k2)) +
  bi_color[[1]] +
  geom_point(size = 3.2, alpha = 0.85, stroke = 0.85) +
  geom_segment(
    data = gene_arrows,
    aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
    color = "gray15",
    linewidth = 0.65
  ) +
  label_repel_geom(gene_arrows, aes(x = .data$label_x, y = .data$label_y, label = .data$feature_label)) +
  bi_color[[2]] +
  scale_k2_shape +
  labs(
    title = paste0("PCA biplot: samples and protein directions (", panel_label, ", ", matrix_meta, ", n = ", nrow(bi_df), ")"),
    x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2], 1), "%)")
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, "pca_biplot_gene_contribution.pdf"), p_biplot, width = 9.5, height = 6.5)

# ---- PCA on ALL proteins of the panel(s) (same samples as above) ----
# Proteins with any missing value or zero variance in these samples are dropped.
expr_all <- npq_wide[rownames(expr), , drop = FALSE]
expr_all <- expr_all[, colSums(is.na(expr_all)) == 0, drop = FALSE]
sd_all <- apply(expr_all, 2, stats::sd)
expr_all <- expr_all[, !is.na(sd_all) & sd_all > 0, drop = FALSE]
message("PCA (all proteins): n = ", nrow(expr_all), ", p = ", ncol(expr_all),
        " (of ", ncol(npq_wide), " proteins in the ", panel_label, " panel)")

pca_all <- stats::prcomp(expr_all, center = TRUE, scale. = TRUE)
var_expl_all <- summary(pca_all)$importance[2, 1:2] * 100
scores_all <- scores
scores_all$PC1 <- pca_all$x[, 1]
scores_all$PC2 <- pca_all$x[, 2]

all_title <- paste0("PCA: all proteins (", panel_label, ", ", matrix_meta, ", ",
                    ncol(expr_all), " proteins)")
all_axes <- labs(
  x = paste0("PC1 (", round(var_expl_all[1], 1), "%)"),
  y = paste0("PC2 (", round(var_expl_all[2], 1), "%)")
)

# progression_rate (binned)
df_all_prog <- scores_all %>% filter(!is.na(progression_rate_bin))
p_all_prog <- ggplot(df_all_prog, aes(x = .data$PC1, y = .data$PC2,
                                      color = .data$progression_rate_bin, shape = .data$k2)) +
  geom_point(size = 3.8, alpha = 0.92, stroke = 0.9) +
  scale_color_manual(values = pr_palette, name = "progression_rate\n(binned)", drop = FALSE) +
  scale_k2_shape +
  labs(title = all_title, subtitle = paste0("n = ", nrow(df_all_prog))) +
  all_axes +
  theme_classic(base_size = 13)
ggsave(file.path(output_dir, "pca_all_proteins_progression_k2.pdf"), p_all_prog, width = 8, height = 6)

# continuous variables: slow_vital_capacity, Nfl
for (v in c("slow_vital_capacity", "Nfl")) {
  if (!v %in% names(scores_all)) next
  df_v <- scores_all[!is.na(scores_all[[v]]), , drop = FALSE]
  if (nrow(df_v) < 3L) next
  p_v <- ggplot(df_v, aes(x = .data$PC1, y = .data$PC2, color = .data[[v]], shape = .data$k2)) +
    geom_point(size = 3.6, alpha = 0.9, stroke = 0.9) +
    viridis::scale_color_viridis(option = "C", name = v) +
    scale_k2_shape +
    labs(title = all_title, subtitle = paste0("n = ", nrow(df_v))) +
    all_axes +
    theme_classic(base_size = 13)
  ggsave(file.path(output_dir, paste0("pca_all_proteins_", v, "_k2.pdf")), p_v, width = 8, height = 6)
}

# ---- Biplot on the all-protein PCA (Nfl plot): directions of the ML-selected proteins ----
# Arrow = correlation of each ML-selected protein with all-protein PC1/PC2.
# Arrows are NOT normalised to equal length: a short arrow means the protein is only
# weakly aligned with the PC1/PC2 plane of the full panel.
ml_in_all <- intersect(colnames(expr), colnames(expr_all))
if ("Nfl" %in% names(scores_all) && length(ml_in_all) >= 1L) {
  bi_all_df <- scores_all %>% filter(!is.na(Nfl))
  if (nrow(bi_all_df) >= 3L) {
    cor_all <- stats::cor(expr_all[, ml_in_all, drop = FALSE], pca_all$x[, 1:2, drop = FALSE])
    len_all <- 0.9 * max(abs(c(bi_all_df$PC1, bi_all_df$PC2)), na.rm = TRUE)
    arrows_all <- data.frame(
      feature_label = ifelse(
        grepl("^CNS_", ml_in_all),
        paste0("CNS:", sub("^CNS_", "", ml_in_all)),
        paste0("Imm:", sub("^Immune_", "", ml_in_all))
      ),
      PC1 = cor_all[, 1] * len_all,
      PC2 = cor_all[, 2] * len_all,
      label_x = cor_all[, 1] * len_all * 1.08,
      label_y = cor_all[, 2] * len_all * 1.08,
      stringsAsFactors = FALSE
    )
    write.csv(
      data.frame(Feature = ml_in_all, cor_PC1 = cor_all[, 1], cor_PC2 = cor_all[, 2],
                 length = sqrt(rowSums(cor_all^2))),
      file.path(output_dir, "pca_all_proteins_ML_protein_correlations.csv"), row.names = FALSE
    )

    p_all_bi <- ggplot(bi_all_df, aes(x = .data$PC1, y = .data$PC2, color = .data$Nfl, shape = .data$k2)) +
      geom_point(size = 3.6, alpha = 0.9, stroke = 0.9) +
      geom_segment(
        data = arrows_all,
        aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
        inherit.aes = FALSE,
        arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
        color = "gray15",
        linewidth = 0.65
      ) +
      label_repel_geom(arrows_all, aes(x = .data$label_x, y = .data$label_y, label = .data$feature_label)) +
      viridis::scale_color_viridis(option = "C", name = "Nfl") +
      scale_k2_shape +
      labs(
        title = all_title,
        subtitle = paste0("Arrows: ML-selected proteins; direction and length = correlation with PC1/PC2. n = ",
                          nrow(bi_all_df))
      ) +
      all_axes +
      theme_classic(base_size = 13)
    ggsave(file.path(output_dir, "pca_all_proteins_biplot_ML_proteins_Nfl_k2.pdf"), p_all_bi,
           width = 9.5, height = 6.5)
  }
}

message("Saved:")
message("  pca_biplot_gene_contribution.pdf")
message("  pca_ML_proteins_{progression,slow_vital_capacity,Nfl}_k2.pdf")
message("  pca_all_proteins_{progression,slow_vital_capacity,Nfl}_k2.pdf")
message("  pca_all_proteins_biplot_ML_proteins_Nfl_k2.pdf, pca_all_proteins_ML_protein_correlations.csv")
