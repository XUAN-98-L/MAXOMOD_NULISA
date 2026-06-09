# PCA + PLS on ML-selected NULISA proteins; MS top5 cross-platform PLS
#
# Rscript Script/09_PCA_ML_Genes.R \
#   -i CNS_immune/Results \
#   -f .../08_Clustering_Feature_Selection/selected_features.txt \
#   -o .../09_PCA_ML_Genes

suppressMessages(library(optparse))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Note: install 'ggrepel' for non-overlapping biplot labels.")
}

k2_colors <- setNames(brewer.pal(8, "Set2")[c(2, 3)], c("alpha", "beta"))

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
  make_option(c("-i", "--input"), type = "character", default = "CNS_immune/Results"),
  make_option(
    c("-f", "--features"),
    type = "character",
    default = "CNS_immune/Results/Combined_CNS_Immune_panels/Without_tears/08_Clustering_Feature_Selection/selected_features.txt"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "CNS_immune/Results/Combined_CNS_Immune_panels/Without_tears/09_PCA_ML_Genes"
  ),
  make_option(c("--require_k2"), type = "logical", default = TRUE),
  make_option(c("--pr_breaks"), type = "character", default = "0,0.5,1,2,4,Inf"),
  make_option(
    c("--MS_expression"),
    type = "character",
    default = "/Users/xliu2942/Documents/Projects/MAXOMOD/14324_0_source_data_169635_t9mdsb/Supp table 12 - intensity_imputed_log2transf_norm_DC.csv"
  ),
  make_option(c("--ms_protein_col"), type = "character", default = "Protein"),
  make_option(c("--ms_top5"), type = "character", default = "PARK7,PTPRS,ATRN,CNTN1,PCSK1N"),
  make_option(c("--pls_ncomp"), type = "integer", default = 5L)
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_panel_paths <- function(input_arg) {
  if (grepl(",", input_arg, fixed = TRUE)) {
    parts <- trimws(strsplit(input_arg, ",", fixed = TRUE)[[1]])
    if (length(parts) != 3L) {
      stop("--input must be a panel base directory or three comma-separated file paths.")
    }
    list(cns = parts[1], immune = parts[2], metadata = parts[3])
  } else {
    subdir <- "Without_tears"
    list(
      cns = file.path(input_arg, "CNS_panel", subdir, "00_Initialization", "npq_counts.xlsx"),
      immune = file.path(input_arg, "Immune_panel", subdir, "00_Initialization", "npq_counts.xlsx"),
      metadata = file.path(input_arg, "CNS_panel", subdir, "00_Initialization", "all_participants_IDs.xlsx")
    )
  }
}

read_ms_matrix <- function(path, protein_col = "Protein") {
  if (!file.exists(path)) {
    stop("MS expression file not found: ", path)
  }
  ms_df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!protein_col %in% colnames(ms_df)) {
    stop("MS file must contain column '", protein_col, "'.")
  }
  prot <- as.character(ms_df[[protein_col]])
  mat <- ms_df[, setdiff(colnames(ms_df), protein_col), drop = FALSE]
  mat <- as.data.frame(lapply(mat, function(x) suppressWarnings(as.numeric(x))))
  rownames(mat) <- prot
  mat
}

npq_to_wide <- function(npq_path, panel_prefix) {
  read_excel(npq_path) %>%
    filter(SampleMatrixType == "CSF", SampleType == "Sample") %>%
    mutate(Feature = paste0(panel_prefix, "_", Target)) %>%
    select(SampleName, Feature, NPQ) %>%
    pivot_wider(names_from = Feature, values_from = NPQ) %>%
    as.data.frame()
}

panel_paths <- parse_panel_paths(opt$input)
output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

selected_features <- trimws(readLines(opt$features, warn = FALSE))
selected_features <- selected_features[nzchar(selected_features)]
if (length(selected_features) < 2L) {
  stop("Need at least 2 features for PCA.")
}

ms_top5 <- trimws(strsplit(opt$ms_top5, ",", fixed = TRUE)[[1]])
ms_top5 <- ms_top5[nzchar(ms_top5)]
if (length(ms_top5) < 2L) {
  stop("--ms_top5 must contain at least 2 comma-separated protein names.")
}

message("Loading CNS + Immune NPQ (CSF samples)...")
npq_wide_cns <- npq_to_wide(panel_paths$cns, "CNS")
npq_wide_immune <- npq_to_wide(panel_paths$immune, "Immune")
common_samples <- intersect(npq_wide_cns$SampleName, npq_wide_immune$SampleName)

npq_wide <- npq_wide_cns %>%
  filter(.data$SampleName %in% common_samples) %>%
  inner_join(
    npq_wide_immune %>% filter(.data$SampleName %in% common_samples),
    by = "SampleName"
  )
rownames(npq_wide) <- npq_wide$SampleName
npq_wide$SampleName <- NULL

missing_feat <- setdiff(selected_features, colnames(npq_wide))
if (length(missing_feat) > 0L) {
  stop("Selected features not in combined NPQ matrix: ", paste(missing_feat, collapse = ", "))
}

meta <- read_excel(panel_paths$metadata) %>%
  filter(Material == "CSF") %>%
  mutate(
    progression_rate = suppressWarnings(as.numeric(progression_rate)),
    k2 = as.character(k2),
    across(any_of(c("slow_vital_capacity", "Nfl")), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  select(Tube_ID, CSF_ID, k2, progression_rate, any_of(c("slow_vital_capacity", "Nfl"))) %>%
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

# ---- PLS: MS top5 vs NULISA (alignment plot + loadings only) ----
message("Running PLS (MS top5 vs NULISA ML features)...")
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  stop("Package 'mixOmics' is required for PLS.")
}
suppressMessages(library(mixOmics))

ms_mat <- read_ms_matrix(opt$MS_expression, protein_col = opt$ms_protein_col)
ms_present <- intersect(ms_top5, rownames(ms_mat))
nulisa_present <- intersect(selected_features, colnames(expr))
if (length(ms_present) < 2L || length(nulisa_present) < 2L) {
  stop("PLS needs >=2 MS and >=2 NULISA proteins present in data.")
}

csf_ids_vec <- as.character(meta$CSF_ID)
keep <- which(!is.na(csf_ids_vec) & csf_ids_vec %in% colnames(ms_mat))
X_pls <- t(ms_mat[ms_present, csf_ids_vec[keep], drop = FALSE])
Y_pls <- as.matrix(expr[keep, nulisa_present, drop = FALSE])
ok <- complete.cases(cbind(as.data.frame(X_pls), as.data.frame(Y_pls)))
X_pls <- X_pls[ok, , drop = FALSE]
Y_pls <- Y_pls[ok, , drop = FALSE]
meta_pls <- meta[keep, , drop = FALSE][ok, , drop = FALSE]
if (nrow(X_pls) < 3L) {
  stop("Fewer than 3 complete paired samples for PLS.")
}

ncomp_use <- min(opt$pls_ncomp, nrow(X_pls) - 1L, ncol(X_pls), ncol(Y_pls))
res_pls <- mixOmics::spls(X = X_pls, Y = Y_pls, ncomp = ncomp_use)

load_X <- as.data.frame(res_pls$loadings$X)
load_Y <- as.data.frame(res_pls$loadings$Y)
load_X$protein <- rownames(load_X)
load_Y$protein <- rownames(load_Y)
write.csv(load_X, file.path(output_dir, "pls_loadings_X_MS_top5.csv"), row.names = FALSE)
write.csv(load_Y, file.path(output_dir, "pls_loadings_Y_NULISA_top6.csv"), row.names = FALSE)

scoresX <- as.data.frame(res_pls$variates$X)
scoresY <- as.data.frame(res_pls$variates$Y)
pls_scores <- data.frame(
  k2 = factor(meta_pls$k2, levels = c("alpha", "beta")),
  X_comp1 = scoresX[, 1],
  Y_comp1 = scoresY[, 1],
  stringsAsFactors = FALSE
)
ct <- suppressWarnings(stats::cor.test(pls_scores$X_comp1, pls_scores$Y_comp1, method = "pearson"))

p_pls <- ggplot(pls_scores, aes(x = .data$X_comp1, y = .data$Y_comp1, color = .data$k2)) +
  geom_point(size = 3.5, alpha = 0.9, stroke = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "gray20") +
  scale_color_manual(values = k2_colors, name = "k2") +
  theme_classic(base_size = 13) +
  labs(
    title = "PLS component alignment (comp1): MS top5 vs NULISA",
    subtitle = paste0(
      "X_comp1 vs Y_comp1; r = ", sprintf("%.3f", ct$estimate),
      ", p = ", signif(ct$p.value, 3), ", n = ", nrow(pls_scores)
    ),
    x = "X_comp1",
    y = "Y_comp1"
  )
ggsave(file.path(output_dir, "pls_component_alignment_comp1.pdf"), p_pls, width = 7.2, height = 6)
ggsave(file.path(output_dir, "pls_component_alignment_comp1.png"), p_pls, width = 7.2, height = 6, dpi = 300)

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
    title = "PCA: ML-selected proteins (CNS + Immune CSF NPQ)",
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
        title = "PCA: ML-selected proteins (CNS+Immune)",
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
        title = "PCA: ML-selected proteins (CNS+Immune)",
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

score_span <- max(abs(c(plot_df$PC1, plot_df$PC2)), na.rm = TRUE)
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

p_biplot <- ggplot(plot_df, aes(x = .data$PC1, y = .data$PC2, color = .data$progression_rate_bin, shape = .data$k2)) +
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
  scale_color_manual(values = pr_palette, name = "progression_rate\n(binned)", drop = FALSE) +
  scale_k2_shape +
  labs(
    title = "PCA biplot: samples and protein directions",
    x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2], 1), "%)")
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, "pca_biplot_gene_contribution.pdf"), p_biplot, width = 9.5, height = 6.5)

message("Saved:")
message("  pls_component_alignment_comp1.pdf/png")
message("  pls_loadings_X_MS_top5.csv, pls_loadings_Y_NULISA_top6.csv")
message("  pca_biplot_gene_contribution.pdf")
message("  pca_ML_proteins_{progression,slow_vital_capacity,Nfl}_k2.pdf")
