# 12_NFL.R
#
# Correlate the clinical NfL measurement (metadata column "Nfl", pg/ml) with the NULISA
# NEFL NPQ in the same patients (CSF by default).
#
# NPQ is already on a log2 scale, so the clinical NfL is log-transformed for Pearson
# correlation / regression; Spearman (rank-based) is reported as the primary statistic.
# Only samples with a non-missing clinical NfL are used (in MAXOMOD: the 43 ALS patients;
# controls have no clinical NfL). All samples are analysed together (no subgroups).
# Outliers = samples whose studentized residual from lm(NPQ ~ log2 NfL) exceeds --outlier_z
# in absolute value; they are labelled in the plot and statistics are also given without them.
#
# Example (run from NULISA_MAXOMOD/):
# Rscript Script/12_NFL.R -i CNS_immune/Results
# Rscript Script/12_NFL.R -i CNS_immune/Results --matrix Serum
#
# The script also redraws the two PCA biplots of 09_PCA_ML_Genes.R (ML-selected proteins and
# all proteins; same samples, same PCA, Nfl colour, alpha/beta shape) and marks the NfL outlier(s),
# so they can be located on the slide "The six proteins separate the subtypes ...".
#   -> CNS_immune/Results/CNS_panel/Without_tears/12_NFL[_<Matrix>]/

suppressMessages(library(optparse))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "CNS_immune/Results",
              help = "Panel base directory (contains CNS_panel/Without_tears/00_Initialization) [default: %default]"),
  make_option(c("--matrix"), type = "character", default = "CSF",
              help = "Sample matrix for the NULISA NEFL value: CSF, Serum or Plasma [default: %default]"),
  make_option(c("--target"), type = "character", default = "NEFL",
              help = "NULISA target to correlate [default: %default]"),
  make_option(c("--clin_col"), type = "character", default = "Nfl",
              help = "Metadata column with the clinical measurement [default: %default]"),
  make_option(c("--outlier_z"), type = "numeric", default = 3,
              help = "|Studentized residual| above which a sample is flagged as an outlier [default: %default]"),
  make_option(c("--features"), type = "character", default = NULL,
              help = paste("selected_features.txt from 08, for the PCA [default:",
                           "<input>/CNS_panel/Without_tears/08_Clustering_Feature_Selection[_<Matrix>]/selected_features.txt]")),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory [default: <input>/CNS_panel/Without_tears/12_NFL[_<Matrix>]]")
)
opt <- parse_args(OptionParser(option_list = option_list))

matrix_key <- toupper(trimws(opt$matrix))
if (!matrix_key %in% c("CSF", "SERUM", "PLASMA")) stop("--matrix must be CSF, Serum or Plasma.")
matrix_label <- switch(matrix_key, CSF = "CSF", SERUM = "Serum", PLASMA = "Plasma")

init_dir <- file.path(opt$input, "CNS_panel", "Without_tears", "00_Initialization")
npq_path <- file.path(init_dir, "npq_counts.xlsx")
meta_path <- file.path(init_dir, "all_participants_IDs.xlsx")
for (f in c(npq_path, meta_path)) if (!file.exists(f)) stop("Not found: ", f)

output_dir <- opt$output
if (is.null(output_dir) || !nzchar(output_dir)) {
  output_dir <- file.path(opt$input, "CNS_panel", "Without_tears",
                          paste0("12_NFL", if (matrix_label == "CSF") "" else paste0("_", matrix_label)))
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("NULISA NPQ: ", npq_path, " (", matrix_label, ", target ", opt$target, ")")
message("Metadata:   ", meta_path, " (column ", opt$clin_col, ")")
message("Output:     ", output_dir)

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------
npq <- read_excel(npq_path) %>%
  filter(toupper(SampleMatrixType) == matrix_key, SampleType == "Sample",
         toupper(Target) == toupper(opt$target)) %>%
  transmute(SampleName = as.character(SampleName), NPQ = as.numeric(NPQ), PlateID, SampleQC)
if (nrow(npq) == 0L) stop("No ", opt$target, " values found for ", matrix_label, ".")

meta <- read_excel(meta_path) %>%
  filter(toupper(Material) == matrix_key) %>%
  mutate(Tube_ID = as.character(Tube_ID))
if (!opt$clin_col %in% names(meta)) stop("Column '", opt$clin_col, "' not in metadata.")

dat <- meta %>%
  transmute(
    Tube_ID, Patient = if ("Patient" %in% names(meta)) Patient else NA,
    type = if ("type" %in% names(meta)) type else NA,
    k2 = if ("k2" %in% names(meta)) k2 else NA,
    sex = if ("sex" %in% names(meta)) sex else NA,
    onset = if ("onset" %in% names(meta)) onset else NA,
    age = if ("age" %in% names(meta)) suppressWarnings(as.numeric(age)) else NA,
    clin = suppressWarnings(as.numeric(as.character(.data[[opt$clin_col]])))
  ) %>%
  inner_join(npq, by = c("Tube_ID" = "SampleName")) %>%
  mutate(group = case_when(
    !is.na(k2) & k2 %in% c("alpha", "beta") ~ as.character(k2),
    !is.na(type) ~ as.character(type),
    TRUE ~ "unknown"
  ))

n_all <- nrow(dat)
dat <- dat %>% filter(!is.na(clin), !is.na(NPQ), clin > 0)
message("Samples with ", opt$target, " NPQ: ", n_all, "; with a clinical ", opt$clin_col, " value: ", nrow(dat))
if (nrow(dat) < 4L) stop("Fewer than 4 samples with both values.")
dat$log10_clin <- log10(dat$clin)
dat$log2_clin <- log2(dat$clin)

write.csv(dat, file.path(output_dir, paste0("NfL_vs_", opt$target, "_", matrix_label, "_data.csv")),
          row.names = FALSE)

# ------------------------------------------------------------
# Outliers: studentized residuals from the log-log regression
# ------------------------------------------------------------
fit_all <- lm(NPQ ~ log2_clin, data = dat)
dat$stud_resid <- rstudent(fit_all)
dat$cooks_d <- cooks.distance(fit_all)
dat$outlier <- abs(dat$stud_resid) > opt$outlier_z
dat$label <- ifelse(is.na(dat$Patient), dat$Tube_ID, as.character(dat$Patient))
# Outlier annotation: ID, sex, onset site, subtype (alpha/beta)
na_dash <- function(x) ifelse(is.na(x) | !nzchar(as.character(x)), "NA", as.character(x))
dat$label_full <- paste0(dat$label, "\n", tolower(na_dash(dat$sex)), ", ",
                         na_dash(dat$onset), " onset, ", na_dash(dat$k2))
n_out <- sum(dat$outlier)
message("Outliers (|studentized residual| > ", opt$outlier_z, "): ", n_out,
        if (n_out > 0) paste0(" -> ", paste(sub("\n", " (", dat$label_full[dat$outlier]), collapse = "), "), ")") else "")
write.csv(dat, file.path(output_dir, paste0("NfL_vs_", opt$target, "_", matrix_label, "_data.csv")),
          row.names = FALSE)
if (n_out > 0) {
  write.csv(dat[dat$outlier, c("Tube_ID", "Patient", "sex", "onset", "k2", "clin", "NPQ", "stud_resid", "cooks_d", "SampleQC")],
            file.path(output_dir, paste0("NfL_vs_", opt$target, "_", matrix_label, "_outliers.csv")),
            row.names = FALSE)
}

# ------------------------------------------------------------
# Correlation statistics: all samples, and without the flagged outliers
# ------------------------------------------------------------
cor_stats <- function(d, label) {
  if (nrow(d) < 4L) return(NULL)
  sp <- suppressWarnings(cor.test(d$clin, d$NPQ, method = "spearman", exact = FALSE))
  pe <- cor.test(d$log2_clin, d$NPQ, method = "pearson")
  fit <- lm(NPQ ~ log2_clin, data = d)
  data.frame(
    samples = label, n = nrow(d),
    spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
    pearson_r_log = unname(pe$estimate), pearson_p = pe$p.value,
    pearson_ci_low = pe$conf.int[1], pearson_ci_high = pe$conf.int[2],
    slope_NPQ_per_log2_NfL = unname(coef(fit)[2]), r_squared = summary(fit)$r.squared
  )
}
stats <- bind_rows(
  cor_stats(dat, "all"),
  if (n_out > 0) cor_stats(dat[!dat$outlier, ], "without outliers")
)
write.csv(stats, file.path(output_dir, paste0("NfL_vs_", opt$target, "_", matrix_label, "_correlation.csv")),
          row.names = FALSE)
cat("\nCorrelation of clinical", opt$clin_col, "with NULISA", opt$target, "NPQ (", matrix_label, "):\n")
print(stats %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)

# ------------------------------------------------------------
# Plot: all samples together, outliers labelled
# ------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
lab_stats <- function(s) sprintf("%s: Spearman rho = %.2f (p = %s), Pearson r (log NfL) = %.2f, n = %d",
                                 s$samples, s$spearman_rho, fmt_p(s$spearman_p), s$pearson_r_log, s$n)
subtitle <- paste(vapply(seq_len(nrow(stats)), function(k) lab_stats(stats[k, ]), character(1)),
                  collapse = "\n")
dat$point_type <- factor(ifelse(dat$outlier, "outlier", "sample"), levels = c("sample", "outlier"))

label_geom <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(data = dat[dat$outlier, ], aes(label = .data$label_full), size = 3.3,
                           color = "#B2182B", fontface = "bold", min.segment.length = 0,
                           box.padding = 0.6, seed = 1)
} else {
  geom_text(data = dat[dat$outlier, ], aes(label = .data$label_full), size = 3.3, color = "#B2182B",
            fontface = "bold", vjust = -1)
}

p_all <- ggplot(dat, aes(x = .data$clin, y = .data$NPQ)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "gray25", fill = "gray80", linewidth = 0.8) +
  geom_point(aes(color = .data$point_type, shape = .data$point_type), size = 3, alpha = 0.9) +
  label_geom +
  scale_x_log10() +
  scale_color_manual(values = c(sample = "#4575B4", outlier = "#B2182B"), drop = FALSE,
                     labels = c(sample = "sample", outlier = paste0("outlier (|stud. resid.| > ", opt$outlier_z, ")")),
                     name = NULL) +
  scale_shape_manual(values = c(sample = 16, outlier = 17), drop = FALSE,
                     labels = c(sample = "sample", outlier = paste0("outlier (|stud. resid.| > ", opt$outlier_z, ")")),
                     name = NULL) +
  labs(#title = paste0("Clinical NfL vs NULISA ", opt$target, " (", matrix_label, ")"),
       subtitle = subtitle,
       x = paste0("Clinical ", opt$clin_col, " (pg/ml, log scale)"),
       y = paste0("NULISA ", opt$target, " NPQ (", matrix_label, ")")) +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(size = 9), legend.position = "bottom")
ggsave(file.path(output_dir, paste0("NfL_vs_", opt$target, "_", matrix_label, ".pdf")), p_all, width = 7, height = 6)
ggsave(file.path(output_dir, paste0("NfL_vs_", opt$target, "_", matrix_label, ".png")), p_all, width = 7, height = 6, dpi = 300)

message("\nSaved in ", normalizePath(output_dir), ":")
message("  NfL_vs_", opt$target, "_", matrix_label, ".pdf/png          (all samples, outliers labelled)")
message("  NfL_vs_", opt$target, "_", matrix_label, "_correlation.csv  (all samples / without outliers)")
message("  NfL_vs_", opt$target, "_", matrix_label, "_outliers.csv     (flagged samples, if any)")
message("  NfL_vs_", opt$target, "_", matrix_label, "_data.csv         (values used, residuals, outlier flag)")

# ============================================================
# PCA as in 09_PCA_ML_Genes.R (--panel cns), with the NfL outlier(s) marked
#   same samples (CSF/matrix samples with k2), same PCA (prcomp, centred + scaled, raw NPQ),
#   Nfl colour, k2 shape, same arrows (ML proteins: equal length; all proteins: length = r)
# ============================================================
features_path <- opt$features
if (is.null(features_path) || !nzchar(features_path)) {
  features_path <- file.path(opt$input, "CNS_panel", "Without_tears",
                             paste0("08_Clustering_Feature_Selection",
                                    if (matrix_label == "CSF") "" else paste0("_", matrix_label)),
                             "selected_features.txt")
}
if (!file.exists(features_path)) {
  message("Note: ", features_path, " not found; PCA with outliers skipped.")
} else {
  selected_features <- trimws(readLines(features_path, warn = FALSE))
  selected_features <- selected_features[nzchar(selected_features)]

  # All CNS-panel proteins for this matrix, prefixed as in 08/09
  npq_wide <- read_excel(npq_path) %>%
    filter(toupper(SampleMatrixType) == matrix_key, SampleType == "Sample") %>%
    mutate(Feature = paste0("CNS_", Target), SampleName = as.character(SampleName)) %>%
    select(SampleName, Feature, NPQ) %>%
    tidyr::pivot_wider(names_from = Feature, values_from = NPQ) %>%
    as.data.frame()
  rownames(npq_wide) <- npq_wide$SampleName
  npq_wide$SampleName <- NULL

  pmeta <- meta %>%
    transmute(Tube_ID, Patient = if ("Patient" %in% names(meta)) as.character(Patient) else Tube_ID,
              k2 = as.character(k2),
              Nfl = suppressWarnings(as.numeric(as.character(.data[[opt$clin_col]])))) %>%
    filter(Tube_ID %in% rownames(npq_wide), !is.na(k2), nzchar(k2))

  expr <- npq_wide[pmeta$Tube_ID, intersect(selected_features, colnames(npq_wide)), drop = FALSE]
  expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  pmeta <- pmeta[match(rownames(expr), pmeta$Tube_ID), , drop = FALSE]
  expr <- expr[, apply(expr, 2, function(x) { s <- stats::sd(x); !is.na(s) && s > 0 }), drop = FALSE]

  expr_all <- npq_wide[rownames(expr), , drop = FALSE]
  expr_all <- expr_all[, colSums(is.na(expr_all)) == 0, drop = FALSE]
  sd_all <- apply(expr_all, 2, stats::sd)
  expr_all <- expr_all[, !is.na(sd_all) & sd_all > 0, drop = FALSE]

  out_ids <- dat$Tube_ID[dat$outlier]
  pmeta$outlier <- pmeta$Tube_ID %in% out_ids
  pmeta$out_label <- dat$label_full[match(pmeta$Tube_ID, dat$Tube_ID)]

  feat_lab <- function(f) ifelse(grepl("^CNS_", f), paste0("CNS:", sub("^CNS_", "", f)),
                                 paste0("Imm:", sub("^Immune_", "", f)))
  repel_text <- function(data, mapping, ...) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = data, mapping = mapping, inherit.aes = FALSE, seed = 1,
                               min.segment.length = 0, max.overlaps = Inf, ...)
    } else {
      geom_text(data = data, mapping = mapping, inherit.aes = FALSE, ...)
    }
  }

  pca_outlier_plot <- function(X, arrow_mode, title, file_stem) {
    pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
    ve <- summary(pca)$importance[2, 1:2] * 100
    sc <- data.frame(Tube_ID = rownames(X), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE) %>%
      left_join(pmeta, by = "Tube_ID")
    sc$k2 <- factor(sc$k2, levels = c("alpha", "beta"))
    bi <- sc[!is.na(sc$Nfl), , drop = FALSE]

    r <- stats::cor(X, pca$x[, 1:2, drop = FALSE])
    span <- max(abs(c(bi$PC1, bi$PC2)), na.rm = TRUE)
    if (arrow_mode == "normalised") {
      af <- colnames(X)
      len <- sqrt(rowSums(r^2)); len[len == 0] <- 1
      d <- r / len; L <- 0.38 * span; lab_mult <- 1.06
    } else {
      af <- intersect(selected_features, colnames(X))
      d <- r[af, , drop = FALSE]; L <- 0.9 * span; lab_mult <- 1.08
    }
    d <- d[af, , drop = FALSE]
    arr <- data.frame(lab = feat_lab(af), PC1 = d[, 1] * L, PC2 = d[, 2] * L,
                      lx = d[, 1] * L * lab_mult, ly = d[, 2] * L * lab_mult)

    p <- ggplot(bi, aes(x = .data$PC1, y = .data$PC2)) +
      geom_point(aes(color = .data$Nfl, shape = .data$k2), size = 3.2, alpha = 0.85, stroke = 0.85) +
      geom_point(data = bi[bi$outlier, ], shape = 21, size = 6.5, stroke = 1.3, color = "#B2182B",
                 fill = NA) +
      geom_segment(data = arr, aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2), inherit.aes = FALSE,
                   arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                   color = "gray15", linewidth = 0.65) +
      repel_text(arr, aes(x = .data$lx, y = .data$ly, label = .data$lab), size = 3.2,
                 fontface = "italic", color = "gray10", segment.color = "gray65", box.padding = 0.45) +
      repel_text(bi[bi$outlier, ], aes(x = .data$PC1, y = .data$PC2, label = .data$out_label),
                 size = 3.3, fontface = "bold", color = "#B2182B", box.padding = 1.2,
                 point.padding = 1.2, segment.color = "#B2182B") +
      scale_color_viridis_c(option = "C", name = "Nfl") +
      scale_shape_manual(values = c(alpha = 16, beta = 17), name = "k2") +
      labs(title = title,
           subtitle = paste0("Same PCA as 09_PCA_ML_Genes.R; red ring = NfL outlier (|stud. resid.| > ",
                             opt$outlier_z, "); n = ", nrow(bi)),
           x = paste0("PC1 (", round(ve[1], 1), "%)"), y = paste0("PC2 (", round(ve[2], 1), "%)")) +
      theme_classic(base_size = 12)
    ggsave(file.path(output_dir, paste0(file_stem, ".pdf")), p, width = 9.5, height = 6.5)
    ggsave(file.path(output_dir, paste0(file_stem, ".png")), p, width = 9.5, height = 6.5, dpi = 300)
    sc$PCA <- file_stem
    sc
  }

  sc_ml <- pca_outlier_plot(expr, "normalised",
                            paste0("PCA biplot: ML-selected proteins (", matrix_label, ", ", ncol(expr), " proteins)"),
                            "pca_ML_proteins_biplot_NfL_outlier")
  sc_all <- pca_outlier_plot(expr_all, "correlation",
                             paste0("PCA: all proteins (", matrix_label, ", ", ncol(expr_all), " proteins)"),
                             "pca_all_proteins_biplot_NfL_outlier")
  pos <- bind_rows(sc_ml, sc_all) %>% filter(.data$outlier) %>%
    select(PCA, Tube_ID, Patient, k2, Nfl, PC1, PC2)
  if (nrow(pos) > 0) {
    write.csv(pos, file.path(output_dir, "NfL_outlier_PCA_position.csv"), row.names = FALSE)
    cat("\nNfL outlier position in the 09 PCAs:\n"); print(pos, row.names = FALSE)
  }
  message("  pca_ML_proteins_biplot_NfL_outlier.pdf/png, pca_all_proteins_biplot_NfL_outlier.pdf/png")
}
