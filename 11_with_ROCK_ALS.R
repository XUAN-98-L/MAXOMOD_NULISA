# 11_with_ROCK_ALS.R
#
# Apply the alpha/beta models trained in 08_Clustering_Feature_Selection.R (MAXOMOD)
# to the CSF samples of the ROCK-ALS trial and predict their k2 group.
#
# Inputs
#   --rock_samples : ROCK_ALS/Results/01_Samples_summary  (samples_all_biofluids.csv: sample info)
#   --rock_expr    : ROCK_ALS CSF NULISA data, long format (Target, NPQ, PatNr, ...)
#                    [default: <rock_samples>/../csf/00_Initialization/Expression.csv]
#   --model_dir    : 08_Clustering_Feature_Selection output folder
#                    (08_final_models.rds + 08_model_results.rds)
#
# Harmonisation (--harmonize), applied per protein before prediction:
#   none        : raw ROCK-ALS NPQ
#   center      : shift ROCK-ALS so its mean equals the MAXOMOD training mean
#   standardize : match ROCK-ALS mean and SD to the MAXOMOD training mean and SD [default]
#   Both "center" and "standardize" assume the two cohorts have a broadly similar mix of
#   alpha/beta; "none" assumes no batch/cohort shift between the NULISA runs.
#   Predictions for ALL three modes are written so their agreement can be checked.
#
# Example (run from NULISA_MAXOMOD/):
# Rscript Script/11_with_ROCK_ALS.R \
#   --rock_samples /Users/xliu2942/Documents/Projects/ROCK_ALS/Results/01_Samples_summary \
#   --model_dir CNS_immune/Results/CNS_panel/Without_tears/08_Clustering_Feature_Selection \
#   -o CNS_immune/Results/CNS_panel/Without_tears/11_with_ROCK_ALS

suppressMessages(library(optparse))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(caret))
# Packages needed to predict from the saved caret models
for (pkg in c("glmnet", "randomForest", "kernlab", "xgboost")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Note: package '", pkg, "' not installed; models that need it will be skipped.")
  } else {
    suppressMessages(library(pkg, character.only = TRUE))
  }
}

option_list <- list(
  make_option(
    c("--rock_samples"), type = "character",
    default = "/Users/xliu2942/Documents/Projects/ROCK_ALS/Results/01_Samples_summary",
    help = "ROCK-ALS 01_Samples_summary folder [default: %default]"
  ),
  make_option(
    c("--rock_expr"), type = "character", default = NULL,
    help = "ROCK-ALS CSF Expression.csv [default: <rock_samples>/../csf/00_Initialization/Expression.csv]"
  ),
  make_option(
    c("--model_dir"), type = "character",
    default = "CNS_immune/Results/CNS_panel/Without_tears/08_Clustering_Feature_Selection",
    help = "Output folder of 08_Clustering_Feature_Selection.R [default: %default]"
  ),
  make_option(
    c("-o", "--output"), type = "character",
    default = "CNS_immune/Results/CNS_panel/Without_tears/11_with_ROCK_ALS",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("--harmonize"), type = "character", default = "standardize",
    help = "Main harmonisation mode: none, center or standardize [default: %default]"
  ),
  make_option(
    c("--threshold"), type = "numeric", default = 0.5,
    help = "P(beta) >= threshold -> beta [default: %default]"
  ),
  make_option(
    c("--ensemble_exclude"), type = "character", default = "ridge",
    help = paste(
      "Comma-separated models left out of the ensemble (still reported individually).",
      "Default 'ridge': its tuned lambda hit the grid maximum in 08, so its P(beta) is squeezed",
      "around the training prevalence. Use 'none' to keep all [default: %default]"
    )
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

harm_modes <- c("none", "center", "standardize")
main_mode <- tolower(opt$harmonize)
if (!main_mode %in% harm_modes) stop("--harmonize must be none, center or standardize.")
thr <- opt$threshold

output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 1. Load MAXOMOD models (from 08)
# ============================================================
models_rds  <- file.path(opt$model_dir, "08_final_models.rds")
results_rds <- file.path(opt$model_dir, "08_model_results.rds")
for (f in c(models_rds, results_rds)) {
  if (!file.exists(f)) stop("Not found: ", f, "\nRun 08_Clustering_Feature_Selection.R (--vis FALSE) first.")
}
final_models <- readRDS(models_rds)
res08 <- readRDS(results_rds)
st <- res08$settings
features <- res08$selected_final

message("Models from: ", opt$model_dir)
message("  08 settings: panel = ", st$panel_choice, ", matrix = ", st$matrix,
        ", splits = ", st$BS_number, ", cutoff = ", st$feature_freq_cutoff)
message("  Models: ", paste(names(final_models), collapse = ", "))
message("  Features: ", paste(features, collapse = ", "))
if (!is.null(st$matrix) && toupper(st$matrix) != "CSF") {
  warning("The 08 models were trained on ", st$matrix, " but ROCK-ALS samples here are CSF.")
}

# MAXOMOD training data (raw NPQ) and labels, taken from the fitted caret models
train_data <- final_models[[1]]$trainingData
train_y <- factor(train_data$.outcome, levels = c("alpha", "beta"))
train_X <- as.data.frame(train_data[, features, drop = FALSE])
ref_mean <- vapply(train_X, mean, numeric(1))
ref_sd   <- vapply(train_X, stats::sd, numeric(1))

# ============================================================
# 2. Load ROCK-ALS CSF samples
# ============================================================
samples_csv <- file.path(opt$rock_samples, "samples_all_biofluids.csv")
if (!file.exists(samples_csv)) stop("Not found: ", samples_csv)
rock_samples <- read_csv(samples_csv, show_col_types = FALSE, na = c("", "NA")) %>%
  filter(toupper(fluid) == "CSF")
message("ROCK-ALS CSF samples in 01_Samples_summary: ", nrow(rock_samples))

rock_expr_path <- opt$rock_expr
if (is.null(rock_expr_path) || !nzchar(rock_expr_path)) {
  rock_expr_path <- file.path(dirname(normalizePath(opt$rock_samples)), "csf", "00_Initialization", "Expression.csv")
}
if (!file.exists(rock_expr_path)) {
  stop("ROCK-ALS CSF expression file not found: ", rock_expr_path, "\nGive it with --rock_expr.")
}
message("ROCK-ALS CSF NPQ: ", rock_expr_path)

rock_long <- read_csv(rock_expr_path, show_col_types = FALSE, na = c("", "NA")) %>%
  mutate(PatNr = gsub("-", "_", PatNr))

# Map 08 feature names (e.g. CNS_FCN2) to ROCK-ALS targets (FCN2)
feat_map <- data.frame(
  Feature = features,
  Panel = sub("^(CNS|Immune)_.*$", "\\1", features),
  Target = sub("^(CNS|Immune)_", "", features),
  stringsAsFactors = FALSE
)
missing_t <- setdiff(feat_map$Target, unique(rock_long$Target))
if (length(missing_t) > 0L) {
  stop("ROCK-ALS data lacks model features: ", paste(missing_t, collapse = ", "),
       "\n(ROCK-ALS has the CNS panel; models using Immune-panel features cannot be applied.)")
}

rock_wide <- rock_long %>%
  filter(Target %in% feat_map$Target, PatNr %in% rock_samples$PatNr) %>%
  group_by(PatNr, Target) %>%
  summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Target, values_from = NPQ) %>%
  as.data.frame()
rownames(rock_wide) <- rock_wide$PatNr
rock_X_raw <- rock_wide[, feat_map$Target, drop = FALSE]
colnames(rock_X_raw) <- feat_map$Feature

n_before <- nrow(rock_X_raw)
rock_X_raw <- rock_X_raw[stats::complete.cases(rock_X_raw), , drop = FALSE]
message("ROCK-ALS CSF samples with all ", length(features), " features: ", nrow(rock_X_raw),
        if (nrow(rock_X_raw) < n_before) paste0(" (", n_before - nrow(rock_X_raw), " dropped: missing values)") else "")
not_in_expr <- setdiff(rock_samples$PatNr, rownames(rock_wide))
if (length(not_in_expr) > 0L) message("  ", length(not_in_expr), " CSF samples in 01_Samples_summary have no NPQ data.")

# ============================================================
# 3. Harmonisation
# ============================================================
rock_mean <- vapply(rock_X_raw, mean, numeric(1))
rock_sd   <- vapply(rock_X_raw, stats::sd, numeric(1))

harmonize <- function(X, mode) {
  out <- X
  for (f in colnames(X)) {
    out[[f]] <- switch(mode,
      none = X[[f]],
      center = X[[f]] - rock_mean[[f]] + ref_mean[[f]],
      standardize = if (rock_sd[[f]] > 0) {
        (X[[f]] - rock_mean[[f]]) / rock_sd[[f]] * ref_sd[[f]] + ref_mean[[f]]
      } else {
        X[[f]] - rock_mean[[f]] + ref_mean[[f]]
      }
    )
  }
  out
}

shift_tbl <- data.frame(
  Feature = features,
  MAXOMOD_mean = ref_mean, MAXOMOD_sd = ref_sd,
  ROCK_mean = rock_mean, ROCK_sd = rock_sd,
  mean_diff = rock_mean - ref_mean,
  mean_diff_in_MAXOMOD_SD = (rock_mean - ref_mean) / ref_sd,
  sd_ratio = rock_sd / ref_sd,
  row.names = NULL
)
write.csv(shift_tbl, file.path(output_dir, "feature_shift_MAXOMOD_vs_ROCK.csv"), row.names = FALSE)
message("Cohort shift (ROCK - MAXOMOD, in MAXOMOD SD units):")
for (i in seq_len(nrow(shift_tbl))) {
  message(sprintf("  %-14s %+6.2f SD (SD ratio %.2f)", shift_tbl$Feature[i],
                  shift_tbl$mean_diff_in_MAXOMOD_SD[i], shift_tbl$sd_ratio[i]))
}

# ============================================================
# 4. Predict with every model under every harmonisation mode
# ============================================================
predict_beta <- function(fit, newdata) {
  p <- predict(fit, newdata = newdata, type = "prob")
  col <- if ("beta" %in% colnames(p)) "beta" else colnames(p)[2]
  as.numeric(p[[col]])
}

pred_long <- list()
for (mode in harm_modes) {
  X_m <- harmonize(rock_X_raw, mode)
  for (m in names(final_models)) {
    p <- tryCatch(predict_beta(final_models[[m]], X_m), error = function(e) {
      message("  Prediction failed for ", m, " (", mode, "): ", conditionMessage(e))
      NULL
    })
    if (is.null(p)) next
    pred_long[[paste(mode, m)]] <- data.frame(
      PatNr = rownames(X_m), harmonize = mode, model = m, P_beta = p, stringsAsFactors = FALSE
    )
  }
}
pred_long <- bind_rows(pred_long)
if (nrow(pred_long) == 0L) stop("No model could produce predictions.")
pred_long$pred_k2 <- ifelse(pred_long$P_beta >= thr, "beta", "alpha")
models_ok <- unique(pred_long$model)
ens_exclude <- trimws(strsplit(opt$ensemble_exclude, ",", fixed = TRUE)[[1]])
ens_exclude <- setdiff(ens_exclude[nzchar(ens_exclude)], "none")
ens_models <- setdiff(models_ok, ens_exclude)
if (length(ens_models) == 0L) stop("--ensemble_exclude removes every model.")
message("Ensemble models: ", paste(ens_models, collapse = ", "))

# Ensemble: mean P(beta) and majority vote across models
ens <- pred_long %>%
  filter(model %in% ens_models) %>%
  group_by(PatNr, harmonize) %>%
  summarise(
    ensemble_P_beta = mean(P_beta),
    n_models = n(),
    n_vote_beta = sum(pred_k2 == "beta"),
    .groups = "drop"
  ) %>%
  mutate(
    ensemble_k2 = ifelse(ensemble_P_beta >= thr, "beta", "alpha"),
    vote_k2 = ifelse(n_vote_beta > n_models / 2, "beta",
                     ifelse(n_vote_beta < n_models / 2, "alpha", "tie")),
    model_agreement = pmax(n_vote_beta, n_models - n_vote_beta) / n_models
  )

meta_cols <- intersect(
  c("PatNr", "treatment", "group", "arm", "stratum", "sex", "age", "C9orf72", "PlateID", "SampleQC",
    "visit_label", "ALSFRS-R Score", "Calculated Vital capacity (%)", "NfL", "GFAP"),
  colnames(rock_samples)
)
rock_meta <- rock_samples[, meta_cols, drop = FALSE] %>% distinct(PatNr, .keep_all = TRUE)

write_mode_table <- function(mode) {
  wide_p <- pred_long %>%
    filter(harmonize == mode) %>%
    select(PatNr, model, P_beta) %>%
    pivot_wider(names_from = model, values_from = P_beta, names_prefix = "P_beta_")
  wide_c <- pred_long %>%
    filter(harmonize == mode) %>%
    select(PatNr, model, pred_k2) %>%
    pivot_wider(names_from = model, values_from = pred_k2, names_prefix = "k2_")
  out <- ens %>%
    filter(harmonize == mode) %>%
    select(-harmonize) %>%
    left_join(wide_p, by = "PatNr") %>%
    left_join(wide_c, by = "PatNr") %>%
    left_join(rock_meta, by = "PatNr") %>%
    relocate(any_of(meta_cols)) %>%
    arrange(PatNr)
  write.csv(out, file.path(output_dir, paste0("ROCK_ALS_CSF_predicted_k2_", mode, ".csv")), row.names = FALSE)
  out
}
tabs <- lapply(harm_modes, write_mode_table)
names(tabs) <- harm_modes
main_tab <- tabs[[main_mode]]

# Summary counts
summ <- pred_long %>%
  count(harmonize, model, pred_k2) %>%
  pivot_wider(names_from = pred_k2, values_from = n, values_fill = 0)
summ_ens <- ens %>%
  count(harmonize, ensemble_k2) %>%
  mutate(model = "ensemble_mean") %>%
  rename(pred_k2 = ensemble_k2) %>%
  pivot_wider(names_from = pred_k2, values_from = n, values_fill = 0)
summ_all <- bind_rows(summ, summ_ens) %>% arrange(harmonize, model)
write.csv(summ_all, file.path(output_dir, "prediction_counts_by_model.csv"), row.names = FALSE)

# Agreement of the ensemble call between harmonisation modes
agree <- ens %>%
  select(PatNr, harmonize, ensemble_k2) %>%
  pivot_wider(names_from = harmonize, values_from = ensemble_k2)
agree_summary <- data.frame(
  comparison = c("none vs standardize", "center vs standardize", "none vs center"),
  agreement = c(mean(agree$none == agree$standardize),
                mean(agree$center == agree$standardize),
                mean(agree$none == agree$center))
)
write.csv(agree_summary, file.path(output_dir, "harmonization_agreement.csv"), row.names = FALSE)

# By treatment arm (main mode)
if ("treatment" %in% colnames(main_tab)) {
  by_trt <- main_tab %>% count(treatment, ensemble_k2) %>%
    group_by(treatment) %>% mutate(fraction = n / sum(n)) %>% ungroup()
  write.csv(by_trt, file.path(output_dir, paste0("predicted_k2_by_treatment_", main_mode, ".csv")), row.names = FALSE)
}

# ============================================================
# 5. Plots
# ============================================================
k2_cols <- c(alpha = "#FC8D62", beta = "#8DA0CB")

# (a) Feature distributions: MAXOMOD (true k2) vs ROCK-ALS (harmonised values, ensemble-predicted k2)
rock_X_main <- harmonize(rock_X_raw, main_mode)
rock_pred_k2 <- main_tab$ensemble_k2[match(rownames(rock_X_main), main_tab$PatNr)]
rock_cohort_lab <- paste0("ROCK-ALS\n(", main_mode, ", predicted k2)")
dist_df <- bind_rows(
  data.frame(train_X, cohort = "MAXOMOD\n(training, k2)", k2 = as.character(train_y), check.names = FALSE),
  data.frame(rock_X_main, cohort = rock_cohort_lab, k2 = rock_pred_k2, check.names = FALSE)
) %>%
  pivot_longer(-c(cohort, k2), names_to = "Feature", values_to = "NPQ")
dist_df$cohort <- factor(dist_df$cohort, levels = c("MAXOMOD\n(training, k2)", rock_cohort_lab))
dist_df$k2 <- factor(dist_df$k2, levels = c("alpha", "beta"))
grp_n <- bind_rows(
  data.frame(cohort = "MAXOMOD\n(training, k2)", k2 = as.character(train_y)),
  data.frame(cohort = rock_cohort_lab, k2 = rock_pred_k2)
) %>% count(cohort, k2)
p_dist <- ggplot(dist_df, aes(x = .data$cohort, y = .data$NPQ, fill = .data$k2)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.9, alpha = 0.6, color = "gray20") +
  facet_wrap(~ Feature, scales = "free_y") +
  scale_fill_manual(values = k2_cols, name = "k2") +
  labs(
    title = "Model features: MAXOMOD training vs ROCK-ALS CSF",
    subtitle = paste0(
      "ROCK-ALS values harmonised (", main_mode, ") to the MAXOMOD scale; ROCK-ALS groups = ensemble prediction.\nn: ",
      paste(sprintf("%s %s = %d", sub("\n.*", "", grp_n$cohort), grp_n$k2, grp_n$n), collapse = ", ")
    ),
    x = NULL, y = "NPQ"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5))
ggsave(file.path(output_dir, "feature_distribution_MAXOMOD_vs_ROCK.pdf"), p_dist,
       width = 10, height = 3 + 2.6 * ceiling(length(features) / 3))

# (b) P(beta) per model (main mode)
pl_main <- pred_long %>% filter(harmonize == main_mode)
pl_main <- bind_rows(
  pl_main,
  ens %>% filter(harmonize == main_mode) %>%
    transmute(PatNr, harmonize, model = "ensemble_mean", P_beta = ensemble_P_beta, pred_k2 = ensemble_k2)
)
pl_main$model <- factor(pl_main$model, levels = c(models_ok, "ensemble_mean"))
p_prob <- ggplot(pl_main, aes(x = .data$model, y = .data$P_beta)) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.shape = NA, fill = "gray92", width = 0.6) +
  geom_jitter(aes(color = .data$pred_k2), width = 0.18, height = 0, size = 1.6, alpha = 0.85) +
  scale_color_manual(values = k2_cols, name = "Predicted") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = paste0("ROCK-ALS CSF: predicted P(beta) per model (harmonize = ", main_mode, ")"),
       x = NULL, y = "P(beta)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(output_dir, paste0("predicted_P_beta_by_model_", main_mode, ".pdf")), p_prob, width = 8, height = 5)

# (d) Predicted group by treatment arm (ensemble, main mode)
if ("treatment" %in% colnames(main_tab)) {
  p_trt <- ggplot(main_tab, aes(x = .data$treatment, fill = .data$ensemble_k2)) +
    geom_bar(position = "fill", width = 0.65) +
    geom_text(stat = "count", aes(label = after_stat(count)),
              position = position_fill(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = k2_cols, name = "Predicted k2") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "ROCK-ALS CSF: predicted k2 by treatment arm (ensemble)",
         subtitle = paste0("harmonize = ", main_mode, ", threshold = ", thr), x = NULL, y = "Fraction") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  ggsave(file.path(output_dir, paste0("predicted_k2_by_treatment_", main_mode, ".pdf")), p_trt, width = 7, height = 5)
}

# ============================================================
# 6. PCA of ROCK-ALS CSF samples (as in 09_PCA_ML_Genes.R)
#    Shape = ensemble-predicted k2 (main harmonisation mode).
#    PCA uses ROCK-ALS NPQ_adj (age + sex + PlateID removed, treatment group protected;
#    made by ROCK_ALS 00_Initialization.R), as in ROCK_ALS 06_PCA.R. Falls back to raw NPQ
#    if NPQ_adj is absent. PC1/PC2 signs are aligned to ROCK_ALS csf/06_PCA/pca_scores.csv
#    when that file exists (PCA signs are arbitrary). Predictions above still use raw NPQ.
# ============================================================
pca_dir <- file.path(output_dir, paste0("PCA_predicted_k2_", main_mode))
if (!dir.exists(pca_dir)) dir.create(pca_dir, recursive = TRUE)

label_repel_geom <- function(data, mapping, ...) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = data, mapping = mapping, inherit.aes = FALSE, size = 3.2,
                             fontface = "italic", color = "gray10", segment.color = "gray65",
                             segment.size = 0.3, min.segment.length = 0, box.padding = 0.45,
                             point.padding = 0.25, force = 2.5, max.overlaps = Inf, ...)
  } else {
    geom_text(data = data, mapping = mapping, inherit.aes = FALSE, size = 3.2,
              fontface = "italic", color = "gray10", ...)
  }
}
feat_label <- function(f) {
  ifelse(grepl("^CNS_", f), paste0("CNS:", sub("^CNS_", "", f)), paste0("Imm:", sub("^Immune_", "", f)))
}
scale_k2_shape <- scale_shape_manual(values = c(alpha = 16, beta = 17), name = "Predicted k2")

# Sample annotation (ROCK-ALS clinical variables used for colouring)
pca_meta <- main_tab %>%
  select(PatNr, ensemble_k2, any_of(c("treatment", "ALSFRS-R Score", "Calculated Vital capacity (%)",
                                      "NfL", "GFAP"))) %>%
  mutate(ensemble_k2 = factor(ensemble_k2, levels = c("alpha", "beta")))
color_vars <- list(
  list(col = "ALSFRS-R Score", file = "ALSFRS_R", log = FALSE),
  list(col = "Calculated Vital capacity (%)", file = "Vital_capacity", log = FALSE),
  list(col = "NfL", file = "NfL", log = TRUE),
  list(col = "GFAP", file = "GFAP", log = TRUE),
  list(col = "treatment", file = "treatment", log = FALSE)
)

# Value used for PCA: NPQ_adj (as ROCK_ALS 06_PCA.R) if available, else raw NPQ
pca_value <- if ("NPQ_adj" %in% colnames(rock_long)) "NPQ_adj" else "NPQ"
if (pca_value == "NPQ") message("Note: NPQ_adj not in ROCK-ALS Expression.csv; PCA uses raw NPQ.")

# Reference PCA scores from ROCK_ALS 06_PCA.R, used only to align PC signs
ref_scores_path <- file.path(dirname(normalizePath(opt$rock_samples)), "csf", "06_PCA", "pca_scores.csv")
ref_scores <- if (file.exists(ref_scores_path)) {
  read_csv(ref_scores_path, show_col_types = FALSE) %>%
    mutate(PatNr = gsub("-", "_", PatNr)) %>%
    select(PatNr, PC1, PC2)
} else {
  message("Note: ", ref_scores_path, " not found; PC signs not aligned to ROCK_ALS 06_PCA.")
  NULL
}

# All ROCK-ALS CNS-panel proteins (prefixed as in 08), same samples as the predictions
all_wide <- rock_long %>%
  filter(PatNr %in% pca_meta$PatNr) %>%
  group_by(PatNr, Target) %>%
  summarise(NPQ = mean(.data[[pca_value]], na.rm = TRUE), .groups = "drop") %>%
  mutate(Feature = paste0("CNS_", Target)) %>%
  select(PatNr, Feature, NPQ) %>%
  pivot_wider(names_from = Feature, values_from = NPQ) %>%
  as.data.frame()
rownames(all_wide) <- all_wide$PatNr
all_wide$PatNr <- NULL
all_wide <- all_wide[pca_meta$PatNr, , drop = FALSE]
all_wide <- all_wide[, colSums(is.na(all_wide)) == 0, drop = FALSE]
sd_all <- apply(all_wide, 2, stats::sd)
all_wide <- all_wide[, !is.na(sd_all) & sd_all > 0, drop = FALSE]

ml_wide <- all_wide[, intersect(features, colnames(all_wide)), drop = FALSE]

# One PCA + its plots. arrows: "normalised" (equal length, as 09 ML biplot) or
# "correlation" (length = correlation, as 09 all-protein biplot)
make_pca_plots <- function(X, set_name, set_title, arrow_features, arrow_mode, ref = NULL) {
  pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  # Align PC1/PC2 signs with a reference PCA (sign of a PC is arbitrary)
  if (!is.null(ref)) {
    ids <- intersect(rownames(pca$x), ref$PatNr)
    for (k in 1:2) {
      r_k <- stats::cor(pca$x[ids, k], ref[[paste0("PC", k)]][match(ids, ref$PatNr)])
      if (!is.na(r_k) && r_k < 0) {
        pca$x[, k] <- -pca$x[, k]
        pca$rotation[, k] <- -pca$rotation[, k]
      }
    }
  }
  ve <- summary(pca)$importance[2, 1:2] * 100
  sc <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PatNr = rownames(X), stringsAsFactors = FALSE) %>%
    left_join(pca_meta, by = "PatNr")
  ttl <- paste0("ROCK-ALS CSF PCA: ", set_title, " (", ncol(X), " proteins, ",
                if (pca_value == "NPQ_adj") "age/sex/PlateID-adjusted NPQ" else "raw NPQ", ")")
  axes <- labs(x = paste0("PC1 (", round(ve[1], 1), "%)"), y = paste0("PC2 (", round(ve[2], 1), "%)"))
  write.csv(sc, file.path(pca_dir, paste0("pca_", set_name, "_scores.csv")), row.names = FALSE)

  colour_scale <- function(cv) {
    if (cv$col == "treatment") {
      scale_color_manual(values = c("Placebo" = "gray45", "Fasudil 15mg (30mg/Day)" = "#66C2A5",
                                    "Fasudil 30mg (60mg/Day)" = "#1B7837"), name = "treatment")
    } else if (isTRUE(cv$log)) {
      scale_color_viridis_c(option = "C", name = paste0(cv$col, "\n(log10 scale)"), trans = "log10",
                            na.value = "gray75")
    } else {
      scale_color_viridis_c(option = "C", name = cv$col, na.value = "gray75")
    }
  }

  for (cv in color_vars) {
    if (!cv$col %in% names(sc)) next
    df <- sc                                   # always all samples; missing values drawn in grey
    n_na <- sum(is.na(df[[cv$col]]))
    p <- ggplot(df, aes(x = .data$PC1, y = .data$PC2, color = .data[[cv$col]], shape = .data$ensemble_k2)) +
      geom_point(size = 3.6, alpha = 0.9, stroke = 0.9) +
      colour_scale(cv) + scale_k2_shape +
      labs(title = ttl, subtitle = paste0("Shape = predicted k2 (ensemble, ", main_mode, "); n = ", nrow(df),
                                          if (n_na > 0) paste0(" (", n_na, " with missing ", cv$col, " in grey)") else "")) +
      axes + theme_classic(base_size = 13)
    ggsave(file.path(pca_dir, paste0("pca_", set_name, "_", cv$file, "_k2.pdf")), p, width = 8, height = 6)
  }

  # Biplot: all samples, coloured (and shaped) by predicted k2 only - no clinical variables
  bdf <- sc
  af <- intersect(arrow_features, colnames(X))
  if (length(af) >= 1L && nrow(bdf) >= 3L) {
    r <- stats::cor(X[, af, drop = FALSE], pca$x[, 1:2, drop = FALSE])
    span <- max(abs(c(bdf$PC1, bdf$PC2)), na.rm = TRUE)
    if (arrow_mode == "normalised") {
      len <- sqrt(rowSums(r^2)); len[len == 0] <- 1
      dir_m <- r / len
      L <- 0.38 * span
      sub_txt <- "Arrows: ML-selected proteins (direction = correlation with PC1/PC2; equal length)"
    } else {
      dir_m <- r
      L <- 0.9 * span
      sub_txt <- "Arrows: ML-selected proteins (direction and length = correlation with PC1/PC2)"
    }
    arr <- data.frame(feature_label = feat_label(af), PC1 = dir_m[, 1] * L, PC2 = dir_m[, 2] * L,
                      label_x = dir_m[, 1] * L * 1.08, label_y = dir_m[, 2] * L * 1.08)
    write.csv(data.frame(Feature = af, cor_PC1 = r[, 1], cor_PC2 = r[, 2], length = sqrt(rowSums(r^2))),
              file.path(pca_dir, paste0("pca_", set_name, "_ML_protein_correlations.csv")), row.names = FALSE)
    p_bi <- ggplot(bdf, aes(x = .data$PC1, y = .data$PC2, color = .data$ensemble_k2, shape = .data$ensemble_k2)) +
      geom_point(size = 3.4, alpha = 0.85, stroke = 0.85) +
      geom_segment(data = arr, aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2), inherit.aes = FALSE,
                   arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                   color = "gray15", linewidth = 0.65) +
      label_repel_geom(arr, aes(x = .data$label_x, y = .data$label_y, label = .data$feature_label)) +
      scale_color_manual(values = k2_cols, name = "Predicted k2") + scale_k2_shape +
      labs(title = ttl, subtitle = paste0(sub_txt, "\nColour/shape = predicted k2 (ensemble, ", main_mode, "); n = ", nrow(bdf))) +
      axes + theme_classic(base_size = 12)
    ggsave(file.path(pca_dir, paste0("pca_", set_name, "_biplot_k2.pdf")), p_bi,
           width = 9.5, height = 6.5)
  }
  invisible(pca)
}

message("PCA (ROCK-ALS CSF): ML-selected proteins n = ", nrow(ml_wide), ", p = ", ncol(ml_wide),
        "; all proteins p = ", ncol(all_wide))
make_pca_plots(ml_wide, "ML_proteins", "ML-selected proteins", features, "normalised")
make_pca_plots(all_wide, "all_proteins", "all proteins", features, "correlation", ref = ref_scores)

# ============================================================
# Console summary
# ============================================================
cat("\n============================================================\n")
cat("ROCK-ALS CSF predicted k2 (ensemble mean P(beta) >= ", thr, "; models: ",
    paste(ens_models, collapse = ", "), ")\n", sep = "")
cat("============================================================\n")
print(as.data.frame(ens %>% count(harmonize, ensemble_k2) %>%
                      pivot_wider(names_from = ensemble_k2, values_from = n, values_fill = 0)))
cat("\nAgreement of ensemble calls between harmonisation modes:\n")
print(agree_summary, row.names = FALSE)

message("\nOutput: ", normalizePath(output_dir))
message("  ROCK_ALS_CSF_predicted_k2_{none,center,standardize}.csv  (per-sample predictions, main = ", main_mode, ")")
message("  prediction_counts_by_model.csv, harmonization_agreement.csv, feature_shift_MAXOMOD_vs_ROCK.csv")
message("  predicted_k2_by_treatment_", main_mode, ".csv/.pdf")
message("  feature_distribution_MAXOMOD_vs_ROCK.pdf, predicted_P_beta_by_model_", main_mode, ".pdf")
message("  PCA_predicted_k2_", main_mode, "/  (PCA of ML-selected and all proteins; shape = predicted k2)")
