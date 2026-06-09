message("[1/6] Loading libraries...")
suppressMessages(library(optparse))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(glmnet))
suppressMessages(library(pROC))
suppressMessages(library(caret))
suppressMessages(library(randomForest))
suppressMessages(library(e1071))
suppressMessages(library(xgboost))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(fastshap))
if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
  message("      Note: install 'ggbeeswarm' for beeswarm layout; using jitter fallback.")
}
message("      Libraries loaded.")


# Rscript Script/08_Clustering_Feature_Selection.R \
#   -i CNS_immune/Results \
#   -o CNS_immune/Results/Combined_CNS_Immune_panels/Without_tears/08_Clustering_Feature_Selection \
#   -b 100 \
#   -f 0.4 \
#   -l min
# ============================================================
# Command-line options
# ============================================================
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "CNS_immune/Results",
    help = paste(
      "Panel base directory (contains CNS_panel/ and Immune_panel/ under Without_tears),",
      "OR comma-separated paths: cns_npq.xlsx,immune_npq.xlsx,metadata.xlsx"
    )
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "CNS_immune/Results/Combined_CNS_Immune_panels/Without_tears/08_Clustering_Feature_Selection",
    help = "Output directory for results and plots"
  ),
  make_option(
    c("-b", "--n_boot"),
    type = "integer",
    default = 100L,
    help = "Number of stratified bootstrap iterations (LASSO FS and ML) [default: %default]"
  ),
  make_option(
    c("-f", "--feature_freq_cutoff"),
    type = "numeric",
    default = 0.4,
    help = "LASSO feature selection frequency cutoff (0-1) [default: %default]"
  ),
  make_option(
    c("-l", "--lambda"),
    type = "character",
    default = "min",
    help = "Lambda for LASSO feature selection: 'min' or '1se' [default: %default]"
  ),
  make_option(
    c("--shap_nsim"),
    type = "integer",
    default = 2000L,
    help = "Monte Carlo samples for SHAP (fastshap) per model [default: %default]"
  ),
  make_option(
    c("--shap_top_n"),
    type = "integer",
    default = 20L,
    help = "Top N features to show in SHAP beeswarm plots [default: %default]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

parse_panel_paths <- function(input_arg) {
  if (grepl(",", input_arg, fixed = TRUE)) {
    parts <- trimws(strsplit(input_arg, ",", fixed = TRUE)[[1]])
    if (length(parts) != 3L) {
      stop(
        "--input must be a panel base directory, or exactly three comma-separated paths:\n",
        "  cns_npq.xlsx,immune_npq.xlsx,all_participants_IDs.xlsx"
      )
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

# ============================================================
# Config
# ============================================================
set.seed(42)
BS_number <- opt$n_boot
feature_freq_cutoff <- opt$feature_freq_cutoff
lambda_value <- opt$lambda
shap_nsim <- opt$shap_nsim
shap_top_n <- opt$shap_top_n
cv_folds  <- 10

if (!lambda_value %in% c("min", "1se")) {
  stop("--lambda must be 'min' or '1se'.")
}
if (feature_freq_cutoff <= 0 || feature_freq_cutoff > 1) {
  stop("--feature_freq_cutoff must be between 0 and 1.")
}

panel_paths <- parse_panel_paths(opt$input)
cns_npq_path <- panel_paths$cns
immune_npq_path <- panel_paths$immune
metadata_path <- panel_paths$metadata

output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("      Input (CNS NPQ): ", cns_npq_path)
message("      Input (Immune NPQ): ", immune_npq_path)
message("      Metadata: ", metadata_path)
message("      Output: ", output_dir)
message("      Bootstrap iterations: ", BS_number)
message("      LASSO feature freq cutoff: ", feature_freq_cutoff)
message("      LASSO lambda rule: ", lambda_value)
message("      SHAP nsim: ", shap_nsim, ", top features in plot: ", shap_top_n)

# Long CSF sample table -> wide matrix (samples x prefixed features)
npq_to_wide <- function(npq_path, panel_prefix) {
  read_excel(npq_path) %>%
    filter(SampleMatrixType == "CSF", SampleType == "Sample") %>%
    mutate(Feature = paste0(panel_prefix, "_", Target)) %>%
    select(SampleName, Feature, NPQ) %>%
    pivot_wider(names_from = Feature, values_from = NPQ) %>%
    as.data.frame()
}

# ============================================================
# Load data: CNS + Immune panels -> combined wide matrix
# ============================================================
message("[2/6] Loading and preparing data (CNS + Immune panels)...")

npq_wide_cns <- npq_to_wide(cns_npq_path, "CNS")
npq_wide_immune <- npq_to_wide(immune_npq_path, "Immune")

common_samples <- intersect(npq_wide_cns$SampleName, npq_wide_immune$SampleName)
message("      CNS features: ", ncol(npq_wide_cns) - 1L, ", samples: ", nrow(npq_wide_cns))
message("      Immune features: ", ncol(npq_wide_immune) - 1L, ", samples: ", nrow(npq_wide_immune))
message("      Samples in both panels: ", length(common_samples))

npq_wide <- npq_wide_cns %>%
  filter(.data$SampleName %in% common_samples) %>%
  inner_join(
    npq_wide_immune %>% filter(.data$SampleName %in% common_samples),
    by = "SampleName"
  )

rownames(npq_wide) <- npq_wide$SampleName
npq_wide$SampleName <- NULL

feature_panel_map <- data.frame(
  Feature = colnames(npq_wide),
  Panel = ifelse(grepl("^CNS_", colnames(npq_wide)), "CNS", "Immune"),
  Target = sub("^(CNS|Immune)_", "", colnames(npq_wide)),
  stringsAsFactors = FALSE
)
write.csv(feature_panel_map, file.path(output_dir, "feature_panel_map.csv"), row.names = FALSE)
message("      Combined matrix: ", nrow(npq_wide), " samples x ", ncol(npq_wide), " features")

all_participants_IDs <- read_excel(metadata_path)
all_participants_IDs <- all_participants_IDs %>% filter(Material == "CSF")

# Merge with k2 labels (alpha vs beta)
labels <- all_participants_IDs %>%
  select(Tube_ID, k2) %>%
  filter(!is.na(k2)) %>%
  filter(Tube_ID %in% rownames(npq_wide))

npq_model <- npq_wide[labels$Tube_ID, , drop = FALSE]
y <- factor(labels$k2, levels = c("alpha", "beta"))

cat("Model data: n =", nrow(npq_model), ", p =", ncol(npq_model), "\n")
cat("Outcome distribution:\n")
print(table(y))

# Remove columns with NA or zero variance
npq_model <- npq_model[, colSums(is.na(npq_model)) == 0, drop = FALSE]
var_filter <- apply(npq_model, 2, var) > 0
npq_model <- npq_model[, var_filter, drop = FALSE]
cat("Features after filtering:", ncol(npq_model), "\n")
message("      Data preparation complete.")

X <- as.matrix(npq_model)

# Min-max scale per feature (used for LASSO feature selection only; ML uses raw NPQ)
scale_manual <- function(df) {
  as.data.frame(apply(df, 2, function(x) (x - min(x)) / diff(range(x))))
}

stratified_boot_idx <- function(y_vec) {
  class_levels <- levels(y_vec)
  max_n <- max(table(y_vec))
  unlist(lapply(class_levels, function(cl) {
    sample(which(y_vec == cl), size = max_n, replace = TRUE)
  }))
}

# ============================================================
# Bootstrap LASSO feature selection (15_ML_multi_models.R)
# ============================================================
message("[3/6] LASSO feature selection (", BS_number, " bootstrap iterations)...")

sample_names_mat <- rownames(X)
feat_names <- colnames(X)
selected_features <- matrix(0, nrow = BS_number, ncol = length(feat_names),
                            dimnames = list(NULL, feat_names))
coef_names <- c("(Intercept)", feat_names)
coef_mat <- matrix(NA_real_, nrow = BS_number, ncol = length(coef_names),
                   dimnames = list(paste0("run_", seq_len(BS_number)), coef_names))

X_lasso <- scale_manual(as.data.frame(X))
rownames(X_lasso) <- sample_names_mat

for (i in seq_len(BS_number)) {
  if (i %% 10 == 0 || i == 1L) {
    message("  [LASSO BS ", i, "/", BS_number, "]")
  }
  boot_idx <- stratified_boot_idx(y)
  X_boot <- as.matrix(X_lasso[boot_idx, , drop = FALSE])
  y_boot <- y[boot_idx]

  nfolds_lasso <- max(3L, min(10L, nrow(X_boot)))
  cv_fit <- tryCatch(
    cv.glmnet(X_boot, y_boot, alpha = 1, family = "binomial", nfolds = nfolds_lasso),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) next

  lambda_pick <- if (lambda_value == "min") cv_fit$lambda.min else cv_fit$lambda.1se
  lasso_fit <- glmnet(X_boot, y_boot, alpha = 1, lambda = lambda_pick, family = "binomial")

  beta <- coef(lasso_fit)
  selected_features[i, ] <- as.numeric(beta[-1, 1] != 0)
  beta_vec <- as.numeric(beta)
  names(beta_vec) <- coef_names
  coef_mat[i, ] <- beta_vec
}

feature_freq <- colMeans(selected_features)
selected_final <- names(feature_freq[feature_freq > feature_freq_cutoff])

message("      Features selected in >", feature_freq_cutoff * 100, "% of bootstraps: ", length(selected_final))
if (length(selected_final) == 0L) {
  stop(
    "No features passed the frequency cutoff. Lower --feature_freq_cutoff ",
    "(current: ", feature_freq_cutoff, ") or increase --n_boot."
  )
}

writeLines(selected_final, file.path(output_dir, "selected_features.txt"))

feature_freq_df <- data.frame(
  Feature = names(feature_freq),
  Frequency = as.numeric(feature_freq),
  stringsAsFactors = FALSE
) %>%
  filter(.data$Frequency > 0)

feature_se <- apply(selected_features, 2, function(x) stats::sd(x) / sqrt(length(x)))
feature_freq_df <- feature_freq_df %>%
  left_join(
    data.frame(Feature = names(feature_se), SE = as.numeric(feature_se)),
    by = "Feature"
  )

write.csv(feature_freq_df, file.path(output_dir, "feature_selection_stability.csv"), row.names = FALSE)

feature_plot <- ggplot(feature_freq_df, aes(x = reorder(.data$Feature, -.data$Frequency), y = .data$Frequency)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = .data$Frequency - .data$SE, ymax = .data$Frequency + .data$SE),
                width = 0.2, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = feature_freq_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "LASSO Feature Selection Stability",
       subtitle = paste0("Dashed line: cutoff = ", feature_freq_cutoff, " (lambda.", lambda_value, ")"),
       x = "Feature", y = "Selection frequency") +
  theme_classic()

ggsave(
  file.path(output_dir, "feature_selection_stability.pdf"),
  feature_plot,
  width = 10,
  height = max(6, nrow(feature_freq_df) / 8)
)

# Subset to selected features (raw NPQ) for downstream ML
X <- X[, selected_final, drop = FALSE]
npq_model <- npq_model[, selected_final, drop = FALSE]
feature_panel_map <- feature_panel_map %>%
  filter(.data$Feature %in% selected_final)
write.csv(feature_panel_map, file.path(output_dir, "feature_panel_map_selected.csv"), row.names = FALSE)

message("      ML will use ", ncol(X), " selected features (raw NPQ, no scaling)")

# ============================================================
# Helper functions (from 15_ML_multi_models.R)
# ============================================================
get_roc_points <- function(y_true, y_score, run_id) {
  roc_i <- pROC::roc(response = y_true, predictor = y_score, quiet = TRUE,
                     levels = c("alpha", "beta"), direction = "<")
  Sp <- roc_i$specificities
  Sn <- roc_i$sensitivities
  eps <- 1e-12
  if (!any(abs(Sp - 1) < eps & abs(Sn - 0) < eps)) { Sp <- c(1, Sp); Sn <- c(0, Sn) }
  if (!any(abs(Sp - 0) < eps & abs(Sn - 1) < eps)) { Sp <- c(Sp, 0); Sn <- c(Sn, 1) }
  df <- data.frame(Sp = Sp, Sn = Sn, n = seq_along(Sn), run = run_id)
  df$auc <- as.numeric(pROC::auc(roc_i))
  df
}

get_probs <- function(fit, newx, y_true) {
  probs <- predict(fit, newx, type = "prob")
  pos <- levels(y_true)[2]
  if (!pos %in% colnames(probs)) pos <- colnames(probs)[2]
  probs[[pos]]
}

caret_prob_beta <- function(fit, newdata, feature_names) {
  nd <- as.data.frame(newdata)
  if (ncol(nd) == length(feature_names)) {
    colnames(nd) <- feature_names
  }
  probs <- predict(fit, newdata = nd, type = "prob")
  pos <- "beta"
  if (!pos %in% colnames(probs)) pos <- colnames(probs)[2]
  out <- probs[, pos, drop = TRUE]
  as.numeric(out)
}

fit_caret_model <- function(model_key, train_df, train_control) {
  switch(
    model_key,
    lasso = {
      lasso_grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 50))
      train(class ~ ., data = train_df, method = "glmnet",
            metric = "ROC", trControl = train_control, tuneGrid = lasso_grid)
    },
    elastic_net = {
      enet_grid <- expand.grid(alpha = 0.5, lambda = 10^seq(-3, 1, length = 50))
      train(class ~ ., data = train_df, method = "glmnet",
            metric = "ROC", trControl = train_control, tuneGrid = enet_grid)
    },
    ridge = {
      ridge_grid <- expand.grid(alpha = 0, lambda = 10^seq(-3, 1, length = 50))
      train(class ~ ., data = train_df, method = "glmnet",
            metric = "ROC", trControl = train_control, tuneGrid = ridge_grid)
    },
    rf = train(class ~ ., data = train_df, method = "rf",
               metric = "ROC", trControl = train_control, tuneLength = 3, ntree = 500),
    svm_radial = train(class ~ ., data = train_df, method = "svmRadial",
                       metric = "ROC", trControl = train_control, tuneLength = 3),
    svm_linear = train(class ~ ., data = train_df, method = "svmLinear",
                       metric = "ROC", trControl = train_control, tuneLength = 3),
    xgb = train(class ~ ., data = train_df, method = "xgbTree",
                metric = "ROC", trControl = train_control, tuneLength = 3, verbosity = 0),
    stop("Unknown model: ", model_key)
  )
}

# fastshap::explain() fails with caret models (foreach/cbind path); use explain_column directly.
compute_shap_matrix <- function(fit, X_df, nsim = 100L) {
  X <- as.data.frame(X_df)
  if (!is.null(rownames(X_df))) rownames(X) <- rownames(X_df)
  feat_cols <- colnames(X)

  pred_wrapper <- function(object, newdata) {
    caret_prob_beta(object, newdata, feat_cols)
  }

  n <- nrow(X)
  p <- ncol(X)
  shap_mat <- matrix(
    0,
    nrow = n,
    ncol = p,
    dimnames = list(rownames(X), feat_cols)
  )

  for (j in seq_len(p)) {
    phi_runs <- matrix(NA_real_, nrow = n, ncol = nsim)
    for (s in seq_len(nsim)) {
      phi_runs[, s] <- fastshap:::explain_column(
        object = fit,
        X = X,
        column = j,
        pred_wrapper = pred_wrapper,
        newdata = X
      )
    }
    shap_mat[, j] <- rowMeans(phi_runs, na.rm = TRUE)
  }

  shap_mat
}

plot_shap_beeswarm <- function(shap_mat, X_df, top_n = 20L, title = NULL) {
  mean_abs <- sort(colMeans(abs(shap_mat)), decreasing = TRUE)
  feats <- names(mean_abs)[seq_len(min(top_n, length(mean_abs)))]
  if (length(feats) == 0L) return(NULL)

  shap_long <- as.data.frame(shap_mat[, feats, drop = FALSE]) %>%
    mutate(Sample = rownames(X_df)) %>%
    pivot_longer(-Sample, names_to = "Feature", values_to = "SHAP")

  feat_long <- as.data.frame(X_df[, feats, drop = FALSE]) %>%
    mutate(Sample = rownames(X_df)) %>%
    pivot_longer(-Sample, names_to = "Feature", values_to = "Feature_value")

  plot_df <- left_join(shap_long, feat_long, by = c("Sample", "Feature"))
  plot_df$Feature <- factor(plot_df$Feature, levels = rev(feats))

  n_samples <- length(unique(plot_df$Sample))
  x_rng <- range(plot_df$SHAP, na.rm = TRUE)
  x_pad <- max(0.02, diff(x_rng) * 0.12)
  x_lo <- x_rng[1] - x_pad
  x_hi <- x_rng[2] + x_pad

  # groupOnX = FALSE: spread all samples vertically within each feature row (standard SHAP beeswarm).
  # groupOnX = TRUE (ggbeeswarm default) only stacks points with similar SHAP → thin "line" when n is small.
  beeswarm_geom <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    # orientation = "y": vertical spread within each feature row (standard SHAP beeswarm).
    ggbeeswarm::geom_quasirandom(
      width = 0.5,
      varwidth = TRUE,
      orientation = "y",
      size = 1.8,
      alpha = 0.85
    )
  } else {
    message("      Install 'ggbeeswarm' for standard SHAP beeswarm layout; using wide jitter fallback.")
    geom_jitter(width = 0, height = 0.45, size = 1.8, alpha = 0.85)
  }

  ggplot(plot_df, aes(x = SHAP, y = Feature, color = Feature_value)) +
    beeswarm_geom +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
    annotate(
      "text",
      x = x_lo,
      y = 0.35,
      label = "SHAP < 0\n\u2190 toward alpha",
      hjust = 0,
      vjust = 1,
      size = 3.2,
      fontface = "italic",
      color = "#2166AC",
      lineheight = 0.95
    ) +
    annotate(
      "text",
      x = x_hi,
      y = 0.35,
      label = "SHAP > 0\ntoward beta \u2192",
      hjust = 1,
      vjust = 1,
      size = 3.2,
      fontface = "italic",
      color = "#B2182B",
      lineheight = 0.95
    ) +
    scale_color_viridis_c(option = "C", name = "NPQ\n(feature value)") +
    scale_x_continuous(limits = c(x_lo, x_hi), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(0.8, 0.8))) +
    coord_cartesian(clip = "off") +
    labs(
      title = title,
      subtitle = paste0(
        "Outcome explained: P(beta). Each point = one sample (n = ", n_samples, "). ",
        "Color = raw NPQ (purple low, yellow high)."
      ),
      x = "SHAP value (contribution to predicted P[beta])",
      y = NULL,
      caption = paste0(
        "Vertical spread within a row shows how many samples share similar SHAP values; ",
        "with n = ", n_samples, " the cloud is thinner than in large-cohort papers (often n > 500)."
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40", margin = ggplot2::margin(t = 6)),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 12, r = 8, b = 8, l = 8)
    )
}

# ============================================================
# Stratified Bootstrap: train on bootstrap, test on OOB
# Models: LASSO, Elastic Net, Ridge, RF, SVM-radial, SVM-linear, XGBoost
# ============================================================
models <- c("lasso", "elastic_net", "ridge", "rf", "svm_radial", "svm_linear", "xgb")

auc_list   <- setNames(lapply(models, function(m) numeric(BS_number)), models)
rocc_store <- setNames(lapply(models, function(m) data.frame()), models)

train_control <- trainControl(
  method = "cv", number = cv_folds,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Prepare full data frame for caret
full_df <- as.data.frame(X)
full_df$class <- y

sample_names <- rownames(X)
split_record <- vector("list", BS_number)

message("[4/6] Running ", BS_number, " bootstrap iterations (7 models each)...")
message("      Models: ", paste(models, collapse = ", "))
t_start <- proc.time()

for (i in seq_len(BS_number)) {
  message("  [BS ", i, "/", BS_number, "] training & evaluating...")

  # Stratified bootstrap (sample with replacement, balanced)
  class_levels <- levels(y)
  max_n <- max(table(y))
  train_idx <- unlist(lapply(class_levels, function(cl) {
    sample(which(y == cl), size = max_n, replace = TRUE)
  }))

  # OOB = samples not in bootstrap

  oob_idx <- setdiff(seq_len(nrow(X)), unique(train_idx))
  if (length(oob_idx) < 3) next

  # Record which samples went to train vs OOB
  split_record[[i]] <- data.frame(
    bootstrap = i,
    SampleName = sample_names,
    role = ifelse(seq_len(nrow(X)) %in% unique(train_idx), "train", "OOB"),
    class = as.character(y),
    stringsAsFactors = FALSE
  )

  train_df <- full_df[train_idx, ]
  X_test   <- full_df[oob_idx, -ncol(full_df), drop = FALSE]
  y_test   <- full_df$class[oob_idx]

  # Train models on bootstrap sample
  # --- LASSO (alpha=1) ---
  lasso_grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 50))
  lasso_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "glmnet",
          metric = "ROC", trControl = train_control, tuneGrid = lasso_grid),
    error = function(e) NULL)


  # --- Elastic Net (alpha=0.5) ---
  enet_grid <- expand.grid(alpha = 0.5, lambda = 10^seq(-3, 1, length = 50))
  enet_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "glmnet",
          metric = "ROC", trControl = train_control, tuneGrid = enet_grid),
    error = function(e) NULL)

  # --- Ridge (alpha=0) ---
  ridge_grid <- expand.grid(alpha = 0, lambda = 10^seq(-3, 1, length = 50))
  ridge_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "glmnet",
          metric = "ROC", trControl = train_control, tuneGrid = ridge_grid),
    error = function(e) NULL)

  # --- Random Forest ---
  rf_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "rf",
          metric = "ROC", trControl = train_control, tuneLength = 3, ntree = 500),
    error = function(e) NULL)

  # --- SVM Radial ---
  svr_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "svmRadial",
          metric = "ROC", trControl = train_control, tuneLength = 3),
    error = function(e) NULL)

  # --- SVM Linear ---
  svl_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "svmLinear",
          metric = "ROC", trControl = train_control, tuneLength = 3),
    error = function(e) NULL)

  # --- XGBoost ---
  xgb_fit <- tryCatch(
    train(class ~ ., data = train_df, method = "xgbTree",
          metric = "ROC", trControl = train_control, tuneLength = 3,
          verbosity = 0),
    error = function(e) NULL)

  # Evaluate on OOB
  fit_list <- list(lasso = lasso_fit, elastic_net = enet_fit, ridge = ridge_fit,
                   rf = rf_fit, svm_radial = svr_fit, svm_linear = svl_fit, xgb = xgb_fit)

  # Need both classes in OOB for ROC
  if (length(unique(y_test)) < 2) next

  for (m in models) {
    fit_m <- fit_list[[m]]
    if (is.null(fit_m)) { auc_list[[m]][i] <- NA; next }
    p_m <- tryCatch(get_probs(fit_m, X_test, y_test), error = function(e) NULL)
    if (is.null(p_m) || length(p_m) != length(y_test)) { auc_list[[m]][i] <- NA; next }
    roc_pts <- tryCatch(get_roc_points(y_test, p_m, i), error = function(e) NULL)
    if (is.null(roc_pts)) { auc_list[[m]][i] <- NA; next }
    rocc_store[[m]] <- rbind(rocc_store[[m]], roc_pts[, c("Sp", "Sn", "n", "run")])
    auc_list[[m]][i] <- roc_pts$auc[1]
  }
}

t_elapsed <- (proc.time() - t_start)["elapsed"]
message("      Bootstrap complete. Elapsed: ", round(t_elapsed / 60, 1), " min")

# Save train/OOB split record
split_df <- bind_rows(split_record)
write.csv(split_df, file.path(output_dir, "bootstrap_sample_splits.csv"), row.names = FALSE)
message("      Sample split table saved: bootstrap_sample_splits.csv")

# Save per-bootstrap AUC for each model
auc_per_bs <- data.frame(bootstrap = seq_len(BS_number))
for (m in models) {
  auc_per_bs[[m]] <- auc_list[[m]]
}
write.csv(auc_per_bs, file.path(output_dir, "bootstrap_AUC_per_iteration.csv"), row.names = FALSE)
message("      Per-iteration AUC table saved: bootstrap_AUC_per_iteration.csv")

# ============================================================
# Summarize AUC
# ============================================================
message("[5/6] Summarizing results...")
auc_summary <- data.frame(
  Model = models,
  Mean_AUC = sapply(auc_list, function(x) mean(x, na.rm = TRUE)),
  SD_AUC   = sapply(auc_list, function(x) sd(x, na.rm = TRUE)),
  N_valid  = sapply(auc_list, function(x) sum(!is.na(x))),
  stringsAsFactors = FALSE
)
auc_summary <- auc_summary %>% arrange(desc(Mean_AUC))

cat("\n============================================================\n")
cat("        Model Performance (OOB Bootstrap AUC)\n")
cat("============================================================\n")
print(auc_summary, row.names = FALSE)

write.csv(auc_summary, file.path(output_dir, "model_AUC_summary.csv"), row.names = FALSE)

# ============================================================
# SHAP values: final model on all samples + beeswarm plots
# ============================================================
message("      Computing SHAP values (final model on all samples)...")
shap_dir <- file.path(output_dir, "shap")
if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)

X_shap <- full_df[, colnames(full_df) != "class", drop = FALSE]
rownames(X_shap) <- rownames(full_df)

for (m in models) {
  message("        SHAP: ", m, " ...")
  fit_final <- tryCatch(
    fit_caret_model(m, full_df, train_control),
    error = function(e) NULL
  )
  if (is.null(fit_final)) {
    message("          skipped (model fit failed)")
    next
  }

  shap_mat <- tryCatch(
    compute_shap_matrix(fit_final, X_shap, nsim = shap_nsim),
    error = function(e) {
      message("          SHAP error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(shap_mat)) {
    message("          skipped (SHAP computation failed)")
    next
  }

  shap_out <- as.data.frame(shap_mat)
  shap_out <- cbind(SampleName = rownames(X_shap), shap_out)
  write.csv(shap_out, file.path(shap_dir, paste0("SHAP_values_", m, ".csv")), row.names = FALSE)

  shap_summary <- data.frame(
    Feature = colnames(shap_mat),
    mean_abs_SHAP = colMeans(abs(shap_mat)),
    mean_SHAP = colMeans(shap_mat),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(.data$mean_abs_SHAP))
  write.csv(shap_summary, file.path(shap_dir, paste0("SHAP_summary_", m, ".csv")), row.names = FALSE)

  p_shap <- plot_shap_beeswarm(
    shap_mat,
    X_shap,
    top_n = shap_top_n,
    title = paste0("SHAP beeswarm: ", m)
  )
  if (!is.null(p_shap)) {
    ggsave(
      file.path(shap_dir, paste0("SHAP_beeswarm_", m, ".pdf")),
      p_shap,
      width = 7,
      height = min(7, shap_top_n * 0.45 + 2.5)
    )
    ggsave(
      file.path(shap_dir, paste0("SHAP_beeswarm_", m, ".png")),
      p_shap,
      width = 7,
      height = min(7, shap_top_n * 0.45 + 2.5),
      dpi = 300
    )
  }

  message("          done (top feature: ", shap_summary$Feature[1], ")")
}

# ============================================================
# Mean ROC plot (like 15_ML_multi_models.R)
# ============================================================
message("[6/6] Generating plots...")
summary_df <- bind_rows(lapply(models, function(m) {
  rocc <- rocc_store[[m]]
  if (nrow(rocc) == 0) return(NULL)
  Sp_mean <- aggregate(Sp ~ n, rocc, mean)$Sp
  Sn_mean <- aggregate(Sn ~ n, rocc, mean)$Sn
  Sp_sd   <- aggregate(Sp ~ n, rocc, sd)$Sp
  Sn_sd   <- aggregate(Sn ~ n, rocc, sd)$Sn

  out <- data.frame(
    method   = m,
    FPR      = 1 - Sp_mean,
    TPR      = Sn_mean,
    FPR_sd   = Sp_sd,
    TPR_sd   = Sn_sd,
    mean_auc = mean(auc_list[[m]], na.rm = TRUE),
    sd_auc   = sd(auc_list[[m]], na.rm = TRUE)
  )

  eps <- 1e-12
  if (!any(out$FPR <= eps & out$TPR <= eps)) {
    out <- rbind(
      data.frame(method = m, FPR = 0, TPR = 0, FPR_sd = 0, TPR_sd = 0,
                 mean_auc = out$mean_auc[1], sd_auc = out$sd_auc[1]),
      out)
  }
  if (!any(abs(out$FPR - 1) <= eps & abs(out$TPR - 1) <= eps)) {
    out <- rbind(out,
      data.frame(method = m, FPR = 1, TPR = 1, FPR_sd = 0, TPR_sd = 0,
                 mean_auc = out$mean_auc[1], sd_auc = out$sd_auc[1]))
  }
  out[order(out$FPR, out$TPR), ]
}))

auc_text <- summary_df %>%
  group_by(method) %>%
  summarise(txt = sprintf("%s: %.3f \u00B1 %.3f",
                          unique(method), unique(mean_auc), unique(sd_auc)),
            .groups = "drop") %>%
  arrange(desc(as.numeric(sub(":.*", "", gsub(".*?(\\d+\\.\\d+).*", "\\1", .$txt))))) %>%
  pull(txt) %>%
  paste(collapse = "\n")

# Reorder for legend (by AUC)
summary_df$method <- factor(summary_df$method,
  levels = auc_summary$Model)

p_roc <- ggplot(summary_df, aes(x = FPR, y = TPR, color = method, fill = method)) +
  geom_ribbon(aes(ymin = pmax(0, TPR - 0.95 * TPR_sd),
                  ymax = pmin(1, TPR + 0.95 * TPR_sd)),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       title = "Multi-Model Comparison: alpha vs beta (CNS + Immune CSF NPQ)",
       subtitle = paste0("OOB Bootstrap (", BS_number, " iterations) AUC \u00B1 SD:\n", auc_text)) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        plot.subtitle = element_text(size = 9))

ggsave(file.path(output_dir, "mean_ROC_multimodel.pdf"), p_roc, width = 8, height = 7)
ggsave(file.path(output_dir, "mean_ROC_multimodel.png"), p_roc, width = 8, height = 7, dpi = 300)
cat("\nMean ROC plot saved to:", output_dir, "\n")

# ============================================================
# AUC boxplot comparison
# ============================================================
auc_long <- bind_rows(lapply(models, function(m) {
  data.frame(Model = m, AUC = auc_list[[m]][!is.na(auc_list[[m]])])
}))
auc_long$Model <- factor(auc_long$Model, levels = auc_summary$Model)

p_box <- ggplot(auc_long, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "AUC Distribution Across Bootstrap Iterations",
       x = NULL, y = "AUC (OOB)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(output_dir, "AUC_boxplot_multimodel.pdf"), p_box, width = 7, height = 5)
ggsave(file.path(output_dir, "AUC_boxplot_multimodel.png"), p_box, width = 7, height = 5, dpi = 300)

message("\n=== All done! ===")
message("Output directory: ", normalizePath(output_dir))
message("Files generated:")
message("  - feature_panel_map.csv")
message("  - feature_panel_map_selected.csv")
message("  - selected_features.txt")
message("  - feature_selection_stability.csv/pdf")
message("  - model_AUC_summary.csv")
message("  - bootstrap_AUC_per_iteration.csv")
message("  - bootstrap_sample_splits.csv")
message("  - shap/SHAP_values_<model>.csv")
message("  - shap/SHAP_summary_<model>.csv")
message("  - shap/SHAP_beeswarm_<model>.pdf/png")
message("  - mean_ROC_multimodel.pdf/png")
message("  - AUC_boxplot_multimodel.pdf/png")
