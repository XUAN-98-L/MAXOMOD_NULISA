message("[1/7] Loading libraries...")
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

# ------------------------------------------------------------
# Workflow (per train/validation split, no information leakage):
#   1. Generate all stratified repeated random splits up front:
#        within each class, --train_frac (default 70%) -> train, the rest -> validation
#        (sampling without replacement; class proportions kept in both sets)
#   2. For each split:
#        a. ONE LASSO fit (cv.glmnet) on the TRAIN samples only
#           (min-max scaling fitted on train); non-zero coefficients = selected
#        b. Train 7 models on train using that split's selected features
#        c. Evaluate AUC on the validation samples
#   3. Final feature set = features selected in more than -f
#      (--feature_freq_cutoff) of the splits -> selected_features.txt,
#      used for SHAP on a model fitted to all samples (interpretation only).
#
# --vis FALSE : run splits/feature selection/models/SHAP values, save to 08_model_results.rds
#               and 08_final_models.rds (no plots)
# --vis TRUE  : only make plots from those .rds files (same -i/--panel/-o so the folder matches)
# (omitted)   : both
# Note: --panel has no "-p" short form (the future/parallelly package claims -p).
# Examples:
# Rscript Script/08_Clustering_Feature_Selection.R -i CNS_immune/Results --panel both
# Rscript Script/08_Clustering_Feature_Selection.R -i CNS_immune/Results --panel cns
# Rscript Script/08_Clustering_Feature_Selection.R -i CNS_immune/Results --panel cns --matrix Serum
# Rscript Script/08_Clustering_Feature_Selection.R -i CNS_immune/Results --panel immune \
#   -b 100 -f 0.4 -l min
# ------------------------------------------------------------

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
      "OR comma-separated paths. --panel both: cns_npq.xlsx,immune_npq.xlsx,metadata.xlsx;",
      "--panel cns/immune: npq.xlsx,metadata.xlsx (three paths also accepted)"
    )
  ),
  make_option(
    c("--panel"),
    type = "character",
    default = "both",
    help = "Which panel(s) to use: 'cns', 'immune' or 'both' [default: %default]"
  ),
  make_option(
    c("--metric"),
    type = "character",
    default = "ROC",
    help = paste(
      "Tuning metric for the penalised models (LASSO, elastic net, ridge): 'ROC' (maximise CV ROC AUC)",
      "or 'logLoss' (minimise CV log loss). RF/SVM/XGBoost always use ROC. With logLoss the output",
      "folder gets a '_logLoss' suffix [default: %default]"
    )
  ),
  make_option(
    c("--matrix"),
    type = "character",
    default = "CSF",
    help = "Sample matrix to use: 'CSF', 'Serum' or 'Plasma' (case-insensitive) [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = paste(
      "Output directory. Default: <input>/Combined_CNS_Immune_panels (both),",
      "<input>/CNS_panel (cns) or <input>/Immune_panel (immune),",
      "followed by /Without_tears/08_Clustering_Feature_Selection"
    )
  ),
  make_option(
    c("-b", "--n_boot"),
    type = "integer",
    default = 100L,
    help = "Number of stratified repeated random splits (train / validation) [default: %default]"
  ),
  make_option(
    c("--train_frac"),
    type = "numeric",
    default = 0.7,
    help = "Fraction of each class assigned to train in every split; the rest is validation [default: %default]"
  ),
  make_option(
    c("-f", "--feature_freq_cutoff"),
    type = "numeric",
    default = 0.4,
    help = paste(
      "Fraction of splits (0-1) in which a feature must be selected by LASSO",
      "to enter the final feature set (selected_features.txt, SHAP) [default: %default]"
    )
  ),
  make_option(
    c("--min_features"),
    type = "integer",
    default = 2L,
    help = "Minimum features per split; if LASSO keeps fewer, the lambda path is followed until N are non-zero [default: %default]"
  ),
  make_option(
    c("-l", "--lambda"),
    type = "character",
    default = "min",
    help = "Lambda for LASSO feature selection: 'min' or '1se' [default: %default]"
  ),
  make_option(
    c("--vis"),
    type = "character",
    default = NULL,
    help = paste(
      "FALSE: run models only and save results (no plots);",
      "TRUE: plots only, from results saved by a previous --vis FALSE run;",
      "omitted: run models then plot"
    )
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 42L,
    help = "Random seed [default: %default]"
  ),
  make_option(
    c("--shap_nsim"),
    type = "integer",
    default = 2000L,
    help = "Monte Carlo samples for SHAP (fastshap) per model [default: %default]"
  ),
  make_option(
    c("--plot_min_freq"),
    type = "numeric",
    default = 0.1,
    help = "Only features with selection frequency >= this are drawn in feature_selection_stability.pdf [default: %default]"
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

# ============================================================
# Config
# ============================================================
panel_choice <- tolower(opt$panel)
if (!panel_choice %in% c("cns", "immune", "both")) {
  stop("--panel must be 'cns', 'immune' or 'both'.")
}
panel_prefixes <- switch(panel_choice,
  cns = "CNS",
  immune = "Immune",
  both = c("CNS", "Immune")
)
panel_label <- switch(panel_choice,
  cns = "CNS",
  immune = "Immune",
  both = "CNS + Immune"
)

# Sample matrix: NPQ files use upper case (CSF/SERUM/PLASMA), metadata uses CSF/Serum/Plasma
matrix_key <- toupper(trimws(opt$matrix))
if (!matrix_key %in% c("CSF", "SERUM", "PLASMA")) {
  stop("--matrix must be 'CSF', 'Serum' or 'Plasma'.")
}
matrix_npq <- matrix_key
matrix_meta <- switch(matrix_key, CSF = "CSF", SERUM = "Serum", PLASMA = "Plasma")
panel_label <- paste0(panel_label, ", ", matrix_meta)   # used in plot titles

# Tuning metric for the penalised (glmnet) models
tune_metric <- switch(tolower(trimws(opt$metric)),
  "roc" = "ROC",
  "logloss" = "logLoss",
  stop("--metric must be 'ROC' or 'logLoss'.")
)
if (tune_metric == "logLoss") panel_label <- paste0(panel_label, ", logLoss-tuned")

vis_arg <- if (is.null(opt$vis)) "BOTH" else toupper(trimws(opt$vis))
vis_mode <- switch(vis_arg,
  "FALSE" = "compute", "F" = "compute",
  "TRUE" = "plot", "T" = "plot",
  "BOTH" = "both",
  stop("--vis must be TRUE or FALSE (or omitted to run both).")
)
run_compute <- vis_mode %in% c("compute", "both")
run_plot <- vis_mode %in% c("plot", "both")

set.seed(opt$seed)
BS_number <- opt$n_boot
feature_freq_cutoff <- opt$feature_freq_cutoff
min_features <- max(2L, opt$min_features)  # caret glmnet needs >= 2 predictors
lambda_value <- opt$lambda
shap_nsim <- opt$shap_nsim
shap_top_n <- opt$shap_top_n
plot_min_freq <- opt$plot_min_freq
cv_folds <- 10L

train_frac <- opt$train_frac
if (train_frac <= 0 || train_frac >= 1) stop("--train_frac must be between 0 and 1 (exclusive).")
if (!lambda_value %in% c("min", "1se")) stop("--lambda must be 'min' or '1se'.")
if (feature_freq_cutoff <= 0 || feature_freq_cutoff > 1) {
  stop("--feature_freq_cutoff must be between 0 and 1.")
}

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
      npq <- setNames(list(parts[1]), prefixes)
      return(list(npq = npq, metadata = parts[2]))
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
  # Metadata: taken from the first panel used (CNS for 'cns'/'both', Immune for 'immune')
  metadata <- file.path(input_arg, panel_dir[[prefixes[1]]], subdir, "00_Initialization",
                        "all_participants_IDs.xlsx")
  list(npq = npq, metadata = metadata)
}

default_output_dir <- function(input_arg, panel_choice, matrix_meta, tune_metric = "ROC") {
  # CSF keeps the original folder name; Serum/Plasma get a suffix; logLoss tuning adds "_logLoss"
  leaf_name <- if (matrix_meta == "CSF") {
    "08_Clustering_Feature_Selection"
  } else {
    paste0("08_Clustering_Feature_Selection_", matrix_meta)
  }
  if (tune_metric == "logLoss") leaf_name <- paste0(leaf_name, "_logLoss")
  leaf <- file.path("Without_tears", leaf_name)
  if (grepl(",", input_arg, fixed = TRUE)) {
    return(paste0(leaf_name, "_", panel_choice))
  }
  sub_dir <- switch(panel_choice,
    both = "Combined_CNS_Immune_panels",
    cns = "CNS_panel",
    immune = "Immune_panel"
  )
  file.path(input_arg, sub_dir, leaf)
}

panel_paths <- parse_panel_paths(opt$input, panel_choice)
metadata_path <- panel_paths$metadata

output_dir <- if (is.null(opt$output) || !nzchar(opt$output)) {
  default_output_dir(opt$input, panel_choice, matrix_meta, tune_metric)
} else {
  opt$output
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
results_rds <- file.path(output_dir, "08_model_results.rds")
models_rds <- file.path(output_dir, "08_final_models.rds")

message("      Mode: ", switch(vis_mode, compute = "models only (--vis FALSE)",
                               plot = "plots only (--vis TRUE)", both = "models + plots"))
message("      Panel(s) and matrix: ", panel_label)
for (p in names(panel_paths$npq)) message("      Input (", p, " NPQ): ", panel_paths$npq[[p]])
message("      Metadata: ", metadata_path)
message("      Output: ", output_dir)
message("      Splits: ", BS_number, " stratified random splits (train ", train_frac * 100,
        "% / validation ", (1 - train_frac) * 100, "% per class)")
message("      Feature selection: one LASSO per train split (lambda.", lambda_value, ")")
message("      Final feature cutoff (fraction of splits): ", feature_freq_cutoff)
message("      Tuning metric (LASSO/elastic net/ridge): ", tune_metric,
        if (tune_metric == "logLoss") " (minimised)" else " (maximised)", "; RF/SVM/XGBoost: ROC")
message("      SHAP nsim: ", shap_nsim, ", top features in plot: ", shap_top_n)

# Long sample table (one matrix) -> wide matrix (samples x prefixed features)
npq_to_wide <- function(npq_path, panel_prefix) {
  read_excel(npq_path) %>%
    filter(toupper(SampleMatrixType) == matrix_npq, SampleType == "Sample") %>%
    mutate(Feature = paste0(panel_prefix, "_", Target)) %>%
    select(SampleName, Feature, NPQ) %>%
    pivot_wider(names_from = Feature, values_from = NPQ) %>%
    as.data.frame()
}

# ============================================================
# Helpers: resampling
# ============================================================
# Stratified random split (without replacement): within each class, round(n * train_frac)
# samples go to train (at least 2) and the rest to validation (at least 1).
make_split <- function(y_vec, train_frac) {
  train_idx <- unlist(lapply(levels(y_vec), function(cl) {
    w <- which(y_vec == cl)
    if (length(w) < 3L) stop("Class '", cl, "' has fewer than 3 samples; cannot split.")
    n_tr <- min(length(w) - 1L, max(2L, round(length(w) * train_frac)))
    w[sample.int(length(w), n_tr)]
  }))
  train_idx <- sort(train_idx)
  list(train = train_idx, val = setdiff(seq_along(y_vec), train_idx))
}

# Stratified CV folds that keep all copies of a bootstrapped sample in the same fold
# (otherwise duplicated samples leak between CV-train and CV-test folds).
# Returns the fold id for every row of sample_id.
grouped_strat_foldid <- function(sample_id, y_vec, k) {
  uid <- unique(sample_id)
  ycl <- y_vec[match(uid, sample_id)]
  fold <- integer(length(uid))
  for (cl in unique(as.character(ycl))) {
    w <- which(ycl == cl)
    fold[w] <- sample(rep_len(seq_len(k), length(w)))
  }
  fold[match(sample_id, uid)]
}

n_cv_folds <- function(sample_id, y_vec, k_max = cv_folds) {
  uid <- unique(sample_id)
  min_class <- min(table(factor(y_vec[match(uid, sample_id)], levels = levels(y_vec))))
  as.integer(max(2L, min(k_max, min_class)))
}

# CV summary: ROC/Sens/Spec, plus log loss when tuning by logLoss
cv_summary <- if (tune_metric == "logLoss") {
  function(data, lev = NULL, model = NULL) {
    c(twoClassSummary(data, lev, model), mnLogLoss(data, lev, model))
  }
} else {
  twoClassSummary
}

# caret trainControl whose CV folds are grouped by original sample
grouped_train_control <- function(sample_id, y_vec) {
  k <- n_cv_folds(sample_id, y_vec)
  foldid <- grouped_strat_foldid(sample_id, y_vec, k)
  index <- lapply(seq_len(k), function(j) which(foldid != j))
  names(index) <- sprintf("Fold%02d", seq_len(k))
  trainControl(
    method = "cv",
    index = index,
    classProbs = TRUE,
    summaryFunction = cv_summary,
    savePredictions = TRUE
  )
}

# ============================================================
# Helper: single LASSO feature selection on one training set
# ============================================================
# train_idx: row indices of the training set of one split.
# Min-max scaling parameters are estimated from the unique training samples only.
# CV folds for lambda are grouped by sample so duplicates never straddle folds.
lasso_fs_once <- function(X_all, y_all, train_idx, lambda_rule, min_feat) {
  feat_all <- colnames(X_all)
  train_u <- unique(train_idx)
  mins <- apply(X_all[train_u, , drop = FALSE], 2, min)
  rngs <- apply(X_all[train_u, , drop = FALSE], 2, function(x) diff(range(x)))
  keep <- rngs > 0
  Xs <- sweep(sweep(X_all[train_idx, keep, drop = FALSE], 2, mins[keep], "-"), 2, rngs[keep], "/")
  y_tr <- y_all[train_idx]

  coefs <- setNames(rep(0, length(feat_all)), feat_all)
  k <- n_cv_folds(train_idx, y_tr)
  foldid <- grouped_strat_foldid(train_idx, y_tr, k)
  cv_fit <- tryCatch(
    suppressWarnings(cv.glmnet(Xs, y_tr, alpha = 1, family = "binomial", foldid = foldid)),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) {
    return(list(coef = coefs, selected = character(0), fallback = TRUE, failed = TRUE))
  }

  s_name <- if (lambda_rule == "min") "lambda.min" else "lambda.1se"
  beta <- coef(cv_fit, s = s_name)[-1, 1]
  fallback <- FALSE

  if (sum(beta != 0) < min_feat) {
    # Follow the lambda path to the largest lambda with >= min_feat non-zero coefficients
    fallback <- TRUE
    path_beta <- as.matrix(cv_fit$glmnet.fit$beta)
    nz <- colSums(path_beta != 0)
    j <- which(nz >= min_feat)
    if (length(j) > 0L) {
      beta <- path_beta[, j[1]]
    } else {
      last <- path_beta[, ncol(path_beta)]
      top <- names(sort(abs(last), decreasing = TRUE))[seq_len(min(min_feat, length(last)))]
      beta <- setNames(rep(0, length(last)), names(last))
      beta[top] <- ifelse(last[top] == 0, 1e-8, last[top])
    }
  }

  coefs[names(beta)] <- as.numeric(beta)
  list(coef = coefs, selected = names(coefs)[coefs != 0], fallback = fallback, failed = FALSE)
}

# ============================================================
# Helpers: models, ROC, SHAP
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
  as.numeric(probs[, pos, drop = TRUE])
}

fit_caret_model <- function(model_key, train_df, train_control) {
  switch(
    model_key,
    lasso = {
      lasso_grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 50))
      train(class ~ ., data = train_df, method = "glmnet",
            metric = tune_metric, maximize = (tune_metric == "ROC"),
            trControl = train_control, tuneGrid = lasso_grid)
    },
    elastic_net = {
      enet_grid <- expand.grid(alpha = 0.5, lambda = 10^seq(-3, 1, length = 50))
      train(class ~ ., data = train_df, method = "glmnet",
            metric = tune_metric, maximize = (tune_metric == "ROC"),
            trControl = train_control, tuneGrid = enet_grid)
    },
    ridge = {
      ridge_grid <- expand.grid(alpha = 0, lambda = 10^seq(-3, 1, length = 50))
      train(class ~ ., data = train_df, method = "glmnet",
            metric = tune_metric, maximize = (tune_metric == "ROC"),
            trControl = train_control, tuneGrid = ridge_grid)
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
  shap_mat <- matrix(0, nrow = n, ncol = p, dimnames = list(rownames(X), feat_cols))

  # Newer fastshap (>= 0.1.x dev) no longer exports explain_column; use explain() there.
  if (!exists("explain_column", envir = asNamespace("fastshap"), inherits = FALSE)) {
    sh <- fastshap::explain(fit, X = X, pred_wrapper = pred_wrapper, nsim = nsim, newdata = X)
    shap_mat[, ] <- as.matrix(sh)[, feat_cols, drop = FALSE]
    return(shap_mat)
  }

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

# ---- SHAP plots in the style of the Python 'shap' package ----
shap_blue <- "#1E88E5"
shap_cmap <- c("#008AFB", "#7B4BD6", "#FF0052")   # low -> high (shap red/blue map)

shap_top_features <- function(shap_mat, top_n) {
  mean_abs <- sort(colMeans(abs(shap_mat)), decreasing = TRUE)
  mean_abs[seq_len(min(top_n, length(mean_abs)))]
}

shap_theme <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 11, color = "gray15"),
      axis.text.x = element_text(color = "gray25"),
      axis.title.x = element_text(size = 11)
    )
}

# Panel C style: mean(|SHAP|) bar chart
plot_shap_bar <- function(shap_mat, top_n = 20L, title = NULL) {
  imp <- shap_top_features(shap_mat, top_n)
  if (length(imp) == 0L) return(NULL)
  df <- data.frame(Feature = factor(names(imp), levels = rev(names(imp))), value = as.numeric(imp))
  ggplot(df, aes(x = .data$value, y = .data$Feature)) +
    geom_col(fill = shap_blue, width = 0.7) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, y = NULL,
         x = "mean(|SHAP value|) (average impact on model output magnitude)") +
    shap_theme()
}

# Panel D style: beeswarm, colour = feature value scaled within each feature (low -> high)
plot_shap_beeswarm <- function(shap_mat, X_df, top_n = 20L, title = NULL) {
  imp <- shap_top_features(shap_mat, top_n)
  feats <- names(imp)
  if (length(feats) == 0L) return(NULL)

  # Per-feature colour scaling as in shap: clip to the 5th-95th percentile, then 0-1
  scale_feat <- function(v) {
    q <- stats::quantile(v, c(0.05, 0.95), na.rm = TRUE)
    if (!is.finite(q[2] - q[1]) || q[2] <= q[1]) q <- range(v, na.rm = TRUE)
    if (q[2] <= q[1]) return(rep(0.5, length(v)))
    pmin(1, pmax(0, (v - q[1]) / (q[2] - q[1])))
  }

  plot_df <- bind_rows(lapply(feats, function(f) {
    data.frame(
      Sample = rownames(X_df),
      Feature = f,
      SHAP = shap_mat[, f],
      Value_scaled = scale_feat(X_df[[f]]),
      stringsAsFactors = FALSE
    )
  }))
  plot_df$Feature <- factor(plot_df$Feature, levels = rev(feats))
  plot_df <- plot_df[sample.int(nrow(plot_df)), , drop = FALSE]   # random draw order, as in shap

  swarm <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    ggbeeswarm::geom_quasirandom(orientation = "y", width = 0.35, varwidth = FALSE,
                                 size = 1.6, alpha = 0.9, stroke = 0)
  } else {
    geom_jitter(width = 0, height = 0.25, size = 1.6, alpha = 0.9, stroke = 0)
  }

  ggplot(plot_df, aes(x = .data$SHAP, y = .data$Feature, color = .data$Value_scaled)) +
    geom_vline(xintercept = 0, color = "gray45", linewidth = 0.5) +
    swarm +
    scale_color_gradientn(
      colours = shap_cmap, limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
      name = "Feature value",
      guide = guide_colorbar(title.position = "right", barwidth = grid::unit(0.3, "cm"),
                             barheight = grid::unit(12, "lines"), ticks = FALSE)
    ) +
    labs(title = title, y = NULL,
         x = "SHAP value (impact on model output: < 0 toward alpha, > 0 toward beta)") +
    shap_theme() +
    theme(panel.grid.major.y = element_line(color = "gray80", linetype = "dotted", linewidth = 0.4),
          legend.title = element_text(angle = 90, hjust = 0.5, size = 10),
          legend.text = element_text(size = 10))
}

# ============================================================
# MODEL PART (--vis FALSE or omitted)
# ============================================================
if (run_compute) {

  # ============================================================
  # Load data
  # ============================================================
  message("[2/7] Loading and preparing data (", panel_label, ")...")

  wide_list <- lapply(names(panel_paths$npq), function(p) npq_to_wide(panel_paths$npq[[p]], p))
  names(wide_list) <- names(panel_paths$npq)
  for (p in names(wide_list)) {
    message("      ", p, " features: ", ncol(wide_list[[p]]) - 1L,
            ", samples: ", nrow(wide_list[[p]]))
  }

  common_samples <- Reduce(intersect, lapply(wide_list, function(w) w$SampleName))
  if (length(wide_list) > 1L) {
    message("      Samples in all selected panels: ", length(common_samples))
  }

  npq_wide <- Reduce(
    function(a, b) inner_join(a, b, by = "SampleName"),
    lapply(wide_list, function(w) w[w$SampleName %in% common_samples, , drop = FALSE])
  )
  rownames(npq_wide) <- npq_wide$SampleName
  npq_wide$SampleName <- NULL

  feature_panel_map <- data.frame(
    Feature = colnames(npq_wide),
    Panel = sub("^(CNS|Immune)_.*$", "\\1", colnames(npq_wide)),
    Target = sub("^(CNS|Immune)_", "", colnames(npq_wide)),
    stringsAsFactors = FALSE
  )
  write.csv(feature_panel_map, file.path(output_dir, "feature_panel_map.csv"), row.names = FALSE)
  message("      Matrix: ", nrow(npq_wide), " samples x ", ncol(npq_wide), " features")

  all_participants_IDs <- read_excel(metadata_path) %>% filter(toupper(Material) == matrix_npq)

  # k2 labels (alpha vs beta)
  labels <- all_participants_IDs %>%
    select(Tube_ID, k2) %>%
    filter(!is.na(k2)) %>%
    filter(Tube_ID %in% rownames(npq_wide))

  npq_model <- npq_wide[labels$Tube_ID, , drop = FALSE]
  y <- factor(labels$k2, levels = c("alpha", "beta"))

  cat("Model data: n =", nrow(npq_model), ", p =", ncol(npq_model), "\n")
  cat("Outcome distribution:\n")
  print(table(y))

  # Remove columns with NA or zero variance (label-free filter, safe before splitting)
  npq_model <- npq_model[, colSums(is.na(npq_model)) == 0, drop = FALSE]
  npq_model <- npq_model[, apply(npq_model, 2, var) > 0, drop = FALSE]
  cat("Features after filtering:", ncol(npq_model), "\n")

  X <- as.matrix(npq_model)
  sample_names <- rownames(X)
  all_features <- colnames(X)

  full_df <- as.data.frame(X)
  full_df$class <- y
  message("      Data preparation complete.")

  # ============================================================
  # Step 3: generate all stratified random splits up front
  # ============================================================
  message("[3/7] Generating ", BS_number, " stratified random splits (",
          train_frac * 100, "% train / ", (1 - train_frac) * 100, "% validation per class)...")
  splits <- lapply(seq_len(BS_number), function(i) make_split(y, train_frac))

  split_df <- bind_rows(lapply(seq_len(BS_number), function(i) {
    sp <- splits[[i]]
    n_in_train <- tabulate(sp$train, nbins = length(y))
    data.frame(
      bootstrap = i,
      SampleName = sample_names,
      class = as.character(y),
      role = ifelse(n_in_train > 0, "train", "validation"),
      n_times_in_train = n_in_train,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(split_df, file.path(output_dir, "bootstrap_sample_splits.csv"), row.names = FALSE)

  split_summary <- split_df %>%
    group_by(.data$bootstrap) %>%
    summarise(
      n_train_unique = sum(.data$role == "train"),
      n_train_rows = sum(.data$n_times_in_train),
      n_validation = sum(.data$role == "validation"),
      n_validation_alpha = sum(.data$role == "validation" & .data$class == "alpha"),
      n_validation_beta = sum(.data$role == "validation" & .data$class == "beta"),
      .groups = "drop"
    )
  message("      Splits saved: bootstrap_sample_splits.csv")
  message("      Mean unique train samples: ", round(mean(split_summary$n_train_unique), 1),
          ", mean validation samples: ", round(mean(split_summary$n_validation), 1))

  # ============================================================
  # Step 4: per split -> feature selection on train -> 7 models -> validation AUC
  # ============================================================
  models <- c("lasso", "elastic_net", "ridge", "rf", "svm_radial", "svm_linear", "xgb")

  auc_list   <- setNames(lapply(models, function(m) rep(NA_real_, BS_number)), models)
  rocc_store <- setNames(lapply(models, function(m) data.frame()), models)
  fs_records <- vector("list", BS_number)
  pred_records <- vector("list", BS_number)
  n_feat_per_split <- integer(BS_number)
  fallback_per_split <- logical(BS_number)

  message("[4/7] Per-split feature selection + model training (", BS_number, " splits x ",
          "1 LASSO + ", length(models), " models)...")
  message("      Models: ", paste(models, collapse = ", "))
  t_start <- proc.time()

  for (i in seq_len(BS_number)) {
    sp <- splits[[i]]

    # --- (a) One LASSO on this split's training samples only ---
    fs <- lasso_fs_once(X, y, sp$train, lambda_rule = lambda_value, min_feat = min_features)
    feats_i <- fs$selected
    n_feat_per_split[i] <- length(feats_i)
    fallback_per_split[i] <- fs$fallback

    fs_records[[i]] <- data.frame(
      bootstrap = i,
      Feature = names(fs$coef),
      lasso_coef = as.numeric(fs$coef),
      selected = as.integer(names(fs$coef) %in% feats_i),
      stringsAsFactors = FALSE
    )

    message("  [Split ", i, "/", BS_number, "] ", length(feats_i), " features selected",
            if (fs$fallback) paste0(" (fallback: lambda path to >= ", min_features, " features)") else "",
            "; training models...")

    # --- (b) Train models on the training set with selected features ---
    train_df <- full_df[sp$train, c(feats_i, "class"), drop = FALSE]
    X_val    <- full_df[sp$val, feats_i, drop = FALSE]
    y_val    <- full_df$class[sp$val]
    tc_i     <- grouped_train_control(sp$train, y[sp$train])

    pred_i <- list()
    for (m in models) {
      fit_m <- tryCatch(fit_caret_model(m, train_df, tc_i), error = function(e) NULL)
      if (is.null(fit_m)) next

      # --- (c) Evaluate on validation samples ---
      p_m <- tryCatch(get_probs(fit_m, X_val, y_val), error = function(e) NULL)
      if (is.null(p_m) || length(p_m) != length(y_val)) next
      roc_pts <- tryCatch(get_roc_points(y_val, p_m, i), error = function(e) NULL)
      if (is.null(roc_pts)) next

      rocc_store[[m]] <- rbind(rocc_store[[m]], roc_pts[, c("Sp", "Sn", "n", "run")])
      auc_list[[m]][i] <- roc_pts$auc[1]
      pred_i[[m]] <- data.frame(
        bootstrap = i,
        SampleName = sample_names[sp$val],
        class = as.character(y_val),
        model = m,
        prob_beta = p_m,
        stringsAsFactors = FALSE
      )
    }
    pred_records[[i]] <- bind_rows(pred_i)
  }

  t_elapsed <- (proc.time() - t_start)["elapsed"]
  message("      Done. Elapsed: ", round(t_elapsed / 60, 1), " min")

  # ============================================================
  # Step 5: feature selection stability across outer splits
  # ============================================================
  message("[5/7] Summarising feature selection across splits...")

  fs_long <- bind_rows(fs_records)
  write.csv(fs_long, file.path(output_dir, "per_split_feature_selection.csv"), row.names = FALSE)

  per_split_features <- data.frame(
    bootstrap = seq_len(BS_number),
    n_features = n_feat_per_split,
    fallback_used = fallback_per_split,
    features = vapply(fs_records, function(d) paste(d$Feature[d$selected == 1L], collapse = ";"),
                      character(1)),
    stringsAsFactors = FALSE
  ) %>%
    left_join(split_summary, by = "bootstrap")
  write.csv(per_split_features, file.path(output_dir, "per_split_selected_features.csv"), row.names = FALSE)

  feature_freq_df <- fs_long %>%
    group_by(.data$Feature) %>%
    summarise(
      Frequency = mean(.data$selected),               # fraction of outer splits selecting it
      SE = stats::sd(.data$selected) / sqrt(n()),
      Mean_coef_when_selected = if (any(.data$selected == 1L)) mean(.data$lasso_coef[.data$selected == 1L]) else NA_real_,
      Sign_consistency = if (any(.data$selected == 1L)) abs(mean(sign(.data$lasso_coef[.data$selected == 1L]))) else NA_real_,
      .groups = "drop"
    ) %>%
    left_join(feature_panel_map, by = "Feature") %>%
    arrange(desc(.data$Frequency))

  selected_final <- feature_freq_df$Feature[feature_freq_df$Frequency > feature_freq_cutoff]
  if (length(selected_final) < min_features) {
    warning("Fewer than ", min_features, " features exceed --feature_freq_cutoff (", feature_freq_cutoff,
            "); using the top ", min_features, " by outer-split frequency.")
    selected_final <- feature_freq_df$Feature[seq_len(min(min_features, nrow(feature_freq_df)))]
  }
  feature_freq_df$In_final_set <- feature_freq_df$Feature %in% selected_final

  write.csv(feature_freq_df, file.path(output_dir, "feature_selection_stability.csv"), row.names = FALSE)
  writeLines(selected_final, file.path(output_dir, "selected_features.txt"))
  write.csv(feature_panel_map %>% filter(.data$Feature %in% selected_final),
            file.path(output_dir, "feature_panel_map_selected.csv"), row.names = FALSE)

  message("      Features per split: median ", stats::median(n_feat_per_split),
          " (range ", min(n_feat_per_split), "-", max(n_feat_per_split), ")",
          if (any(fallback_per_split)) paste0("; fallback used in ", sum(fallback_per_split), " splits") else "")
  message("      Final features (selected in >", feature_freq_cutoff * 100, "% of splits): ",
          length(selected_final))

  # ============================================================
  # Step 6: AUC summary, per-split tables, SHAP
  # ============================================================
  message("[6/7] Summarising model performance and computing SHAP...")

  auc_per_bs <- data.frame(bootstrap = seq_len(BS_number), n_features = n_feat_per_split)
  for (m in models) auc_per_bs[[m]] <- auc_list[[m]]
  write.csv(auc_per_bs, file.path(output_dir, "bootstrap_AUC_per_iteration.csv"), row.names = FALSE)

  pred_df <- bind_rows(pred_records)
  write.csv(pred_df, file.path(output_dir, "bootstrap_validation_predictions.csv"), row.names = FALSE)

  auc_summary <- data.frame(
    Model = models,
    Mean_AUC = sapply(auc_list, function(x) mean(x, na.rm = TRUE)),
    SD_AUC   = sapply(auc_list, function(x) sd(x, na.rm = TRUE)),
    N_valid  = sapply(auc_list, function(x) sum(!is.na(x))),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(.data$Mean_AUC))

  cat("\n============================================================\n")
  cat("  Model Performance (validation AUC, FS inside each split)\n")
  cat("  Panel(s): ", panel_label, "\n", sep = "")
  cat("============================================================\n")
  print(auc_summary, row.names = FALSE)
  write.csv(auc_summary, file.path(output_dir, "model_AUC_summary.csv"), row.names = FALSE)

  # SHAP: final models on all samples with the final feature set (interpretation only)
  shap_dir <- file.path(output_dir, "shap")
  if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)

  final_df <- full_df[, c(selected_final, "class"), drop = FALSE]
  X_shap <- final_df[, selected_final, drop = FALSE]
  tc_final <- grouped_train_control(seq_len(nrow(final_df)), y)

  # SHAP: final models on all samples with the final feature set (interpretation only)
  shap_dir <- file.path(output_dir, "shap")
  if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)

  final_df <- full_df[, c(selected_final, "class"), drop = FALSE]
  X_shap <- final_df[, selected_final, drop = FALSE]
  tc_final <- grouped_train_control(seq_len(nrow(final_df)), y)
  final_models <- list()
  shap_list <- list()

  for (m in models) {
    message("        SHAP: ", m, " ...")
    fit_final <- tryCatch(fit_caret_model(m, final_df, tc_final), error = function(e) NULL)
    if (is.null(fit_final)) {
      message("          skipped (model fit failed)")
      next
    }
    final_models[[m]] <- fit_final

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
    shap_list[[m]] <- shap_mat

    shap_out <- cbind(SampleName = rownames(X_shap), as.data.frame(shap_mat))
    write.csv(shap_out, file.path(shap_dir, paste0("SHAP_values_", m, ".csv")), row.names = FALSE)

    shap_summary <- data.frame(
      Feature = colnames(shap_mat),
      mean_abs_SHAP = colMeans(abs(shap_mat)),
      mean_SHAP = colMeans(shap_mat),
      stringsAsFactors = FALSE
    ) %>% arrange(desc(.data$mean_abs_SHAP))
    write.csv(shap_summary, file.path(shap_dir, paste0("SHAP_summary_", m, ".csv")), row.names = FALSE)
    message("          done (top feature: ", shap_summary$Feature[1], ")")
  }

  # ------------------------------------------------------------
  # Save everything the plotting step needs
  # ------------------------------------------------------------
  saveRDS(
    list(
      settings = list(
        panel_choice = panel_choice,
        matrix = matrix_meta,
        tune_metric = tune_metric,
        panel_label = panel_label,
        BS_number = BS_number,
        train_frac = train_frac,
        feature_freq_cutoff = feature_freq_cutoff,
        lambda_value = lambda_value,
        min_features = min_features,
        seed = opt$seed,
        shap_nsim = shap_nsim,
        created = Sys.time()
      ),
      models = models,
      auc_list = auc_list,
      auc_summary = auc_summary,
      rocc_store = rocc_store,
      feature_freq_df = feature_freq_df,
      selected_final = selected_final,
      splits = splits,
      shap_list = shap_list,
      X_shap = X_shap
    ),
    results_rds
  )
  saveRDS(final_models, models_rds)
  message("      Results saved: ", basename(results_rds), ", ", basename(models_rds))
}
# ============================================================
# PLOT PART (--vis TRUE or omitted)
# ============================================================
if (run_plot) {
  if (!run_compute) {
    if (!file.exists(results_rds)) {
      stop("No saved results found at ", results_rds,
           "\nRun first with --vis FALSE (same -i/--panel/-o).")
    }
    message("[plots] Loading saved results: ", results_rds)
    res <- readRDS(results_rds)
    st <- res$settings
    for (nm in c("BS_number", "train_frac", "feature_freq_cutoff", "lambda_value")) {
      if (!identical(as.character(st[[nm]]), as.character(get(nm)))) {
        message("      Note: saved ", nm, " = ", st[[nm]], " (command line: ", get(nm),
                "); plots use the saved value.")
      }
    }
    panel_label <- st$panel_label
    BS_number <- st$BS_number
    if (!is.null(st$train_frac)) train_frac <- st$train_frac
    feature_freq_cutoff <- st$feature_freq_cutoff
    lambda_value <- st$lambda_value
    models <- res$models
    auc_list <- res$auc_list
    auc_summary <- res$auc_summary
    rocc_store <- res$rocc_store
    feature_freq_df <- res$feature_freq_df
    shap_list <- res$shap_list
    X_shap <- res$X_shap
  }

  message("[plots] Feature selection stability...")
  plot_freq_df <- feature_freq_df %>% filter(.data$Frequency >= plot_min_freq - 1e-9, .data$Frequency > 0)
  feature_plot <- ggplot(plot_freq_df,
                         aes(x = reorder(.data$Feature, .data$Frequency), y = .data$Frequency,
                             fill = .data$In_final_set)) +
    geom_col() +
    geom_errorbar(aes(ymin = pmax(0, .data$Frequency - .data$SE),
                      ymax = pmin(1, .data$Frequency + .data$SE)),
                  width = 0.2, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = feature_freq_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
    scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "gray70"),
                      name = "Final set") +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = paste0("Feature selection stability across ", BS_number, " train splits (", panel_label, ")"),
         subtitle = paste0("One LASSO (lambda.", lambda_value, ") per train split. ",
                           "Dashed line: cutoff = ", feature_freq_cutoff,
                           ". Showing features with frequency >= ", plot_min_freq,
                           " (", nrow(plot_freq_df), " of ", nrow(feature_freq_df), ")"),
         x = "Feature", y = "Fraction of splits selecting feature") +
    theme_classic()

  ggsave(file.path(output_dir, "feature_selection_stability.pdf"), feature_plot,
         width = 10, height = max(6, nrow(plot_freq_df) / 6))


  message("[plots] SHAP bar + beeswarm...")
  shap_dir <- file.path(output_dir, "shap")
  if (!dir.exists(shap_dir)) dir.create(shap_dir, recursive = TRUE)
  for (m in names(shap_list)) {
    shap_mat <- shap_list[[m]]
    n_show <- min(shap_top_n, ncol(shap_mat))
    h <- max(3.5, n_show * 0.42 + 1.6)
    ttl <- paste0("SHAP Feature Importance: ", m, " (", panel_label, ")")

    p_bar <- plot_shap_bar(shap_mat, top_n = shap_top_n, title = ttl)
    if (!is.null(p_bar)) {
      ggsave(file.path(shap_dir, paste0("SHAP_bar_", m, ".pdf")), p_bar, width = 7, height = h)
      ggsave(file.path(shap_dir, paste0("SHAP_bar_", m, ".png")), p_bar, width = 7, height = h, dpi = 300)
    }
    p_bee <- plot_shap_beeswarm(shap_mat, X_shap, top_n = shap_top_n, title = ttl)
    if (!is.null(p_bee)) {
      ggsave(file.path(shap_dir, paste0("SHAP_beeswarm_", m, ".pdf")), p_bee, width = 7.5, height = h)
      ggsave(file.path(shap_dir, paste0("SHAP_beeswarm_", m, ".png")), p_bee, width = 7.5, height = h, dpi = 300)
    }
  }

  # ============================================================
  # Step 7: plots
  # ============================================================
  message("[plots] ROC and AUC boxplot...")
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
    out$FPR_sd[is.na(out$FPR_sd)] <- 0
    out$TPR_sd[is.na(out$TPR_sd)] <- 0

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

  auc_text <- paste(
    sprintf("%s: %.3f ± %.3f", auc_summary$Model, auc_summary$Mean_AUC, auc_summary$SD_AUC),
    collapse = "\n"
  )

  summary_df$method <- factor(summary_df$method, levels = auc_summary$Model)

  p_roc <- ggplot(summary_df, aes(x = FPR, y = TPR, color = method, fill = method)) +
    geom_ribbon(aes(ymin = pmax(0, TPR - 0.95 * TPR_sd),
                    ymax = pmin(1, TPR + 0.95 * TPR_sd)),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         title = paste0("Multi-Model Comparison: alpha vs beta (", panel_label, " NPQ)"),
         subtitle = paste0("Validation set, ", BS_number, " stratified random splits (",
                           round(train_frac * 100), "/", round((1 - train_frac) * 100),
                           "), feature selection within each train set.\nAUC ± SD:\n",
                           auc_text)) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          plot.subtitle = element_text(size = 9))

  ggsave(file.path(output_dir, "mean_ROC_multimodel.pdf"), p_roc, width = 8, height = 7)
  ggsave(file.path(output_dir, "mean_ROC_multimodel.png"), p_roc, width = 8, height = 7, dpi = 300)

  auc_long <- bind_rows(lapply(models, function(m) {
    data.frame(Model = m, AUC = auc_list[[m]][!is.na(auc_list[[m]])])
  }))
  auc_long$Model <- factor(auc_long$Model, levels = auc_summary$Model)

  p_box <- ggplot(auc_long, aes(x = Model, y = AUC, fill = Model)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste0("Validation AUC across ", BS_number, " splits (", panel_label, ")"),
         x = NULL, y = "AUC (validation)") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(output_dir, "AUC_boxplot_multimodel.pdf"), p_box, width = 7, height = 5)
  ggsave(file.path(output_dir, "AUC_boxplot_multimodel.png"), p_box, width = 7, height = 5, dpi = 300)

}

message("\n=== All done! ===")
message("Output directory: ", normalizePath(output_dir))
if (run_compute) {
  message("Model results:")
  message("  - 08_model_results.rds / 08_final_models.rds (used by --vis TRUE)")
  message("  - bootstrap_sample_splits.csv            (train/validation role per sample per split)")
  message("  - per_split_feature_selection.csv        (LASSO coefficient per feature per split)")
  message("  - per_split_selected_features.csv        (features used in each split)")
  message("  - feature_selection_stability.csv        (fraction of splits selecting each feature)")
  message("  - selected_features.txt                  (final feature set, input for 09)")
  message("  - feature_panel_map.csv / feature_panel_map_selected.csv")
  message("  - model_AUC_summary.csv, bootstrap_AUC_per_iteration.csv")
  message("  - bootstrap_validation_predictions.csv   (validation predicted P(beta))")
  message("  - shap/SHAP_values_<model>.csv, shap/SHAP_summary_<model>.csv")
}
if (run_plot) {
  message("Plots:")
  message("  - feature_selection_stability.pdf")
  message("  - shap/SHAP_bar_<model>.pdf/png, shap/SHAP_beeswarm_<model>.pdf/png")
  message("  - mean_ROC_multimodel.pdf/png, AUC_boxplot_multimodel.pdf/png")
}
