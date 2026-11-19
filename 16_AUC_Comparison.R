#!/usr/bin/env Rscript
# ============================================================
# 16_AUC_Comparison.R
#
# Held-out AUC vs the number of proteins LASSO selected in each split, per fluid.
# Uses bootstrap_AUC_per_iteration.csv from 08_Clustering_Feature_Selection.R:
#   one row per split (n_features = number of proteins LASSO selected on that split's
#   training samples; one AUC column per model, computed on that split's held-out samples).
# Nothing is refitted, so there is no leakage: each AUC is the one 08 already reported.
#
# Note: n_features is chosen by LASSO, not set by us. Splits where LASSO keeps more proteins can
# also be "harder" or "easier" splits, so this is a descriptive association, not a dose-response.
#
# Outputs (-o, default <base>/16_AUC_Comparison[_logLoss])
#   AUC_vs_nfeatures_<Fluid>.pdf            top: mean AUC +/- SD per n_features bin and model
#                                           bottom: number of splits per n_features
#   AUC_vs_nfeatures_all_fluids.pdf         the binned plots of all fluids side by side
#   AUC_by_nfeatures_bin.csv                n splits, mean, SD, median AUC per fluid/model/bin
#   AUC_nfeatures_spearman.csv              Spearman rho and p per fluid/model (BH within fluid)
#
# Part 2 (--topn TRUE): AUC when proteins are added one at a time in order of selection frequency
#   (like MAXOMOD Fig. d). Uses the same 100 train/validation splits as 08
#   (bootstrap_sample_splits.csv) and the same 7 models / tuning / CV folds as 08. For every split
#   and every n = --topn_min..--topn_max the models are trained on the top-n proteins (train
#   samples only) and the AUC is computed on that split's validation samples.
#   --rank global : order = LASSO selection frequency in 08 (feature_selection_stability.csv).
#                   That frequency was computed from all samples, so every validation sample has
#                   already helped to rank the proteins -> AUCs are optimistic (state this).
#   --rank nested : order re-estimated inside each split from --inner_B LASSO fits on random
#                   subsamples (--inner_frac) of that split's TRAINING samples only -> no leakage.
#   Slow (splits x n x 7 models): results are cached in topn_AUC_per_split_<Fluid>_<rank>.csv and
#   reused on the next run unless --topn_recompute TRUE.
#   topn_AUC_<Fluid>_<rank>.pdf             mean AUC +/- SD vs number of top proteins
#   topn_AUC_all_fluids_<rank>.pdf          all fluids side by side
#   topn_AUC_summary_<rank>.csv             mean, SD, median AUC per fluid/model/n
#   topn_protein_order_<Fluid>_<rank>.csv   protein order used (nested: per split + overall)
#
# Model colours: as in 08 (RColorBrewer Set2 in the order of that fluid's model_AUC_summary.csv).
#
# Example (run from NULISA_MAXOMOD/):
# Rscript Script/16_AUC_Comparison.R --metric logLoss --fluids CSF,Plasma,Serum
# Rscript Script/16_AUC_Comparison.R --metric logLoss --topn TRUE --rank global --topn_max 15 --cores 6
# Rscript Script/16_AUC_Comparison.R --metric logLoss --topn TRUE --rank nested --topn_max 15 --cores 6
# ============================================================

suppressMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(parallel)
})

option_list <- list(
  make_option(c("-b", "--base"), type = "character",
              default = "CNS_immune/Results/CNS_panel/Without_tears",
              help = "Folder that contains the 08 output folders [default: %default]"),
  make_option(c("--fluids"), type = "character", default = "CSF,Plasma,Serum",
              help = "Comma-separated fluids: CSF, Plasma, Serum [default: %default]"),
  make_option(c("--metric"), type = "character", default = "logLoss",
              help = "Tuning metric used in 08: 'logLoss' (folders end in _logLoss) or 'ROC' [default: %default]"),
  make_option(c("--bins"), type = "character", default = "2,3,5,8,12,Inf",
              help = paste("Upper limits of the n_features bins, e.g. '2,3,5,8,12,Inf' gives",
                           "2 | 3 | 4-5 | 6-8 | 9-12 | >12 [default: %default]")),
  make_option(c("--min_splits"), type = "integer", default = 3,
              help = "Bins with fewer splits are drawn as points only (no error bar) [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder [default: <base>/16_AUC_Comparison[_logLoss]]"),
  # ---- part 2: top-n proteins ----
  make_option(c("--topn"), type = "logical", default = FALSE,
              help = "Also run the top-n analysis (proteins added in order of frequency) [default: %default]"),
  make_option(c("--rank"), type = "character", default = "global",
              help = "'global' (08 frequency, all samples; optimistic) or 'nested' (per split, train only) [default: %default]"),
  make_option(c("--topn_min"), type = "integer", default = 2,
              help = "Smallest n (caret glmnet needs >= 2 proteins) [default: %default]"),
  make_option(c("--topn_max"), type = "integer", default = 15,
              help = "Largest n [default: %default]"),
  make_option(c("--topn_models"), type = "character",
              default = "lasso,elastic_net,ridge,rf,svm_radial,svm_linear,xgb",
              help = "Models to train [default: %default]"),
  make_option(c("--inner_B"), type = "integer", default = 50,
              help = "--rank nested: LASSO fits per split on subsamples of the train set [default: %default]"),
  make_option(c("--inner_frac"), type = "double", default = 0.8,
              help = "--rank nested: fraction of train samples (per class) in each subsample [default: %default]"),
  make_option(c("--lambda"), type = "character", default = "min",
              help = "--rank nested: cv.glmnet lambda rule, 'min' or '1se' (as 08 -l) [default: %default]"),
  make_option(c("--cores"), type = "integer", default = max(1L, parallel::detectCores() - 1L),
              help = "Splits run in parallel (forked; use 1 on Windows) [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 123,
              help = "Seed; split i uses seed + i [default: %default]"),
  make_option(c("--topn_splits"), type = "integer", default = NA,
              help = "Use only the first N splits (quick test) [default: all]"),
  make_option(c("--topn_recompute"), type = "logical", default = FALSE,
              help = "Recompute even if the cached per-split csv exists [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

tune_metric <- match.arg(opt$metric, c("logLoss", "ROC"))
metric_suffix <- if (tune_metric == "logLoss") "_logLoss" else ""
fluid_name <- function(x) {
  x <- tolower(trimws(x))
  switch(x, csf = "CSF", plasma = "Plasma", serum = "Serum",
         stop("Unknown fluid '", x, "' (use CSF, Plasma or Serum)."))
}
fluids <- vapply(strsplit(opt$fluids, ",")[[1]], fluid_name, character(1), USE.NAMES = FALSE)
dir_08 <- function(fl) {
  leaf <- if (fl == "CSF") "08_Clustering_Feature_Selection" else paste0("08_Clustering_Feature_Selection_", fl)
  file.path(opt$base, paste0(leaf, metric_suffix))
}
output_dir <- if (is.null(opt$output) || !nzchar(opt$output)) {
  file.path(opt$base, paste0("16_AUC_Comparison", metric_suffix))
} else opt$output

# bins
ub <- suppressWarnings(as.numeric(strsplit(opt$bins, ",")[[1]]))
if (anyNA(ub) || is.unsorted(ub, strictly = TRUE)) stop("--bins must be increasing numbers, e.g. '2,3,5,8,12,Inf'.")
if (is.finite(tail(ub, 1))) ub <- c(ub, Inf)
breaks <- c(-Inf, ub)
bin_labels <- vapply(seq_along(ub), function(i) {
  lo <- if (i == 1) NA else floor(ub[i - 1]) + 1
  hi <- ub[i]
  if (i == 1) { if (hi <= 2) as.character(hi) else paste0("<=", hi) }
  else if (!is.finite(hi)) paste0(">", ub[i - 1])
  else if (lo == hi) as.character(hi)
  else paste0(lo, "-", hi)
}, character(1))

model_label <- c(lasso = "LASSO", elastic_net = "Elastic net", ridge = "Ridge", svm_linear = "SVM linear",
                 svm_radial = "SVM radial", rf = "RF", xgb = "XGBoost")
lab_of <- function(m) ifelse(m %in% names(model_label), model_label[m], m)
set2 <- RColorBrewer::brewer.pal(8, "Set2")

message("[1/3] Reading 08 results (", tune_metric, " tuning)")
dat <- list(); cols <- list()
for (fl in fluids) {
  d8 <- dir_08(fl)
  f_auc <- file.path(d8, "bootstrap_AUC_per_iteration.csv")
  if (!file.exists(f_auc)) { message("      Skipped ", fl, ": not found ", f_auc); next }
  a <- read.csv(f_auc, check.names = FALSE)
  if (!all(c("bootstrap", "n_features") %in% names(a))) stop(f_auc, " needs columns bootstrap and n_features.")
  models <- setdiff(names(a), c("bootstrap", "n_features"))
  f_sum <- file.path(d8, "model_AUC_summary.csv")
  col_order <- if (file.exists(f_sum)) as.character(read.csv(f_sum)$Model) else character(0)
  col_order <- c(intersect(col_order, models), setdiff(models, col_order))
  if (length(col_order) > 8L) stop("More than 8 models: Set2 has only 8 colours.")
  cols[[fl]] <- setNames(set2[seq_along(col_order)], lab_of(col_order))
  dat[[fl]] <- a %>%
    pivot_longer(all_of(models), names_to = "model", values_to = "AUC") %>%
    filter(!is.na(.data$AUC)) %>%
    mutate(fluid = fl, model = factor(lab_of(.data$model), levels = lab_of(col_order)),
           bin = cut(.data$n_features, breaks = breaks, labels = bin_labels, right = TRUE))
  message("      ", fl, ": ", nrow(a), " splits, n_features ", min(a$n_features), "-", max(a$n_features),
          " (median ", median(a$n_features), "); ", d8)
}
if (length(dat) == 0L) stop("No 08 results found under ", opt$base)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
all_long <- bind_rows(dat)

message("[2/3] Summaries and Spearman correlations")
bin_sum <- all_long %>%
  group_by(.data$fluid, .data$model, .data$bin) %>%
  summarise(n_splits = n(), mean_AUC = mean(.data$AUC), SD_AUC = sd(.data$AUC),
            median_AUC = median(.data$AUC), .groups = "drop")
write.csv(bin_sum, file.path(output_dir, "AUC_by_nfeatures_bin.csv"), row.names = FALSE)

spear <- all_long %>%
  group_by(.data$fluid, .data$model) %>%
  summarise(n_splits = n(),
            rho = suppressWarnings(cor.test(.data$n_features, .data$AUC, method = "spearman", exact = FALSE))$estimate,
            p = suppressWarnings(cor.test(.data$n_features, .data$AUC, method = "spearman", exact = FALSE))$p.value,
            .groups = "drop") %>%
  group_by(.data$fluid) %>% mutate(p_BH = p.adjust(.data$p, "BH")) %>% ungroup()
write.csv(spear, file.path(output_dir, "AUC_nfeatures_spearman.csv"), row.names = FALSE)
print(as.data.frame(spear %>% mutate(across(c(rho, p, p_BH), ~ signif(.x, 3)))), row.names = FALSE)

message("[3/3] Plots")
binned_plot <- function(fl) {
  bs <- bin_sum %>% filter(.data$fluid == fl, .data$n_splits > 0)
  n_bin <- all_long %>% filter(.data$fluid == fl) %>% distinct(.data$bootstrap, .data$bin) %>%
    count(.data$bin, name = "n")
  xlab_map <- setNames(paste0(n_bin$bin, "\n(n=", n_bin$n, ")"), n_bin$bin)
  bs <- bs %>% mutate(SD_plot = ifelse(.data$n_splits >= opt$min_splits, .data$SD_AUC, NA_real_))
  pd <- position_dodge(width = 0.35)
  ggplot(bs, aes(x = .data$bin, y = .data$mean_AUC, colour = .data$model, group = .data$model)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey55") +
    geom_errorbar(aes(ymin = .data$mean_AUC - .data$SD_plot, ymax = .data$mean_AUC + .data$SD_plot),
                  width = 0.25, position = pd, na.rm = TRUE) +
    geom_line(linewidth = 0.8, position = pd) +
    geom_point(size = 2.2, position = pd) +
    scale_colour_manual(values = cols[[fl]], name = NULL) +
    scale_x_discrete(labels = xlab_map, drop = TRUE) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste0(fl, ": held-out AUC by number of LASSO-selected proteins"),
         subtitle = paste0("Mean +/- SD over the splits in each bin (", tune_metric,
                           " tuning); no error bar if < ", opt$min_splits, " splits"),
         x = "Proteins selected by LASSO in the split (n = number of splits)",
         y = "Mean AUC +/- SD") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
}
count_plot <- function(fl) {
  cnt <- dat[[fl]] %>% distinct(.data$bootstrap, .data$n_features) %>% count(.data$n_features)
  ggplot(cnt, aes(x = .data$n_features, y = .data$n)) +
    geom_col(fill = "grey60", width = 0.8) +
    geom_text(aes(label = .data$n), vjust = -0.3, size = 3) +
    scale_x_continuous(breaks = seq(min(cnt$n_features), max(cnt$n_features), by = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Proteins selected by LASSO in the split", y = "Number of splits",
         title = paste0(fl, ": proteins selected per split (", sum(cnt$n), " splits)")) +
    theme_classic(base_size = 11)
}

binned <- list()
for (fl in names(dat)) {
  pb <- binned_plot(fl); binned[[fl]] <- pb
  g <- arrangeGrob(pb, count_plot(fl), ncol = 1, heights = c(2.2, 1))
  ggsave(file.path(output_dir, paste0("AUC_vs_nfeatures_", fl, ".pdf")), g, width = 8.5, height = 8)
}
if (length(binned) > 1L) {
  g_all <- arrangeGrob(grobs = lapply(binned, function(p) p + labs(subtitle = NULL) +
                                        theme(legend.position = "bottom") +
                                        guides(colour = guide_legend(nrow = 2))),
                       nrow = 1)
  ggsave(file.path(output_dir, "AUC_vs_nfeatures_all_fluids.pdf"), g_all,
         width = 6.5 * length(binned), height = 5.5)
}

# ============================================================
# Part 2: top-n proteins (in order of selection frequency), like MAXOMOD Fig. d
# ============================================================
if (isTRUE(opt$topn)) {
  suppressMessages({
    library(readxl); library(caret); library(glmnet); library(pROC)
    library(randomForest); library(kernlab); library(xgboost)
  })
  rank_mode <- match.arg(opt$rank, c("global", "nested"))
  lambda_rule <- match.arg(opt$lambda, c("min", "1se"))
  topn_models <- trimws(strsplit(opt$topn_models, ",")[[1]])
  if (opt$topn_min < 2L) stop("--topn_min must be >= 2 (caret glmnet needs >= 2 predictors).")
  init_dir <- file.path(opt$base, "00_Initialization")
  npq_file <- file.path(init_dir, "npq_counts.xlsx")
  meta_file <- file.path(init_dir, "all_participants_IDs.xlsx")
  if (!file.exists(npq_file) || !file.exists(meta_file)) stop("Need ", npq_file, " and ", meta_file)
  message("[topn] Rank: ", rank_mode,
          if (rank_mode == "global") " (08 frequency from all samples -> AUC optimistic)"
          else paste0(" (", opt$inner_B, " LASSO fits on ", opt$inner_frac * 100,
                      "% subsamples of each split's train set)"),
          "; n = ", opt$topn_min, "-", opt$topn_max, "; models: ", paste(topn_models, collapse = ", "),
          "; cores: ", opt$cores)

  # ---- helpers copied from 08 (same CV folds, tuning, models) ----
  cv_folds <- 10L
  grouped_strat_foldid <- function(sample_id, y_vec, k) {
    uid <- unique(sample_id); ycl <- y_vec[match(uid, sample_id)]
    fold <- integer(length(uid))
    for (cl in unique(as.character(ycl))) { w <- which(ycl == cl); fold[w] <- sample(rep_len(seq_len(k), length(w))) }
    fold[match(sample_id, uid)]
  }
  n_cv_folds <- function(sample_id, y_vec, k_max = cv_folds) {
    uid <- unique(sample_id)
    min_class <- min(table(factor(y_vec[match(uid, sample_id)], levels = levels(y_vec))))
    as.integer(max(2L, min(k_max, min_class)))
  }
  cv_summary <- if (tune_metric == "logLoss") {
    function(data, lev = NULL, model = NULL) c(twoClassSummary(data, lev, model), mnLogLoss(data, lev, model))
  } else twoClassSummary
  grouped_train_control <- function(sample_id, y_vec) {
    k <- n_cv_folds(sample_id, y_vec)
    foldid <- grouped_strat_foldid(sample_id, y_vec, k)
    index <- lapply(seq_len(k), function(j) which(foldid != j)); names(index) <- sprintf("Fold%02d", seq_len(k))
    trainControl(method = "cv", index = index, classProbs = TRUE, summaryFunction = cv_summary,
                 savePredictions = FALSE)
  }
  fit_caret_model <- function(model_key, train_df, train_control) {
    mx <- (tune_metric == "ROC")
    grid <- function(a) expand.grid(alpha = a, lambda = 10^seq(-3, 1, length = 50))
    switch(model_key,
      lasso       = train(class ~ ., data = train_df, method = "glmnet", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneGrid = grid(1)),
      elastic_net = train(class ~ ., data = train_df, method = "glmnet", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneGrid = grid(0.5)),
      ridge       = train(class ~ ., data = train_df, method = "glmnet", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneGrid = grid(0)),
      rf          = train(class ~ ., data = train_df, method = "rf", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneLength = 3, ntree = 500),
      svm_radial  = train(class ~ ., data = train_df, method = "svmRadial", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneLength = 3),
      svm_linear  = train(class ~ ., data = train_df, method = "svmLinear", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneLength = 3),
      xgb         = train(class ~ ., data = train_df, method = "xgbTree", metric = tune_metric, maximize = mx,
                          trControl = train_control, tuneLength = 3, verbosity = 0, nthread = 1),
      stop("Unknown model: ", model_key))
  }
  # LASSO (min-max scaling from the samples used) -> selected features
  lasso_select <- function(X, y, idx) {
    mins <- apply(X[idx, , drop = FALSE], 2, min)
    rngs <- apply(X[idx, , drop = FALSE], 2, function(v) diff(range(v)))
    keep <- rngs > 0
    Xs <- sweep(sweep(X[idx, keep, drop = FALSE], 2, mins[keep], "-"), 2, rngs[keep], "/")
    yy <- y[idx]
    foldid <- grouped_strat_foldid(idx, yy, n_cv_folds(idx, yy))
    fit <- tryCatch(suppressWarnings(cv.glmnet(Xs, yy, alpha = 1, family = "binomial", foldid = foldid)),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b <- coef(fit, s = if (lambda_rule == "min") "lambda.min" else "lambda.1se")[-1, 1]
    out <- setNames(rep(0, ncol(X)), colnames(X)); out[names(b)] <- abs(b); out
  }
  strat_subsample <- function(y, idx, frac) {
    sort(unlist(lapply(levels(y), function(cl) {
      w <- idx[y[idx] == cl]; w[sample.int(length(w), max(2L, round(length(w) * frac)))]
    })))
  }

  load_fluid <- function(fl) {
    key <- toupper(fl)
    w <- read_excel(npq_file) %>%
      filter(toupper(.data$SampleMatrixType) == key, .data$SampleType == "Sample") %>%
      mutate(Feature = paste0("CNS_", .data$Target)) %>%
      select("SampleName", "Feature", "NPQ") %>%
      pivot_wider(names_from = "Feature", values_from = "NPQ") %>% as.data.frame()
    rownames(w) <- w$SampleName; w$SampleName <- NULL
    lab <- read_excel(meta_file) %>% filter(toupper(.data$Material) == key, !is.na(.data$k2),
                                            .data$Tube_ID %in% rownames(w))
    X <- as.matrix(w[lab$Tube_ID, , drop = FALSE])
    X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
    X <- X[, apply(X, 2, var) > 0, drop = FALSE]
    list(X = X, y = factor(lab$k2, levels = c("alpha", "beta")))
  }

  topn_res <- list(); topn_labels <- list(); panel_n <- list()
  for (fl in fluids) {
    d8 <- dir_08(fl)
    f_split <- file.path(d8, "bootstrap_sample_splits.csv")
    f_stab <- file.path(d8, "feature_selection_stability.csv")
    if (!file.exists(f_split) || !file.exists(f_stab)) { message("[topn] Skipped ", fl, ": 08 files missing"); next }
    stab <- read.csv(f_stab, stringsAsFactors = FALSE) %>% arrange(desc(.data$Frequency))
    f_sel <- file.path(d8, "selected_features.txt")
    panel_n[[fl]] <- if (file.exists(f_sel)) length(readLines(f_sel)) else NA_integer_
    cache <- file.path(output_dir, paste0("topn_AUC_per_split_", fl, "_", rank_mode, ".csv"))
    f_order <- file.path(output_dir, paste0("topn_protein_order_", fl, "_", rank_mode, ".csv"))

    if (file.exists(cache) && !isTRUE(opt$topn_recompute)) {
      message("[topn] ", fl, ": using cached ", cache)
      res <- read.csv(cache, stringsAsFactors = FALSE)
    } else {
      dd <- load_fluid(fl); X <- dd$X; y <- dd$y
      sp <- read.csv(f_split, colClasses = c(SampleName = "character"))
      if (!setequal(unique(sp$SampleName), rownames(X)))
        stop(fl, ": samples in ", f_split, " differ from the NPQ/metadata samples.")
      n_max <- if (rank_mode == "global") min(opt$topn_max, sum(stab$Frequency > 0), ncol(X)) else min(opt$topn_max, ncol(X))
      if (n_max < opt$topn_min) stop(fl, ": fewer than ", opt$topn_min, " proteins available.")
      ns <- seq(opt$topn_min, n_max)
      global_order <- intersect(stab$Feature, colnames(X))
      splits <- sort(unique(sp$bootstrap))
      if (!is.na(opt$topn_splits)) splits <- head(splits, opt$topn_splits)
      message("[topn] ", fl, ": ", nrow(X), " samples, ", length(splits), " splits, n = ", min(ns), "-", max(ns),
              " (", length(splits) * length(ns) * length(topn_models), " model fits)")
      t0 <- proc.time()
      one_split <- function(i) {
        set.seed(opt$seed + i)
        s <- sp[sp$bootstrap == i, ]
        tr <- match(s$SampleName[s$role == "train"], rownames(X))
        va <- match(s$SampleName[s$role == "validation"], rownames(X))
        ord <- if (rank_mode == "global") global_order else {
          sel <- matrix(0, 0, ncol(X)); coefsum <- setNames(rep(0, ncol(X)), colnames(X))
          cnt <- setNames(rep(0, ncol(X)), colnames(X)); nb <- 0
          for (b in seq_len(opt$inner_B)) {
            cb <- lasso_select(X, y, strat_subsample(y, tr, opt$inner_frac))
            if (is.null(cb)) next
            nb <- nb + 1; cnt <- cnt + (cb > 0); coefsum <- coefsum + cb
          }
          freq <- cnt / max(nb, 1)
          names(freq)[order(-freq, -coefsum)]
        }
        tc <- grouped_train_control(tr, y[tr])   # same CV folds for every n in this split
        out <- list()
        for (n in ns) {
          f <- ord[seq_len(n)]
          trdf <- data.frame(X[tr, f, drop = FALSE], class = y[tr], check.names = FALSE)
          vadf <- data.frame(X[va, f, drop = FALSE], check.names = FALSE)
          for (m in topn_models) {
            fit <- tryCatch(suppressWarnings(fit_caret_model(m, trdf, tc)), error = function(e) NULL)
            auc <- NA_real_
            if (!is.null(fit)) {
              p <- tryCatch(predict(fit, newdata = vadf, type = "prob")[, "beta"], error = function(e) NULL)
              if (!is.null(p)) auc <- as.numeric(pROC::auc(pROC::roc(y[va], p, levels = c("alpha", "beta"),
                                                                     direction = "<", quiet = TRUE)))
            }
            out[[length(out) + 1]] <- data.frame(bootstrap = i, n_proteins = n, model = m, AUC = auc,
                                                 added = f[n], stringsAsFactors = FALSE)
          }
        }
        list(auc = do.call(rbind, out),
             order = data.frame(bootstrap = i, rank = seq_len(max(ns)),
                                Feature = ord[seq_len(max(ns))], stringsAsFactors = FALSE))
      }
      rl <- if (opt$cores > 1L) mclapply(splits, one_split, mc.cores = opt$cores, mc.preschedule = FALSE)
            else lapply(splits, one_split)
      bad <- vapply(rl, function(r) inherits(r, "try-error") || is.null(r), logical(1))
      if (any(bad)) stop(fl, ": ", sum(bad), " splits failed: ", paste(unlist(rl[bad])[1], collapse = " "))
      res <- bind_rows(lapply(rl, `[[`, "auc"))
      write.csv(res, cache, row.names = FALSE)
      ord_df <- bind_rows(lapply(rl, `[[`, "order"))
      if (rank_mode == "global") ord_df <- ord_df %>% filter(.data$bootstrap == min(.data$bootstrap)) %>% select(-"bootstrap")
      write.csv(ord_df, f_order, row.names = FALSE)
      message("[topn] ", fl, ": done in ", round((proc.time() - t0)["elapsed"] / 60, 1), " min; ",
              sum(is.na(res$AUC)), " failed fits")
    }
    res$fluid <- fl
    topn_res[[fl]] <- res
    # x-axis labels: global -> "n\n+PROTEIN" (the protein added at that step)
    topn_labels[[fl]] <- if (rank_mode == "global") {
      a <- res %>% distinct(.data$n_proteins, .data$added) %>% arrange(.data$n_proteins)
      ord_g <- read.csv(f_order, stringsAsFactors = FALSE) %>% arrange(.data$rank)
      short <- function(v) sub("^(CNS|Immune)_", "", v)
      first <- paste(short(ord_g$Feature[seq_len(min(a$n_proteins))]), collapse = "\n+")
      setNames(ifelse(a$n_proteins == min(a$n_proteins), paste0(a$n_proteins, "\n", first),
                      paste0(a$n_proteins, "\n+", short(a$added))), a$n_proteins)
    } else NULL
  }

  if (length(topn_res) > 0L) {
    topn_sum <- bind_rows(topn_res) %>%
      filter(!is.na(.data$AUC)) %>%
      group_by(.data$fluid, .data$model, .data$n_proteins) %>%
      summarise(n_splits = n(), mean_AUC = mean(.data$AUC), SD_AUC = sd(.data$AUC),
                median_AUC = median(.data$AUC), .groups = "drop") %>%
      mutate(rank = rank_mode)
    write.csv(topn_sum, file.path(output_dir, paste0("topn_AUC_summary_", rank_mode, ".csv")), row.names = FALSE)

    topn_plot <- function(fl) {
      d <- topn_sum %>% filter(.data$fluid == fl) %>%
        mutate(model = factor(lab_of(.data$model), levels = names(cols[[fl]])))
      pn <- panel_n[[fl]]
      pd <- position_dodge(width = 0.4)
      p <- ggplot(d, aes(x = .data$n_proteins, y = .data$mean_AUC, colour = .data$model, group = .data$model)) +
        geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey55")
      if (!is.na(pn) && pn >= min(d$n_proteins) && pn <= max(d$n_proteins)) {
        p <- p + geom_vline(xintercept = pn, linetype = "dashed", colour = "red") +
          annotate("text", x = pn, y = 0.02, label = paste0(pn, " proteins (08 panel)"),
                   hjust = -0.05, size = 3, colour = "red")
      }
      brk <- sort(unique(d$n_proteins))
      lab <- if (!is.null(topn_labels[[fl]])) topn_labels[[fl]][as.character(brk)] else brk
      p + geom_errorbar(aes(ymin = .data$mean_AUC - .data$SD_AUC, ymax = .data$mean_AUC + .data$SD_AUC),
                        width = 0.3, position = pd) +
        geom_line(linewidth = 0.8, position = pd) +
        geom_point(size = 2, position = pd) +
        scale_colour_manual(values = cols[[fl]], name = NULL, drop = TRUE) +
        scale_x_continuous(breaks = brk, labels = lab) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(title = paste0(fl, ": held-out AUC with the top-n proteins"),
             subtitle = if (rank_mode == "global")
               "Order = 08 LASSO selection frequency (all samples; optimistic). Mean +/- SD over splits"
             else "Order re-estimated within each split's training samples (nested). Mean +/- SD over splits",
             x = if (rank_mode == "global") "Number of top proteins (protein added at this step)"
                 else "Number of top proteins (order chosen within each split)",
             y = "Mean AUC +/- SD") +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold"),
              axis.text.x = element_text(size = if (rank_mode == "global") 7 else 9))
    }
    tp <- list()
    for (fl in names(topn_res)) {
      tp[[fl]] <- topn_plot(fl)
      ggsave(file.path(output_dir, paste0("topn_AUC_", fl, "_", rank_mode, ".pdf")), tp[[fl]],
             width = 9, height = 5.5)
    }
    if (length(tp) > 1L) {
      g_all <- arrangeGrob(grobs = lapply(tp, function(p) p + labs(subtitle = NULL) +
                                            theme(legend.position = "bottom") +
                                            guides(colour = guide_legend(nrow = 2))), nrow = 1)
      ggsave(file.path(output_dir, paste0("topn_AUC_all_fluids_", rank_mode, ".pdf")), g_all,
             width = 7 * length(tp), height = 5.5)
    }
  }
}
message("Done. Output: ", output_dir)
