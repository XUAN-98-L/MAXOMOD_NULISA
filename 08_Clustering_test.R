message("[1/5] Loading libraries...")
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
message("      Libraries loaded.")

# ============================================================
# Config
# ============================================================
set.seed(42)
BS_number <- 100
cv_folds  <- 10
output_dir <- "CNS_immune/Results/CNS_panel/Without_tears/08_Clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# Load data
# ============================================================
message("[2/5] Loading and preparing data...")
npq_counts <- read_excel("CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx")
npq_counts <- npq_counts %>%
  filter(SampleMatrixType == "CSF", SampleType == "Sample")

all_participants_IDs <- read_excel("CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx")
all_participants_IDs <- all_participants_IDs %>% filter(Material == "CSF")

# Build wide matrix: samples x proteins
npq_wide <- npq_counts %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = Target, values_from = NPQ) %>%
  as.data.frame()

rownames(npq_wide) <- npq_wide$SampleName
npq_wide$SampleName <- NULL

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

# ============================================================
# Stratified Bootstrap: train on bootstrap, test on OOB
# Models: LASSO, Elastic Net, Ridge, RF, SVM-radial, SVM-linear, XGBoost
# ============================================================
models <- c("lasso", "elastic_net", "ridge", "rf", "svm_radial", "svm_linear", "xgb")

auc_list   <- setNames(lapply(models, function(m) numeric(BS_number)), models)
rocc_store <- setNames(lapply(models, function(m) data.frame()), models)
imp_store  <- setNames(lapply(models, function(m) list()), models)

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

message("[3/5] Running ", BS_number, " bootstrap iterations (7 models each)...")
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

    # Extract variable importance
    vi <- tryCatch({
      imp <- caret::varImp(fit_m, scale = FALSE)$importance
      data.frame(Protein = rownames(imp),
                 Importance = rowMeans(imp, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(vi)) imp_store[[m]][[length(imp_store[[m]]) + 1]] <- vi
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
message("[4/5] Summarizing results...")
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
# Variable Importance: aggregate across bootstraps & plot
# ============================================================
message("      Computing variable importance per model...")
imp_dir <- file.path(output_dir, "variable_importance")
if (!dir.exists(imp_dir)) dir.create(imp_dir, recursive = TRUE)

top_n_imp <- 20  # show top N proteins per model

for (m in models) {
  imp_list_m <- imp_store[[m]]
  if (length(imp_list_m) == 0) {
    message("        ", m, ": no importance data collected, skipping.")
    next
  }

  imp_all <- bind_rows(imp_list_m)
  imp_summary_m <- imp_all %>%
    group_by(Protein) %>%
    summarise(
      mean_importance = mean(Importance, na.rm = TRUE),
      sd_importance   = sd(Importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_importance))

  # Save full table
  write.csv(imp_summary_m, file.path(imp_dir, paste0("VarImp_", m, ".csv")), row.names = FALSE)

  # Plot top N
  plot_df <- head(imp_summary_m, top_n_imp)
  plot_df$Protein <- factor(plot_df$Protein, levels = rev(plot_df$Protein))

  p_imp <- ggplot(plot_df, aes(x = Protein, y = mean_importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = pmax(0, mean_importance - sd_importance),
                      ymax = mean_importance + sd_importance),
                  width = 0.4, color = "black") +
    coord_flip() +
    labs(title = paste0("Top ", top_n_imp, " Proteins - ", m),
         subtitle = paste0("Mean \u00B1 SD importance across ", length(imp_list_m), " bootstrap runs"),
         x = NULL, y = "Importance") +
    theme_classic(base_size = 12)

  ggsave(file.path(imp_dir, paste0("VarImp_", m, ".pdf")), p_imp,
         width = 8, height = max(5, top_n_imp * 0.3 + 1))

  message("        ", m, ": done (", nrow(imp_summary_m), " proteins ranked)")
}

# ============================================================
# Mean ROC plot (like 15_ML_multi_models.R)
# ============================================================
message("[5/5] Generating plots...")
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
       title = "Multi-Model Comparison: alpha vs beta (CSF NPQ)",
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
message("  - model_AUC_summary.csv")
message("  - bootstrap_AUC_per_iteration.csv")
message("  - bootstrap_sample_splits.csv")
message("  - variable_importance/VarImp_<model>.csv  (full ranking)")
message("  - variable_importance/VarImp_<model>.pdf  (top ", top_n_imp, " bar plot)")
message("  - mean_ROC_multimodel.pdf/png")
message("  - AUC_boxplot_multimodel.pdf/png")
