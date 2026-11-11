#!/usr/bin/env Rscript
# ============================================================
# 14_ROCK_ALS_allAgree.R
#
# Confidence filter for the ROCK-ALS subtype predictions of 11_with_ROCK_ALS.R:
# keep only the samples on which ALL ML models agree (every model calls the same
# subtype at the 11 threshold), then compare clinical variables between the predicted
# alpha and beta groups with the same tests and plots as 11.
#
# The 11 tertile filter (--Tertile TRUE, files *_tertile*) is NOT used: all 69 ROCK-ALS samples
# start from the plain ensemble call, and the only filter is "all models agree".
#
# Input (--input): an 11_with_ROCK_ALS output folder, read-only (nothing is written there). Used files:
#   ROCK_ALS_CSF_predicted_k2_<mode>.csv          per-sample ensemble call, per-model calls (k2_<model>)
#   PCA_predicted_k2_<mode>/pca_all_proteins_scores.csv   raw NULISA NEFL NPQ per sample
# Model colours (--model_dir): same as 08_Clustering_Feature_Selection.R, i.e. RColorBrewer "Set2"
#   assigned to the models in the order of 08's model_AUC_summary.csv (highest mean AUC first).
#     (fallback: NfL_vs_NEFL/NfL_vs_NEFL_CSF_data.csv)
#
# Tests (same as 11):
#   categorical: Fisher's exact test (predicted k2 x category, missing excluded)
#   continuous : both groups Shapiro p > 0.05 -> t-test (Student's if Levene p > 0.05, else
#                Welch's); otherwise, or n < 3 in a group, Wilcoxon rank-sum
#
# Outputs go to their own folder (default:
#   .../CNS_panel/Without_tears/14_ROCK_ALS_allAgree_logLoss, next to the 11 folder, not inside it):
#   allAgree_samples_<mode>.csv                       all samples: votes, all_agree flag
#   model_votes_<mode>.pdf                            ensemble P(beta), kept vs excluded
#   disagree_samples_dotplot_<mode>.pdf/.png          samples on which the models disagree: one row
#                                                    per sample, one point per model's P(beta)
#                                                    (circle = linear, triangle = non-linear model),
#                                                    black bar = ensemble P(beta), dashed = 0.5
#   disagree_samples_<mode>.csv                       the same samples: per-model P(beta) and call
#   predicted_k2_by_treatment_<mode>_allAgree.csv/.pdf
#   clinical_by_predicted_k2_<mode>_allAgree.csv/.pdf
#   clinical_continuous_by_predicted_k2_<mode>_allAgree/  tests csv + violin plots
#   comparison_all_vs_allAgree_<mode>.csv            p-values: all samples vs all-agree samples
#
# Example (run from NULISA_MAXOMOD/):
# Rscript Script/14_ROCK_ALS_allAgree.R
# Rscript Script/14_ROCK_ALS_allAgree.R \
#   -i CNS_immune/Results/CNS_panel/Without_tears/11_with_ROCK_ALS_logLoss \
#   -o CNS_immune/Results/CNS_panel/Without_tears/14_ROCK_ALS_allAgree_logLoss --mode standardize
# ============================================================

suppressMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/CNS_immune/Results/CNS_panel/Without_tears/11_with_ROCK_ALS_logLoss",
              help = "11_with_ROCK_ALS output folder [default: %default]"),
  make_option(c("-o", "--output"), type = "character",
              default = "/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/CNS_immune/Results/CNS_panel/Without_tears/14_ROCK_ALS_allAgree_logLoss",
              help = "Output folder, separate from the 11 folder [default: %default]"),
  make_option(c("--model_dir"), type = "character",
              default = "/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/CNS_immune/Results/CNS_panel/Without_tears/08_Clustering_Feature_Selection_logLoss",
              help = paste("08 output folder whose model_AUC_summary.csv sets the model colours, so they match",
                           "the 08 plots (Set2 in order of mean AUC) [default: %default]")),
  make_option(c("--mode"), type = "character", default = "standardize",
              help = "Harmonisation mode of the 11 predictions: standardize, center or none [default: %default]"),
  make_option(c("--pr_cutoff"), type = "numeric", default = 0.58,
              help = "Progression-rate cutoff used in 11 for fast/slow (labels only) [default: %default]"),
  make_option(c("--cont_vars"), type = "character",
              default = "age,ALSFRS-R Score,Calculated Vital capacity (%),NfL,NEFL_NPQ,GFAP,progression_rate,months_onset_to_BL",
              help = "Comma-separated continuous variables [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

mode <- tolower(trimws(opt$mode))
if (!mode %in% c("standardize", "center", "none")) stop("--mode must be standardize, center or none.")
in_dir <- opt$input
pred_file <- file.path(in_dir, paste0("ROCK_ALS_CSF_predicted_k2_", mode, ".csv"))
if (!file.exists(pred_file)) stop("Not found: ", pred_file, "\nRun 11_with_ROCK_ALS.R first (or check --input/--mode).")

output_dir <- opt$output
if (is.null(output_dir) || !nzchar(output_dir)) stop("-o/--output must not be empty.")
if (normalizePath(output_dir, mustWork = FALSE) == normalizePath(in_dir)) {
  stop("--output must be a separate folder, not the 11 input folder.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("      Output: ", output_dir)
pr_cutoff <- opt$pr_cutoff
k2_cols <- c(alpha = "#FC8D62", beta = "#8DA0CB")
# Shared colour palette (model colours in the dot plot follow 08 instead; see --model_dir)
my_palette <- c("#66C2A5", "#E78AC3", "#ffc773",
                "#1B9E77", "#7570B3", "#FEE090",
                "#A6D854", "#8DA0CB", "#FC8D62",
                "#8FBC94", "#b0a4e3", "#ffa631",
                "#0aa344", "#e4c6d0", "#ffa400",
                "#519a73", "#4b5cc4", "#eedeb0",
                "#549688", "#ffb3a7", "#b35c44",
                "#7fecad", "#a1afc9", "#a78e44",
                "#519a73", "#2e4e7e", "#955539")
sfx_out <- "_allAgree"

message("[1/6] Reading ", pred_file)
d <- read.csv(pred_file, check.names = FALSE, stringsAsFactors = FALSE)

# ---- per-model calls -> all-agree flag ----
model_cols <- grep("^k2_", names(d), value = TRUE)
if (length(model_cols) == 0L) stop("No per-model call columns (k2_<model>) in ", basename(pred_file))
models <- sub("^k2_", "", model_cols)
calls <- as.matrix(d[, model_cols, drop = FALSE])
n_mod <- rowSums(!is.na(calls))
n_beta <- rowSums(calls == "beta", na.rm = TRUE)
d$n_models_used <- n_mod
d$n_vote_beta_models <- n_beta
d$all_agree <- n_mod == length(models) & (n_beta == 0 | n_beta == n_mod)
# the unanimous call must match the ensemble call (sanity check)
unan_call <- ifelse(n_beta == n_mod, "beta", "alpha")
bad <- d$all_agree & unan_call != d$ensemble_k2
if (any(bad)) warning(sum(bad), " all-agree samples differ from the ensemble call: ", paste(d$PatNr[bad], collapse = ", "))

message("      Models (", length(models), "): ", paste(models, collapse = ", "))
message("      All samples: ", nrow(d), " (", sum(d$ensemble_k2 == "alpha"), " alpha / ",
        sum(d$ensemble_k2 == "beta"), " beta)")
message("      All models agree: ", sum(d$all_agree), " (", sum(d$all_agree & d$ensemble_k2 == "alpha"), " alpha / ",
        sum(d$all_agree & d$ensemble_k2 == "beta"), " beta); excluded: ", sum(!d$all_agree))
cat("Votes for beta per sample (of ", length(models), " models):\n", sep = "")
print(table(n_vote_beta = factor(n_beta, levels = 0:length(models))))

write.csv(d %>% select(PatNr, ensemble_P_beta, ensemble_k2, n_models_used, n_vote_beta_models, all_agree,
                       all_of(model_cols), starts_with("P_beta_")) %>% arrange(ensemble_P_beta),
          file.path(output_dir, paste0("allAgree_samples_", mode, ".csv")), row.names = FALSE)

# ---- NEFL NPQ (raw NULISA, from the 11 outputs) ----
nefl_file <- file.path(in_dir, paste0("PCA_predicted_k2_", mode), "pca_all_proteins_scores.csv")
nefl_alt <- file.path(in_dir, "NfL_vs_NEFL", "NfL_vs_NEFL_CSF_data.csv")
if (file.exists(nefl_file)) {
  s <- read.csv(nefl_file, check.names = FALSE)
  if ("NULISA NEFL NPQ" %in% names(s)) d$NEFL_NPQ <- s[["NULISA NEFL NPQ"]][match(d$PatNr, s$PatNr)]
}
if (!"NEFL_NPQ" %in% names(d) && file.exists(nefl_alt)) {
  s <- read.csv(nefl_alt, check.names = FALSE)
  if ("NEFL_NPQ" %in% names(s)) d$NEFL_NPQ <- s$NEFL_NPQ[match(d$PatNr, s$PatNr)]
}
if (!"NEFL_NPQ" %in% names(d)) message("Note: NULISA NEFL NPQ not found in the 11 outputs; NEFL_NPQ skipped.")

# ---- figure: ensemble P(beta), kept vs excluded ----
message("[2/6] Vote plot")
d$agree_class <- factor(ifelse(d$all_agree, d$ensemble_k2, paste0("excluded (", d$ensemble_k2, ")")),
                        levels = c("alpha", "excluded (alpha)", "excluded (beta)", "beta"))
p_votes <- ggplot(d, aes(x = .data$ensemble_P_beta, fill = .data$agree_class)) +
  geom_histogram(binwidth = 0.025, boundary = 0, color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray55") +
  scale_fill_manual(values = c(alpha = "#FC8D62", "excluded (alpha)" = "#FDD0B8",
                               "excluded (beta)" = "#D3DBEC", beta = "#8DA0CB"), drop = FALSE, name = "Class") +
  labs(title = "ROCK-ALS CSF: samples on which all models agree",
       subtitle = paste0("harmonize = ", mode, "; kept if all ", length(models), " models give the same call (",
                         sum(d$all_agree & d$ensemble_k2 == "alpha"), " alpha, ",
                         sum(d$all_agree & d$ensemble_k2 == "beta"), " beta)\n",
                         sum(!d$all_agree), " samples excluded (models disagree); dotted = 0.5"),
       x = "Ensemble P(beta)", y = "Samples") +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, paste0("model_votes_", mode, ".pdf")), p_votes, width = 7.5, height = 4.5)

# ---- dot plot: samples on which the models disagree ----
message("[3/6] Dot plot of the samples on which the models disagree")
linear_models <- c("lasso", "elastic_net", "ridge", "svm_linear")
model_label <- c(lasso = "LASSO", elastic_net = "Elastic net", ridge = "Ridge", svm_linear = "SVM linear",
                 svm_radial = "SVM radial", rf = "RF", xgb = "XGBoost")
model_order <- c(intersect(c(linear_models, "svm_radial", "rf", "xgb"), models),
                 setdiff(models, c(linear_models, "svm_radial", "rf", "xgb")))   # linear first
lab_of <- function(m) ifelse(m %in% names(model_label), model_label[m], m)
# model colours: as in 08_Clustering_Feature_Selection.R (RColorBrewer Set2 in the order of
# 08's model_AUC_summary.csv, highest mean AUC first); models missing there are appended
auc_file <- file.path(opt$model_dir, "model_AUC_summary.csv")
col_order <- if (file.exists(auc_file)) {
  as.character(read.csv(auc_file, stringsAsFactors = FALSE)$Model)
} else {
  message("      Note: ", auc_file, " not found; model colours follow the model order instead of 08.")
  character(0)
}
col_order <- c(intersect(col_order, model_order), setdiff(model_order, col_order))
if (length(col_order) > 8L) stop("More than 8 models: RColorBrewer Set2 has only 8 colours.")
set2 <- RColorBrewer::brewer.pal(8, "Set2")
dot_cols <- setNames(set2[seq_along(col_order)], lab_of(col_order))
message("      Model colours (as in 08): ", paste0(names(dot_cols), " ", dot_cols, collapse = ", "))

dis <- d[!d$all_agree, , drop = FALSE]
dis <- dis[order(dis$ensemble_P_beta), , drop = FALSE]
if (nrow(dis) > 0L) {
  write.csv(dis %>% select(PatNr, ensemble_P_beta, ensemble_k2, n_vote_beta_models,
                           all_of(paste0("P_beta_", model_order)), all_of(paste0("k2_", model_order))),
            file.path(output_dir, paste0("disagree_samples_", mode, ".csv")), row.names = FALSE)
  row_lv <- dis$PatNr
  dot_long <- bind_rows(lapply(model_order, function(m) {
    data.frame(PatNr = dis$PatNr, model = lab_of(m),
               type = ifelse(m %in% linear_models, "linear", "non-linear"),
               P = dis[[paste0("P_beta_", m)]], stringsAsFactors = FALSE)
  })) %>% filter(!is.na(P))
  dot_long$PatNr <- factor(dot_long$PatNr, levels = row_lv)
  dot_long$model <- factor(dot_long$model, levels = lab_of(model_order))
  ens <- data.frame(PatNr = factor(dis$PatNr, levels = row_lv), P = dis$ensemble_P_beta,
                    k2 = dis$ensemble_k2, votes = paste0(dis$n_vote_beta_models, "/", dis$n_models_used, " beta"))

  p_dot <- ggplot(dot_long, aes(x = .data$P, y = .data$PatNr)) +
    annotate("rect", xmin = -Inf, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = k2_cols[["alpha"]], alpha = 0.08) +
    annotate("rect", xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = k2_cols[["beta"]], alpha = 0.10) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40") +
    geom_line(aes(group = .data$PatNr), color = "gray75") +
    geom_point(aes(color = .data$model, shape = .data$type), size = 2.8, alpha = 0.9) +
    geom_point(data = ens, aes(x = .data$P, y = .data$PatNr), inherit.aes = FALSE,
               shape = 124, size = 6, color = "black") +
    geom_text(data = ens, aes(x = 1.03, y = .data$PatNr, label = .data$votes, color = NULL), inherit.aes = FALSE,
              hjust = 0, size = 3.2, color = k2_cols[ens$k2]) +
    scale_color_manual(values = dot_cols, name = "Model") +
    scale_shape_manual(values = c(linear = 16, "non-linear" = 17), name = "Model type") +
    scale_x_continuous(limits = c(0, 1.17), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0))) +
    labs(title = "ROCK-ALS CSF: samples where models disagree",
         subtitle = paste0(
           nrow(dis), " of ", nrow(d), " samples  |  harmonize = ", mode, "\n",
           "Ordered by ensemble P(beta); black bar = ensemble; dashed = 0.5\n",
           "Right: n voting beta (colour = ensemble call)"
         ),
         x = "P(beta)", y = NULL) +
    theme_classic(base_size = 12)
  h_dot <- max(3.5, 1.8 + 0.32 * nrow(dis))
  ggsave(file.path(output_dir, paste0("disagree_samples_dotplot_", mode, ".pdf")), p_dot, width = 9, height = h_dot)
  ggsave(file.path(output_dir, paste0("disagree_samples_dotplot_", mode, ".png")), p_dot, width = 9, height = h_dot, dpi = 300)
} else {
  message("      All models agree on every sample; no dot plot.")
}

# ============================================================
# Clinical comparisons (same as 11), on a given sample set
# ============================================================
clin_vars <- c("treatment" = "Treatment arm", "sex" = "Sex", "stratum" = "Onset (stratum)",
               "C9orf72" = "C9orf72 status", "progression_group" = "Progression")

categorical_tests <- function(src) {
  vars <- clin_vars[names(clin_vars) %in% colnames(src)]
  long <- bind_rows(lapply(names(vars), function(v) {
    cat_v <- as.character(src[[v]])
    if (v == "progression_group") {
      cat_v <- ifelse(is.na(cat_v), NA, paste0(cat_v, ifelse(cat_v == "fast", paste0(" (>", pr_cutoff, ")"),
                                                                       paste0(" (<=", pr_cutoff, ")"))))
    }
    data.frame(variable = vars[[v]], k2 = src$ensemble_k2,
               category = ifelse(is.na(cat_v), "missing", sub("^C9orf72_", "C9orf72 ", cat_v)),
               stringsAsFactors = FALSE)
  }))
  tab <- long %>% count(variable, k2, category, name = "n") %>%
    group_by(variable, k2) %>% mutate(k2_total = sum(n), percent = 100 * n / k2_total) %>% ungroup()
  p <- long %>% filter(category != "missing") %>% group_by(variable) %>%
    summarise(fisher_p = tryCatch(fisher.test(table(k2, category))$p.value, error = function(e) NA_real_),
              .groups = "drop")
  list(tab = tab %>% left_join(p, by = "variable") %>% arrange(variable, k2, category), p = p, vars = vars)
}

choose_test_two_groups <- function(y, g) {
  g <- droplevels(factor(g))
  grp <- split(y, g)
  if (min(lengths(grp)) < 3) return(list(name = "Wilcoxon rank-sum", test = wilcox.test(y ~ g, exact = FALSE)))
  sw <- vapply(grp, function(v) tryCatch(shapiro.test(v)$p.value, error = function(e) NA_real_), numeric(1))
  if (all(!is.na(sw) & sw > 0.05)) {
    lev <- if (requireNamespace("car", quietly = TRUE)) {
      tryCatch(car::leveneTest(y ~ g)$`Pr(>F)`[1], error = function(e) NA_real_)
    } else NA_real_
    if (!is.na(lev) && lev > 0.05) list(name = "Student's t-test", test = t.test(y ~ g, var.equal = TRUE))
    else list(name = "Welch's t-test", test = t.test(y ~ g))
  } else {
    list(name = "Wilcoxon rank-sum", test = wilcox.test(y ~ g, exact = FALSE))
  }
}
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 2))
  sprintf("%.3f", p)
}
cont_label <- c("progression_rate" = "Progression rate (points/month)",
                "months_onset_to_BL" = "Months from onset to baseline",
                "NfL" = "NfL (pg/ml)", "NEFL_NPQ" = "NULISA NEFL NPQ (CSF)", "GFAP" = "GFAP (pg/ml)",
                "age" = "Age (years)")
cont_vars <- trimws(strsplit(opt$cont_vars, ",", fixed = TRUE)[[1]])
cont_vars <- intersect(cont_vars[nzchar(cont_vars)], names(d))

continuous_tests <- function(src, xlab = "Predicted k2", make_plots = FALSE) {
  res <- list(); plots <- list()
  for (v in cont_vars) {
    dd <- data.frame(k2 = src$ensemble_k2, value = suppressWarnings(as.numeric(src[[v]])))
    dd <- dd[!is.na(dd$value) & !is.na(dd$k2), ]
    dd$k2 <- droplevels(factor(dd$k2, levels = c("alpha", "beta")))
    if (nlevels(dd$k2) != 2L || nrow(dd) < 4L) next
    tst <- choose_test_two_groups(dd$value, dd$k2)
    n_by <- table(dd$k2)
    res[[v]] <- data.frame(variable = v, n_alpha = as.integer(n_by["alpha"]), n_beta = as.integer(n_by["beta"]),
                           median_alpha = median(dd$value[dd$k2 == "alpha"]),
                           median_beta = median(dd$value[dd$k2 == "beta"]),
                           mean_alpha = mean(dd$value[dd$k2 == "alpha"]),
                           mean_beta = mean(dd$value[dd$k2 == "beta"]),
                           test = tst$name, p_value = tst$test$p.value, check.names = FALSE)
    if (make_plots) {
      lab <- if (v %in% names(cont_label)) cont_label[[v]] else v
      plots[[v]] <- ggplot(dd, aes(x = .data$k2, y = .data$value, fill = .data$k2)) +
        geom_violin(trim = FALSE, color = "black", fill = "white") +
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8, color = "black") +
        scale_fill_manual(values = k2_cols, name = "Predicted k2") +
        scale_x_discrete(labels = paste0(levels(dd$k2), "\n(n = ", as.integer(n_by[levels(dd$k2)]), ")")) +
        labs(x = xlab, y = lab, subtitle = paste0(tst$name, ", p = ", format_p(tst$test$p.value))) +
        theme_classic(base_size = 14) +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
              legend.text = element_text(size = 12), legend.title = element_text(size = 14),
              plot.subtitle = element_text(size = 11))
    }
  }
  list(stats = bind_rows(res), plots = plots)
}

src_all <- d
src_agree <- d[d$all_agree, ]
note <- paste0("All ", length(models), " models agree: ", sum(src_agree$ensemble_k2 == "alpha"), " alpha, ",
               sum(src_agree$ensemble_k2 == "beta"), " beta (", sum(!d$all_agree), " samples excluded)")

# ---- treatment ----
message("[4/6] Treatment arm and categorical variables")
if ("treatment" %in% names(src_agree)) {
  by_trt <- src_agree %>% count(treatment, ensemble_k2) %>% group_by(treatment) %>%
    mutate(fraction = n / sum(n)) %>% ungroup()
  write.csv(by_trt, file.path(output_dir, paste0("predicted_k2_by_treatment_", mode, sfx_out, ".csv")), row.names = FALSE)
  p_trt <- ggplot(src_agree, aes(x = .data$treatment, fill = .data$ensemble_k2)) +
    geom_bar(position = "fill", width = 0.65) +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_fill(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = k2_cols, name = "Predicted k2") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "ROCK-ALS CSF: predicted k2 by treatment arm (all models agree)",
         subtitle = paste0("harmonize = ", mode, "\n", note), x = NULL, y = "Fraction") +
    theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 15, hjust = 1))
  ggsave(file.path(output_dir, paste0("predicted_k2_by_treatment_", mode, sfx_out, ".pdf")), p_trt, width = 7, height = 5)
}

# ---- categorical ----
cat_agree <- categorical_tests(src_agree)
write.csv(cat_agree$tab, file.path(output_dir, paste0("clinical_by_predicted_k2_", mode, sfx_out, ".csv")),
          row.names = FALSE)
if (nrow(cat_agree$tab) > 0L) {
  tab <- cat_agree$tab
  facet_lab <- setNames(paste0(cat_agree$p$variable, "\n(Fisher p = ", signif(cat_agree$p$fisher_p, 2), ")"),
                        cat_agree$p$variable)
  tab$variable <- factor(tab$variable, levels = unname(cat_agree$vars))
  tab$k2_lab <- paste0(tab$k2, "\n(n = ", tab$k2_total, ")")
  tab$label <- paste0(tab$category, "\n", tab$n, " (", round(tab$percent), "%)")
  cats <- sort(unique(tab$category))
  cat_cols <- setNames(rep(c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                             "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F", "#D9D9D9"),
                           length.out = length(cats)), cats)
  if ("missing" %in% names(cat_cols)) cat_cols[["missing"]] <- "gray85"
  p_clin <- ggplot(tab, aes(x = .data$k2_lab, y = .data$percent, fill = .data$category)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = .data$label), position = position_stack(vjust = 0.5), size = 2.6, lineheight = 0.85) +
    facet_wrap(~ variable, nrow = 1, scales = "free_x", labeller = labeller(variable = facet_lab)) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "ROCK-ALS CSF: clinical characteristics by predicted k2 (all models agree)",
         subtitle = paste0("Ensemble prediction (harmonize = ", mode, "); labels = n (% within predicted group)\n", note),
         x = "Predicted k2", y = "Patients (%)") +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(file.path(output_dir, paste0("clinical_by_predicted_k2_", mode, sfx_out, ".pdf")), p_clin,
         width = 3 * length(cat_agree$vars) + 1, height = 5.5)
}

# ---- continuous ----
message("[5/6] Continuous variables")
cont_dir <- file.path(output_dir, paste0("clinical_continuous_by_predicted_k2_", mode, sfx_out))
dir.create(cont_dir, recursive = TRUE, showWarnings = FALSE)
cont_agree <- continuous_tests(src_agree, xlab = "Predicted k2 (all models agree)", make_plots = TRUE)
write.csv(cont_agree$stats, file.path(cont_dir, "continuous_by_predicted_k2_tests.csv"), row.names = FALSE)
safe_name <- function(x) gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
for (v in names(cont_agree$plots)) {
  suppressMessages(ggsave(file.path(cont_dir, paste0("k2_vs_", safe_name(v), ".pdf")), cont_agree$plots[[v]],
                          width = 4.8, height = 4.2))
}
if (length(cont_agree$plots) > 0L && requireNamespace("gridExtra", quietly = TRUE)) {
  n_p <- length(cont_agree$plots)
  ncol_m <- if (n_p <= 2) n_p else if (n_p <= 4) 2 else 3
  nrow_m <- ceiling(n_p / ncol_m)
  merged <- suppressMessages(do.call(gridExtra::arrangeGrob,
                    c(lapply(cont_agree$plots, function(p) p + theme(legend.position = "none")),
                      ncol = ncol_m, nrow = nrow_m)))
  suppressMessages(ggsave(file.path(cont_dir, "k2_vs_merged.pdf"), merged, width = 4.8 * ncol_m, height = 4.2 * nrow_m))
}

# ---- comparison: all samples vs all-agree samples ----
message("[6/6] Comparison with all samples")
cat_all <- categorical_tests(src_all)
cont_all <- continuous_tests(src_all)
cmp_cont <- full_join(
  cont_all$stats %>% transmute(variable, type = "continuous",
                               n_all = paste0(n_alpha, "/", n_beta),
                               median_alpha_all = median_alpha, median_beta_all = median_beta,
                               test_all = test, p_all = p_value),
  cont_agree$stats %>% transmute(variable, n_allAgree = paste0(n_alpha, "/", n_beta),
                                 median_alpha_allAgree = median_alpha, median_beta_allAgree = median_beta,
                                 test_allAgree = test, p_allAgree = p_value),
  by = "variable")
cmp_cat <- full_join(cat_all$p %>% rename(p_all = fisher_p), cat_agree$p %>% rename(p_allAgree = fisher_p),
                     by = "variable") %>%
  mutate(type = "categorical", test_all = "Fisher's exact", test_allAgree = "Fisher's exact",
         n_all = paste0(sum(src_all$ensemble_k2 == "alpha"), "/", sum(src_all$ensemble_k2 == "beta")),
         n_allAgree = paste0(sum(src_agree$ensemble_k2 == "alpha"), "/", sum(src_agree$ensemble_k2 == "beta")))
cmp <- bind_rows(cmp_cont, cmp_cat) %>%
  select(variable, type, n_all, n_allAgree, any_of(c("median_alpha_all", "median_beta_all",
         "median_alpha_allAgree", "median_beta_allAgree")), test_all, p_all, test_allAgree, p_allAgree)
write.csv(cmp, file.path(output_dir, paste0("comparison_all_vs_allAgree_", mode, ".csv")), row.names = FALSE)

cat("\nAll samples vs all-agree samples (", mode, "):\n", sep = "")
print(cmp %>% select(variable, n_all, n_allAgree, p_all, p_allAgree) %>%
        mutate(across(c(p_all, p_allAgree), ~ signif(.x, 3))), row.names = FALSE)
message("Done. Output: ", output_dir)
