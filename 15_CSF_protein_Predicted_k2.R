#!/usr/bin/env Rscript
# ============================================================
# 15_CSF_protein_Predicted_k2.R
#
# CSF levels (raw NPQ) of the ML-selected proteins (selected_features.txt from 08) by subtype.
# One figure per protein:
#   left  = MAXOMOD CSF, true k2 (alpha vs beta)
#   right = ROCK-ALS CSF, ensemble-predicted k2 from 11_with_ROCK_ALS.R
# Each panel: violin + box + points, Wilcoxon rank-sum p-value (and BH over the proteins of that
# cohort). The two panels have their own y axis (the cohorts were run separately on NULISA).
#
# Inputs
#   --features     selected_features.txt (08)
#   --npq / --metadata   MAXOMOD 00_Initialization/npq_counts.xlsx and all_participants_IDs.xlsx
#   --rock_pred    11 output ROCK_ALS_CSF_predicted_k2_<mode>.csv (PatNr, ensemble_k2)
#   --rock_expr    ROCK-ALS CSF NULISA data, long format (Target, NPQ, PatNr), as used by 11
#
# Outputs (-o, default .../Without_tears/15_CSF_protein_Predicted_k2)
#   <Protein>_MAXOMOD_vs_ROCK_ALS.pdf       one figure per protein
#   all_proteins_MAXOMOD_vs_ROCK_ALS.pdf    all proteins, one page each
#   wilcox_alpha_vs_beta.csv                n, medians, Wilcoxon p, BH p, direction per cohort and protein
#   expression_values_long.csv              the raw NPQ values that were plotted
#
# Example (run from NULISA_MAXOMOD/):
# Rscript Script/15_CSF_protein_Predicted_k2.R
# ============================================================

suppressMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base <- "CNS_immune/Results/CNS_panel/Without_tears"
option_list <- list(
  make_option(c("--features"), type = "character",
              default = file.path(base, "08_Clustering_Feature_Selection_logLoss", "selected_features.txt"),
              help = "selected_features.txt from 08 [default: %default]"),
  make_option(c("--npq"), type = "character", default = file.path(base, "00_Initialization", "npq_counts.xlsx"),
              help = "MAXOMOD NPQ (long format) [default: %default]"),
  make_option(c("--metadata"), type = "character",
              default = file.path(base, "00_Initialization", "all_participants_IDs.xlsx"),
              help = "MAXOMOD metadata (Tube_ID, Material, k2) [default: %default]"),
  make_option(c("--rock_pred"), type = "character",
              default = file.path(base, "11_with_ROCK_ALS_logLoss", "ROCK_ALS_CSF_predicted_k2_standardize.csv"),
              help = "11 ROCK-ALS predictions (PatNr, ensemble_k2) [default: %default]"),
  make_option(c("--rock_expr"), type = "character",
              default = "/Users/xliu2942/Documents/Projects/ROCK_ALS/Results/csf/00_Initialization/Expression.csv",
              help = "ROCK-ALS CSF NULISA data, long format (Target, NPQ, PatNr) [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = file.path(base, "15_CSF_protein_Predicted_k2"),
              help = "Output folder [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(opt$features, opt$npq, opt$metadata, opt$rock_pred)) if (!file.exists(f)) stop("Not found: ", f)
if (!file.exists(opt$rock_expr)) stop("ROCK-ALS expression file not found: ", opt$rock_expr, "\nGive it with --rock_expr.")

k2_cols <- c(alpha = "#FC8D62", beta = "#8DA0CB")
cohort_lv <- c("MAXOMOD (true k2)", "ROCK-ALS (predicted k2)")

features <- trimws(readLines(opt$features, warn = FALSE))
features <- features[nzchar(features)]
targets <- sub("^(CNS|Immune)_", "", features)
message("[1/4] Proteins (", length(targets), "): ", paste(targets, collapse = ", "))

# ------------------------------------------------------------
# 1. MAXOMOD CSF, true k2 (raw NPQ)
# ------------------------------------------------------------
message("[2/4] MAXOMOD CSF")
meta <- read_excel(opt$metadata) %>%
  filter(toupper(Material) == "CSF", k2 %in% c("alpha", "beta")) %>%
  transmute(Sample = as.character(Tube_ID), k2 = as.character(k2))
maxo <- read_excel(opt$npq) %>%
  filter(toupper(SampleMatrixType) == "CSF", SampleType == "Sample", Target %in% targets) %>%
  mutate(Sample = as.character(SampleName)) %>%
  inner_join(meta, by = "Sample") %>%
  group_by(Sample, k2, Target) %>% summarise(value = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
  mutate(cohort = cohort_lv[1])
miss <- setdiff(targets, unique(maxo$Target))
if (length(miss) > 0L) stop("Proteins not in the MAXOMOD NPQ data: ", paste(miss, collapse = ", "))
message("      alpha / beta: ", n_distinct(maxo$Sample[maxo$k2 == "alpha"]), " / ",
        n_distinct(maxo$Sample[maxo$k2 == "beta"]))

# ------------------------------------------------------------
# 2. ROCK-ALS CSF, ensemble-predicted k2 (raw NPQ)
# ------------------------------------------------------------
message("[3/4] ROCK-ALS CSF (predicted k2)")
pred <- read.csv(opt$rock_pred, check.names = FALSE, stringsAsFactors = FALSE) %>%
  transmute(Sample = gsub("-", "_", PatNr), k2 = ensemble_k2) %>%
  filter(k2 %in% c("alpha", "beta"))
rock <- read.csv(opt$rock_expr, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(Sample = gsub("-", "_", PatNr)) %>%
  filter(Target %in% targets) %>%
  inner_join(pred, by = "Sample") %>%
  group_by(Sample, k2, Target) %>% summarise(value = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
  mutate(cohort = cohort_lv[2])
miss_r <- setdiff(targets, unique(rock$Target))
if (length(miss_r) > 0L) stop("Proteins not in the ROCK-ALS data: ", paste(miss_r, collapse = ", "))
message("      predicted alpha / beta: ", n_distinct(rock$Sample[rock$k2 == "alpha"]), " / ",
        n_distinct(rock$Sample[rock$k2 == "beta"]),
        " (of ", nrow(pred), " predicted samples)")

# ------------------------------------------------------------
# 3. Wilcoxon tests
# ------------------------------------------------------------
message("[4/4] Wilcoxon tests and plots")
all_long <- bind_rows(maxo, rock) %>%
  filter(!is.na(value)) %>%
  rename(Protein = Target) %>%
  mutate(Protein = factor(Protein, levels = targets), k2 = factor(k2, levels = c("alpha", "beta")),
         cohort = factor(cohort, levels = cohort_lv))
write.csv(all_long %>% select(cohort, Sample, k2, Protein, NPQ = value),
          file.path(output_dir, "expression_values_long.csv"), row.names = FALSE)

stats_tab <- all_long %>% group_by(cohort, Protein) %>%
  summarise(n_alpha = sum(k2 == "alpha"), n_beta = sum(k2 == "beta"),
            median_alpha = median(value[k2 == "alpha"]), median_beta = median(value[k2 == "beta"]),
            p_wilcox = tryCatch(wilcox.test(value[k2 == "alpha"], value[k2 == "beta"], exact = FALSE)$p.value,
                                error = function(e) NA_real_),
            .groups = "drop") %>%
  group_by(cohort) %>% mutate(p_BH = p.adjust(p_wilcox, method = "BH")) %>% ungroup() %>%
  mutate(higher_in = ifelse(median_alpha > median_beta, "alpha", "beta"))
dir_wide <- stats_tab %>% select(cohort, Protein, higher_in) %>%
  mutate(cohort = ifelse(grepl("^MAXOMOD", cohort), "MAXOMOD", "ROCK")) %>%
  pivot_wider(names_from = cohort, values_from = higher_in)
stats_tab <- stats_tab %>% left_join(dir_wide %>% transmute(Protein, same_direction = MAXOMOD == ROCK), by = "Protein")
write.csv(stats_tab, file.path(output_dir, "wilcox_alpha_vs_beta.csv"), row.names = FALSE)
print(as.data.frame(stats_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# ------------------------------------------------------------
# 4. One figure per protein: MAXOMOD (left) | ROCK-ALS (right)
# ------------------------------------------------------------
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
star <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
lab_tab <- stats_tab %>%
  left_join(all_long %>% group_by(cohort, Protein) %>%
              summarise(ymax = max(value), yr = diff(range(value)), .groups = "drop"),
            by = c("cohort", "Protein")) %>%
  mutate(label = paste0("Wilcoxon p = ", fmt_p(p_wilcox), " ", star(p_wilcox), "\nBH p = ", fmt_p(p_BH)),
         y = ymax + 0.12 * yr, y_bar = ymax + 0.06 * yr)

# x labels with n per group and cohort
n_tab <- all_long %>% group_by(cohort, Protein, k2) %>% summarise(n = n_distinct(Sample), .groups = "drop")

protein_plot <- function(prot) {
  dd <- all_long %>% filter(Protein == prot)
  ll <- lab_tab %>% filter(Protein == prot)
  nn <- n_tab %>% filter(Protein == prot)
  dd <- dd %>% left_join(nn, by = c("cohort", "Protein", "k2")) %>%
    mutate(x = factor(paste0(k2, "\n(n = ", n, ")")))
  # keep alpha left of beta within each panel
  x_lv <- unique(dd$x[order(dd$k2)])
  dd$x <- factor(dd$x, levels = x_lv)
  ggplot(dd, aes(x = .data$x, y = .data$value)) +
    geom_violin(aes(fill = .data$k2), trim = FALSE, alpha = 0.35, color = "gray30", linewidth = 0.4) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.4) +
    geom_jitter(aes(color = .data$k2), width = 0.12, height = 0, size = 1.6, alpha = 0.8) +
    geom_segment(data = ll, aes(x = 1, xend = 2, y = .data$y_bar, yend = .data$y_bar), inherit.aes = FALSE,
                 linewidth = 0.4) +
    geom_text(data = ll, aes(x = 1.5, y = .data$y, label = .data$label), inherit.aes = FALSE, size = 3.3,
              vjust = 0, lineheight = 0.9) +
    facet_wrap(~ cohort, nrow = 1, scales = "free") +
    scale_fill_manual(values = k2_cols, guide = "none") +
    scale_color_manual(values = k2_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = paste0(prot, ": CSF level by subtype"),
         subtitle = "Raw NPQ; Wilcoxon rank-sum test (BH over the proteins within each cohort)",
         x = NULL, y = paste0(prot, " NPQ (CSF)")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12))
}

plots <- lapply(targets, protein_plot)
names(plots) <- targets
safe <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
for (p in targets) {
  ggsave(file.path(output_dir, paste0(safe(p), "_MAXOMOD_vs_ROCK_ALS.pdf")), plots[[p]], width = 7.5, height = 5)
}
pdf(file.path(output_dir, "all_proteins_MAXOMOD_vs_ROCK_ALS.pdf"), width = 7.5, height = 5)
for (p in targets) print(plots[[p]])
invisible(dev.off())

message("Done. Output: ", output_dir)
