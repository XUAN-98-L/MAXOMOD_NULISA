#!/usr/bin/env Rscript
# ============================================================
# 13_BetweenFluid_ML_Proteins.R
#
# For the ML-selected proteins (selected_features.txt from 08), measure how well
# blood levels track CSF levels: CSF is the baseline, and for every protein the
# Spearman correlation of Plasma vs CSF and Serum vs CSF is computed across
# participants (samples paired by participant ID).
#
# Heatmap: fill = Spearman rho, text = "rho (p value)" (raw by default; --p_adjust BH for FDR).
#
# Usage (from the project root):
# Rscript Script/13_BetweenFluid_ML_Proteins.R -i CNS_immune/Results --panel cns
# Rscript Script/13_BetweenFluid_ML_Proteins.R -i CNS_immune/Results --panel cns \
#   --features CNS_immune/Results/CNS_panel/Without_tears/08_Clustering_Feature_Selection_logLoss/selected_features.txt
# Rscript Script/13_BetweenFluid_ML_Proteins.R -i CNS_immune/Results --panel cns --group ALS
#
# Outputs (default: <panel>/Without_tears/13_BetweenFluid_ML_Proteins/):
#   BetweenFluid_spearman_heatmap.pdf / .png
#   BetweenFluid_spearman_results.csv   (rho, p, BH-adjusted p, n per protein x fluid)
#   BetweenFluid_paired_values.csv      (the paired NPQ values used)
# ============================================================

suppressMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "CNS_immune/Results",
              help = "Results root containing CNS_panel/ and Immune_panel/ [default: %default]"),
  make_option(c("--panel"), type = "character", default = "cns",
              help = "'cns', 'immune' or 'both' (which panels the selected features come from) [default: %default]"),
  make_option(c("-f", "--features"), type = "character", default = NULL,
              help = paste("Path to selected_features.txt. Default:",
                           "<panel dir>/Without_tears/08_Clustering_Feature_Selection_logLoss/selected_features.txt")),
  make_option(c("--id_col"), type = "character", default = NULL,
              help = paste("Metadata column identifying the participant (same value for CSF/Plasma/Serum",
                           "of one person). Default: auto-detect from common names")),
  make_option(c("--group"), type = "character", default = "all",
              help = paste("Participants to include: 'all', 'ALS' (k2 alpha/beta only),",
                           "'alpha', 'beta' or 'control' (no k2 label) [default: %default]")),
  make_option(c("--p_adjust"), type = "character", default = "none",
              help = "p value shown in the heatmap: 'none' (raw p) or 'BH' (FDR across all cells) [default: %default]"),
  make_option(c("--include_ps"), type = "logical", default = FALSE,
              help = "Also add a Plasma vs Serum column [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder [default: <panel dir>/Without_tears/13_BetweenFluid_ML_Proteins]")
)
opt <- parse_args(OptionParser(option_list = option_list))

panel_choice <- tolower(trimws(opt$panel))
if (!panel_choice %in% c("cns", "immune", "both")) stop("--panel must be 'cns', 'immune' or 'both'.")
group_choice <- tolower(trimws(opt$group))
if (!group_choice %in% c("all", "als", "alpha", "beta", "control")) {
  stop("--group must be 'all', 'ALS', 'alpha', 'beta' or 'control'.")
}
p_adjust <- toupper(trimws(opt$p_adjust))
if (!p_adjust %in% c("NONE", "BH")) stop("--p_adjust must be 'none' or 'BH'.")

panel_dirs <- c(CNS = "CNS_panel", Immune = "Immune_panel")
panels_used <- switch(panel_choice, cns = "CNS", immune = "Immune", both = c("CNS", "Immune"))
panel_root <- function(p) file.path(opt$input, panel_dirs[[p]], "Without_tears")

# Default locations follow the folder layout of 08/09
base_dir <- if (panel_choice == "both") {
  file.path(opt$input, "Combined_CNS_Immune_panels", "Without_tears")
} else {
  panel_root(panels_used)
}
features_path <- if (!is.null(opt$features)) opt$features else
  file.path(base_dir, "08_Clustering_Feature_Selection_logLoss", "selected_features.txt")
output_dir <- if (!is.null(opt$output)) opt$output else
  file.path(base_dir, "13_BetweenFluid_ML_Proteins")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- if (group_choice == "all") "" else paste0("_", group_choice)   # file suffix for --group

if (!file.exists(features_path)) stop("Feature file not found: ", features_path)

message("[1/4] Settings")
message("      Features: ", features_path)
message("      Panel(s): ", paste(panels_used, collapse = " + "), "; group: ", opt$group)
message("      Output:   ", output_dir)

# ============================================================
# Load selected features
# ============================================================
features <- trimws(readLines(features_path, warn = FALSE))
features <- unique(features[nzchar(features)])
if (length(features) == 0L) stop("No features in ", features_path)

has_prefix <- grepl("^(CNS|Immune)_", features)
feat_map <- data.frame(
  Feature = features,
  Panel = ifelse(has_prefix, sub("^(CNS|Immune)_.*$", "\\1", features),
                 if (length(panels_used) == 1L) panels_used else NA_character_),
  Target = sub("^(CNS|Immune)_", "", features),
  stringsAsFactors = FALSE
)
if (anyNA(feat_map$Panel)) {
  stop("Features without a CNS_/Immune_ prefix are ambiguous with --panel both: ",
       paste(feat_map$Feature[is.na(feat_map$Panel)], collapse = ", "))
}
feat_map <- feat_map[feat_map$Panel %in% panels_used, , drop = FALSE]
if (nrow(feat_map) == 0L) stop("None of the selected features belong to panel(s) ", paste(panels_used, collapse = ", "))
message("      Selected proteins: ", nrow(feat_map))

# ============================================================
# Load NPQ + metadata, pair samples across fluids by participant
# ============================================================
message("[2/4] Loading NPQ and metadata")

id_candidates <- c("Participant_ID", "participant_id", "ParticipantID", "Participant",
                   "Patient_ID", "PatientID", "Subject_ID", "SubjectID", "Subject",
                   "MAXOMOD_ID", "Proband_ID", "Individual_ID", "Donor_ID", "ID")

pick_id_col <- function(meta) {
  if (!is.null(opt$id_col)) {
    if (!opt$id_col %in% names(meta)) {
      stop("--id_col '", opt$id_col, "' not in metadata. Columns: ", paste(names(meta), collapse = ", "))
    }
    return(opt$id_col)
  }
  hit <- intersect(id_candidates, names(meta))
  # Must give one sample per fluid for (most) participants, i.e. not be the tube ID
  for (h in hit) {
    per_id <- table(meta[[h]])
    if (h != "Tube_ID" && stats::median(per_id) > 1) return(h)
  }
  stop("Could not auto-detect the participant ID column. Use --id_col. Metadata columns: ",
       paste(names(meta), collapse = ", "))
}

load_panel <- function(p) {
  init_dir <- file.path(panel_root(p), "00_Initialization")
  npq <- read_excel(file.path(init_dir, "npq_counts.xlsx"))
  meta <- read_excel(file.path(init_dir, "all_participants_IDs.xlsx"))
  id_col <- pick_id_col(meta)

  meta <- meta %>%
    mutate(Fluid = dplyr::recode(toupper(.data$Material), CSF = "CSF", PLASMA = "Plasma", SERUM = "Serum",
                                 .default = NA_character_),
           Participant = as.character(.data[[id_col]]),
           Tube_ID = as.character(.data$Tube_ID)) %>%
    filter(!is.na(.data$Fluid), !is.na(.data$Participant))

  # Group filter uses k2 (alpha/beta = ALS subtypes); k2 is taken per participant
  if ("k2" %in% names(meta)) {
    k2_pp <- meta %>% filter(!is.na(.data$k2)) %>% distinct(.data$Participant, .data$k2)
    keep_ids <- switch(group_choice,
      all = unique(meta$Participant),
      als = k2_pp$Participant,
      alpha = k2_pp$Participant[k2_pp$k2 == "alpha"],
      beta = k2_pp$Participant[k2_pp$k2 == "beta"],
      control = setdiff(unique(meta$Participant), k2_pp$Participant))
    meta <- meta[meta$Participant %in% keep_ids, , drop = FALSE]
  } else if (group_choice != "all") {
    stop("Metadata has no k2 column; --group must be 'all'.")
  }

  targets <- feat_map$Target[feat_map$Panel == p]
  out <- npq %>%
    filter(.data$SampleType == "Sample", .data$Target %in% targets) %>%
    mutate(SampleName = as.character(.data$SampleName)) %>%
    inner_join(meta %>% select("Tube_ID", "Participant", "Fluid"), by = c("SampleName" = "Tube_ID")) %>%
    group_by(.data$Participant, .data$Fluid, .data$Target) %>%
    summarise(NPQ = mean(.data$NPQ, na.rm = TRUE), .groups = "drop") %>%
    mutate(Panel = p)

  miss <- setdiff(targets, unique(out$Target))
  if (length(miss) > 0L) warning(p, " targets not found in NPQ data: ", paste(miss, collapse = ", "))
  message("      ", p, ": participant ID column '", id_col, "', ",
          n_distinct(out$Participant), " participants, ",
          paste(names(table(out$Fluid)), collapse = "/"), " samples matched")
  out
}

long <- bind_rows(lapply(panels_used, load_panel)) %>%
  left_join(feat_map, by = c("Panel", "Target"))

wide <- long %>%
  select("Feature", "Panel", "Target", "Participant", "Fluid", "NPQ") %>%
  pivot_wider(names_from = "Fluid", values_from = "NPQ")
for (fl in c("CSF", "Plasma", "Serum")) if (!fl %in% names(wide)) wide[[fl]] <- NA_real_
write.csv(wide, file.path(output_dir, paste0("BetweenFluid_paired_values", suffix, ".csv")), row.names = FALSE)

# ============================================================
# Spearman correlations (CSF = baseline)
# ============================================================
message("[3/4] Spearman correlations")

comparisons <- list(c("Plasma", "CSF"), c("Serum", "CSF"))
if (isTRUE(opt$include_ps)) comparisons <- c(comparisons, list(c("Plasma", "Serum")))

cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(c(rho = NA_real_, p = NA_real_, n = n))
  }
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = n)
}

res <- bind_rows(lapply(comparisons, function(cmp) {
  wide %>%
    group_by(.data$Feature, .data$Panel, .data$Target) %>%
    summarise(stat = list(cor_one(.data[[cmp[1]]], .data[[cmp[2]]])), .groups = "drop") %>%
    mutate(Comparison = paste0(cmp[1], " vs ", cmp[2]),
           rho = vapply(.data$stat, `[[`, numeric(1), "rho"),
           p = vapply(.data$stat, `[[`, numeric(1), "p"),
           n = vapply(.data$stat, `[[`, numeric(1), "n")) %>%
    select(-"stat")
})) %>%
  mutate(p_BH = stats::p.adjust(.data$p, method = "BH"))

write.csv(res %>% arrange(.data$Comparison, desc(.data$rho)),
          file.path(output_dir, paste0("BetweenFluid_spearman_results", suffix, ".csv")), row.names = FALSE)
print(as.data.frame(res %>% select("Target", "Comparison", "rho", "p", "p_BH", "n")), digits = 3)

# ============================================================
# Heatmap
# ============================================================
message("[4/4] Heatmap")

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
}
p_show <- if (p_adjust == "BH") res$p_BH else res$p
p_lab <- if (p_adjust == "BH") "FDR (BH)" else "p"

# Row label: target only when one panel, else keep the panel prefix
res$Label <- if (length(panels_used) == 1L) res$Target else res$Feature
# Rows ordered by mean rho across comparisons (highest at top)
row_order <- res %>% group_by(.data$Label) %>%
  summarise(m = mean(.data$rho, na.rm = TRUE), .groups = "drop") %>%
  arrange(.data$m) %>% pull(.data$Label)

n_lab <- res %>% group_by(.data$Comparison) %>%
  summarise(n = paste0(min(.data$n), if (min(.data$n) != max(.data$n)) paste0("-", max(.data$n)) else ""),
            .groups = "drop")
col_levels <- vapply(comparisons, function(cmp) paste0(cmp[1], " vs ", cmp[2]), character(1))
col_labels <- setNames(paste0(col_levels, "\n(n = ", n_lab$n[match(col_levels, n_lab$Comparison)], ")"),
                       col_levels)

plot_df <- res %>%
  mutate(Label = factor(.data$Label, levels = row_order),
         Comparison = factor(.data$Comparison, levels = col_levels),
         cell_text = ifelse(is.na(.data$rho), "NA",
                            paste0(sprintf("%.2f", .data$rho), " (", fmt_p(p_show), ")")),
         txt_col = ifelse(!is.na(.data$rho) & abs(.data$rho) > 0.6, "white", "black"))

group_txt <- switch(group_choice, all = "all participants", als = "ALS (alpha + beta)",
                    alpha = "ALS alpha", beta = "ALS beta", control = "controls")

hm <- ggplot(plot_df, aes(x = .data$Comparison, y = .data$Label, fill = .data$rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = .data$cell_text, color = .data$txt_col), size = 3.2) +
  scale_color_identity() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), na.value = "grey85", name = "Spearman\nrho") +
  scale_x_discrete(labels = col_labels, position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Blood vs CSF correlation of ML-selected proteins",
       subtitle = paste0(paste(panels_used, collapse = " + "), " panel, NPQ, ", group_txt,
                         "; CSF = baseline. Cell text: rho (", p_lab, ")"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.title = element_text(size = 9))

h <- max(3.5, 1.6 + 0.32 * length(row_order))
w <- 3.2 + 1.9 * length(col_levels)
ggsave(file.path(output_dir, paste0("BetweenFluid_spearman_heatmap", suffix, ".pdf")), hm, width = w, height = h)
ggsave(file.path(output_dir, paste0("BetweenFluid_spearman_heatmap", suffix, ".png")), hm, width = w, height = h, dpi = 300)

message("Done. Results in ", output_dir)
