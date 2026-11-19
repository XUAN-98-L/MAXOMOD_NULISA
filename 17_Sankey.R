#!/usr/bin/env Rscript
# ============================================================
# 17_Sankey.R
#
# Stability of the ROCK-ALS subtype calls between two visits (default: serum V1 -> V20).
#
#   1. For each visit separately, find the patients on which ALL ML models agree
#      (every k2_<model> column of the 11b output gives the same call).
#      Per visit each patient is then "alpha" (all models alpha), "beta" (all models beta)
#      or "models disagree".
#   2. Sankey (alluvial) plot V1 -> V20 of the patients that are all-agree at BOTH visits:
#      did any patient switch alpha <-> beta?
#      Exact McNemar test (binomial test on the switchers alpha->beta vs beta->alpha).
#   3. Context plots: all patients, with "models disagree" as a third group, and the
#      both-visits all-agree Sankey split by treatment arm.
#
# Input (read-only): two 11b_with_ROCK_ALS_serum.R output folders
#   --v1  : .../11b_with_ROCK_ALS_Serum_V1_logLoss    (ROCK_ALS_Serum_V1_predicted_k2_<mode>.csv)
#   --v20 : .../11b_with_ROCK_ALS_Serum_V20_logLoss   (ROCK_ALS_Serum_V20_predicted_k2_<mode>.csv)
#   --mode: harmonisation mode of the 11b file to use (none / center / standardize)
#
# Outputs (-o, default .../Without_tears/17_Sankey_Serum_logLoss):
#   per_patient_calls_<mode>.csv            every patient: per-visit status, ensemble P(beta), votes, switched flag
#   allAgree_counts_per_visit_<mode>.csv    alpha / beta / models disagree per visit
#   transitions_allAgree_both_<mode>.csv    V1 -> V20 transition counts (all-agree at both visits)
#   switched_patients_<mode>.csv            patients all-agree at both visits whose call changed
#   sankey_allAgree_both_<mode>.pdf/.png    main plot
#   sankey_all_patients_<mode>.pdf/.png     all patients, incl. "models disagree"
#   sankey_allAgree_both_by_arm_<mode>.pdf/.png
#
# Example (run from NULISA_MAXOMOD/):
# Rscript Script/17_Sankey.R
# Rscript Script/17_Sankey.R \
#   --v1  CNS_immune/Results/CNS_panel/Without_tears/11b_with_ROCK_ALS_Serum_V1_logLoss \
#   --v20 CNS_immune/Results/CNS_panel/Without_tears/11b_with_ROCK_ALS_Serum_V20_logLoss \
#   -o CNS_immune/Results/CNS_panel/Without_tears/17_Sankey_Serum_logLoss --mode standardize
# ============================================================

suppressMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
})

base_dir <- "CNS_immune/Results/CNS_panel/Without_tears"
option_list <- list(
  make_option(c("--v1"), type = "character",
              default = file.path(base_dir, "11b_with_ROCK_ALS_Serum_V1_logLoss"),
              help = "11b output folder of the first visit [default: %default]"),
  make_option(c("--v20"), type = "character",
              default = file.path(base_dir, "11b_with_ROCK_ALS_Serum_V20_logLoss"),
              help = "11b output folder of the second visit [default: %default]"),
  make_option(c("--labels"), type = "character", default = "V1,V20",
              help = "Axis labels of the two visits, comma-separated [default: %default]"),
  make_option(c("-o", "--output"), type = "character",
              default = file.path(base_dir, "17_Sankey_Serum_logLoss"),
              help = "Output folder [default: %default]"),
  make_option(c("--mode"), type = "character", default = "standardize",
              help = "Harmonisation mode of the 11b predictions: none, center or standardize [default: %default]"),
  make_option(c("--exclude_models"), type = "character", default = "none",
              help = paste("Comma-separated models NOT used for the all-agree check",
                           "(e.g. 'ridge'); 'none' = use every k2_<model> column [default: %default]"))
)
opt <- parse_args(OptionParser(option_list = option_list))

mode <- tolower(trimws(opt$mode))
if (!mode %in% c("none", "center", "standardize")) stop("--mode must be none, center or standardize.")
visit_labels <- trimws(strsplit(opt$labels, ",", fixed = TRUE)[[1]])
if (length(visit_labels) != 2L) stop("--labels needs exactly two labels, e.g. 'V1,V20'.")
excl <- trimws(strsplit(opt$exclude_models, ",", fixed = TRUE)[[1]])
excl <- setdiff(excl[nzchar(excl)], "none")

output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

status_levels <- c("alpha", "beta", "models disagree")
status_cols <- c(alpha = "#FC8D62", beta = "#8DA0CB", "models disagree" = "grey75")

# ------------------------------------------------------------
# 1. Read one visit and flag the all-agree patients
# ------------------------------------------------------------
read_visit <- function(dir, label) {
  f <- list.files(dir, pattern = paste0("^ROCK_ALS_.+_predicted_k2_", mode, "\\.csv$"), full.names = TRUE)
  if (length(f) != 1L) {
    stop("Expected one *_predicted_k2_", mode, ".csv in ", dir, ", found ", length(f), ".")
  }
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  k2_cols <- grep("^k2_", colnames(d), value = TRUE)
  k2_cols <- k2_cols[!sub("^k2_", "", k2_cols) %in% excl]
  if (length(k2_cols) < 2L) stop("Fewer than 2 model columns (k2_<model>) left in ", f)
  calls <- as.matrix(d[, k2_cols, drop = FALSE])
  n_beta <- rowSums(calls == "beta", na.rm = TRUE)
  n_alpha <- rowSums(calls == "alpha", na.rm = TRUE)
  n_ok <- rowSums(!is.na(calls))
  status <- ifelse(n_ok == length(k2_cols) & n_beta == length(k2_cols), "beta",
            ifelse(n_ok == length(k2_cols) & n_alpha == length(k2_cols), "alpha", "models disagree"))
  fluid_label <<- sub("^ROCK_ALS_([^_]+)_.*$", "\\1", basename(f))
  message(label, ": ", basename(f), " | ", nrow(d), " patients, ", length(k2_cols), " models (",
          paste(sub("^k2_", "", k2_cols), collapse = ", "), ")")
  message("   all agree: ", sum(status != "models disagree"), " (alpha ", sum(status == "alpha"),
          ", beta ", sum(status == "beta"), "); models disagree: ", sum(status == "models disagree"))
  calls_df <- data.frame(
    PatNr = d$PatNr,
    status = status,
    n_vote_beta = n_beta,
    n_models = length(k2_cols),
    ensemble_P_beta = if ("ensemble_P_beta" %in% colnames(d)) d$ensemble_P_beta else NA_real_,
    ensemble_k2 = if ("ensemble_k2" %in% colnames(d)) d$ensemble_k2 else NA_character_,
    stringsAsFactors = FALSE
  ) %>%
    rename_with(~ paste0(label, "_", .x), -PatNr)
  list(calls = calls_df,
       meta = d[, intersect(c("PatNr", "treatment", "arm", "group", "stratum", "sex"), colnames(d)), drop = FALSE])
}

v1 <- read_visit(opt$v1, visit_labels[1])
v2 <- read_visit(opt$v20, visit_labels[2])
s1 <- paste0(visit_labels[1], "_status")
s2 <- paste0(visit_labels[2], "_status")

only1 <- setdiff(v1$calls$PatNr, v2$calls$PatNr)
only2 <- setdiff(v2$calls$PatNr, v1$calls$PatNr)
if (length(only1) || length(only2)) {
  message("Patients only at ", visit_labels[1], ": ", length(only1), "; only at ", visit_labels[2], ": ",
          length(only2), " (shown as 'not measured' in the all-patients plot)")
}

meta <- bind_rows(v1$meta, v2$meta) %>% distinct(PatNr, .keep_all = TRUE)
pp <- full_join(v1$calls, v2$calls, by = "PatNr") %>%
  left_join(meta, by = "PatNr") %>%
  mutate(
    across(all_of(c(s1, s2)), ~ ifelse(is.na(.x), "not measured", .x)),
    allAgree_both = .data[[s1]] %in% c("alpha", "beta") & .data[[s2]] %in% c("alpha", "beta"),
    switched = allAgree_both & .data[[s1]] != .data[[s2]],
    transition = paste(.data[[s1]], "->", .data[[s2]])
  ) %>%
  relocate(PatNr, any_of(c("treatment", "arm", "stratum", "sex")), all_of(c(s1, s2)),
           allAgree_both, switched, transition)
write.csv(pp, file.path(output_dir, paste0("per_patient_calls_", mode, ".csv")), row.names = FALSE)

counts_visit <- bind_rows(
  count(pp, status = .data[[s1]]) %>% mutate(visit = visit_labels[1]),
  count(pp, status = .data[[s2]]) %>% mutate(visit = visit_labels[2])
) %>% select(visit, status, n)
write.csv(counts_visit, file.path(output_dir, paste0("allAgree_counts_per_visit_", mode, ".csv")), row.names = FALSE)

# ------------------------------------------------------------
# 2. Patients all-agree at both visits: transitions + McNemar
# ------------------------------------------------------------
both <- pp %>% filter(allAgree_both)
trans <- both %>%
  count(from = .data[[s1]], to = .data[[s2]]) %>%
  complete(from = c("alpha", "beta"), to = c("alpha", "beta"), fill = list(n = 0L))
write.csv(trans, file.path(output_dir, paste0("transitions_allAgree_both_", mode, ".csv")), row.names = FALSE)

n_ab <- trans$n[trans$from == "alpha" & trans$to == "beta"]
n_ba <- trans$n[trans$from == "beta" & trans$to == "alpha"]
n_sw <- n_ab + n_ba
mcnemar_p <- if (n_sw > 0) stats::binom.test(n_ab, n_sw, 0.5)$p.value else NA_real_

switched_tab <- both %>% filter(switched) %>%
  select(PatNr, any_of(c("treatment", "arm", "stratum", "sex")), transition,
         starts_with(paste0(visit_labels[1], "_")), starts_with(paste0(visit_labels[2], "_")))
write.csv(switched_tab, file.path(output_dir, paste0("switched_patients_", mode, ".csv")), row.names = FALSE)

message("\nAll-agree at both visits: ", nrow(both), " patients")
print(as.data.frame(trans), row.names = FALSE)
message("Switched: ", n_sw, " (alpha->beta ", n_ab, ", beta->alpha ", n_ba, "); exact McNemar p = ",
        if (is.na(mcnemar_p)) "NA" else signif(mcnemar_p, 3))

# ------------------------------------------------------------
# 3. Plots
# ------------------------------------------------------------
save_plot <- function(p, name, w = 7, h = 6) {
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("strata appear at multiple axes", conditionMessage(w))) invokeRestart("muffleWarning")
  })
  quiet(ggsave(file.path(output_dir, paste0(name, "_", mode, ".pdf")), p, width = w, height = h))
  quiet(ggsave(file.path(output_dir, paste0(name, "_", mode, ".png")), p, width = w, height = h, dpi = 300))
}

# Number of patients on each flow, written just right of the first axis where the flow leaves it.
# ggrepel (if installed) spreads the labels of thin, adjacent flows vertically so they do not overlap.
flow_label_layer <- function() {
  # anchor at the right edge of the first stratum (stratum width = 1/4)
  lab_aes <- aes(label = after_stat(ifelse(flow == "from", count, NA)), x = after_stat(x + 1 / 8))
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(stat = "flow", mapping = lab_aes, nudge_x = 0.08, direction = "y",
                             hjust = 0, size = 3.2, fontface = "bold", segment.size = 0.3,
                             segment.color = "grey30", min.segment.length = 0, box.padding = 0.15,
                             point.padding = 0, seed = 1, na.rm = TRUE)
  } else {
    geom_text(stat = "flow", mapping = lab_aes, nudge_x = 0.03, hjust = 0,
              size = 3.2, fontface = "bold", na.rm = TRUE)
  }
}

sankey <- function(df, from_col, to_col, levels_use) {
  agg <- df %>%
    count(from = .data[[from_col]], to = .data[[to_col]], across(any_of("facet"))) %>%
    mutate(from = factor(from, levels = levels_use), to = factor(to, levels = levels_use))
  ggplot(agg, aes(axis1 = from, axis2 = to, y = n)) +
    geom_alluvium(aes(fill = from), width = 1 / 4, alpha = 0.55, knot.pos = 0.4, color = "grey40",
                  linewidth = 0.2) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1 / 4, color = "grey20") +
    geom_text(stat = "stratum", aes(label = paste0(after_stat(stratum), "\nn = ", after_stat(count))),
              size = 3.3, lineheight = 0.9) +
    # number of patients on each flow, written where the flow leaves the first axis
    flow_label_layer() +
    scale_x_discrete(limits = visit_labels, expand = c(0.08, 0.08)) +
    scale_fill_manual(values = c(status_cols, "not measured" = "white"), name = paste0("Call at ", visit_labels[1]),
                      breaks = intersect(levels_use, unique(as.character(agg$from)))) +
    labs(x = NULL, y = "Patients") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}

# (a) main: all-agree at both visits
p_main <- sankey(both, s1, s2, c("alpha", "beta")) +
  labs(
    title = paste0("ROCK-ALS ", fluid_label, ": subtype call ", visit_labels[1], " -> ", visit_labels[2],
                   " (all models agree at both visits)"),
    subtitle = paste0(nrow(both), " of ", nrow(pp), " patients; ", n_sw, " switched (alpha->beta ", n_ab,
                      ", beta->alpha ", n_ba, ")\nexact McNemar p = ",
                      if (is.na(mcnemar_p)) "NA" else format(signif(mcnemar_p, 3)))
  ) +
  theme(plot.title = element_text(size = 12, face = "bold"))
save_plot(p_main, "sankey_allAgree_both")

# (b) all patients incl. "models disagree"
lv_all <- c(status_levels, if (length(only1) || length(only2)) "not measured")
trans_all <- pp %>% count(from = .data[[s1]], to = .data[[s2]])
write.csv(trans_all, file.path(output_dir, paste0("transitions_all_patients_", mode, ".csv")), row.names = FALSE)
p_all <- sankey(pp, s1, s2, lv_all) +
  labs(
    title = paste0("ROCK-ALS ", fluid_label, ": subtype call ", visit_labels[1], " -> ", visit_labels[2], " (all patients)"),
    subtitle = paste0(nrow(pp), " patients; alpha / beta = all models agree, grey = models disagree")
  ) +
  theme(plot.title = element_text(size = 12, face = "bold"))
save_plot(p_all, "sankey_all_patients", w = 7.5, h = 6.5)

# (c) all-agree at both visits, by treatment arm
arm_col <- intersect(c("arm", "treatment"), colnames(both))[1]
if (!is.na(arm_col)) {
  both_arm <- both %>% mutate(facet = .data[[arm_col]])
  arm_lab <- both_arm %>% group_by(facet) %>%
    summarise(lab = paste0(first(facet), " (n = ", n(), "; switched ", sum(switched), ")"), .groups = "drop")
  both_arm$facet <- factor(both_arm$facet, levels = arm_lab$facet, labels = arm_lab$lab)
  p_arm <- sankey(both_arm, s1, s2, c("alpha", "beta")) +
    facet_wrap(~ facet, scales = "free_y") +
    labs(title = paste0("All-agree at both visits, by ", arm_col)) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  save_plot(p_arm, "sankey_allAgree_both_by_arm", w = 10, h = 6)
  arm_tab <- both_arm %>% count(arm = .data[[arm_col]], transition)
  write.csv(arm_tab, file.path(output_dir, paste0("transitions_allAgree_both_by_arm_", mode, ".csv")),
            row.names = FALSE)
}

message("\nOutputs written to ", output_dir)
