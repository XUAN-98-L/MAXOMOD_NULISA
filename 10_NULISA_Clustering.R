# Unsupervised clustering on combined CNS + Immune NULISA CSF panels
# Clustering workflow adapted from MAXOMOD_CSF/Script/08_Clustering_subclusters.R
#
# Rscript Script/10_NULISA_Clustering.R \
#   -i CNS_immune/Results \
#   -o CNS_immune/Results/Combined_CNS_Immune_panels/Without_tears/10_NULISA_Clustering \
#   -e 9

suppressMessages(library(optparse))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(cluster))
suppressMessages(library(mclust))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggalluvial))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(RColorBrewer))

final_colours <- list(
  clustering = brewer.pal(n = 8, "Set2")[c(2, 3, 5)]
)
names(final_colours$clustering) <- c("alpha", "beta", "theta")

ha_colors <- list(
  sex = c(Female = "#E78AC3", Male = "#66C2A5"),
  onset = c(bulbar = "#FEE090", spinal = "#FDAE61"),
  progression_group = c(SP = "#1B9E77", IP = "#D95F02", FP = "#7570B3"),
  cluster = c(alpha = "#FC8D62", beta = "#8DA0CB", theta = "#A6D854")
)

.choose_test_two_groups <- function(y, g, mode = c("auto", "parametric", "nonparametric")) {
  mode <- match.arg(mode)
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) stop("choose_test_two_groups: g must have exactly 2 groups.")
  grp <- split(y, g)
  n1 <- sum(!is.na(grp[[1]]))
  n2 <- sum(!is.na(grp[[2]]))

  if (mode == "parametric") {
    lev <- tryCatch(car::leveneTest(y ~ g), error = function(e) NULL)
    if (!is.null(lev) && lev$`Pr(>F)`[1] > 0.05) {
      return(list(name = "Student's t-test", test = t.test(y ~ g, var.equal = TRUE)))
    }
    return(list(name = "Welch's t-test", test = t.test(y ~ g)))
  }
  if (mode == "nonparametric") {
    return(list(name = "Wilcoxon rank-sum", test = wilcox.test(y ~ g)))
  }

  if (n1 < 3 || n2 < 3) {
    return(list(name = "Wilcoxon rank-sum", test = wilcox.test(y ~ g)))
  }

  sw1 <- tryCatch(shapiro.test(grp[[1]]), error = function(e) NULL)
  sw2 <- tryCatch(shapiro.test(grp[[2]]), error = function(e) NULL)
  normal1 <- !is.null(sw1) && sw1$p.value > 0.05
  normal2 <- !is.null(sw2) && sw2$p.value > 0.05

  if (normal1 && normal2) {
    lev <- tryCatch(car::leveneTest(y ~ g), error = function(e) NULL)
    if (!is.null(lev) && lev$`Pr(>F)`[1] > 0.05) {
      list(name = "Student's t-test", test = t.test(y ~ g, var.equal = TRUE))
    } else {
      list(name = "Welch's t-test", test = t.test(y ~ g))
    }
  } else {
    list(name = "Wilcoxon rank-sum", test = wilcox.test(y ~ g))
  }
}

.format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 2))
  sprintf("%.3f", p)
}

plot_cat_by_k2 <- function(
    df, k2_col = "k2", var_col = "onset",
    ha_colors = NULL, output_dir = ".", file_name = NULL,
    linesize = 0.5, table_base_size = 12, return_plot_only = FALSE) {
  stopifnot(all(c(k2_col, var_col) %in% names(df)))
  d <- subset(df, !is.na(df[[k2_col]]) & !is.na(df[[var_col]]))
  d[[k2_col]] <- factor(d[[k2_col]])
  d[[var_col]] <- factor(d[[var_col]])

  tab <- table(k2 = d[[k2_col]], var = d[[var_col]])
  chi_try <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (any(chi_try$expected < 5)) {
    test_res <- fisher.test(tab)
    test_type <- "Fisher's exact test"
  } else {
    test_res <- chi_try
    test_type <- "Chi-squared test"
  }

  df_plot <- as.data.frame(tab)
  colnames(df_plot) <- c("cluster", "variable", "Freq")
  coeff <- data.frame(
    Method = test_type,
    `p-value` = sprintf("%.3f", test_res$p.value),
    check.names = FALSE
  )
  tbl_grob <- tableGrob(coeff, rows = NULL, theme = ttheme_minimal(base_size = table_base_size))

  if (!is.null(ha_colors) && var_col %in% names(ha_colors)) {
    pal <- ha_colors[[var_col]]
    pal <- pal[levels(d[[var_col]])]
  } else {
    nlev <- nlevels(d[[var_col]])
    if (nlev > 8) stop("Too many levels for default palette; please supply ha_colors.")
    pal <- brewer.pal(8, "Set2")[seq_len(nlev)]
  }

  ylab <- sprintf("proportions of %s", var_col)
  p <- ggplot(df_plot, aes(fill = variable, y = Freq, x = cluster)) +
    geom_bar(position = "fill", stat = "identity", color = "white", width = 0.7) +
    annotation_custom(tbl_grob, ymin = 0, ymax = 0.2) +
    scale_fill_manual(values = pal, name = var_col) +
    scale_y_continuous(labels = percent_format(accuracy = 10), expand = expansion(mult = c(0, 0.08))) +
    labs(x = k2_col, y = ylab) +
    theme_classic(base_size = 14) +
    theme(
      panel.background = element_rect(linewidth = linesize),
      plot.title = element_text(hjust = 0),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )

  if (!return_plot_only) {
    if (is.null(file_name)) file_name <- sprintf("k2_vs_%s.svg", var_col)
    save_path <- file.path(output_dir, file_name)
    ggsave(save_path, plot = p, width = 4.5, height = 4.2, dpi = 300)
  } else {
    save_path <- NULL
  }

  result <- list(
    plot = p, table = tab, proportions = prop.table(tab, margin = 1),
    test_type = test_type, test = test_res
  )
  if (!return_plot_only) result$save_path <- save_path
  result
}

plot_cont_by_k2 <- function(
    df,
    k2_col = "k2",
    var_col = "age",
    ha_colors = NULL,
    output_dir = ".",
    file_name = NULL,
    linesize = 0.5,
    table_base_size = 12,
    test_pref = c("auto", "parametric", "nonparametric"),
    return_plot_only = FALSE) {
  test_pref <- match.arg(test_pref)
  stopifnot(all(c(k2_col, var_col) %in% names(df)))
  d <- subset(df, !is.na(df[[k2_col]]) & !is.na(df[[var_col]]))
  d[[k2_col]] <- droplevels(factor(d[[k2_col]]))
  if (nlevels(d[[k2_col]]) != 2) stop("plot_cont_by_k2: k2 must have exactly 2 groups for this auto test.")
  d[[var_col]] <- suppressWarnings(as.numeric(d[[var_col]]))
  if (all(is.na(d[[var_col]]))) stop(sprintf("Variable '%s' is not numeric.", var_col))

  tst <- .choose_test_two_groups(d[[var_col]], d[[k2_col]], mode = test_pref)
  test_name <- tst$name
  pval <- tst$test$p.value

  coeff <- data.frame(
    Method = test_name,
    `p-value` = .format_p(pval),
    check.names = FALSE
  )
  tbl_grob <- tableGrob(coeff, rows = NULL, theme = ttheme_minimal(base_size = table_base_size))

  if (!is.null(ha_colors) && "cluster" %in% names(ha_colors)) {
    pal <- ha_colors[["cluster"]]
    pal <- pal[levels(d[[k2_col]])]
  } else {
    nlev <- nlevels(d[[k2_col]])
    pal <- brewer.pal(8, "Set2")[seq_len(nlev)]
  }

  ylab <- var_col
  p <- ggplot(d, aes(y = .data[[var_col]], x = .data[[k2_col]], fill = .data[[k2_col]])) +
    geom_violin(trim = FALSE, color = "black", fill = "white") +
    geom_dotplot(
      binaxis = "y", stackdir = "center", dotsize = 0.8,
      color = "black", aes(fill = .data[[k2_col]])
    ) +
    scale_fill_manual(values = pal, name = k2_col) +
    labs(x = k2_col, y = ylab) +
    theme_classic(base_size = 14) +
    theme(
      panel.background = element_rect(linewidth = linesize),
      plot.title = element_text(hjust = 0),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )

  y_min <- min(d[[var_col]], na.rm = TRUE)
  y_span <- diff(range(d[[var_col]], na.rm = TRUE))
  p <- p + annotation_custom(tbl_grob, ymin = y_min, ymax = y_min + 0.25 * y_span)

  if (!return_plot_only) {
    if (is.null(file_name)) file_name <- sprintf("k2_vs_%s.svg", var_col)
    save_path <- file.path(output_dir, file_name)
    ggsave(save_path, plot = p, width = 4.8, height = 4.2, dpi = 300)
  } else {
    save_path <- NULL
  }

  result <- list(plot = p, test_name = test_name, p_value = pval)
  if (!return_plot_only) result$save_path <- save_path
  result
}

plot_clinical_features_merged <- function(df, var_cols, ha_colors, output_dir, file_name, nrow = 1) {
  plot_list <- list()
  for (var_col in var_cols) {
    if (!var_col %in% names(df)) {
      warning(sprintf("Variable '%s' not found in data, skipping.", var_col))
      next
    }
    result <- tryCatch(
      if (is.character(df[[var_col]]) || is.factor(df[[var_col]])) {
        plot_cat_by_k2(
          df, k2_col = "k2", var_col = var_col, ha_colors = ha_colors,
          output_dir = output_dir, return_plot_only = TRUE
        )
      } else {
        plot_cont_by_k2(
          df, k2_col = "k2", var_col = var_col, test_pref = "auto",
          ha_colors = ha_colors, output_dir = output_dir, return_plot_only = TRUE
        )
      },
      error = function(e) {
        warning(sprintf("Could not plot '%s': %s", var_col, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(result)) plot_list[[var_col]] <- result$plot
  }
  if (length(plot_list) == 0L) return(invisible(NULL))

  ncol <- ceiling(length(plot_list) / nrow)
  merged_plot <- do.call(arrangeGrob, c(plot_list, ncol = ncol, nrow = nrow))
  save_path <- file.path(output_dir, file_name)
  ggsave(save_path, plot = merged_plot, width = 4.8 * ncol, height = 4.2 * nrow, dpi = 300)
  message("Merged clinical plot saved to: ", save_path)
}

compare_adjusted_ari <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 2L) return(NA_real_)
  adjustedRandIndex(factor(x[keep]), factor(y[keep]))
}

read_ms_cluster_labels <- function(path, disease = "als") {
  if (!file.exists(path)) stop("MS cluster file not found: ", path)
  ms_df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("Patient.ID" %in% names(ms_df)) "Patient.ID" else "Patient ID"
  if (!id_col %in% names(ms_df)) stop("MS cluster file must contain Patient.ID.")
  if (!all(c("k2", "k3") %in% names(ms_df))) stop("MS cluster file must contain k2 and k3 columns.")
  ms_df %>%
    mutate(
      ID = as.character(.data[[id_col]]),
      disease = as.character(.data$disease),
      k2_ms = as.character(k2),
      k3_ms = as.character(k3)
    ) %>%
    filter(.data$disease == disease) %>%
    select(ID, k2_ms, k3_ms)
}

plot_ms_crosstab <- function(tab_df, title, output_path) {
  p <- ggplot(tab_df, aes(x = MS, y = NULISA, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), color = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "#FC8D62") +
    labs(title = title, x = "MS-based cluster", y = "NULISA cluster") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0))
  ggsave(output_path, plot = p, width = 5, height = 4.2, dpi = 300)
  p
}

compare_nulisa_to_ms_clusters <- function(
    cluster_assignments,
    cluster_assignments_2,
    clin_df,
    ms_cluster_path,
    output_dir,
    disease_label = "als") {
  compare_dir <- file.path(output_dir, "MS_cluster_comparison")
  if (!dir.exists(compare_dir)) dir.create(compare_dir, recursive = TRUE)

  ms_labels <- read_ms_cluster_labels(ms_cluster_path, disease = disease_label)
  id_map <- clin_df %>%
    transmute(
      patid = as.character(label),
      ID = as.character(ID)
    )

  compare_df <- cluster_assignments %>%
    left_join(id_map, by = "patid") %>%
    left_join(ms_labels, by = "ID") %>%
    left_join(
      cluster_assignments_2 %>%
        transmute(
          patid = as.character(patid),
          kmeans_k2_named = as.character(`kmeans_k=2`),
          kmeans_k3_named = as.character(`kmeans_k=3`)
        ),
      by = "patid"
    )

  write.csv(compare_df, file.path(compare_dir, "nulisa_vs_ms_cluster_assignments.csv"), row.names = FALSE)

  k2_cols <- grep("_k=2$", colnames(cluster_assignments), value = TRUE)
  k3_cols <- grep("_k=3$", colnames(cluster_assignments), value = TRUE)

  ari_k2 <- vapply(
    c(k2_cols, "kmeans_k2_named"),
    function(col) compare_adjusted_ari(compare_df[[col]], compare_df$k2_ms),
    numeric(1)
  )
  names(ari_k2) <- c(k2_cols, "kmeans_k2_named")

  ari_k3 <- vapply(
    c(k3_cols, "kmeans_k3_named"),
    function(col) compare_adjusted_ari(compare_df[[col]], compare_df$k3_ms),
    numeric(1)
  )
  names(ari_k3) <- c(k3_cols, "kmeans_k3_named")

  ari_summary <- data.frame(
    comparison = c(
      paste0(names(ari_k2), "_vs_MS_k2"),
      paste0(names(ari_k3), "_vs_MS_k3")
    ),
    adjusted_rand_index = c(unname(ari_k2), unname(ari_k3)),
    n_samples = c(
      rep(sum(!is.na(compare_df$k2_ms)), length(ari_k2)),
      rep(sum(!is.na(compare_df$k3_ms)), length(ari_k3))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(ari_summary, file.path(compare_dir, "adjusted_rand_index_summary.csv"), row.names = FALSE)

  message("Adjusted Rand index vs MS-based clusters (n = ", sum(!is.na(compare_df$k2_ms)), " ALS samples):")
  print(ari_summary)

  tab_k2 <- as.data.frame(table(
    NULISA = compare_df$kmeans_k2_named,
    MS = compare_df$k2_ms,
    useNA = "no"
  ))
  colnames(tab_k2) <- c("NULISA", "MS", "Freq")
  write.csv(tab_k2, file.path(compare_dir, "crosstab_kmeans_k2_vs_MS_k2.csv"), row.names = FALSE)
  plot_ms_crosstab(
    tab_k2,
    sprintf("NULISA k-means k=2 vs MS k2 (ARI = %.3f)", ari_k2["kmeans_k2_named"]),
    file.path(compare_dir, "crosstab_kmeans_k2_vs_MS_k2.svg")
  )

  tab_k3 <- as.data.frame(table(
    NULISA = compare_df$kmeans_k3_named,
    MS = compare_df$k3_ms,
    useNA = "no"
  ))
  colnames(tab_k3) <- c("NULISA", "MS", "Freq")
  write.csv(tab_k3, file.path(compare_dir, "crosstab_kmeans_k3_vs_MS_k3.csv"), row.names = FALSE)
  plot_ms_crosstab(
    tab_k3,
    sprintf("NULISA k-means k=3 vs MS k3 (ARI = %.3f)", ari_k3["kmeans_k3_named"]),
    file.path(compare_dir, "crosstab_kmeans_k3_vs_MS_k3.svg")
  )

  ari_plot_df <- rbind(
    data.frame(k = 2L, method = names(ari_k2), ARI = unname(ari_k2), stringsAsFactors = FALSE),
    data.frame(k = 3L, method = names(ari_k3), ARI = unname(ari_k3), stringsAsFactors = FALSE)
  )
  ari_plot_df$method <- gsub("_named$", "", ari_plot_df$method)

  p_ari <- ggplot(ari_plot_df, aes(x = reorder(method, ARI), y = ARI, fill = factor(k))) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("2" = "#8DA0CB", "3" = "#FC8D62"), name = "k") +
    labs(
      title = "Adjusted Rand index: NULISA vs MS-based clusters",
      x = NULL,
      y = "Adjusted Rand index"
    ) +
    theme_classic(base_size = 14)
  ggsave(file.path(compare_dir, "adjusted_rand_index_barplot.svg"), plot = p_ari, width = 7, height = 5, dpi = 300)

  invisible(list(compare_df = compare_df, ari_summary = ari_summary))
}

calc_SS <- function(df) sum(as.matrix(dist(df)^2)) / (2 * nrow(df))

calc_TWSS <- function(df, clusters) {
  number_clusters <- length(levels(as.factor(clusters)))
  sum_of_squares <- vapply(seq_len(number_clusters), function(i) {
    calc_SS(df[clusters == i, , drop = FALSE])
  }, numeric(1))
  sum(sum_of_squares)
}

BIC2 <- function(df, clusters) {
  m <- ncol(df)
  n <- nrow(df)
  k <- length(levels(as.factor(clusters)))
  D <- calc_TWSS(df, clusters)
  data.frame(AIC = D + 2 * m * k, BIC = D + log(n) * m * k)
}

perform_clustering <- function(assay_data, seed, clin_labels) {
  set.seed(seed)

  silhouette_scores <- TWSS_scores <- AIC_scores <- BIC_scores <- as.data.frame(matrix())
  rownames(silhouette_scores) <- rownames(TWSS_scores) <- rownames(AIC_scores) <- rownames(BIC_scores) <- "hclust"

  cluster_assignments <- data.frame(patid = clin_labels, stringsAsFactors = FALSE)
  dist_mat <- dist(assay_data, method = "euclidean")

  for (i in 2:10) {
    title <- paste0("hclust_k=", i)
    cl <- hclust(dist_mat, method = "ward.D")
    cluster_assignments[[title]] <- cutree(cl, k = i)
    ss <- silhouette(cluster_assignments[[title]], dist_mat)
    silhouette_scores["hclust", as.character(i)] <- mean(ss[, 3])
    TWSS_scores["hclust", as.character(i)] <- calc_TWSS(assay_data, cluster_assignments[[title]])
    bic_vals <- BIC2(assay_data, cluster_assignments[[title]])
    AIC_scores["hclust", as.character(i)] <- bic_vals$AIC
    BIC_scores["hclust", as.character(i)] <- bic_vals$BIC

    title <- paste0("mclust_k=", i)
    cl <- Mclust(assay_data, G = i)
    cluster_assignments[[title]] <- cl$classification
    ss <- silhouette(cluster_assignments[[title]], dist_mat)
    silhouette_scores["mclust", as.character(i)] <- mean(ss[, 3])
    TWSS_scores["mclust", as.character(i)] <- calc_TWSS(assay_data, cluster_assignments[[title]])
    bic_vals <- BIC2(assay_data, cluster_assignments[[title]])
    AIC_scores["mclust", as.character(i)] <- bic_vals$AIC
    BIC_scores["mclust", as.character(i)] <- bic_vals$BIC

    title <- paste0("kmeans_k=", i)
    cl <- kmeans(assay_data, centers = i, nstart = 25, iter.max = 50)
    cluster_assignments[[title]] <- cl$cluster
    ss <- silhouette(cluster_assignments[[title]], dist_mat)
    silhouette_scores["kmeans", as.character(i)] <- mean(ss[, 3])
    TWSS_scores["kmeans", as.character(i)] <- calc_TWSS(assay_data, cluster_assignments[[title]])
    bic_vals <- BIC2(assay_data, cluster_assignments[[title]])
    AIC_scores["kmeans", as.character(i)] <- bic_vals$AIC
    BIC_scores["kmeans", as.character(i)] <- bic_vals$BIC

    title <- paste0("pam_k=", i)
    cl <- pam(assay_data, k = i)
    cluster_assignments[[title]] <- cl$clustering
    ss <- silhouette(cluster_assignments[[title]], dist_mat)
    silhouette_scores["pam", as.character(i)] <- mean(ss[, 3])
    TWSS_scores["pam", as.character(i)] <- calc_TWSS(assay_data, cluster_assignments[[title]])
    bic_vals <- BIC2(assay_data, cluster_assignments[[title]])
    AIC_scores["pam", as.character(i)] <- bic_vals$AIC
    BIC_scores["pam", as.character(i)] <- bic_vals$BIC
  }

  cluster_cols <- 2:ncol(cluster_assignments)
  cluster_assignments[, cluster_cols] <- cluster_assignments[, cluster_cols] - 1L

  list(
    cluster_assignments = cluster_assignments,
    AIC_scores = AIC_scores,
    BIC_scores = BIC_scores,
    silhouette_scores = silhouette_scores,
    TWSS_scores = TWSS_scores
  )
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "CNS_immune/Results"),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "CNS_immune/Results/Combined_CNS_Immune_panels/Without_tears/10_NULISA_Clustering"
  ),
  make_option(c("-e", "--seed"), type = "integer", default = 9L),
  make_option(
    c("-d", "--disease"),
    type = "logical",
    default = TRUE,
    help = "Cluster ALS patients (TRUE) or controls (FALSE)."
  ),
  make_option(c("-r", "--reverse"), type = "logical", default = FALSE),
  make_option(
    c("-c", "--ms_cluster"),
    type = "character",
    default = "Data/clinical_data_with_cluster_VC.csv",
    help = "MS-based cluster labels (Patient.ID, k2, k3) for adjusted Rand index comparison."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_panel_paths <- function(input_arg) {
  if (grepl(",", input_arg, fixed = TRUE)) {
    parts <- trimws(strsplit(input_arg, ",", fixed = TRUE)[[1]])
    if (length(parts) != 3L) {
      stop("--input must be a panel base directory or three comma-separated file paths.")
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

npq_to_wide <- function(npq_path, panel_prefix) {
  read_excel(npq_path) %>%
    filter(SampleMatrixType == "CSF", SampleType == "Sample") %>%
    mutate(Feature = paste0(panel_prefix, "_", Target)) %>%
    select(SampleName, Feature, NPQ) %>%
    pivot_wider(names_from = Feature, values_from = NPQ) %>%
    as.data.frame()
}

scale_manual <- function(df) {
  as.data.frame(apply(df, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) 0 else (x - rng[1]) / diff(rng)
  }))
}

panel_paths <- parse_panel_paths(opt$input)
output_dir <- opt$output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Loading CNS + Immune NPQ (CSF samples)...")
npq_wide_cns <- npq_to_wide(panel_paths$cns, "CNS")
npq_wide_immune <- npq_to_wide(panel_paths$immune, "Immune")
common_samples <- intersect(npq_wide_cns$SampleName, npq_wide_immune$SampleName)

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

use_features <- colnames(npq_wide)
message("Using all combined NULISA proteins: ", length(use_features))

clin <- read_excel(panel_paths$metadata) %>%
  filter(Material == "CSF") %>%
  mutate(
    label = as.character(Tube_ID),
    ID = as.character(ID),
    disease = as.character(disease),
    sex = as.factor(sex),
    onset = as.factor(onset),
    limb = as.factor(limb),
    progression_group = as.factor(progression_group),
    across(
      any_of(c("age", "Nfl", "pNFh", "progression_rate", "slow_vital_capacity", "age_at_onset", "disease_duration")),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  filter(Tube_ID %in% rownames(npq_wide))

if (isTRUE(opt$disease)) {
  clin <- clin %>% filter(disease == "als")
  message("Clustering ALS patients only: n = ", nrow(clin))
} else {
  clin <- clin %>% filter(disease == "ctrl")
  message("Clustering controls only: n = ", nrow(clin))
}

if (nrow(clin) < 4L) {
  stop("Need at least 4 samples for clustering after filtering.")
}

assay_df <- npq_wide[clin$Tube_ID, use_features, drop = FALSE]
rownames(assay_df) <- clin$label
assay_df <- assay_df[, colSums(is.na(assay_df)) == 0, drop = FALSE]
var_ok <- apply(assay_df, 2, function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  !is.na(s) && s > 0
})
assay_df <- assay_df[, var_ok, drop = FALSE]

if (ncol(assay_df) < 2L) {
  stop("Need at least 2 non-constant features for clustering.")
}

assay_scaled <- scale_manual(assay_df)
assay <- as.data.frame(assay_scaled)
rownames(assay) <- rownames(assay_df)

message("Clustering matrix: ", nrow(assay), " samples x ", ncol(assay), " features")
write.csv(assay_df, file.path(output_dir, "clustering_input_unscaled.csv"))
write.csv(assay, file.path(output_dir, "clustering_input_scaled.csv"))

set.seed(opt$seed)
cluster_assignments_result <- perform_clustering(assay, seed = opt$seed, clin_labels = clin$label)
cluster_assignments <- cluster_assignments_result$cluster_assignments
AIC_scores <- cluster_assignments_result$AIC_scores
BIC_scores <- cluster_assignments_result$BIC_scores
silhouette_scores <- cluster_assignments_result$silhouette_scores
TWSS_scores <- cluster_assignments_result$TWSS_scores

saveRDS(cluster_assignments, file.path(output_dir, "cluster_assignments.rds"))
write.csv(cluster_assignments, file.path(output_dir, "cluster_assignments.csv"), row.names = FALSE)

scores <- list(AIC = AIC_scores, BIC = BIC_scores, silhouette = silhouette_scores, TWSS = TWSS_scores)
plots <- vector("list", length(scores))

for (i in seq_along(scores)) {
  score_df <- scores[[i]]
  colnames(score_df) <- as.character(2:10)
  score_df$method <- rownames(score_df)
  melt_df <- reshape2::melt(score_df, id.vars = "method")
  plots[[i]] <- ggplot(melt_df, aes(x = variable, y = value, group = method, colour = method)) +
    geom_line() +
    geom_point() +
    ggtitle(names(scores)[i]) +
    theme_few()
}
names(plots) <- names(scores)

ggarrange(plotlist = plots, labels = LETTERS[seq_along(plots)], ncol = 2, nrow = 2)
ggsave(file.path(output_dir, "fit_scores.pdf"), width = 11, height = 8, units = "in")

number_of_clusters <- as.character(2:5)
plots <- vector("list", length(number_of_clusters))
for (i in seq_along(number_of_clusters)) {
  data <- cluster_assignments[, c(1, grep(number_of_clusters[i], colnames(cluster_assignments)))]
  data <- reshape2::melt(data)
  colnames(data) <- c("patid", "method", "cluster")
  data <- as.data.frame(lapply(data, as.factor))
  plots[[i]] <- ggplot(
    data,
    aes(x = method, stratum = cluster, alluvium = patid, fill = cluster, label = cluster)
  ) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom") +
    ggtitle(paste0("Cluster agreement across methods (k=", number_of_clusters[i], ")")) +
    theme_few()
}
ggarrange(plotlist = plots, labels = seq_along(plots), ncol = 2, nrow = 2)
ggsave(file.path(output_dir, "Sankey_plots_across_methods.pdf"), width = 11 * 2, height = 8 * 2, units = "in")

methods <- c("kmeans", "hclust", "mclust", "pam")
plots <- vector("list", length(methods))
for (i in seq_along(methods)) {
  data <- cluster_assignments[, c("patid", paste0(methods[i], "_k=", number_of_clusters))]
  data <- reshape2::melt(data, id = "patid")
  colnames(data) <- c("patid", "method", "cluster")
  data <- as.data.frame(lapply(data, as.factor))
  plots[[i]] <- ggplot(
    data,
    aes(x = method, stratum = cluster, alluvium = patid, fill = cluster, label = cluster)
  ) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
    geom_stratum() +
    theme(legend.position = "bottom") +
    ggtitle(paste0("Cluster stability within ", methods[i])) +
    theme_few() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}
ggarrange(plotlist = plots, labels = seq_along(plots), ncol = 2, nrow = 2)
ggsave(file.path(output_dir, "Sankey_plots_within_method.pdf"), width = 11 * 2, height = 8 * 2, units = "in")

cluster_assignments_2 <- cluster_assignments[, c("patid", "kmeans_k=2", "kmeans_k=3")]
cluster_counts <- table(cluster_assignments_2$`kmeans_k=2`)
larger_cluster <- names(cluster_counts)[which.max(cluster_counts)]
cluster_assignments_2$`kmeans_k=2` <- as.factor(cluster_assignments_2$`kmeans_k=2`)
if (isTRUE(opt$reverse)) {
  levels(cluster_assignments_2$`kmeans_k=2`) <- ifelse(
    levels(cluster_assignments_2$`kmeans_k=2`) == larger_cluster, "beta", "alpha"
  )
} else {
  levels(cluster_assignments_2$`kmeans_k=2`) <- ifelse(
    levels(cluster_assignments_2$`kmeans_k=2`) == larger_cluster, "alpha", "beta"
  )
}
cluster_assignments_2$`kmeans_k=2` <- ordered(cluster_assignments_2$`kmeans_k=2`, levels = c("alpha", "beta", "theta"))

message("kmeans k=2:")
print(table(cluster_assignments_2$`kmeans_k=2`))

adjusted_clusters <- cluster_assignments_2 %>%
  group_by(`kmeans_k=3`) %>%
  summarize(
    majority_alpha = sum(`kmeans_k=2` == "alpha"),
    majority_beta = sum(`kmeans_k=2` == "beta"),
    total = n(),
    .groups = "drop"
  ) %>%
  mutate(
    alpha_beta = majority_alpha / total - majority_beta / total,
    kmeans_k3_adjusted = case_when(
      alpha_beta == max(alpha_beta) ~ "alpha",
      alpha_beta == min(alpha_beta) ~ "beta",
      TRUE ~ "theta"
    )
  ) %>%
  select(`kmeans_k=3`, kmeans_k3_adjusted)

cluster_assignments_2 <- cluster_assignments_2 %>%
  left_join(adjusted_clusters, by = "kmeans_k=3") %>%
  select(-`kmeans_k=3`) %>%
  rename(`kmeans_k=3` = kmeans_k3_adjusted)
cluster_assignments_2$`kmeans_k=3` <- ordered(
  as.factor(cluster_assignments_2$`kmeans_k=3`),
  levels = c("alpha", "beta", "theta")
)

write.csv(cluster_assignments_2, file.path(output_dir, "cluster_assignments_2.csv"), row.names = FALSE)

message("kmeans k=3:")
print(table(cluster_assignments_2$`kmeans_k=3`))

if (isTRUE(opt$disease)) {
  message("Comparing NULISA clusters with MS-based clusters (adjusted Rand index)...")
  compare_nulisa_to_ms_clusters(
    cluster_assignments = cluster_assignments,
    cluster_assignments_2 = cluster_assignments_2,
    clin_df = clin,
    ms_cluster_path = opt$ms_cluster,
    output_dir = output_dir,
    disease_label = "als"
  )
}

data23 <- cluster_assignments_2 %>%
  pivot_longer(cols = c(`kmeans_k=2`, `kmeans_k=3`), names_to = "k", values_to = "cluster")
ggplot(
  data23,
  aes(x = k, stratum = cluster, alluvium = patid, fill = cluster, label = cluster)
) +
  scale_fill_manual(values = final_colours$clustering) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("K-means Sankey: k=2 to k=3") +
  theme_few()
ggsave(file.path(output_dir, "Sankey_plot_kmeans_2_3_patients.pdf"), width = 5.5, height = 4, units = "in")

changed_cluster <- cluster_assignments_2 %>%
  filter(`kmeans_k=2` == "alpha", `kmeans_k=3` == "beta")
if (nrow(changed_cluster) > 0L) {
  message("Patients switching alpha -> beta between k=2 and k=3:")
  print(changed_cluster)
}

clinical_dir <- file.path(output_dir, "Clinical_features")
if (!dir.exists(clinical_dir)) dir.create(clinical_dir, recursive = TRUE)

clin_df <- as.data.frame(clin[match(cluster_assignments_2$patid, clin$label), , drop = FALSE])
rownames(clin_df) <- clin_df$label
clin_df$k2 <- cluster_assignments_2$`kmeans_k=2`

if (isTRUE(opt$disease)) {
  clin_vars <- c("sex", "progression_group", "onset", "age", "age_at_onset", "Nfl", "pNFh")
} else {
  clin_vars <- c("sex", "age", "Nfl")
}
clin_vars <- intersect(clin_vars, names(clin_df))

message("Clinical features vs k-means k=2 clusters...")
for (var_col in clin_vars) {
  tryCatch(
    if (is.character(clin_df[[var_col]]) || is.factor(clin_df[[var_col]])) {
      plot_cat_by_k2(
        clin_df, k2_col = "k2", var_col = var_col,
        ha_colors = ha_colors, output_dir = clinical_dir
      )
    } else {
      plot_cont_by_k2(
        clin_df, k2_col = "k2", var_col = var_col, test_pref = "auto",
        ha_colors = ha_colors, output_dir = clinical_dir
      )
    },
    error = function(e) warning(sprintf("Could not plot '%s': %s", var_col, conditionMessage(e)))
  )
}

if (isTRUE(opt$disease)) {
  cont_vars <- intersect(c("age", "age_at_onset", "Nfl", "pNFh"), clin_vars)
  cat_vars <- intersect(c("sex", "progression_group", "onset"), clin_vars)
  if (length(cont_vars) > 1L) {
    plot_clinical_features_merged(
      clin_df, cont_vars, ha_colors, clinical_dir,
      "k2_vs_merged_age_age_at_onset_Nfl_pNFh.svg", nrow = 1
    )
  }
  if (length(cat_vars) > 1L) {
    plot_clinical_features_merged(
      clin_df, cat_vars, ha_colors, clinical_dir,
      "k2_vs_merged_progression_group_sex_onset.svg", nrow = 1
    )
  }
}

if (isTRUE(opt$disease) && "age_at_onset" %in% names(clin_df)) {
  mean_age_at_onset <- as.data.frame(
    tapply(clin_df$age_at_onset, cluster_assignments_2$`kmeans_k=3`, mean, na.rm = TRUE)
  )
  colnames(mean_age_at_onset)[1] <- "mean_age_at_onset"
  write.table(
    mean_age_at_onset,
    file.path(output_dir, "mean_age_at_onset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    quote = FALSE
  )
}

message("Done. Results written to: ", output_dir)
