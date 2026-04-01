# Trajectory of protein levels across progression rate
######### Functions#########################
# set colors
my_palette = c("#66C2A5","#E78AC3", "#ffc773",
                "#1B9E77","#7570B3","#FEE090",
                "#A6D854","#8DA0CB","#FC8D62",
                "#8FBC94","#b0a4e3","#ffa631",
                "#0aa344","#e4c6d0","#ffa400",
                "#519a73","#4b5cc4","#eedeb0",
                "#549688","#ffb3a7","#b35c44",
                "#7fecad","#a1afc9","#a78e44",
                "#519a73","#2e4e7e","#955539")
# Analogue of Monocle::genSmoothCurves + plot_pseudotime_heatmap: fit a smooth
# NPQ ~ progression curve per protein, evaluated on an evenly spaced progression grid.
gen_smooth_curves_npq <- function(df, prog_col, value_col, id_col, n_grid = 100L) {
  prog_min <- min(df[[prog_col]], na.rm = TRUE)
  prog_max <- max(df[[prog_col]], na.rm = TRUE)
  if (!is.finite(prog_min) || !is.finite(prog_max) || prog_min >= prog_max) {
    stop("Progression column must span at least two distinct finite values.")
  }
  grid_x <- seq(prog_min, prog_max, length.out = n_grid)
  targets <- unique(as.character(df[[id_col]]))
  out <- matrix(NA_real_, nrow = length(targets), ncol = n_grid,
                dimnames = list(targets, NULL))
  for (i in seq_along(targets)) {
    tg <- targets[i]
    sub <- df[as.character(df[[id_col]]) == tg, c(prog_col, value_col), drop = FALSE]
    names(sub) <- c("px", "py")
    sub <- sub[is.finite(sub$px) & is.finite(sub$py), , drop = FALSE]
    sub <- sub[order(sub$px), , drop = FALSE]
    sub <- stats::aggregate(py ~ px, sub, mean)
    px <- sub$px
    py <- sub$py
    n <- length(px)
    if (n < 1L) {
      out[i, ] <- NA_real_
    } else if (n == 1L) {
      out[i, ] <- py[1]
    } else if (n == 2L) {
      out[i, ] <- stats::approx(px, py, xout = grid_x, rule = 2)$y
    } else {
      # GCV can fail ("smoothing parameter value too small") with few unique x
      # or nearly degenerate layouts; fall back to linear interpolation.
      sp <- tryCatch(
        stats::smooth.spline(px, py, cv = FALSE),
        error = function(e) NULL
      )
      if (is.null(sp)) {
        out[i, ] <- stats::approx(px, py, xout = grid_x, rule = 2)$y
      } else {
        out[i, ] <- stats::predict(sp, x = grid_x)$y
      }
    }
  }
  attr(out, "progression_grid") <- grid_x
  out
}

###############################################
### Libraries
suppressMessages(library("optparse"))
suppressMessages(library("readxl"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("grid"))
suppressMessages(library("ggplot2"))
###############################################
#=================Parse command line arguments=================
option_list = list(
  make_option(c("-i", "--input"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx", help="Input protein_data_IDs.xlsx file"),
  make_option(c("-o", "--output"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/06_Trajectory_progressionRate", help="Output trajectory progression rate path"),
  # column name in the all_participants_IDs.xlsx file
  make_option(c("-c", "--column_name"), type="character", default="progression_rate", help="Column name in the all_participants_IDs.xlsx file, should be a numeric column"),
  make_option(c("-m","--matrix_type"), type="character", default="CSF", help="Choose from: CSF, PLASMA, SERUM, TEARS"),
  make_option(c("--n_row_clusters"), type="integer", default=4, help="Number of row blocks after clustering (0 = no split; k must be <= number of proteins)"),
  make_option(c("--n_smooth_grid"), type="integer", default=20, help="Number of evenly spaced progression points for smooth curves (Monocle-style; 0 = use raw per-sample columns, no smoothing)"),
  make_option(c("--heatmap_width_cm"), type="double", default=5, help="Target width of the heatmap body in cm (maps to pheatmap cellwidth; narrower column like reference figures)"),
  make_option(c("--pdf_width"), type="double", default=6, help="PDF page width in inches"),
  make_option(c("--pdf_height"), type="double", default=8, help="PDF page height in inches"),
  make_option(
    c("--progression_sd_k"),
    type = "double",
    default = 0,
    help = "If >0, exclude samples whose mean progression is outside cohort mean ± k*SD (0 = off; e.g. 3 for 3 SD). Applied before matrices and plots."
  )
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#=================options setting=================
if (is.null(opt$input)) {
  print("NO INPUT PROTEIN DATA IDS FILE SUPPLIED, EXITING!")
  stop("Please provide the input protein data IDs file path!")
} else {
  input_file <- opt$input
  protein_data_IDs <- read_excel(input_file)
}

if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

if (is.null(opt$column_name)) {
  print("NO COLUMN NAME SUPPLIED, EXITING!")
  stop("Please provide the column name!")
} else {
  column_name <- opt$column_name
  # check if it's a numeric column
  protein_data_IDs[[column_name]] = as.numeric(protein_data_IDs[[column_name]])
}

if (is.null(opt$matrix_type)) {
  print("NO MATRIX TYPE SUPPLIED, EXITING!")
  stop("Please provide the matrix type!")
} else {
  matrix_type <- toupper(opt$matrix_type)
  if (!matrix_type %in% c("CSF", "PLASMA", "SERUM", "TEARS")) {
    print("INVALID MATRIX TYPE, EXITING!")
    stop("Please provide a valid matrix type!")
  }
}

n_row_clusters <- as.integer(opt$n_row_clusters)
if (is.na(n_row_clusters) || n_row_clusters < 0L) {
  stop("--n_row_clusters must be a non-negative integer.")
}
n_smooth_grid <- as.integer(opt$n_smooth_grid)
if (is.na(n_smooth_grid) || n_smooth_grid < 0L) {
  stop("--n_smooth_grid must be a non-negative integer.")
}
if (n_smooth_grid == 1L) {
  stop("--n_smooth_grid must be 0 (raw per-sample columns) or >= 2 (smooth grid).")
}
heatmap_width_cm <- as.numeric(opt$heatmap_width_cm)
if (is.na(heatmap_width_cm) || heatmap_width_cm <= 0) {
  stop("--heatmap_width_cm must be a positive number.")
}
pdf_width_in <- as.numeric(opt$pdf_width)
pdf_height_in <- as.numeric(opt$pdf_height)
if (is.na(pdf_width_in) || pdf_width_in <= 0 || is.na(pdf_height_in) || pdf_height_in <= 0) {
  stop("--pdf_width and --pdf_height must be positive numbers (inches).")
}
progression_sd_k <- as.numeric(opt$progression_sd_k)
if (is.na(progression_sd_k) || progression_sd_k < 0) {
  stop("--progression_sd_k must be a non-negative number (0 = no SD-based sample exclusion).")
}
###############################################
# filter the protein_data_IDs dataframe by the matrix type
protein_data_IDs <- protein_data_IDs %>%
  filter(SampleMatrixType == matrix_type)

###############################################
# filter protein_data_IDs when column_name is not NA
if (!column_name %in% names(protein_data_IDs)) {
  stop("Column `", column_name, "` not found in input file. Available: ",
       paste(names(protein_data_IDs), collapse = ", "))
}

protein_data_IDs <- protein_data_IDs %>%
  filter(!is.na(!!sym(column_name)))

if (!"NPQ" %in% names(protein_data_IDs)) {
  stop("Input file must contain an `NPQ` column.")
}

protein_data_IDs <- protein_data_IDs %>%
  filter(!is.na(NPQ))

if (nrow(protein_data_IDs) == 0L) {
  stop("No rows left after filtering for non-NA `", column_name, "` and NPQ.")
}

# --- Samples ordered by progression (columns = increasing progression_rate) ---
sample_order_df <- protein_data_IDs %>%
  group_by(SampleName) %>%
  summarise(!!sym(column_name) := mean(!!sym(column_name), na.rm = TRUE), .groups = "drop") %>%
  arrange(!!sym(column_name))

sample_order <- sample_order_df$SampleName
prog_per_sample <- stats::setNames(sample_order_df[[column_name]], sample_order_df$SampleName)

# Optional: exclude samples whose mean progression lies outside cohort mean ± k*SD (before all plots)
if (progression_sd_k > 0) {
  pv <- sample_order_df[[column_name]]
  mu_pr <- mean(pv, na.rm = TRUE)
  sigma_pr <- stats::sd(pv, na.rm = TRUE)
  if (!is.finite(mu_pr) || !is.finite(sigma_pr) || sigma_pr < 1e-12) {
    stop(
      "Cannot apply --progression_sd_k: need finite mean and non-zero SD of per-sample ",
      column_name, " (n=", nrow(sample_order_df), ")."
    )
  }
  z_pr <- (pv - mu_pr) / sigma_pr
  keep_sd <- is.finite(z_pr) & abs(z_pr) <= progression_sd_k
  if (sum(keep_sd) < nrow(sample_order_df)) {
    excluded_df <- sample_order_df[!keep_sd, , drop = FALSE]
    excluded_df$z_vs_cohort <- z_pr[!keep_sd]
    excluded_df$cohort_mean <- mu_pr
    excluded_df$cohort_sd <- sigma_pr
    excl_path <- file.path(
      output_dir,
      paste0("excluded_samples_", column_name, "_outside_", progression_sd_k, "sd.csv")
    )
    utils::write.csv(excluded_df, file = excl_path, row.names = FALSE)
    message(
      "Excluded ", nrow(excluded_df), " sample(s) outside mean ± ", progression_sd_k,
      " SD of ", column_name, ". Wrote ", normalizePath(excl_path, winslash = "/", mustWork = FALSE)
    )
    sample_order_df <- sample_order_df[keep_sd, , drop = FALSE]
    protein_data_IDs <- protein_data_IDs %>%
      dplyr::filter(SampleName %in% sample_order_df$SampleName)
    sample_order <- sample_order_df$SampleName
    prog_per_sample <- stats::setNames(sample_order_df[[column_name]], sample_order_df$SampleName)
  }
  if (nrow(sample_order_df) < 2L) {
    stop("After --progression_sd_k filter need at least 2 samples; left ", nrow(sample_order_df), ".")
  }
  if (all(keep_sd)) {
    message(
      "Progression SD filter (k=", progression_sd_k, "): all ", nrow(sample_order_df),
      " sample(s) within mean ± ", progression_sd_k, " SD of ", column_name, "."
    )
  }
}

# --- Matrix: rows = proteins (Target), columns = samples ---
mat_wide <- protein_data_IDs %>%
  select(Target, SampleName, NPQ) %>%
  group_by(Target, SampleName) %>%
  summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = SampleName, values_from = NPQ)
mat <- as.matrix(mat_wide[, -1, drop = FALSE])
rownames(mat) <- mat_wide$Target

cols_present <- intersect(sample_order, colnames(mat))
if (length(cols_present) < 2L) {
  stop("Need at least 2 samples with NPQ after filtering; found ", length(cols_present), ".")
}
mat <- mat[, cols_present, drop = FALSE]

if (any(is.na(mat))) {
  mat <- mat[rowSums(is.na(mat)) == 0L, , drop = FALSE]
}
if (nrow(mat) < 1L) {
  stop("No proteins with complete NPQ across all samples after filtering.")
}

# Smooth NPQ along progression (cf. Monocle::genSmoothCurves in plot_pseudotime_heatmap)
if (n_smooth_grid >= 2L) {
  df_smooth <- protein_data_IDs %>%
    filter(Target %in% rownames(mat)) %>%
    select(Target, NPQ, !!sym(column_name))
  mat_smooth <- gen_smooth_curves_npq(
    df_smooth,
    prog_col = column_name,
    value_col = "NPQ",
    id_col = "Target",
    n_grid = n_smooth_grid
  )
  mat_plot <- mat_smooth[rownames(mat), , drop = FALSE]
  prog_grid <- attr(mat_smooth, "progression_grid")
} else {
  mat_plot <- mat
  prog_grid <- unname(prog_per_sample[colnames(mat)])
}

# Per-sample NPQ (columns = samples in progression order). Used for raw heatmap only:
# x-axis = raw progression_rate per sample (no smooth grid / binning).
mat_raw <- mat[rownames(mat_plot), , drop = FALSE]

# Row-wise z-score along trajectory columns (below/above average along progression)
mat_z <- t(scale(t(mat_plot)))
mat_z[is.nan(mat_z)] <- 0

# Diverging color scale: limits from min/max of scaled NPQ (no fixed clip)
lo <- ceiling(min(mat_z, na.rm = TRUE))
hi <- floor(max(mat_z, na.rm = TRUE))
if (lo < 0 && hi > 0) {
  col_fun <- circlize::colorRamp2(c(lo, 0, hi), c("#79B6E4", "white", "#E7352E"))
  leg_at <- c(lo, 0, hi)
} else {
  lim <- max(abs(c(lo, hi)), na.rm = TRUE)
  if (!is.finite(lim) || lim < 1e-15) lim <- 1
  col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#79B6E4", "white", "#E7352E"))
  leg_at <- c(-lim, 0, lim)
}
leg_labs <- format(signif(leg_at, 4), trim = TRUE)
leg_labs[abs(leg_at) < 1e-12] <- "0"

# Raw NPQ heatmap colors: sequential scale from actual per-sample NPQ (not spline grid)
lo_r <- min(mat_raw, na.rm = TRUE)
hi_r <- max(mat_raw, na.rm = TRUE)
if (!is.finite(lo_r) || !is.finite(hi_r)) {
  lo_r <- 0
  hi_r <- 1
}
if (hi_r <= lo_r) {
  hi_r <- lo_r + 1e-6
}
mid_r <- (lo_r + hi_r) / 2
col_fun_raw <- circlize::colorRamp2(
  c(lo_r, mid_r, hi_r),
  c("#79B6E4", "white", "#E7352E")
)
leg_at_raw <- c(lo_r, mid_r, hi_r)
leg_labs_raw <- format(signif(leg_at_raw, 4), trim = TRUE)

# Row-wise z-scores on per-sample matrix (same columns as raw heatmap); differs from mat_z when mat_plot is smoothed
mat_raw_z <- t(scale(t(mat_raw)))
mat_raw_z[is.nan(mat_raw_z)] <- 0
lo_rz <- ceiling(min(mat_raw_z, na.rm = TRUE))
hi_rz <- floor(max(mat_raw_z, na.rm = TRUE))
if (lo_rz < 0 && hi_rz > 0) {
  col_fun_raw_z <- circlize::colorRamp2(c(lo_rz, 0, hi_rz), c("#79B6E4", "white", "#E7352E"))
  leg_at_rz <- c(lo_rz, 0, hi_rz)
} else {
  lim_rz <- max(abs(c(lo_rz, hi_rz)), na.rm = TRUE)
  if (!is.finite(lim_rz) || lim_rz < 1e-15) lim_rz <- 1
  col_fun_raw_z <- circlize::colorRamp2(c(-lim_rz, 0, lim_rz), c("#79B6E4", "white", "#E7352E"))
  leg_at_rz <- c(-lim_rz, 0, lim_rz)
}
leg_labs_rz <- format(signif(leg_at_rz, 4), trim = TRUE)
leg_labs_rz[abs(leg_at_rz) < 1e-12] <- "0"

# Optional row splits (e.g. 4 blocks) — pheatmap uses cutree_rows (integer k)
nr <- nrow(mat_z)
cutree_k <- NA
if (n_row_clusters > 0L && nr >= 2L) {
  k_use <- min(as.integer(n_row_clusters), nr)
  if (k_use >= 2L) {
    cutree_k <- k_use
  }
}

# Row clustering (must match ComplexHeatmap::pheatmap: euclidean + complete)
k_assign <- if (!is.na(cutree_k) && cutree_k >= 2L) {
  min(as.integer(cutree_k), nr)
} else {
  1L
}
if (nr > 1L) {
  d_mat <- stats::dist(mat_z, method = "euclidean")
  hc_rows <- stats::hclust(d_mat, method = "complete")
  cluster_vec <- stats::cutree(hc_rows, k = k_assign)
  row_order_names <- rownames(mat_z)[hc_rows$order]
} else {
  cluster_vec <- stats::setNames(1L, rownames(mat_z)[1])
  row_order_names <- rownames(mat_z)
}

cluster_tbl <- data.frame(
  Target = names(cluster_vec),
  cluster = as.integer(cluster_vec),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(cluster, Target)

out_csv <- file.path(output_dir, paste0("NPQ_clusters_", matrix_type, "_", column_name, ".csv"))
utils::write.csv(cluster_tbl, file = out_csv, row.names = FALSE)
message("Wrote ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))

col_labels <- if (n_smooth_grid >= 2L) {
  format(round(prog_grid, digits = 2), trim = TRUE)
} else {
  format(round(unname(prog_per_sample[colnames(mat_z)]), digits = 2), trim = TRUE)
}
show_col_names <- ncol(mat_z) <= 60L

# Raw heatmap: one column per sample; labels = that sample's progression_rate (continuous, unbinned)
col_labels_raw <- format(round(unname(prog_per_sample[colnames(mat_raw)]), digits = 2), trim = TRUE)
show_col_names_raw <- ncol(mat_raw) <= 60L

# Narrow heatmap body (reference-style): pheatmap sets heatmap width = ncol * cellwidth (pt).
# 1 in = 72 pt; 1 cm ≈ 28.35 pt.
cellwidth_pt <- (heatmap_width_cm * 28.3464567) / ncol(mat_z)
cellwidth_pt <- max(cellwidth_pt, 0.35)
cellwidth_pt_raw <- (heatmap_width_cm * 28.3464567) / ncol(mat_raw)
cellwidth_pt_raw <- max(cellwidth_pt_raw, 0.35)

# Row labels = Target (matrix rownames); shrink font when many proteins
fontsize_row <- if (nr > 120L) 4 else if (nr > 70L) 5 else 6

# Row annotation: cluster strip using my_palette (cycles if k > length(my_palette))
cl_per_row <- as.integer(cluster_vec[rownames(mat_z)])
k_levels <- sort(unique(cl_per_row))
cluster_level_names <- as.character(k_levels)
cluster_colors <- stats::setNames(
  my_palette[((seq_along(k_levels) - 1L) %% length(my_palette)) + 1L],
  cluster_level_names
)
annotation_row <- data.frame(
  Cluster = factor(cl_per_row, levels = k_levels),
  row.names = rownames(mat_z),
  stringsAsFactors = FALSE
)
annotation_colors <- list(Cluster = cluster_colors)

# Keep raw heatmaps in the exact same row order as the binned heatmap clustering,
# but draw raw panels without a left dendrogram.
mat_raw_ord <- mat_raw[row_order_names, , drop = FALSE]
mat_raw_z_ord <- mat_raw_z[row_order_names, , drop = FALSE]
annotation_row_ord <- annotation_row[row_order_names, , drop = FALSE]

# Use pheatmap-style names only: cluster_cols / color / show_* / labels_*.
# Do not pass cluster_columns, col, column_labels, row_split — those duplicate
# ht_param fields (see ComplexHeatmap::pheatmap) and trigger "matched multiple times".

# (1) Raw NPQ heatmap first — per-sample NPQ; x = raw progression_rate (column labels), no grid smoothing
ht_raw <- ComplexHeatmap::pheatmap(
  mat_raw_ord,
  name = "NPQ",
  color = col_fun_raw,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = show_col_names_raw,
  labels_col = col_labels_raw,
  fontsize_col = if (show_col_names_raw) 7 else 6,
  fontsize_row = fontsize_row,
  heatmap_legend_param = list(at = leg_at_raw, labels = leg_labs_raw),
  use_raster = ncol(mat_raw_ord) > 80L || nr > 200L,
  raster_quality = 2,
  column_title = paste0(column_name, " (per sample)"),
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = grid::gpar(fontsize = 11),
  border_color = NA,
  cellwidth = cellwidth_pt_raw,
  annotation_row = annotation_row_ord,
  annotation_colors = annotation_colors
)

out_pdf_raw <- file.path(output_dir, paste0("heatmap_rawNPQ_", matrix_type, "_", column_name, ".pdf"))
grDevices::pdf(out_pdf_raw, width = pdf_width_in, height = pdf_height_in)
ComplexHeatmap::draw(ht_raw, heatmap_legend_side = "left", padding = grid::unit(c(3, 3, 3, 3), "mm"))
grDevices::dev.off()

message("Wrote ", normalizePath(out_pdf_raw, winslash = "/", mustWork = FALSE))

# (1b) Same per-sample layout as raw, row z-scores of NPQ (only needed when smooth grid differs from mat_raw)
if (n_smooth_grid >= 2L) {
  ht_raw_z <- ComplexHeatmap::pheatmap(
    mat_raw_z_ord,
    name = "Estimated protein level",
    color = col_fun_raw_z,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = show_col_names_raw,
    labels_col = col_labels_raw,
    fontsize_col = if (show_col_names_raw) 7 else 6,
    fontsize_row = fontsize_row,
    heatmap_legend_param = list(at = leg_at_rz, labels = leg_labs_rz),
    use_raster = ncol(mat_raw_z_ord) > 80L || nr > 200L,
    raster_quality = 2,
    column_title = paste0(column_name, " (per sample; row z-score)"),
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 11),
    border_color = NA,
    cellwidth = cellwidth_pt_raw,
    annotation_row = annotation_row_ord,
    annotation_colors = annotation_colors
  )
  out_pdf_raw_z <- file.path(
    output_dir,
    paste0("heatmap_rawNPQ_rowZscore_", matrix_type, "_", column_name, ".pdf")
  )
  grDevices::pdf(out_pdf_raw_z, width = pdf_width_in, height = pdf_height_in)
  ComplexHeatmap::draw(ht_raw_z, heatmap_legend_side = "left", padding = grid::unit(c(3, 3, 3, 3), "mm"))
  grDevices::dev.off()
  message("Wrote ", normalizePath(out_pdf_raw_z, winslash = "/", mustWork = FALSE))
}

# (2) Binned / scaled heatmap — row z-scores (Estimated protein level) on same progression columns
ht <- ComplexHeatmap::pheatmap(
  mat_z,
  name = "Estimated protein level",
  color = col_fun,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  cutree_rows = cutree_k,
  show_rownames = TRUE,
  show_colnames = show_col_names,
  labels_col = col_labels,
  fontsize_col = if (show_col_names) 7 else 6,
  fontsize_row = fontsize_row,
  heatmap_legend_param = list(at = leg_at, labels = leg_labs),
  use_raster = ncol(mat_z) > 80L || nr > 200L,
  raster_quality = 2,
  column_title = column_name,
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = grid::gpar(fontsize = 11),
  border_color = NA,
  cellwidth = cellwidth_pt,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)

out_pdf <- file.path(output_dir, paste0("heatmap_NPQ_", matrix_type, "_", column_name, ".pdf"))
grDevices::pdf(out_pdf, width = pdf_width_in, height = pdf_height_in)
ComplexHeatmap::draw(ht, heatmap_legend_side = "left", padding = grid::unit(c(3, 3, 3, 3), "mm"))
grDevices::dev.off()

message("Wrote ", normalizePath(out_pdf, winslash = "/", mustWork = FALSE))

# LOESS trajectory panels: raw NPQ vs progression, one facet per cluster
loess_df <- protein_data_IDs %>%
  dplyr::filter(Target %in% cluster_tbl$Target) %>%
  dplyr::mutate(
    cluster = cluster_vec[as.character(Target)],
    cluster = factor(cluster, levels = sort(unique(as.integer(cluster))))
  )

p_loess <- ggplot2::ggplot(loess_df, ggplot2::aes(x = .data[[column_name]], y = NPQ)) +
  ggplot2::geom_line(ggplot2::aes(group = Target), alpha = 0.2, linewidth = 0.3, color = "grey45") +
  ggplot2::geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.9,
    color = "#1a5276",
    fill = "#79B6E4",
    alpha = 0.25
  ) +
  ggplot2::facet_wrap(~cluster, ncol = 1L, scales = "free_y",
                      labeller = ggplot2::labeller(cluster = function(x) paste("Cluster", x))) +
  ggplot2::labs(
    title = paste0(matrix_type, ": protein NPQ trajectories by heatmap cluster"),
    x = column_name,
    y = "NPQ"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
  )

loess_pdf <- file.path(output_dir, paste0("trajectory_LOESS_", matrix_type, "_", column_name, ".pdf"))
loess_h_in <- max(5, min(22, 1.9 * length(levels(loess_df$cluster))))
ggplot2::ggsave(
  filename = loess_pdf,
  plot = p_loess,
  width = 2,
  height = loess_h_in,
  units = "in",
  limitsize = FALSE
)
message("Wrote ", normalizePath(loess_pdf, winslash = "/", mustWork = FALSE))

# Per-protein scatter plots: NPQ vs progression rate, one file per Target
plot_dir <- file.path(output_dir, "Plot")
if (!file.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (nchar(x) == 0L) x <- "unknown_target"
  x
}

target_ids <- sort(unique(as.character(loess_df$Target)))
for (tg in target_ids) {
  df_tg <- loess_df %>% dplyr::filter(Target == tg)
  if (nrow(df_tg) < 1L) next

  p_scatter <- ggplot2::ggplot(df_tg, ggplot2::aes(x = .data[[column_name]], y = NPQ)) +
    ggplot2::geom_point(color = "#2C7FB8", alpha = 0.75, size = 1.4) +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.8,
      color = "#1a5276",
      fill = "#79B6E4",
      alpha = 0.25
    ) +
    ggplot2::labs(
      title = tg,
      subtitle = paste0(matrix_type, ": NPQ vs ", column_name),
      x = column_name,
      y = "NPQ"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35")
    )

  out_scatter <- file.path(plot_dir, paste0(safe_filename(tg), ".pdf"))
  ggplot2::ggsave(
    filename = out_scatter,
    plot = p_scatter,
    width = 3.4,
    height = 2.8,
    units = "in",
    limitsize = FALSE
  )
}
message("Wrote per-protein scatter plots to ", normalizePath(plot_dir, winslash = "/", mustWork = FALSE))

# LOESS on heatmap-aligned data: same values as heatmap body (row z-scores, mat_z) vs progression axis
prog_heatmap <- if (n_smooth_grid >= 2L) {
  as.numeric(prog_grid)
} else {
  as.numeric(unname(prog_per_sample[colnames(mat_plot)]))
}
stopifnot(length(prog_heatmap) == ncol(mat_plot))
stopifnot(all(dim(mat_z) == dim(mat_plot)))
loess_df_binned <- data.frame(
  Target = rep(rownames(mat_z), times = ncol(mat_z)),
  progression = rep(prog_heatmap, each = nrow(mat_z)),
  y_heatmap = as.vector(mat_z),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    cluster = cluster_vec[as.character(Target)],
    cluster = factor(cluster, levels = sort(unique(as.integer(cluster))))
  )

y_rng <- range(mat_z, na.rm = TRUE)

p_loess_binned <- ggplot2::ggplot(loess_df_binned, ggplot2::aes(x = progression, y = y_heatmap)) +
  ggplot2::geom_line(ggplot2::aes(group = Target), alpha = 0.2, linewidth = 0.3, color = "grey45") +
  ggplot2::geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.9,
    color = "#1a5276",
    fill = "#79B6E4",
    alpha = 0.25
  ) +
  ggplot2::facet_wrap(~cluster, ncol = 1L, scales = "fixed",
                      labeller = ggplot2::labeller(cluster = function(x) paste("Cluster", x))) +
  ggplot2::scale_y_continuous(limits = y_rng, expand = ggplot2::expansion(mult = 0.02)) +
  ggplot2::labs(
    title = paste0(matrix_type, ": trajectories (heatmap columns)"),
    subtitle = if (n_smooth_grid >= 2L) {
      paste0("Smooth grid (n=", n_smooth_grid, "); y = row z-score (same matrix as heatmap)")
    } else {
      "Per-sample columns; y = row z-score (same matrix as heatmap)"
    },
    x = column_name,
    y = "Estimated protein level"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
  )

loess_binned_pdf <- file.path(output_dir, paste0("trajectory_LOESS_binned_", matrix_type, "_", column_name, ".pdf"))
ggplot2::ggsave(
  filename = loess_binned_pdf,
  plot = p_loess_binned,
  width = 2,
  height = loess_h_in,
  units = "in",
  limitsize = FALSE
)
message("Wrote ", normalizePath(loess_binned_pdf, winslash = "/", mustWork = FALSE))
