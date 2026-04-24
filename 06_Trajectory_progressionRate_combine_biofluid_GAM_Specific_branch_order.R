# Trajectory NPQ across progression — combine CSF / PLASMA / SERUM on one progression grid,
# per-protein feature = [z-scored smooth CSF | z-scored smooth PLASMA | z-scored smooth SERUM],
# then cluster rows (proteins) on the concatenated matrix.
# Heatmap rows (K>=2): ordered cluster 1..K top to bottom; hclust leaf order within cluster; row dendrogram omitted.
# Between-cluster spacing: gaps_row (pheatmap sets row_gap to unit(4, "bigpts") internally; do not pass row_gap again).
# NPQ ~ progression: Gamma GAM linear vs s(progression); lmrob linear. shape_by_gamma_AIC: nonlinear if smooth p<0.05,
# Edf>2, AIC_smooth<AIC_linear; else linear if lmrob progression p<0.05; else neither. Covariates age_z + sex_f per fluid.
# Run log: trajectory_curve_source_*.csv + table() counts — mgcv vs fallback (see curve_method column).

######### Functions ##########################
my_palette <- c(
  "#66C2A5", "#E78AC3", "#ffc773",
  "#1B9E77", "#7570B3", "#FEE090",
  "#A6D854", "#8DA0CB", "#FC8D62",
  "#8FBC94", "#b0a4e3", "#ffa631",
  "#0aa344", "#e4c6d0", "#ffa400",
  "#519a73", "#4b5cc4", "#eedeb0",
  "#549688", "#ffb3a7", "#b35c44",
  "#7fecad", "#a1afc9", "#a78e44",
  "#519a73", "#2e4e7e", "#955539"
)

# Fallback: unadjusted smoothing spline (same as non-GAM script).
gen_smooth_curves_npq_fixed_grid <- function(df, prog_col, value_col, id_col, grid_x) {
  if (length(grid_x) < 2L) {
    stop("grid_x must have length >= 2.")
  }
  targets <- unique(as.character(df[[id_col]]))
  n_grid <- length(grid_x)
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

# GAM: NPQ ~ s(progression) + age_z + sex, Gamma(log); predict response on grid_x at reference covariates.
# Attributes: curve_method (per Target), gam_gamma_chosen, n_clean_rows — see trajectory_curve_source_*.csv.
gen_gam_curves_npq_fixed_grid <- function(
  df,
  prog_col,
  value_col,
  id_col,
  grid_x,
  min_n_gam = 15L,
  gamma_seq = seq(0.5, 5, by = 1L),
  k_max = 10L
) {
  if (length(grid_x) < 2L) {
    stop("grid_x must have length >= 2.")
  }
  if (!all(c("age_z", "sex_f") %in% names(df))) {
    stop("df must contain columns age_z and sex_f for GAM.")
  }
  targets <- unique(as.character(df[[id_col]]))
  n_grid <- length(grid_x)
  out <- matrix(NA_real_, nrow = length(targets), ncol = n_grid,
                dimnames = list(targets, NULL))
  curve_method <- rep(NA_character_, length(targets))
  names(curve_method) <- targets
  gam_gamma_chosen <- rep(NA_real_, length(targets))
  names(gam_gamma_chosen) <- targets
  n_clean_rows <- rep(NA_integer_, length(targets))
  names(n_clean_rows) <- targets

  for (i in seq_along(targets)) {
    tg <- targets[i]
    sub <- df[as.character(df[[id_col]]) == tg, , drop = FALSE]
    ok <- is.finite(sub[[value_col]]) & is.finite(sub[[prog_col]]) &
      is.finite(sub$age_z) & !is.na(sub$sex_f)
    sub <- sub[ok, , drop = FALSE]
    n_clean_rows[tg] <- as.integer(nrow(sub))
    if (nrow(sub) < 4L) {
      curve_method[tg] <- "skip_clean_rows_lt_4"
      next
    }
    sub$NPQ_pos <- pmax(sub[[value_col]], 1e-9)
    sex_f <- factor(sub$sex_f)
    k_prog <- max(4L, min(as.integer(k_max), as.integer(nrow(sub) - 3L)))
    use_gam <- nrow(sub) >= min_n_gam
    fml <- if (length(levels(sex_f)) >= 2L) {
      stats::as.formula(paste0(
        "NPQ_pos ~ s(", prog_col, ", k=", k_prog, ") + age_z + sex_f"
      ))
    } else {
      stats::as.formula(paste0(
        "NPQ_pos ~ s(", prog_col, ", k=", k_prog, ") + age_z"
      ))
    }
    sub$sex_f <- sex_f
    best_fit <- NULL
    best_aic <- Inf
    best_gamma <- NA_real_
    if (use_gam) {
      for (ga in gamma_seq) {
        fit <- tryCatch(
          mgcv::gam(
            fml,
            family = Gamma(link = "log"),
            data = sub,
            method = "REML",
            select = TRUE,
            gamma = ga
          ),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          aic <- stats::AIC(fit)
          if (is.finite(aic) && aic < best_aic) {
            best_aic <- aic
            best_fit <- fit
            best_gamma <- ga
          }
        }
      }
    }
    nd <- data.frame(grid_x)
    names(nd)[1] <- prog_col
    nd$age_z <- mean(sub$age_z, na.rm = TRUE)
    if (length(levels(sex_f)) >= 2L) {
      nd$sex_f <- factor(levels(sex_f)[1], levels = levels(sex_f))
    }
    if (use_gam && !is.null(best_fit)) {
      pr <- tryCatch(
        as.numeric(stats::predict(best_fit, newdata = nd, type = "response")),
        error = function(e) NULL
      )
      if (!is.null(pr) && length(pr) == n_grid && all(is.finite(pr))) {
        out[i, ] <- pr
        curve_method[tg] <- "mgcv_gam_Gamma_log"
        gam_gamma_chosen[tg] <- best_gamma
        next
      }
    }

    fb_prefix <- if (!use_gam) {
      paste0("fallback_n_lt_", min_n_gam, "_skipped_gam")
    } else if (is.null(best_fit)) {
      "fallback_gam_no_valid_fit"
    } else {
      "fallback_gam_predict_failed"
    }

    sub2 <- sub[, c(prog_col, value_col), drop = FALSE]
    names(sub2) <- c("px", "py")
    sub2 <- sub2[order(sub2$px), , drop = FALSE]
    sub2 <- stats::aggregate(py ~ px, sub2, mean)
    px <- sub2$px
    py <- sub2$py
    if (length(px) < 2L) {
      curve_method[tg] <- paste0(fb_prefix, "__no_curve_lt2_unique_prog")
      next
    }
    if (length(px) == 2L) {
      out[i, ] <- stats::approx(px, py, xout = grid_x, rule = 2)$y
      curve_method[tg] <- paste0(fb_prefix, "__unadjusted_linear_interp_2pts")
    } else {
      sp <- tryCatch(
        stats::smooth.spline(px, py, cv = FALSE),
        error = function(e) NULL
      )
      if (is.null(sp)) {
        out[i, ] <- stats::approx(px, py, xout = grid_x, rule = 2)$y
        curve_method[tg] <- paste0(fb_prefix, "__unadjusted_approx")
      } else {
        out[i, ] <- as.numeric(stats::predict(sp, x = grid_x)$y)
        curve_method[tg] <- paste0(fb_prefix, "__unadjusted_smooth_spline")
      }
    }
  }
  attr(out, "progression_grid") <- grid_x
  attr(out, "curve_method") <- curve_method
  attr(out, "gam_gamma_chosen") <- gam_gamma_chosen
  attr(out, "n_clean_rows") <- n_clean_rows
  attr(out, "min_n_gam") <- min_n_gam
  out
}

# s() row of summary(fit)$s.table for prog_col smooth.
extract_gam_smooth_term_edf_p <- function(fit, prog_col) {
  edf <- NA_real_
  pv <- NA_real_
  if (is.null(fit)) {
    return(list(edf = edf, p = pv))
  }
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm) || is.null(sm$s.table) || nrow(sm$s.table) == 0L) {
    return(list(edf = edf, p = pv))
  }
  st <- sm$s.table
  rn <- rownames(st)
  idx <- grep(paste0("^s\\(", prog_col), rn)[1L]
  if (is.na(idx)) {
    idx <- 1L
  }
  cn <- colnames(st)
  edf_col <- grep("^edf$", cn, ignore.case = TRUE)[1L]
  pcol <- grep("p[-. ]?value|Pr\\(", cn, ignore.case = TRUE)[1L]
  if (!is.na(edf_col)) {
    edf <- suppressWarnings(as.numeric(st[idx, edf_col]))
  }
  if (!is.na(pcol)) {
    pv <- suppressWarnings(as.numeric(st[idx, pcol]))
  }
  list(edf = edf, p = pv)
}

# Per protein×fluid: classify nonlinear vs linear vs neither (see file header); lmrob + Gamma GAM terms; delta AIC kept for reference.
build_trajectory_model_comparison <- function(
  df,
  prog_col,
  value_col,
  id_col,
  min_n_gam = 15L,
  gamma_seq = seq(0.5, 5, by = 1L),
  delta_aic_cut = 2,
  k_max = 10L
) {
  if (!all(c("age_z", "sex_f") %in% names(df))) {
    stop("df must contain age_z and sex_f.")
  }
  targets <- unique(as.character(df[[id_col]]))
  rows <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    tg <- targets[i]
    sub <- df[as.character(df[[id_col]]) == tg, , drop = FALSE]
    ok <- is.finite(sub[[value_col]]) & is.finite(sub[[prog_col]]) &
      is.finite(sub$age_z) & !is.na(sub$sex_f)
    sub <- sub[ok, , drop = FALSE]
    n <- nrow(sub)
    if (n < 4L) {
      rows[[i]] <- data.frame(
        Target = tg,
        n = n,
        AIC_gam_linear = NA_real_,
        AIC_gam_smooth = NA_real_,
        delta_AIC_smooth_minus_linear = NA_real_,
        gam_smooth_term_edf = NA_real_,
        gam_smooth_term_pvalue = NA_real_,
        gam_AIC_prefers_smooth = NA,
        shape_by_gamma_AIC = "insufficient_n",
        lmrob_NPQ_progression_coef = NA_real_,
        lmrob_progression_SE = NA_real_,
        lmrob_progression_pvalue = NA_real_,
        lmrob_r_squared = NA_real_,
        lmrob_adj_r_squared = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    sub$NPQ_pos <- pmax(sub[[value_col]], 1e-9)
    sex_f <- factor(sub$sex_f)
    sub$sex_f <- sex_f
    k_prog <- max(4L, min(as.integer(k_max), as.integer(nrow(sub) - 3L)))

    fml_lin <- if (length(levels(sex_f)) >= 2L) {
      stats::as.formula(paste0("NPQ_pos ~ ", prog_col, " + age_z + sex_f"))
    } else {
      stats::as.formula(paste0("NPQ_pos ~ ", prog_col, " + age_z"))
    }
    fit_lin <- tryCatch(
      mgcv::gam(fml_lin, family = Gamma(link = "log"), data = sub, method = "REML"),
      error = function(e) NULL
    )
    aic_lin <- if (!is.null(fit_lin)) stats::AIC(fit_lin) else NA_real_

    fml_smooth <- if (length(levels(sex_f)) >= 2L) {
      stats::as.formula(paste0(
        "NPQ_pos ~ s(", prog_col, ", k=", k_prog, ") + age_z + sex_f"
      ))
    } else {
      stats::as.formula(paste0(
        "NPQ_pos ~ s(", prog_col, ", k=", k_prog, ") + age_z"
      ))
    }
    best_fit <- NULL
    best_aic <- Inf
    smooth_attempted <- n >= min_n_gam
    if (smooth_attempted) {
      for (ga in gamma_seq) {
        fit <- tryCatch(
          mgcv::gam(
            fml_smooth,
            family = Gamma(link = "log"),
            data = sub,
            method = "REML",
            select = TRUE,
            gamma = ga
          ),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          aic <- stats::AIC(fit)
          if (is.finite(aic) && aic < best_aic) {
            best_aic <- aic
            best_fit <- fit
          }
        }
      }
    }
    aic_smooth <- if (!is.null(best_fit) && is.finite(best_aic)) best_aic else NA_real_

    delta <- if (is.finite(aic_lin) && is.finite(aic_smooth)) {
      aic_smooth - aic_lin
    } else {
      NA_real_
    }

    sp_smooth <- extract_gam_smooth_term_edf_p(best_fit, prog_col)
    gam_smooth_edf <- sp_smooth$edf
    gam_smooth_p <- sp_smooth$p
    aic_prefers_smooth <- if (is.finite(aic_lin) && is.finite(aic_smooth)) {
      aic_smooth < aic_lin
    } else {
      NA
    }
    nonlinear_ok <- isTRUE(smooth_attempted) &&
      !is.null(best_fit) &&
      is.finite(gam_smooth_p) && gam_smooth_p < 0.05 &&
      is.finite(gam_smooth_edf) && gam_smooth_edf > 2 &&
      isTRUE(aic_prefers_smooth)

    f_lm <- if (length(levels(sex_f)) >= 2L) {
      stats::as.formula(paste0(value_col, " ~ ", prog_col, " + age_z + sex_f"))
    } else {
      stats::as.formula(paste0(value_col, " ~ ", prog_col, " + age_z"))
    }
    fit_lmrob <- tryCatch(
      robustbase::lmrob(f_lm, data = sub, setting = "KS2014"),
      error = function(e) NULL
    )
    coef_prog <- NA_real_
    se_prog <- NA_real_
    p_prog <- NA_real_
    r2_lmrob <- NA_real_
    adjr2_lmrob <- NA_real_
    if (!is.null(fit_lmrob)) {
      sm <- tryCatch(summary(fit_lmrob), error = function(e) NULL)
      if (!is.null(sm) && is.matrix(sm$coefficients)) {
        ct <- sm$coefficients
        idx <- match(prog_col, rownames(ct))
        if (!is.na(idx)) {
          coef_prog <- as.numeric(ct[idx, 1L])
          se_prog <- as.numeric(ct[idx, 2L])
          pcol <- grep("^Pr\\(", colnames(ct), ignore.case = TRUE)
          if (length(pcol) >= 1L) {
            p_prog <- as.numeric(ct[idx, pcol[1L]])
          }
        }
      }
      if (!is.null(sm)) {
        r2_lmrob <- as.numeric(sm$r.squared)[1]
        adjr2_lmrob <- as.numeric(sm$adj.r.squared)[1]
      }
    }

    lmrob_linear_ok <- is.finite(p_prog) && p_prog < 0.05
    shape <- if (isTRUE(nonlinear_ok)) {
      "nonlinear"
    } else if (lmrob_linear_ok) {
      "linear"
    } else {
      "neither"
    }

    rows[[i]] <- data.frame(
      Target = tg,
      n = n,
      AIC_gam_linear = aic_lin,
      AIC_gam_smooth = aic_smooth,
      delta_AIC_smooth_minus_linear = delta,
      gam_smooth_term_edf = gam_smooth_edf,
      gam_smooth_term_pvalue = gam_smooth_p,
      gam_AIC_prefers_smooth = aic_prefers_smooth,
      shape_by_gamma_AIC = shape,
      lmrob_NPQ_progression_coef = coef_prog,
      lmrob_progression_SE = se_prog,
      lmrob_progression_pvalue = p_prog,
      lmrob_r_squared = r2_lmrob,
      lmrob_adj_r_squared = adjr2_lmrob,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

###############################################
suppressMessages(library(optparse))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(grid))
suppressMessages(library(ggplot2))
suppressMessages(library(mgcv))
suppressMessages(library(robustbase))

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx",
    help = "Input protein_data_IDs.xlsx"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "CNS_immune/Results/CNS_panel/Without_tears/06_Trajectory_progressionRate_combine_biofluid_GAM",
    help = "Output directory"
  ),
  make_option(
    c("-c", "--column_name"),
    type = "character",
    default = "progression_rate",
    help = "Numeric progression column"
  ),
  make_option(
    c("--fluids"),
    type = "character",
    default = "CSF,PLASMA,SERUM",
    help = "Comma-separated SampleMatrixType values (e.g. CSF,PLASMA,SERUM)"
  ),
  make_option(
    c("--n_row_clusters"),
    type = "integer",
    default = 4L,
    help = "Row cluster blocks (0 = no cutree split)"
  ),
  make_option(
    c("--n_smooth_grid"),
    type = "integer",
    default = 20L,
    help = "Grid points on unified progression axis (must be >= 2)"
  ),
  make_option(
    c("--heatmap_width_cm"),
    type = "double",
    default = 14,
    help = "Heatmap body width in cm (wider default: 3 fluid blocks)"
  ),
  make_option(c("--pdf_width"), type = "double", default = 10),
  make_option(c("--pdf_height"), type = "double", default = 8),
  make_option(
    c("--progression_sd_k"),
    type = "double",
    default = 0,
    help = "Exclude samples with mean progression outside cohort mean ± k*SD (pooled across fluids)"
  ),
  make_option(
    c("--min_n_gam"),
    type = "integer",
    default = 15L,
    help = "Minimum samples per protein×fluid to attempt mgcv::gam; else smooth.spline fallback"
  ),
  make_option(
    c("--gam_gamma_step"),
    type = "double",
    default = 1,
    help = "Step for gamma grid (0.5–5) when selecting GAM by AIC; smaller = more models, slower"
  ),
  make_option(
    c("--delta_aic"),
    type = "double",
    default = 2,
    help = "Legacy: passed to model comparison (shape uses GAM smooth p/Edf/AIC + lmrob p; see CSV columns)"
  ),
  make_option(
    c("--k-max"),
    type = "integer",
    default = 10L,
    dest = "k_max",
    help = "Cap for s(progression) basis size k: k_prog = max(4, min(k_max, n_clean-3)); larger allows wigglier smooths"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

input_file <- opt$input
output_dir <- opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

column_name <- opt$column_name
fluids <- toupper(trimws(strsplit(opt$fluids, ",", fixed = TRUE)[[1]]))
fluids <- fluids[nzchar(fluids)]
allowed <- c("CSF", "PLASMA", "SERUM", "TEARS")
if (length(fluids) < 2L) {
  stop("--fluids must list at least two matrix types.")
}
if (!all(fluids %in% allowed)) {
  stop("--fluids must be subset of: ", paste(allowed, collapse = ", "))
}

n_row_clusters <- as.integer(opt$n_row_clusters)
n_smooth_grid <- as.integer(opt$n_smooth_grid)
if (is.na(n_smooth_grid) || n_smooth_grid < 2L) {
  stop("--n_smooth_grid must be >= 2 (unified grid + smoothing).")
}
heatmap_width_cm <- as.numeric(opt$heatmap_width_cm)
pdf_width_in <- as.numeric(opt$pdf_width)
pdf_height_in <- as.numeric(opt$pdf_height)
progression_sd_k <- as.numeric(opt$progression_sd_k)
min_n_gam <- as.integer(opt$min_n_gam)
if (is.na(min_n_gam) || min_n_gam < 4L) {
  stop("--min_n_gam must be an integer >= 4.")
}
gam_gamma_step <- as.numeric(opt$gam_gamma_step)
if (is.na(gam_gamma_step) || gam_gamma_step <= 0) {
  stop("--gam_gamma_step must be positive.")
}
gamma_seq <- seq(0.5, 5, by = gam_gamma_step)
delta_aic_cut <- as.numeric(opt$delta_aic)
if (is.na(delta_aic_cut) || delta_aic_cut < 0) {
  stop("--delta_aic must be non-negative.")
}
k_max <- as.integer(opt$k_max)
if (is.na(k_max) || k_max < 4L) {
  stop("--k-max must be an integer >= 4.")
}

if (is.na(progression_sd_k) || progression_sd_k < 0) {
  stop("--progression_sd_k must be non-negative (0 = no exclusion).")
}

fluid_tag <- paste(fluids, collapse = "_")

protein_data_IDs <- read_excel(input_file)
protein_data_IDs[[column_name]] <- as.numeric(protein_data_IDs[[column_name]])

protein_data_IDs <- protein_data_IDs %>%
  filter(SampleMatrixType %in% fluids)

if (!column_name %in% names(protein_data_IDs)) {
  stop("Column `", column_name, "` not found.")
}
protein_data_IDs <- protein_data_IDs %>% filter(!is.na(!!sym(column_name)))
if (!"NPQ" %in% names(protein_data_IDs)) {
  stop("Input must contain NPQ.")
}
protein_data_IDs <- protein_data_IDs %>% filter(!is.na(NPQ))

if (!all(c("age", "sex") %in% names(protein_data_IDs))) {
  stop("GAM mode requires `age` and `sex` columns in the input (merge clinical metadata into protein_data_IDs).")
}
protein_data_IDs$age <- suppressWarnings(as.numeric(protein_data_IDs$age))
protein_data_IDs$sex <- trimws(as.character(protein_data_IDs$sex))
protein_data_IDs <- protein_data_IDs %>%
  dplyr::filter(is.finite(age), !is.na(sex), nzchar(sex))

if (nrow(protein_data_IDs) == 0L) {
  stop("No rows after filtering fluids and NA checks.")
}

# Pooled per-sample progression (all fluids) for ordering and SD filter
sample_order_df <- protein_data_IDs %>%
  group_by(SampleName) %>%
  summarise(!!sym(column_name) := mean(!!sym(column_name), na.rm = TRUE), .groups = "drop") %>%
  arrange(!!sym(column_name))

if (progression_sd_k > 0) {
  pv <- sample_order_df[[column_name]]
  mu_pr <- mean(pv, na.rm = TRUE)
  sigma_pr <- stats::sd(pv, na.rm = TRUE)
  if (!is.finite(mu_pr) || !is.finite(sigma_pr) || sigma_pr < 1e-12) {
    stop("--progression_sd_k: need finite mean and non-zero SD of per-sample ", column_name, ".")
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
      "Excluded ", nrow(excluded_df), " sample(s). Wrote ",
      normalizePath(excl_path, winslash = "/", mustWork = FALSE)
    )
    sample_order_df <- sample_order_df[keep_sd, , drop = FALSE]
    protein_data_IDs <- protein_data_IDs %>%
      dplyr::filter(SampleName %in% sample_order_df$SampleName)
  }
  if (nrow(sample_order_df) < 2L) {
    stop("After SD filter need at least 2 samples.")
  }
}

# Unified progression range across all retained rows (all fluids)
prog_min <- min(protein_data_IDs[[column_name]], na.rm = TRUE)
prog_max <- max(protein_data_IDs[[column_name]], na.rm = TRUE)
if (!is.finite(prog_min) || !is.finite(prog_max) || prog_min >= prog_max) {
  stop("Progression must span at least two finite values across combined fluids.")
}
grid_x <- seq(prog_min, prog_max, length.out = n_smooth_grid)

# Per fluid: complete-case matrix then smooth on grid_x; collect targets present
mats_smooth <- vector("list", length(fluids))
names(mats_smooth) <- fluids
targets_per_fluid <- vector("list", length(fluids))
model_cmp_rows <- vector("list", length(fluids))
curve_trace_rows <- vector("list", length(fluids))

for (fi in seq_along(fluids)) {
  fl <- fluids[fi]
  pd <- protein_data_IDs %>% filter(SampleMatrixType == fl)
  if (nrow(pd) == 0L) {
    stop("No rows for fluid ", fl, ".")
  }
  pd$age_z <- as.numeric(scale(pd$age))
  sex_levels <- sort(unique(as.character(pd$sex)))
  pd$sex_f <- factor(pd$sex, levels = sex_levels)

  mat_wide <- pd %>%
    dplyr::select(Target, SampleName, NPQ) %>%
    group_by(Target, SampleName) %>%
    summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = SampleName, values_from = NPQ)
  mat <- as.matrix(mat_wide[, -1, drop = FALSE])
  rownames(mat) <- mat_wide$Target
  if (any(is.na(mat))) {
    mat <- mat[rowSums(is.na(mat)) == 0L, , drop = FALSE]
  }
  if (nrow(mat) < 1L) {
    stop("No proteins with complete NPQ across all samples for fluid ", fl, ".")
  }
  df_s <- pd %>%
    dplyr::filter(Target %in% rownames(mat)) %>%
    group_by(Target, SampleName) %>%
    dplyr::summarise(
      NPQ = mean(NPQ, na.rm = TRUE),
      !!sym(column_name) := mean(!!sym(column_name), na.rm = TRUE),
      age_z = dplyr::first(age_z),
      sex_f = dplyr::first(sex_f),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(!!sym(column_name)), is.finite(NPQ), is.finite(age_z))

  model_cmp_rows[[fi]] <- build_trajectory_model_comparison(
    df_s,
    prog_col = column_name,
    value_col = "NPQ",
    id_col = "Target",
    min_n_gam = min_n_gam,
    gamma_seq = gamma_seq,
    delta_aic_cut = delta_aic_cut,
    k_max = k_max
  )
  model_cmp_rows[[fi]]$SampleMatrixType <- fl

  ms_full <- gen_gam_curves_npq_fixed_grid(
    df_s,
    prog_col = column_name,
    value_col = "NPQ",
    id_col = "Target",
    grid_x = grid_x,
    min_n_gam = min_n_gam,
    gamma_seq = gamma_seq,
    k_max = k_max
  )
  rn_mat <- rownames(mat)
  ms <- ms_full[rn_mat, , drop = FALSE]
  cm_full <- attr(ms_full, "curve_method")
  ggm_full <- attr(ms_full, "gam_gamma_chosen")
  nc_full <- attr(ms_full, "n_clean_rows")
  mn_gam_attr <- attr(ms_full, "min_n_gam")
  if (!is.null(cm_full)) {
    attr(ms, "curve_method") <- cm_full[rn_mat]
    attr(ms, "gam_gamma_chosen") <- ggm_full[rn_mat]
    attr(ms, "n_clean_rows") <- nc_full[rn_mat]
    attr(ms, "min_n_gam") <- mn_gam_attr
  }
  curve_trace_rows[[fi]] <- data.frame(
    Target = rn_mat,
    SampleMatrixType = fl,
    n_clean_rows = as.integer(nc_full[rn_mat]),
    min_n_gam = as.integer(mn_gam_attr),
    curve_method = as.character(cm_full[rn_mat]),
    gam_gamma_chosen = as.numeric(ggm_full[rn_mat]),
    stringsAsFactors = FALSE
  )
  mats_smooth[[fi]] <- ms
  targets_per_fluid[[fi]] <- rownames(ms)
}

model_cmp_all <- dplyr::bind_rows(model_cmp_rows)
path_model_cmp <- file.path(
  output_dir,
  paste0("trajectory_linear_vs_nonlinear_lmrob_", fluid_tag, "_", column_name, ".csv")
)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
utils::write.csv(model_cmp_all, path_model_cmp, row.names = FALSE)
message(
  "Linear vs nonlinear (Gamma GAM AIC) + lmrob(NPQ, no extra log): ",
  normalizePath(path_model_cmp, winslash = "/", mustWork = FALSE)
)

curve_trace_all <- dplyr::bind_rows(curve_trace_rows)
path_curve_trace <- file.path(
  output_dir,
  paste0("trajectory_curve_source_", fluid_tag, "_", column_name, ".csv")
)
utils::write.csv(curve_trace_all, path_curve_trace, row.names = FALSE)
message(
  "Per-target heatmap curve source (mgcv::gam vs unadjusted fallback): ",
  normalizePath(path_curve_trace, winslash = "/", mustWork = FALSE)
)
message("\n=== Counts: curve_method (each protein × fluid row) ===")
print(table(curve_trace_all$curve_method, useNA = "ifany"))
message(
  "Tip: rows with curve_method == \"mgcv_gam_Gamma_log\" used mgcv; ",
  "others used fallback (see prefix fallback_n_lt_* vs fallback_gam_*)."
)

common_targets <- Reduce(intersect, targets_per_fluid)
if (length(common_targets) < 2L) {
  stop(
    "Need at least 2 proteins present with complete data in every fluid. ",
    "Per-fluid counts: ",
    paste(names(mats_smooth), lengths(targets_per_fluid), sep = "=", collapse = "; ")
  )
}

# Per-fluid row z-score along grid, then column-bind
blocks_z <- lapply(seq_along(fluids), function(fi) {
  m <- mats_smooth[[fi]][common_targets, , drop = FALSE]
  z <- t(scale(t(m)))
  z[is.nan(z)] <- 0
  z
})
mat_concat <- do.call(cbind, blocks_z)
rownames(mat_concat) <- common_targets
colnames(mat_concat) <- paste(
  rep(fluids, each = n_smooth_grid),
  rep(seq_len(n_smooth_grid), times = length(fluids)),
  sep = "_"
)

# Clustering on concatenated z-scores (same metric as visualized)
nr <- nrow(mat_concat)
mat_z <- mat_concat

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

cutree_k <- NA
if (n_row_clusters > 0L && nr >= 2L) {
  k_use <- min(as.integer(n_row_clusters), nr)
  if (k_use >= 2L) cutree_k <- k_use
}
k_assign <- if (!is.na(cutree_k) && cutree_k >= 2L) {
  min(as.integer(cutree_k), nr)
} else {
  1L
}

if (nr > 1L) {
  d_mat <- stats::dist(mat_z, method = "euclidean")
  hc_rows <- stats::hclust(d_mat, method = "complete")
  cluster_vec <- stats::cutree(hc_rows, k = k_assign)
} else {
  hc_rows <- NULL
  cluster_vec <- stats::setNames(1L, rownames(mat_z)[1])
}

# Heatmap row order: cluster 1, 2, …, K from top to bottom (within each cluster, keep hclust leaf order).
heatmap_rows_by_cluster_id <- (nr > 1L && k_assign >= 2L)
if (heatmap_rows_by_cluster_id) {
  leaf_order <- rownames(mat_z)[hc_rows$order]
  cl_leaf <- as.integer(cluster_vec[leaf_order])
  row_order_by_cluster <- unlist(
    lapply(seq_len(k_assign), function(cid) leaf_order[cl_leaf == cid]),
    use.names = FALSE
  )
  mat_z <- mat_z[row_order_by_cluster, , drop = FALSE]
  message(
    "Heatmap rows ordered top-to-bottom: cluster 1 .. ", k_assign,
    " (within-cluster order = hclust leaf order)."
  )
}

cluster_tbl <- data.frame(
  Target = names(cluster_vec),
  cluster = as.integer(cluster_vec),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(cluster, Target)

out_csv <- file.path(
  output_dir,
  paste0("NPQ_clusters_", fluid_tag, "_", column_name, ".csv")
)
utils::write.csv(cluster_tbl, out_csv, row.names = FALSE)
message("Wrote ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))

# --- Tables for manuscript / replotting (same numbers as heatmaps & trajectory plots) ---
meta_grid <- data.frame(
  grid_index = seq_len(n_smooth_grid),
  progression = grid_x,
  column_name = column_name,
  prog_min = prog_min,
  prog_max = prog_max,
  n_smooth_grid = n_smooth_grid,
  fluids = paste(fluids, collapse = ","),
  stringsAsFactors = FALSE
)
path_meta <- file.path(
  output_dir,
  paste0("unified_progression_grid_", fluid_tag, "_", column_name, ".csv")
)
utils::write.csv(meta_grid, path_meta, row.names = FALSE)
message("Wrote ", normalizePath(path_meta, winslash = "/", mustWork = FALSE))

mat_z_df <- data.frame(Target = rownames(mat_z), as.data.frame(mat_z, check.names = FALSE), stringsAsFactors = FALSE)
path_wide_z <- file.path(
  output_dir,
  paste0("heatmap_matrix_z_wide_", fluid_tag, "_", column_name, ".csv")
)
utils::write.csv(mat_z_df, path_wide_z, row.names = FALSE)
message("Wrote ", normalizePath(path_wide_z, winslash = "/", mustWork = FALSE))

tr_long_list <- vector("list", length(common_targets) * length(fluids) * n_smooth_grid)
idx_l <- 0L
for (tg in common_targets) {
  for (fi in seq_along(fluids)) {
    fl <- fluids[fi]
    for (j in seq_len(n_smooth_grid)) {
      col_idx <- (fi - 1L) * n_smooth_grid + j
      idx_l <- idx_l + 1L
      tr_long_list[[idx_l]] <- data.frame(
        Target = tg,
        fluid = fl,
        grid_index = j,
        progression = grid_x[j],
        NPQ_smooth = as.numeric(mats_smooth[[fi]][tg, j]),
        z_heatmap = as.numeric(mat_z[tg, col_idx]),
        cluster = as.integer(unname(cluster_vec[tg])),
        stringsAsFactors = FALSE
      )
    }
  }
}
tr_long_df <- dplyr::bind_rows(tr_long_list)
path_long <- file.path(
  output_dir,
  paste0("trajectory_long_smooth_and_z_", fluid_tag, "_", column_name, ".csv")
)
utils::write.csv(tr_long_df, path_long, row.names = FALSE)
message("Wrote ", normalizePath(path_long, winslash = "/", mustWork = FALSE))

loess_export <- protein_data_IDs %>%
  dplyr::filter(Target %in% common_targets, SampleMatrixType %in% fluids) %>%
  dplyr::mutate(cluster = as.integer(cluster_vec[as.character(Target)])) %>%
  dplyr::mutate(SampleMatrixType = factor(SampleMatrixType, levels = fluids)) %>%
  dplyr::arrange(SampleMatrixType, .data[[column_name]], Target, SampleName)
path_loess <- file.path(
  output_dir,
  paste0("sample_level_NPQ_for_LOESS_", fluid_tag, "_", column_name, ".csv")
)
utils::write.csv(loess_export, path_loess, row.names = FALSE)
message("Wrote ", normalizePath(path_loess, winslash = "/", mustWork = FALSE))

# Column split for heatmap: one block per fluid
col_split <- rep(factor(fluids, levels = fluids), each = n_smooth_grid)

col_labels <- format(round(rep(grid_x, times = length(fluids)), digits = 2), trim = TRUE)
show_col_names <- ncol(mat_z) <= 90L

cellwidth_pt <- (heatmap_width_cm * 28.3464567) / ncol(mat_z)
cellwidth_pt <- max(cellwidth_pt, 0.25)

fontsize_row <- if (nr > 120L) 4 else if (nr > 70L) 5 else 6

cl_per_row <- as.integer(cluster_vec[rownames(mat_z)])
k_levels <- sort(unique(cl_per_row))
cluster_colors <- stats::setNames(
  my_palette[((seq_along(k_levels) - 1L) %% length(my_palette)) + 1L],
  as.character(k_levels)
)
annotation_row <- data.frame(
  Cluster = factor(cl_per_row, levels = k_levels),
  row.names = rownames(mat_z),
  stringsAsFactors = FALSE
)
annotation_colors <- list(Cluster = cluster_colors)

# Visual gaps between cluster blocks (only when rows are ordered 1..K; gaps_row = last row index per block except last)
cluster_row_gaps <- if (heatmap_rows_by_cluster_id && k_assign >= 2L) {
  sz <- vapply(seq_len(k_assign), function(cid) {
    sum(as.integer(cluster_vec[rownames(mat_z)]) == cid)
  }, integer(1))
  cumsum(sz)[-k_assign]
} else {
  NULL
}

# Row dendrogram: if rows were reordered by cluster id, row tree no longer matches matrix — omit row tree.
ht <- ComplexHeatmap::pheatmap(
  mat_z,
  name = "Row z (per fluid)\nconcat",
  color = col_fun,
  cluster_cols = FALSE,
  cluster_rows = if (nr > 1L && !heatmap_rows_by_cluster_id) hc_rows else FALSE,
  cutree_rows = if (heatmap_rows_by_cluster_id) NA else cutree_k,
  gaps_row = cluster_row_gaps,
  column_split = col_split,
  column_gap = grid::unit(2, "mm"),
  show_rownames = TRUE,
  show_colnames = show_col_names,
  labels_col = col_labels,
  fontsize_col = if (show_col_names) 5 else 4,
  fontsize_row = fontsize_row,
  heatmap_legend_param = list(at = leg_at, labels = leg_labs),
  use_raster = ncol(mat_z) > 120L || nr > 200L,
  raster_quality = 2,
  column_title = fluids,
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = grid::gpar(fontsize = 11),
  border_color = NA,
  cellwidth = cellwidth_pt,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)

out_pdf <- file.path(
  output_dir,
  paste0("heatmap_NPQ_combine_", fluid_tag, "_", column_name, ".pdf")
)
grDevices::pdf(out_pdf, width = pdf_width_in, height = pdf_height_in)
ComplexHeatmap::draw(ht, heatmap_legend_side = "left", padding = grid::unit(c(3, 3, 3, 3), "mm"))
grDevices::dev.off()
message("Wrote ", normalizePath(out_pdf, winslash = "/", mustWork = FALSE))

# LOESS: NPQ vs progression, facet cluster, colour = fluid
loess_df <- protein_data_IDs %>%
  dplyr::filter(
    Target %in% cluster_tbl$Target,
    SampleMatrixType %in% fluids
  ) %>%
  dplyr::mutate(
    cluster = cluster_vec[as.character(Target)],
    cluster = factor(cluster, levels = sort(unique(as.integer(cluster))))
  )

n_cls <- length(levels(loess_df$cluster))
p_loess <- ggplot2::ggplot(loess_df, ggplot2::aes(x = .data[[column_name]], y = NPQ)) +
  ggplot2::geom_line(
    ggplot2::aes(group = interaction(Target, SampleMatrixType)),
    alpha = 0.12,
    linewidth = 0.25,
    color = "grey45"
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(color = SampleMatrixType, fill = SampleMatrixType),
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.15
  ) +
  ggplot2::facet_wrap(~cluster, ncol = 1L, scales = "free_y",
                      labeller = ggplot2::labeller(cluster = function(x) paste("Cluster", x))) +
  ggplot2::labs(
    title = paste0(fluid_tag, ": NPQ vs ", column_name, " by cluster (colours = fluid)"),
    x = column_name,
    y = "NPQ",
    color = "Fluid",
    fill = "Fluid"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
    legend.position = "bottom"
  )

loess_pdf <- file.path(
  output_dir,
  paste0("trajectory_LOESS_combine_", fluid_tag, "_", column_name, ".pdf")
)
loess_h_in <- max(12, min(24, 2 * n_cls))
ggplot2::ggsave(
  filename = loess_pdf,
  plot = p_loess,
  width = 6,
  height = loess_h_in,
  units = "in",
  limitsize = FALSE
)
message("Wrote ", normalizePath(loess_pdf, winslash = "/", mustWork = FALSE))

# Per-protein plots: NPQ vs progression, one LOESS curve per fluid, saved under Plot/
plot_dir <- file.path(output_dir, "Plot")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (nchar(x) == 0L) {
    x <- "unknown_target"
  }
  x
}

plot_df <- protein_data_IDs %>%
  dplyr::filter(
    Target %in% common_targets,
    SampleMatrixType %in% fluids
  ) %>%
  dplyr::mutate(
    SampleMatrixType = factor(SampleMatrixType, levels = fluids),
    cluster = as.integer(cluster_vec[as.character(Target)])
  )

target_plot <- sort(unique(as.character(plot_df$Target)))
for (tg in target_plot) {
  df_tg <- plot_df %>% dplyr::filter(Target == tg)
  if (nrow(df_tg) < 2L) {
    next
  }
  cl <- unique(df_tg$cluster)
  cl_lab <- if (length(cl) == 1L) paste0("cluster ", cl[1L]) else paste0("clusters ", paste(sort(cl), collapse = ","))

  p_tg <- ggplot2::ggplot(df_tg, ggplot2::aes(x = .data[[column_name]], y = NPQ)) +
    ggplot2::geom_point(
      ggplot2::aes(color = SampleMatrixType),
      alpha = 0.78,
      size = 1.5
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(color = SampleMatrixType, fill = SampleMatrixType),
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.85,
      alpha = 0.18
    ) +
    ggplot2::labs(
      title = tg,
      subtitle = paste0(fluid_tag, ": NPQ vs ", column_name, " | ", cl_lab, " (one LOESS per fluid)"),
      x = column_name,
      y = "NPQ",
      color = "Fluid",
      fill = "Fluid"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey35"),
      legend.position = "bottom"
    )

  out_tg <- file.path(plot_dir, paste0(safe_filename(tg), ".pdf"))
  ggplot2::ggsave(
    filename = out_tg,
    plot = p_tg,
    width = 4.5,
    height = 3.6,
    units = "in",
    limitsize = FALSE
  )
}
message("Wrote per-protein LOESS (one curve per fluid) to ", normalizePath(plot_dir, winslash = "/", mustWork = FALSE))

# Per-protein trajectories on the **same numeric scale as the heatmap** (mat_z: row z-score
# within each fluid block on the unified progression grid; not raw NPQ).
plot_z_dir <- file.path(output_dir, "Plot_heatmap_scale")
if (!dir.exists(plot_z_dir)) {
  dir.create(plot_z_dir, recursive = TRUE)
}

for (tg in rownames(mat_z)) {
  df_z <- dplyr::bind_rows(lapply(seq_along(fluids), function(fi) {
    idx <- (fi - 1L) * n_smooth_grid + seq_len(n_smooth_grid)
    data.frame(
      progression = grid_x,
      z_heatmap = as.numeric(mat_z[tg, idx, drop = TRUE]),
      fluid = fluids[fi],
      stringsAsFactors = FALSE
    )
  }))
  df_z$fluid <- factor(df_z$fluid, levels = fluids)
  cl <- unname(cluster_vec[as.character(tg)])
  zr <- range(df_z$z_heatmap, na.rm = TRUE)
  pad <- 0.06 * max(1e-6, diff(zr))
  if (!is.finite(pad)) {
    pad <- 0.1
  }

  p_z <- ggplot2::ggplot(df_z, ggplot2::aes(x = progression, y = z_heatmap, color = fluid, group = fluid)) +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, color = "grey55") +
    ggplot2::coord_cartesian(ylim = c(zr[1] - pad, zr[2] + pad)) +
    ggplot2::labs(
      title = tg,
      subtitle = paste0(
        fluid_tag, ": same values as heatmap | cluster ", cl,
        " (row z within each fluid; unified grid)"
      ),
      x = column_name,
      y = "Estimated protein level (heatmap scale)",
      color = "Fluid"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey35"),
      legend.position = "bottom"
    )

  out_z <- file.path(plot_z_dir, paste0(safe_filename(tg), ".pdf"))
  ggplot2::ggsave(
    filename = out_z,
    plot = p_z,
    width = 4.5,
    height = 3.6,
    units = "in",
    limitsize = FALSE
  )
}
message(
  "Wrote per-protein trajectories (heatmap z-scale) to ",
  normalizePath(plot_z_dir, winslash = "/", mustWork = FALSE)
)

# Binned LOESS: same values as heatmap body (row z-scores within each fluid block on unified grid)
prog_col_rep <- rep(grid_x, times = length(fluids))
fluid_col_rep <- rep(fluids, each = n_smooth_grid)
stopifnot(length(prog_col_rep) == ncol(mat_z), length(fluid_col_rep) == ncol(mat_z))

loess_df_binned <- data.frame(
  Target = rep(rownames(mat_z), times = ncol(mat_z)),
  progression = rep(prog_col_rep, each = nrow(mat_z)),
  fluid = rep(fluid_col_rep, each = nrow(mat_z)),
  y_heatmap = as.vector(mat_z),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    cluster = cluster_vec[as.character(Target)],
    cluster = factor(cluster, levels = sort(unique(as.integer(cluster_vec))))
  )

y_rng_b <- range(mat_z, na.rm = TRUE)

p_loess_binned <- ggplot2::ggplot(loess_df_binned, ggplot2::aes(x = progression, y = y_heatmap)) +
  ggplot2::geom_line(
    ggplot2::aes(
      group = interaction(Target, fluid),
      color = fluid
    ),
    alpha = 0.12,
    linewidth = 0.22
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(color = fluid, fill = fluid),
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.75,
    alpha = 0.12
  ) +
  ggplot2::facet_wrap(~cluster, ncol = 1L, scales = "fixed",
                      labeller = ggplot2::labeller(cluster = function(x) paste("Cluster", x))) +
  ggplot2::scale_y_continuous(limits = y_rng_b, expand = ggplot2::expansion(mult = 0.02)) +
  ggplot2::labs(
    title = paste0(fluid_tag, ": trajectories (heatmap matrix)"),
    subtitle = paste0(
      "Unified grid n=", n_smooth_grid, "; y = row z-score per fluid (same as heatmap); one LOESS per fluid"
    ),
    x = column_name,
    y = "Estimated protein level",
    color = "Fluid",
    fill = "Fluid"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
    legend.position = "bottom"
  )

loess_binned_pdf <- file.path(
  output_dir,
  paste0("trajectory_LOESS_binned_combine_", fluid_tag, "_", column_name, ".pdf")
)
ggplot2::ggsave(
  filename = loess_binned_pdf,
  plot = p_loess_binned,
  width = 6,
  height = loess_h_in,
  units = "in",
  limitsize = FALSE
)
message("Wrote ", normalizePath(loess_binned_pdf, winslash = "/", mustWork = FALSE))

message(
  "Unified grid [", signif(prog_min, 4), ", ", signif(prog_max, 4), "], n=",
  n_smooth_grid, "; proteins (intersection) = ", length(common_targets)
)
