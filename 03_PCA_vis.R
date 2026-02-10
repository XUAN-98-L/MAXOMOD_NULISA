
### Libraries
suppressMessages(library("optparse"))
suppressMessages(library("readxl"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stats"))
suppressMessages(library("ggrepel"))
suppressMessages(library("umap"))
#=================Functions=================
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

count_subtypes_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid, !is.na(subtype)) %>%
    distinct(SampleName, subtype) %>%
    group_by(subtype) %>%
    summarise(nr_samples = n(), .groups = "drop") %>%
    mutate(biofluid = fluid)
}

# PCA per fluid labelled based on status
# value_col: "NPQ" or "NPQ_adj" (use NPQ_adj when adjust = TRUE)
run_pca <- function(df, matrix_type, value_col = "NPQ") {
  df_matrix <- df %>%
    filter(SampleMatrixType == matrix_type) %>%
    dplyr::select(SampleName, Patient, Target, all_of(value_col), type, subtype, age, sex, any_of("cohort"), any_of("PlateID"))
  
  # Reshape wide
  df_wide <- df_matrix %>%
    pivot_wider(
      names_from = Target,
      values_from = all_of(value_col),
      values_fill = 0
    )
  
  meta <- df_wide %>% dplyr::select(SampleName, Patient, type, subtype, age, sex, any_of("cohort"), any_of("PlateID"))
  X <- df_wide %>% dplyr::select(-SampleName, -Patient, -type, -subtype, -age, -sex, -any_of("cohort"), -any_of("PlateID"))
  
  X <- X %>% select(where(~ {
    s <- sd(.x, na.rm = TRUE)
    !is.na(s) && s > 0
  }))
  
  n_samp <- nrow(X)
  n_feat <- ncol(X)
  if (n_samp < 2 || n_feat < 2) {
    stop(
      "PCA failed for fluid '", matrix_type, "': need at least 2 samples (have ", n_samp,
      ") and 2 features with non-zero variance (have ", n_feat, "). ",
      "With adjust=TRUE and many covariates (e.g. age,sex,PlateID), NPQ_adj can be mostly NA or constant. ",
      "Try --covariates \"age,sex\" without PlateID, or check missing values."
    )
  }
  
  pca <- prcomp(X, scale. = TRUE)
  
  scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(meta)
  
  list(scores = scores, pca = pca)
}

# PCA on all fluids combined; color by fluid (SampleMatrixType)
# value_col: "NPQ" or "NPQ_adj". exclude_targets: optional character vector of protein names to drop (e.g. c("APOE4","CRP")).
run_pca_all <- function(df, value_col = "NPQ", exclude_targets = c("APOE4","CRP")) {
  df <- df %>%
    mutate(across(all_of(value_col), ~ suppressWarnings(as.numeric(.x))))
  df_matrix <- df %>%
    dplyr::select(SampleName, SampleMatrixType, Patient, Target, all_of(value_col), type, subtype, age, sex, any_of("cohort"), any_of("PlateID"))
  df_wide <- df_matrix %>%
    tidyr::pivot_wider(
      names_from = Target,
      values_from = all_of(value_col),
      values_fill = 0
    )
  meta <- df_wide %>%
    dplyr::select(SampleName, SampleMatrixType, Patient, type, subtype, age, sex, any_of("cohort"), any_of("PlateID"))
  X <- df_wide %>%
    dplyr::select(-SampleName, -SampleMatrixType, -Patient, -type, -subtype, -age, -sex, -any_of("cohort"), -any_of("PlateID"))
  X <- X %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  #  Drop columns with all NA or zero variance
  X <- X %>% select(where(~ {
    all_na <- all(is.na(.x))
    s <- sd(.x, na.rm = TRUE)
    !all_na && !is.na(s) && s > 0
  }))

  #  Replace remaining NA with 0 
  X[is.na(X)] <- 0
  if (!is.null(exclude_targets) && length(exclude_targets) > 0) {
    X <- X %>% dplyr::select(-any_of(exclude_targets))
  }
  n_samp <- nrow(X)
  n_feat <- ncol(X)
  if (n_samp < 2 || n_feat < 2) {
    stop(
      "PCA (all fluids) failed: need at least 2 samples (have ", n_samp,
      ") and 2 features with non-zero variance (have ", n_feat, ")."
    )
  }

  # PCA
  pca <- prcomp(X, scale. = TRUE)
  scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(meta)
  list(scores = scores, pca = pca)
}

# Fluid colors for PCA all-fluids plot (user-specified + TEARS)
fluid_palette <- c(
  "CSF"    = "#1B9E77",
  "PLASMA" = "#D95F02",
  "SERUM"  = "#7570B3",
  "TEARS"  = "#E78AC3"
)

# Plot PCA (all fluids) colored by SampleMatrixType
plot_pca_all <- function(pca_res, title = "All fluids", label = FALSE, fluid_palette = fluid_palette) {
  pca <- pca_res$pca
  scores <- pca_res$scores
  pal <- fluid_palette[names(fluid_palette) %in% unique(scores$SampleMatrixType)]
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = SampleMatrixType)) +
    geom_point(size = 4, alpha = 0.8, shape = 16) +
    scale_color_manual(values = pal, name = "Fluid", na.value = "gray70") +
    theme_minimal(base_size = 16) +
    labs(
      title = paste("PCA -", title),
      x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
    ) +
    theme(
      text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 16)
    )
  if (label) p <- p + geom_text_repel(aes(label = Patient), size = 4)
  p
}

# # UMAP per fluid (same wide matrix as PCA)
# # value_col: "NPQ" or "NPQ_adj" (use NPQ_adj when adjust = TRUE)
# run_umap <- function(df, matrix_type, value_col = "NPQ") {
#   df_matrix <- df %>%
#     filter(SampleMatrixType == matrix_type) %>%
#     dplyr::select(SampleName, Patient, Target, all_of(value_col), type, subtype, age, sex, any_of("cohort"), any_of("PlateID"))
#   df_wide <- df_matrix %>%
#     pivot_wider(names_from = Target, values_from = all_of(value_col), values_fill = 0)
#   meta <- df_wide %>% dplyr::select(SampleName, Patient, type, subtype, age, sex, any_of("cohort"), any_of("PlateID"))
#   X <- df_wide %>% dplyr::select(-SampleName, -Patient, -type, -subtype, -age, -sex, -any_of("cohort"), -any_of("PlateID"))
#   X <- X %>% select(where(~ { s <- sd(.x, na.rm = TRUE); !is.na(s) && s > 0 }))
#   umap_out <- umap::umap(as.matrix(X))
#   scores <- as_tibble(umap_out$layout) %>%
#     setNames(c("UMAP1", "UMAP2")) %>%
#     bind_cols(meta)
#   list(scores = scores, umap = umap_out)
# }

# # Plot function for UMAP (same style and legend logic as plot_pca)
# # group_by: "type", "subtype", "sex", "cohort", or "plateid" (must exist in scores)
# plot_umap <- function(umap_res, title, subtype = FALSE, label = FALSE, group_by = "type") {
#   use_subtype <- isTRUE(subtype[1])
#   scores <- umap_res$scores %>%
#     mutate(
#       subtype = if (use_subtype) factor(subtype, levels = c("CTRL", "alpha", "beta")) else factor(subtype, levels = c("CTRL", "ALS"))
#     )
#   col_var <- if (group_by == "plateid") "PlateID" else group_by
#   if (group_by %in% c("sex", "cohort", "plateid") && col_var %in% names(scores)) {
#     scores <- scores %>% filter(!is.na(.data[[col_var]]))
#     nlev <- n_distinct(scores[[col_var]])
#     p <- ggplot(scores, aes(x = UMAP1, y = UMAP2, color = .data[[col_var]])) +
#       geom_point(size = 5, alpha = 0.8, shape = 16) +
#       scale_color_brewer(palette = if (nlev <= 8) "Set2" else "Set3", name = col_var, na.value = "gray70")
#   } else if (group_by == "subtype" || (group_by == "type" && use_subtype)) {
#     p <- ggplot(scores, aes(x = UMAP1, y = UMAP2, color = subtype, shape = type)) +
#       geom_point(size = 5, alpha = 0.8) +
#       scale_shape_manual(values = c(CTRL = 16, ALS = 17), name = "type") +
#       scale_color_manual(values = c('CTRL' = 'black', 'ALS' = '#D73027', 'alpha' = '#FC8D62', 'beta' = '#8DA0CB'), name = "subtype") +
#       guides(
#         shape = guide_legend(order = 1, title = "type", override.aes = list(color = c("black", "#D73027"))),
#         color  = guide_legend(order = 2, title = "subtype")
#       )
#   } else {
#     p <- ggplot(scores, aes(x = UMAP1, y = UMAP2, color = type, shape = subtype)) +
#       geom_point(size = 5, alpha = 0.8) +
#       scale_color_manual(values = c('CTRL' = 'black', 'ALS' = '#D73027', 'alpha' = '#FC8D62', 'beta' = '#8DA0CB'))
#   }
#   p <- p +
#     theme_minimal(base_size = 16) +
#     labs(title = paste("UMAP -", title), x = "UMAP1", y = "UMAP2") +
#     theme(
#       text = element_text(size = 16),
#       axis.title = element_text(size = 18),
#       axis.text = element_text(size = 16),
#       plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
#       legend.title = element_text(size = 17),
#       legend.text = element_text(size = 16)
#     )
#   if (label) p <- p + geom_text_repel(aes(label = Patient), size = 4)
#   p
# }

# Plot function for PCA
# group_by: "type", "subtype", "sex", "cohort", or "plateid"
# plateid_palette: when group_by == "plateid", use this named vector (PlateID -> color) so same PlateID has same color across fluids
plot_pca <- function(pca_res, title, subtype = FALSE, label = FALSE, group_by = "type", plateid_palette = NULL) {
  use_subtype <- isTRUE(subtype[1])
  pca <- pca_res$pca
  scores <- pca_res$scores %>%
    mutate(
      subtype = if (use_subtype) factor(subtype, levels = c("CTRL","alpha","beta")) else factor(subtype, levels = c("CTRL","ALS"))
    )
  col_var <- if (group_by == "plateid") "PlateID" else group_by
  if (group_by %in% c("sex", "cohort", "plateid") && col_var %in% names(scores)) {
    scores <- scores %>% filter(!is.na(.data[[col_var]]))
    if (group_by == "sex") {
      p <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[col_var]])) +
        geom_point(size = 5, alpha = 0.8, shape = 16) +
        scale_color_manual(values = c(Female = "#E78AC3", Male = "#66C2A5"), name = col_var, na.value = "gray70")
    } else if (group_by == "cohort") {
      p <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[col_var]])) +
        geom_point(size = 5, alpha = 0.8, shape = 16) +
        scale_color_manual(values = c(DC = "salmon", VC = "#26b3ff"), name = col_var, na.value = "gray70")
    } else {
      # plateid: use global plateid_palette when provided (same PlateID -> same color across all fluids)
      if (!is.null(plateid_palette)) {
        pal <- plateid_palette
      } else {
        plate_levels <- sort(unique(as.character(scores[[col_var]])))
        pal <- setNames(my_palette[((seq_along(plate_levels) - 1) %% length(my_palette)) + 1], plate_levels)
      }
      p <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[col_var]])) +
        geom_point(size = 5, alpha = 0.8, shape = 16) +
        scale_color_manual(values = pal, name = col_var, na.value = "gray70")
    }
  } else if (group_by == "subtype" || (group_by == "type" && use_subtype)) {
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = subtype, shape = type)) +
      geom_point(size = 5, alpha = 0.8) +
      scale_shape_manual(values = c(CTRL = 16, ALS = 17), name = "type") +
      scale_color_manual(values = c('CTRL' = 'black', 'ALS' = '#D73027', 'alpha' = '#FC8D62', 'beta' = '#8DA0CB'), name = "subtype") +
      guides(
        shape = guide_legend(order = 1, title = "type", override.aes = list(color = c("black", "black"))),
        color  = guide_legend(order = 2, title = "subtype")
      )
  } else {
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = type, shape = subtype)) +
      geom_point(size = 5, alpha = 0.8) +
      scale_color_manual(values = c('CTRL' = 'black', 'ALS' = '#D73027', 'alpha' = '#FC8D62', 'beta' = '#8DA0CB'))
  }
  p <- p +
    theme_minimal(base_size = 16) +
    labs(
      title = paste("PCA -", title),
      x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)")
    ) +
    theme(
      text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 16)
    )
  if (label) p <- p + geom_text_repel(aes(label = Patient), size = 4)
  p
}



# Remove effect of age, sex in data (for PCA and UMAP)
remove_effect_covariates <- function(
    df,
    response = "NPQ",
    feature = "Target",
    sample = "SampleID",
    keep = "SampleMatrixType",
    remove = c("age", "sex")
) {
  
  df %>%
    group_by(.data[[feature]]) %>%
    mutate(
      NPQ_adj = local({
        
        sub_df <- pick(everything())
        vars <- if (is.null(keep)) c(response, remove) else c(response, keep, remove)
        ok <- complete.cases(sub_df[, vars])
        out <- rep(NA_real_, nrow(sub_df))
        if (sum(ok) >= length(vars) + 1) {
          rhs <- if (is.null(keep)) paste(remove, collapse = " + ") else paste(keep, "+", paste(remove, collapse = " + "))
          f <- as.formula(paste(response, "~", rhs))
          fit <- lm(f, data = sub_df[ok, , drop = FALSE])
          X <- model.matrix(fit)
          beta <- coef(fit)
          drop_cols <- colnames(X) %in%
            unlist(lapply(remove, function(v)
              grep(paste0("^", v), colnames(X), value = TRUE)))
          valid <- drop_cols & !is.na(beta[colnames(X)])
          if (any(valid))
            out[ok] <- sub_df[[response]][ok] -
              as.numeric(X[, valid, drop = FALSE] %*% beta[valid])
        }
        
        out
      })
    ) %>%
    ungroup()
}

#=================Parse command line arguments=================
option_list = list(
  make_option(c("-i", "--input"), type="character", default="Results/01_Data_Mining/protein_data_IDs.xlsx", help="Input SampleName with clinical info file"),
  make_option(c("-o", "--output"), type="character", default="Results/03_PCA_vis", help="Output PCA vis path"),
  make_option(c("-m","--metadata"), type="character", default="Results/00_Initialization/all_participants_IDs.xlsx", help="Input all participants IDs file"),
  make_option(c("-s","--subtype"), type="logical", default=FALSE, help="If TRUE, analyze subtype (alpha, beta, CTRL) instead of ALS and CTRL"),
 # make_option(c("-c","--sample_counts"), type="character", default="Results/01_Data_Mining/samples_biofluid_overview.xlsx", help="Input samples counts file"),
  # set seed for reproducibility
  make_option(c("-r","--seed"), type="integer", default=123, help="Set seed for reproducibility"),
  make_option(c("-l","--label"), type="logical", default=FALSE, help="If TRUE, add labels to the plot"),
  # fluids
  make_option(c("-f","--fluids"), type="character", default="CSF,SERUM,PLASMA,TEARS", help="Input fluids to analyze, separated by commas"),
  make_option(c("-a","--adjust"), type="logical", default=FALSE, help="If TRUE, adjust for age and sex in PCA and UMAP"),
  make_option(c("-g","--group_by"), type="character", default="type", help="Variable to display in figure: type, subtype, sex, cohort, plateid"),
  make_option(c("-n","--npq_counts"), type="character", default="Results/00_Initialization/npq_counts.xlsx", help="Input NPQ counts file"),
  make_option("--covariates", type="character", default="age,sex", help="Comma-separated covariates to remove when adjust=TRUE (e.g. age,sex or age,sex,center)")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#=================options setting=================
set.seed(opt$seed)

if (is.null(opt$input)) {
  print("NO INPUT NPQ COUNTS FILE SUPPLIED, EXITING!")
  stop("Please provide the input NPQ counts file path!")
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

if (is.null(opt$metadata)) {
  print("NO METADATA FILE PATH SUPPLIED, EXITING!")
  stop("Please provide the metadata file path!")
} else {
  metadata_file <- opt$metadata
  samples_ID_type <- read_excel(metadata_file)
}

if (is.null(opt$subtype)) {
  print("NO SUBTYPE OPTION SUPPLIED, EXITING!")
  stop("Please provide the subtype option!")
} else {
  subtype <- opt$subtype
}

# if (is.null(opt$sample_counts)) {
#   print("NO SAMPLE COUNTS FILE PATH SUPPLIED, EXITING!")
#   stop("Please provide the sample counts file path!")
# } else {
#   sample_counts <- read_excel(opt$sample_counts)
# }

if (is.null(opt$label)) {
  print("NO LABEL OPTION SUPPLIED, EXITING!")
  stop("Please provide the label option!")
} else {
  label <- opt$label
}

if (is.null(opt$fluids)) {
  print("NO FLUIDS OPTION SUPPLIED, EXITING!")
  stop("Please provide the fluids option!")
} else {
  # split by comma, trim whitespace, drop empty (e.g. trailing comma)
  fluids <- trimws(unlist(strsplit(as.character(opt$fluids[1]), ",")))
  fluids <- fluids[nzchar(fluids)]
  if (length(fluids) == 0) stop("No valid fluid names after parsing --fluids.")
}

if (is.null(opt$adjust)) {
  print("NO ADJUST OPTION SUPPLIED, EXITING!")
  stop("Please provide the adjust option!")
} else {
  adjust <- opt$adjust
}

if (is.null(opt$group_by)) {
  group_by_var <- "type"
} else {
  group_by_var <- tolower(trimws(as.character(opt$group_by[1])))
  if (!group_by_var %in% c("type", "subtype", "sex", "cohort", "plateid")) {
    stop("--group_by must be one of: type, subtype, sex, cohort, plateid")
  }
}

if (is.null(opt$npq_counts)) {
  print("NO INPUT NPQ COUNTS FILE SUPPLIED, EXITING!")
  stop("Please provide the input NPQ counts file path!")
} else {
  input_npq_file <- opt$npq_counts
  npq_counts <- read_excel(input_npq_file)
}

if (is.null(opt$covariates)) {
  covariates <- c("age", "sex")
} else {
  covariates <- trimws(unlist(strsplit(as.character(opt$covariates[1]), ",")))
  covariates <- covariates[nzchar(covariates)]
  if (length(covariates) == 0) covariates <- c("age", "sex")
}

###############################################################################
# Run pipeline
###############################################################################
#PCA by fluid and subtypes 
if (subtype) {
samples_ID_type <- samples_ID_type %>%
  mutate(type = case_when(
    type == "ALS" & k2 == "alpha" ~ "alpha",
    type == "ALS" & k2 == "beta"  ~ "beta",
    TRUE ~ type
  ))
}

# filter npq_counts based one SampleName & only keep one row for each SampleName
npq_counts <- npq_counts %>%
  filter(SampleName %in% protein_data_IDs$SampleName) %>%
  group_by(SampleName) %>%
  slice_head(n = 1) %>%
  ungroup()

protein_data_IDs <- protein_data_IDs %>%
  left_join(npq_counts %>% select(SampleName, PlateID), by = "SampleName")

protein_data_PCA <- protein_data_IDs %>%
  left_join(samples_ID_type %>% rename(subtype = type,
                                       SampleName = `Tube_ID`))


subtypes_counts <- bind_rows(
  lapply(fluids, function(f) count_subtypes_per_fluid(protein_data_PCA, f))
)

# sample_counts = rbind(sample_counts,subtypes_counts %>% rename(type = subtype)) %>%
#   arrange(biofluid) 

protein_data_clean <- protein_data_PCA %>% filter(!is.na(type))

if (adjust) {
  if ("PlateID" %in% covariates && "PlateID" %in% names(protein_data_clean)) {
    fluids_plate <- protein_data_clean %>%
      group_by(SampleMatrixType) %>%
      summarise(n_plate = n_distinct(PlateID, na.rm = TRUE), .groups = "drop")
    protein_data_clean <- bind_rows(lapply(fluids, function(f) {
      sub <- protein_data_clean %>% filter(SampleMatrixType == f)
      remove_f <- if (fluids_plate$n_plate[fluids_plate$SampleMatrixType == f] > 1) covariates else setdiff(covariates, "PlateID")
      remove_effect_covariates(sub, remove = remove_f, keep = NULL)
    }))
  } else {
    protein_data_clean <- remove_effect_covariates(protein_data_clean, remove = covariates)
  }
}

value_col <- if (adjust) "NPQ_adj" else "NPQ"
matrices <- fluids
pca_results_subtype <- lapply(matrices, function(m) run_pca(protein_data_clean, m, value_col = value_col))
names(pca_results_subtype) <- matrices

# When group_by is plateid: assign one color per unique PlateID across ALL fluids (so same PlateID = same color in every figure)
if (group_by_var == "plateid" && "PlateID" %in% names(protein_data_clean)) {
  all_plateids <- sort(unique(as.character(na.omit(protein_data_clean$PlateID))))
  plateid_palette <- setNames(my_palette[((seq_along(all_plateids) - 1) %% length(my_palette)) + 1], all_plateids)
} else {
  plateid_palette <- NULL
}

plots_subtype <- lapply(names(pca_results_subtype), 
                        function(m) plot_pca(pca_results_subtype[[m]], m, subtype = subtype, label = label, group_by = group_by_var, plateid_palette = plateid_palette))


# save PCA figures
for (f in fluids) {
  idx <- which(names(pca_results_subtype) == f)
  if (length(idx) == 0) { warning("No PCA result for fluid: ", f); next }
  pdf(paste0(output_dir, "/PCA_", f, if (group_by_var != "type") paste0("_", group_by_var) else "", if(isTRUE(subtype[1])) "_subtype" else "", if(isTRUE(label[1])) "_label" else "", ".pdf"), width = 8, height = 6.5)
  print(plots_subtype[[idx[1]]])
  dev.off()
}

# Run PCA on all fluids together (one combined plot, colored by fluid)
# Uses value_col: NPQ_adj when adjust=TRUE, NPQ otherwise
pca_all <- tryCatch(
  run_pca_all(protein_data_clean, value_col = value_col),
  error = function(e) { warning("PCA all fluids skipped: ", conditionMessage(e)); NULL }
)
if (!is.null(pca_all)) {
  all_title <- if (adjust) "All fluids (adjusted)" else "All fluids"
  all_suffix <- if (adjust) "_adjusted" else ""
  all_fname <- paste0(output_dir, "/PCA_all_fluids", all_suffix, if(isTRUE(label[1])) "_label" else "", ".pdf")
  pdf(all_fname, width = 8, height = 6.5)
  print(plot_pca_all(pca_all, title = all_title, label = label, fluid_palette = fluid_palette))
  dev.off()
  message("Saved PCA all fluids: ", all_fname)
}

# # UMAP per fluid (same data as PCA)
# umap_results <- lapply(matrices, function(m) run_umap(protein_data_clean, m, value_col = value_col))
# names(umap_results) <- matrices
# plots_umap <- lapply(names(umap_results), function(m) plot_umap(umap_results[[m]], m, subtype = subtype, label = label, group_by = group_by_var))
# for (f in fluids) {
#   idx <- which(names(umap_results) == f)
#   if (length(idx) == 0) { warning("No UMAP result for fluid: ", f); next }
#   pdf(paste0(output_dir, "/UMAP_", f, if (group_by_var != "type") paste0("_", group_by_var) else "", if(isTRUE(subtype[1])) "_subtype" else "", if(isTRUE(label[1])) "_label" else "", ".pdf"), width = 8, height = 6.5)
#   print(plots_umap[[idx[1]]])
#   dev.off()
# }




