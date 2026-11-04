### Libraries
suppressMessages(library("dplyr"))
suppressMessages(library("readxl"))
suppressMessages(library("optparse"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
#=================Functions=================
## 1. Mean NPQ
summarise_mean_NPQ <- function(data, matrix_type,adjusted = FALSE) {
  if(adjusted){
    data %>%
      filter(SampleMatrixType == matrix_type, !is.na(NPQ)) %>%
      group_by(type, Target, UniProtID) %>%
      summarise(mean_NPQ = mean(NPQ_adj), .groups = "drop") %>%
      pivot_wider(
        names_from = type,
        values_from = mean_NPQ,
        names_prefix = "mean_NPQ_"
      )
  } else{
  data %>%
    filter(SampleMatrixType == matrix_type, !is.na(NPQ)) %>%
    group_by(type, Target, UniProtID) %>%
    summarise(mean_NPQ = mean(NPQ), .groups = "drop") %>%
    pivot_wider(
      names_from = type,
      values_from = mean_NPQ,
      names_prefix = "mean_NPQ_"
    )
    }
}


## 2. add log2FC 
add_log2fc <- function(df, subtype = FALSE) {
  if (subtype) {
    df %>%
      mutate(
        log2FC_alpha_beta = mean_NPQ_alpha - mean_NPQ_beta,
        log2FC_alpha_CTRL = mean_NPQ_alpha - mean_NPQ_CTRL,
        log2FC_beta_CTRL = mean_NPQ_beta - mean_NPQ_CTRL
      )
  } else {
    df %>%
      mutate(
        log2FC_ALS_CTR = mean_NPQ_ALS - mean_NPQ_CTRL,
      )
}
}

# 3. add DEx info in data
add_DE_flag <- function(df, lfc, p_col, up, down, alpha) {
  df %>%
    mutate(
      log10_p = -log10(.data[[p_col]]),
      group = case_when(
        .data[[p_col]] < alpha & .data[[lfc]] > 0 ~ up,
        .data[[p_col]] < alpha & .data[[lfc]] < 0 ~ down,
        TRUE ~ "ns"
      )
    )
}

# 4. Volcano plot
volcano_plot <- function(
    df, log2fc, title,
    p_cutoff = 0.05,
    p_label = "p-value",
    colors,
    label_targets = NULL,
    add_x = 0,
    add_y = 0
) {
  to_label <- if (is.null(label_targets) || length(label_targets) == 0) {
    df %>% filter(group != "ns")
  } else {
    df %>% filter(group != "ns" | Target %in% label_targets)
  }
  y_lab <- bquote(-log[10] * "(" * .(p_label) * ")")
  ggplot(df, aes_string(log2fc, "log10_p")) +
    geom_point(aes(color = group, size = log10_p)) +
    geom_text_repel(
      data = to_label,
      aes(label = Target),
      max.overlaps = 40,
      size = 5,
      min.segment.length = 0
    ) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    scale_color_manual(values = colors) +
    annotate(
      "text",
      x = add_x,
      y = -log10(p_cutoff) + add_y,
      label = paste0(p_label, " = ", p_cutoff),
      color = "darkgrey",
      size = 5
    ) +
    theme_minimal() +
    ggtitle(title) +
    scale_size_continuous(range = c(2, 7)) +
    ylab(y_lab) +
    xlab(expression("log"[2]*"(fold-change)")) +
    labs(size = y_lab) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 15),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15))
}

prepare_volcano_df <- function(df, subtype, lfc_col, pvalue_col, padj_col, up, down, cutoff, use_fdr = FALSE) {
  p_col <- if (use_fdr) padj_col else pvalue_col
  if (subtype) {
    df %>%
      select(Target, UniProtID, all_of(c(lfc_col, pvalue_col, padj_col))) %>%
      add_DE_flag(lfc = lfc_col, p_col = p_col, up = up, down = down, alpha = cutoff)
  } else {
    df %>%
      dplyr::select(Target, UniProtID, all_of(c(lfc_col, pvalue_col, padj_col))) %>%
      add_DE_flag(lfc = lfc_col, p_col = p_col, up = up, down = down, alpha = cutoff)
  }
}

save_volcano <- function(df, output_dir, matrix_type, subtype, suffix, title_suffix, colors, label_genes, cutoff, use_fdr) {
  lfc_col <- paste0("log2FC_", ifelse(subtype, "alpha_beta", "ALS_CTR"))
  pvalue_col <- ifelse(subtype, "pvalue_alpha_beta", "pvalue_ALS_CTRL")
  padj_col <- ifelse(subtype, "padj_alpha_beta", "padj_ALS_CTRL")
  up <- ifelse(subtype, "alpha", "ALS")
  down <- ifelse(subtype, "beta", "CTR")
  p_label <- if (use_fdr) "FDR-adjusted p-value" else "raw p-value"

  plot_df <- prepare_volcano_df(df, subtype, lfc_col, pvalue_col, padj_col, up, down, cutoff, use_fdr)
  comparison <- ifelse(subtype, "alpha vs beta", "ALS vs CTR")
  p <- volcano_plot(
    plot_df,
    lfc_col,
    title = paste0(comparison, " in ", matrix_type, title_suffix, " (", p_label, ")"),
    p_cutoff = cutoff,
    p_label = p_label,
    colors = colors,
    label_targets = label_genes,
    add_x = 2,
    add_y = -0.6
  )
  out_file <- paste0(
    output_dir, "/volcano_plot_", matrix_type, "_",
    ifelse(subtype, "subtypes", "groups"), suffix,
    ifelse(use_fdr, "_FDR", "_raw_p"), ".pdf"
  )
  ggsave(out_file, p, width = 10, height = 7)
  invisible(out_file)
}

#=================Parse command line arguments=================
option_list = list(
  make_option(c("-i", "--input"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx", help="Input protein_data_IDs.xlsx file"),
  make_option(c("-o", "--output"), type="character", default="Results/04_Differential_Expression", help="Output differential expression path"),
  make_option(c("-m","--matrix_type"), type="character", default="CSF", help="Choose from: CSF, PLASMA, SERUM, TEARS"),
  make_option(c("-s","--subtype"), type="logical", default=FALSE, help="If TRUE, analyze subtype (alpha, beta, CTRL) instead of ALS and CTRL"),
  make_option(c("-p","--pvalues"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/02_Analysis_Groups/results_ALL.rds", help="Input pvalues file"),
  make_option(c("-c","--cutoff"), type="numeric", default=0.05, help="p-value cutoff for DEx (default is 0.05)"),
  make_option(c("-g","--label_genes"), type="character", default=NULL, help="Custom genes to label on volcano plot: comma-separated list or path to .txt/.csv file (one gene per line or column GeneID). If not set, only significant genes are labeled.")
  # detectability_summary
  #make_option(c("-d","--detectability_summary"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/detectability_summary.xlsx", help="Input detectability summary file")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#=================options setting=================
if (is.null(opt$input)) {
  print("NO INPUT PROTEIN DATA IDS FILE SUPPLIED, EXITING!")
  stop("Please provide the input protein_data_IDs.xlsx file path!")
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

if (is.null(opt$subtype)) {
  print("NO SUBTYPE OPTION SUPPLIED, EXITING!")
  stop("Please provide the subtype option!")
} else {
  subtype <- opt$subtype
}

if (is.null(opt$pvalues)) {
  print("NO PVALUES FILE SUPPLIED, EXITING!")
  stop("Please provide the pvalues file path!")
} else {
  pvalues_file <- opt$pvalues
  # if subtype is TRUE,need to check if "subtypes" is in the pvalues file path, if not, stop and print error message
  if (subtype && !grepl("subtypes", pvalues_file)) {
    print("SUBTYPES PVALUES FILE NOT FOUND, EXITING!")
    stop("Please provide the subtype pvalues file path!")
  } 
  results_ALL <- readRDS(pvalues_file)
}

if (is.null(opt$cutoff)) {
  print("NO CUTOFF OPTION SUPPLIED, EXITING!")
  stop("Please provide the cutoff value!")
} else {
  cutoff <- opt$cutoff
}

# Optional custom gene list for volcano plot labels (if provided, these genes are always annotated)
label_genes <- NULL
if (!is.null(opt$label_genes) && nzchar(trimws(opt$label_genes))) {
  s <- trimws(opt$label_genes)
  if (grepl("\\.(txt|csv)$", s) && file.exists(s)) {
    d <- tryCatch(
      read.table(s, header = TRUE, sep = "\t", comment.char = "", quote = "", stringsAsFactors = FALSE),
      error = function(e) read.csv(s, header = TRUE, stringsAsFactors = FALSE)
    )
    label_genes <- if ("GeneID" %in% names(d)) d$GeneID else d[[1L]]
  } else {
    label_genes <- trimws(strsplit(s, ",")[[1]])
  }
  label_genes <- label_genes[nzchar(label_genes)]
  if (length(label_genes) == 0) label_genes <- NULL
}
# if (is.null(opt$detectability_summary)) {
#   print("NO DETECTABILITY SUMMARY FILE SUPPLIED, EXITING!")
#   stop("Please provide the detectability summary file path!")
# } else {
#   detectability_summary_file <- opt$detectability_summary
#   detectability_summary <- read_excel(detectability_summary_file)
# }
#############Main Script#####################
# If subtype = TRUE, set type to alpha or beta from k2; otherwise keep as is
if (subtype) {
protein_data_IDs <- protein_data_IDs %>%
  mutate(type = case_when(
    type == "ALS" & k2 == "alpha" ~ "alpha",
    type == "ALS" & k2 == "beta"  ~ "beta",
    TRUE ~ type
  ))
}

# Without adjust
DE_subgroups = summarise_mean_NPQ(protein_data_IDs,matrix_type = matrix_type) 

DE_subgroups = DE_subgroups %>%
  add_log2fc(subtype = subtype) %>% 
  left_join(results_ALL[[matrix_type]]$pvals, by = c("Target" = "Target")) %>%
  mutate(Fluid = rep(matrix_type))

writexl::write_xlsx(DE_subgroups,paste0(output_dir, "/DE_subgroups_", matrix_type, "_", ifelse(subtype, "subtypes", "groups"), ".xlsx"))

# With sex and age adjusted
DE_subgroups_adjusted = summarise_mean_NPQ(protein_data_IDs,matrix_type = matrix_type )

DE_subgroups_adjusted = DE_subgroups_adjusted %>%
  add_log2fc(subtype = subtype) %>% 
  left_join(results_ALL[[matrix_type]]$pvals_adjusted, by = c("Target" = "Target")) %>%
  mutate(Fluid = rep(matrix_type))

writexl::write_xlsx(DE_subgroups_adjusted,paste0(output_dir, "/DE_subgroups_", matrix_type, "_", ifelse(subtype, "subtypes", "groups"), "_covariate_adjusted.xlsx"))

####### VOLCANO PLOTS
colors <- if (subtype) {
  c(alpha = "#FC8D62", beta = "#8DA0CB", ns = "lightgrey")
} else {
  c(ALS = "#D73027", CTR = "black", ns = "lightgrey")
}

# Unadjusted: raw p-value and FDR-adjusted p-value
save_volcano(DE_subgroups, output_dir, matrix_type, subtype, "", "", colors, label_genes, cutoff, use_fdr = FALSE)
save_volcano(DE_subgroups, output_dir, matrix_type, subtype, "", "", colors, label_genes, cutoff, use_fdr = TRUE)

# Covariate adjusted (age, sex): raw p-value and FDR-adjusted p-value
save_volcano(
  DE_subgroups_adjusted, output_dir, matrix_type, subtype,
  "_covariate_adjusted", " (covariate adjusted)",
  colors, label_genes, cutoff, use_fdr = FALSE
)
save_volcano(
  DE_subgroups_adjusted, output_dir, matrix_type, subtype,
  "_covariate_adjusted", " (covariate adjusted)",
  colors, label_genes, cutoff, use_fdr = TRUE
)