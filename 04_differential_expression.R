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
add_DE_flag <- function(df, lfc, padj, up, down, alpha) {
  df %>%
    mutate(
      log10_padj = -log10(.data[[padj]]),
      group = case_when(
        .data[[padj]] < alpha & .data[[lfc]] > 0 ~ up,
        .data[[padj]] < alpha & .data[[lfc]] < 0 ~ down,
        TRUE ~ "ns"
      )
    )
}

# 4. Volcano plot
volcano_plot <- function(
    df, log2fc, title,
    fdr = 0.05,
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
  ggplot(df, aes_string(log2fc, "log10_padj")) +
    geom_point(aes(color = group, size = log10_padj)) +
    geom_text_repel(
      data = to_label,
      aes(label = Target),
      max.overlaps = 40,
      size = 5,
      min.segment.length = 0
    ) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    scale_color_manual(values = colors) +
    annotate(
      "text",
      x = add_x,
      y = -log10(fdr) + add_y,
      label = paste0("FDR = ", fdr * 100, "%"),
      color = "darkgrey",
      size = 5
    ) +
    theme_minimal() +
    ggtitle(title) +
    scale_size_continuous(range = c(2, 7)) +
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    labs(size = expression("-log"[10]*"(adjusted p-value)")) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 15),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15))
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
if (subtype) {
  DE_subgroups = DE_subgroups %>%
    select(Target,UniProtID,log2FC_alpha_beta,pvalue_alpha_beta,padj_alpha_beta) %>%
    add_DE_flag(lfc = "log2FC_alpha_beta",padj = "padj_alpha_beta",up = "alpha",down = "beta",alpha = cutoff)

  DE_subgroups_adjusted = DE_subgroups_adjusted %>%
    select(Target,UniProtID,log2FC_alpha_beta,pvalue_alpha_beta,padj_alpha_beta) %>%
    add_DE_flag(lfc = "log2FC_alpha_beta",padj = "padj_alpha_beta",up = "alpha",down = "beta",alpha = cutoff)
} else {
  DE_subgroups = DE_subgroups %>%
    dplyr::select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTRL,padj_ALS_CTRL) %>%
    add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTRL",up = "ALS",down = "CTR",alpha = cutoff)

  DE_subgroups_adjusted = DE_subgroups_adjusted %>%
    dplyr::select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTRL,padj_ALS_CTRL) %>%
    add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTRL",up = "ALS",down = "CTR",alpha = cutoff)
}

# set colors for the volcano plot
colors = case_when(
  subtype ~ c('alpha' = "#FC8D62", 'beta' = "#8DA0CB", 'ns' = 'lightgrey'),
  TRUE ~ c('ALS' = "#D73027", 'CTR' = "black", 'ns' = 'lightgrey')
)

volcano = volcano_plot(DE_subgroups, paste0("log2FC_", ifelse(subtype, "alpha_beta", "ALS_CTR")),
                                      title = paste0(ifelse(subtype, "alpha vs beta", "ALS vs CTR"), " in ", matrix_type),
                                      colors = colors,
                                      label_targets = label_genes,
                                      add_x = 2, add_y = -0.6, fdr = cutoff)

ggsave(paste0(output_dir, "/volcano_plot_", matrix_type, "_", ifelse(subtype, "subtypes", "groups"), ".pdf"), volcano, width = 10, height = 7)

volcano_adjusted = volcano_plot(DE_subgroups_adjusted, paste0("log2FC_", ifelse(subtype, "alpha_beta", "ALS_CTR")),
                                      title = paste0(ifelse(subtype, "alpha vs beta", "ALS vs CTR"), " in ", matrix_type, " (covariate adjusted)"),
                                      colors = colors,
                                      label_targets = label_genes,
                                      add_x = 2, add_y = -0.6, fdr = cutoff)

ggsave(paste0(output_dir, "/volcano_plot_", matrix_type, "_", ifelse(subtype, "subtypes", "groups"), "_covariate_adjusted.pdf"), volcano_adjusted, width = 10, height = 7)