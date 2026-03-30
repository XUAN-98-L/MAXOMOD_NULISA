# Signed pvalue plot
###############################################
### Helper Functions
add_DE_flag <- function(df, lfc, padj, up, down, alpha) {
  df %>%
    mutate(
      log10_padj = -log10(.data[[padj]]),
      DE = case_when(
        .data[[padj]] < alpha & .data[[lfc]] > 0 ~ up,
        .data[[padj]] < alpha & .data[[lfc]] < 0 ~ down,
        TRUE ~ "ns"
      )
    )
}

compute_signed <- function(df, logfc_col, comp_name, up, down, alpha_threshold = 0.1, use_fdr = TRUE){
  
  # Explicitly choose padj or pvalue
  sig_col <- if(use_fdr) paste0("padj_", comp_name) else paste0("pvalue_", comp_name)
  
  #if(!sig_col %in% colnames(df)) sig_col <- paste0("pvalue_", comp_name)
  
  df %>%
    add_DE_flag(lfc = logfc_col, padj = sig_col, up = up, down = down, alpha = alpha_threshold) %>%
    mutate(signed = sign(.data[[logfc_col]]) * log10_padj) %>%
    #mutate(signed = ifelse(is.na(signed), 0, signed)) %>%
    select(Target, UniProtID, signed, DE, any_of("pvalue_anova")) %>% 
    rename(!!paste0("signed_", comp_name) := signed,
           !!paste0("DE_", comp_name) := DE)
}

# Reference: MAXOMOD_CSF/Script/12_Scatterplot_FDR.R — signed -log10(FDR) scatter with quadrants + LM
# Labels: all proteins significant in at least one comparison (|x| or |y| >= cut_off), not quantile top.
scatterplot_signed_FDR <- function(
    data,
    cut_off = -log10(0.05),
    main_title,
    max.overlaps = 10,
    lab_x,
    lab_y,
    text_x = NULL,
    text_y = NULL,
    labels_T_F = TRUE,
    annotate_YN = TRUE
) {
  if (is.null(text_x)) text_x <- paste0("significant in ", lab_x)
  if (is.null(text_y)) text_y <- paste0("significant in ", lab_y)
  data <- data %>%
    mutate(
      omic_type = case_when(
        abs(.data$y) >= cut_off & abs(.data$x) >= cut_off ~ "significant in both",
        abs(.data$y) >= cut_off ~ text_y,
        abs(.data$x) >= cut_off ~ text_x,
        TRUE ~ "ns"
      )
    )
  cols <- c("salmon", "#26b3ff", "grey", "mediumpurple1")
  names(cols) <- c(text_x, text_y, "ns", "significant in both")
  data_labeled <- dplyr::filter(data, abs(.data$x) >= cut_off | abs(.data$y) >= cut_off)
  plot <- ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(colour = omic_type), alpha = 0.5, shape = 16, size = 2) +
    geom_point(
      data = dplyr::filter(data, abs(.data$y) >= cut_off | abs(.data$x) >= cut_off),
      aes(colour = omic_type),
      alpha = 0.5,
      shape = 16,
      size = 3
    ) +
    geom_smooth(method = "lm", color = "#2C3E50", se = TRUE) +
    geom_hline(yintercept = cut_off, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -cut_off, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = cut_off, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = -cut_off, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(
      title = main_title,
      x = lab_x,
      y = lab_y,
      colour = "Differential\nExpression"
    ) +
    theme_classic() +
    theme(
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 15, hjust = 0.5),
      text = element_text(size = 14)
    ) +
    coord_cartesian(clip = "off")
  if (isTRUE(labels_T_F) && nrow(data_labeled) > 0) {
    plot <- plot +
      geom_text_repel(
        data = data_labeled,
        aes(label = name),
        force = 1,
        hjust = 1,
        max.overlaps = max.overlaps,
        segment.size = 0.2,
        min.segment.length = 0,
        size = 2
      )
  }
  if (isTRUE(annotate_YN)) {
    n_y <- sum(data$omic_type == text_y, na.rm = TRUE)
    n_x <- sum(data$omic_type == text_x, na.rm = TRUE)
    n_b <- sum(data$omic_type == "significant in both", na.rm = TRUE)
    xr <- range(data$x, na.rm = TRUE)
    yr <- range(data$y, na.rm = TRUE)
    ax <- xr[1] + 0.02 * diff(xr)
    ay <- yr[2] - 0.02 * diff(yr)
    plot <- plot +
      annotate(
        "text",
        x = ax,
        y = ay,
        hjust = 0,
        vjust = 1,
        label = paste0(n_y, " ", text_y, "\n", n_x, " ", text_x, "\n", n_b, " significant in both"),
        size = 3.2
      )
  }
  plot
}

# ## signed p-values (arrow-style panel)
# signed_plot <- function(data, xvar, yvar, title, comparison, alpha = 0.1,
#                         point_colour = "#2c3e50", label_all = FALSE){
  
#   line_pos <- -log10(alpha)
  
#   get_groups <- function(v) { strsplit(gsub("signed_", "", v), "_")[[1]] }
#   groups_x <- get_groups(xvar)
#   groups_y <- get_groups(yvar)
  
#   x_title <- paste0("signed -log10 p-value (",comparison, ", ", groups_x[1], ")")
#   y_title <- paste0("signed -log10 p-value (", comparison, ", ", groups_y[1], ")")
  
#   limit_val <- max(abs(c(data[[xvar]], data[[yvar]])), na.rm = TRUE) * 1.1
  
#   # Calibrated for a tighter "Interpretation Zone"
#   x_arrow_pos <- -limit_val * 1.40 
#   y_arrow_pos <- -limit_val * 1.42 
  
#   p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
#     geom_point(colour = point_colour, size = 2.25, alpha = 0.55)
  
#   if (isTRUE(label_all) && nrow(data) > 0) {
#     p <- p + geom_text_repel(aes(label = Target), colour = "black",
#                              max.overlaps = Inf, size = 2.8, fontface = "bold", box.padding = 0.35)
#   }
  
#   p <- p + 
#     geom_vline(xintercept = c(line_pos, -line_pos), linetype = "dashed", colour = "grey", alpha = 0.6) + 
#     geom_hline(yintercept = c(line_pos, -line_pos), linetype = "dashed", colour = "grey", alpha = 0.6) +
#     theme_classic(base_size = 12) + 
#     labs(title = title, x = x_title, y = y_title) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 10)),
#       axis.title = element_text(size = 11, face = "plain"),
#       axis.title.x = element_text(margin = margin(t = 8, b = 32)), 
#       axis.title.y = element_text(margin = margin(r = 8, l = 38)), 
#       legend.title = element_text(face = "bold", size = 10),
#       plot.margin = margin(t = 5, r = 5, b = 10, l = 10) 
#     ) + 
#     coord_cartesian(clip = "off", xlim = c(-limit_val, limit_val), ylim = c(-limit_val, limit_val)) 
  
#   p <- p +
#     annotate("segment", x = line_pos, xend = limit_val*0.8, y = x_arrow_pos, yend = x_arrow_pos, 
#              arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
#     annotate("text", x = (line_pos + limit_val*0.8)/2, y = x_arrow_pos, label = groups_x[2], 
#              size = 3.8, fontface = "bold.italic", vjust = 1.4) + 
    
#     annotate("segment", x = -line_pos, xend = -limit_val*0.8, y = x_arrow_pos, yend = x_arrow_pos, 
#              arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
#     annotate("text", x = (-line_pos - limit_val*0.8)/2, y = x_arrow_pos, label = groups_x[3], 
#              size = 3.8, fontface = "bold.italic", vjust = 1.4) + 
    
#     annotate("segment", x = y_arrow_pos, xend = y_arrow_pos, y = line_pos, yend = limit_val*0.8, 
#              arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
#     annotate("text", x = y_arrow_pos, y = (line_pos + limit_val*0.8)/2, label = groups_y[2], 
#              angle = 90, size = 3.8, fontface = "bold.italic", vjust = -1.1) + 
    
#     annotate("segment", x = y_arrow_pos, xend = y_arrow_pos, y = -line_pos, yend = -limit_val*0.8, 
#              arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
#     annotate("text", x = y_arrow_pos, y = (-line_pos - limit_val*0.8)/2, label = groups_y[3], 
#              angle = 90, size = 3.8, fontface = "bold.italic", vjust = -1.1)
  
#   return(p)
# }

###############################################
### Libraries
suppressMessages(library("optparse"))
suppressMessages(library("readxl"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))


# =======================Parse command line arguments======================
option_list = list(
  make_option(c( "--input1"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/04_Differential_Expression_FDR_005/DE_subgroups_CSF_groups_covariate_adjusted.xlsx", help="Input p-values file"),
  make_option(c( "--input2"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/04_Differential_Expression_FDR_005/DE_subgroups_PLASMA_groups_covariate_adjusted.xlsx", help="Input p-values file"),
  make_option(c("-o", "--output"), type="character", default="CNS_immune/Results/CNS_panel/Without_tears/05_Signed_p", help="Output signed p-value plot path"),
  # prefix for the input files
  make_option(c( "--prefix1"), type="character", default="CSF", help="Prefix for the first input files"),
  make_option(c( "--prefix2"), type="character", default="PLASMA", help="Prefix for the second input files"),
  make_option(c( "--comparison"), type="character", default="ALS_CTRL",
                help="Contrast name matching DE columns (e.g. ALS_CTRL for log2FC_ALS_CTR/padj_ALS_CTRL; alpha_beta for subtype)"),
  make_option(c("-a", "--alpha"), type="numeric", default=0.05,
                help="FDR / p-value threshold for DE flags and plot cut-off (default 0.05)")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# =======================options setting======================
if (is.null(opt$input1)) {
  print("NO INPUT P-VALUES FILE SUPPLIED, EXITING!")
  stop("Please provide the input p-values file path!")
} else {
  input_file1 <- opt$input1
  data1 <- read_excel(input_file1)
}

if (is.null(opt$input2)) {
  print("NO INPUT P-VALUES FILE SUPPLIED, EXITING!")
  stop("Please provide the input p-values file path!")
} else {
  input_file2 <- opt$input2
  data2 <- read_excel(input_file2)
}

if (is.null(opt$prefix1)) {
  print("NO PREFIX FOR THE FIRST INPUT FILE SUPPLIED, EXITING!")
  stop("Please provide the prefix for the first input file!")
} else {
  prefix1 <- opt$prefix1
}

if (is.null(opt$prefix2)) {
  print("NO PREFIX FOR THE SECOND INPUT FILE SUPPLIED, EXITING!")
  stop("Please provide the prefix for the second input file!")
} else {
  prefix2 <- opt$prefix2
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

comparison <- trimws(as.character(opt$comparison[1]))
if (!nzchar(comparison)) {
  stop("Please provide a non-empty --comparison (e.g. ALS_CTRL or alpha_beta).")
}

# log2FC / group labels must match 04_differential_expression / DE xlsx columns
log2fc_col <- switch(
  comparison,
  "ALS_CTRL" = "log2FC_ALS_CTR",
  "alpha_beta" = "log2FC_alpha_beta",
  paste0("log2FC_", comparison)
)
grp <- switch(
  comparison,
  "ALS_CTRL" = list(up = "ALS", down = "CTR"),
  "alpha_beta" = list(up = "alpha", down = "beta"),
  list(up = "ALS", down = "CTR")
)

alpha <- as.numeric(opt$alpha)

# =======================data processing======================
use_fdr_flag = TRUE
df1 <- compute_signed(
  data1,
  log2fc_col,
  comparison,
  up = grp$up,
  down = grp$down,
  alpha_threshold = alpha,
  use_fdr = use_fdr_flag
)
df2 <- compute_signed(
  data2,
  log2fc_col,
  comparison,
  up = grp$up,
  down = grp$down,
  alpha_threshold = alpha,
  use_fdr = use_fdr_flag
)

# for data1, add prefix1 to each column name, unless the column name is Target or UniProtID
df1 <- df1 %>%
  rename_with(~ paste0(prefix1, "_", .x), .cols = setdiff(colnames(df1), c("Target", "UniProtID")))

# for data2, expect Target and UniProtID columns, add prefix2 to each column name
df2 <- df2 %>%
  rename_with(~ paste0(prefix2, "_", .x), .cols = setdiff(colnames(df2), c("Target", "UniProtID")))

merged <- inner_join(df1, df2, by = c("Target", "UniProtID"))

# p <- signed_plot(
#     merged,
#     paste0(prefix1,"_","signed_",  comparison),
#     paste0(prefix2,"_","signed_", comparison),
#     paste0(prefix1, " vs ",prefix2," (",comparison,")"),
#     alpha = alpha,
#     comparison = comparison)

# ggsave(
#     paste0(output_dir, "/signed_", prefix1, "vs", prefix2, " (",comparison,").pdf"),
#     p,
#     width = 10,
#     height = 8
# )

# --- Reference 12_Scatterplot_FDR.R style: quadrant colours, LM line, quantile labels, counts ---
x_col <- paste0(prefix1, "_signed_", comparison)
y_col <- paste0(prefix2, "_signed_", comparison)
plot_df <- merged %>%
  transmute(
    name = Target,
    x = .data[[x_col]],
    y = .data[[y_col]]
  )

p_fdr <- scatterplot_signed_FDR(
  plot_df,
  cut_off = -log10(alpha),
  main_title = paste0(prefix1, " vs ", prefix2, " (", comparison, ")"),
  max.overlaps = Inf,
  lab_x = paste0("signed -log10(FDR) (", comparison, ", ", prefix1, ")"),
  lab_y = paste0("signed -log10(FDR) (", comparison, ", ", prefix2, ")"),
  text_x = paste0("significant in ", prefix1),
  text_y = paste0("significant in ", prefix2),
  labels_T_F = TRUE,
  annotate_YN = TRUE
)

ggsave(
  paste0(output_dir, "/signed_FDR_scatter_", prefix1, "vs", prefix2, "_", comparison, ".pdf"),
  p_fdr,
  width = 6,
  height = 4,
  dpi = 300
)

sink(paste0(output_dir, "/cor_pearson_signed_", prefix1, "vs", prefix2, "_", comparison, ".txt"))
print(stats::cor.test(plot_df$x, plot_df$y, method = "pearson"))
sink()