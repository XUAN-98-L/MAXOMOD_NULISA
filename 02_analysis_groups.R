### Libraries
suppressMessages(library("optparse"))
suppressMessages(library("readxl"))
suppressMessages(library("dplyr"))
suppressMessages(library("writexl"))
suppressMessages(library("rstatix"))
suppressMessages(library("tidyr"))
suppressMessages(library("writexl"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
#=================Functions=================
# Full pipeline for one dataset (all samples)
run_full_pipeline <- function(protein_data, sample_map, td, prefix = "ALL",
                              covariates = c("age","sex"),
                              high_detectability = FALSE,
                              output_dir = "Results/02_Analysis_Groups",
                              subtype = FALSE) {
  
  message("Preparing dataset...")
  df <- prepare_dataset(protein_data, sample_map)
  
  fluids <- c("CSF", "SERUM", "PLASMA","TEARS")
  results <- list()
  
  for (fluid in fluids) {
    message("Running stats for ", fluid, " ...")
    
    stats <- run_stats_for_fluid(df, fluid,value_col = "NPQ", adjust = FALSE)
    pvals <- extract_pvalues(stats)
    
    # Save p-values
    if(high_detectability) {
      write_xlsx(pvals,
                 paste0(output_dir, "/", prefix, "_pvalues_", fluid, "_high_detectability.xlsx"))
    } else{
      write_xlsx(pvals,
                 paste0(output_dir, "/", prefix, "_pvalues_", fluid, ".xlsx"))
    }
  
    # Adjusted results
    df_fluid <- df %>% filter(SampleMatrixType == fluid)
    stats_adj <- run_stats_for_fluid(df_fluid,fluid,"NPQ",adjust = TRUE, covariates = covariates)
    pvals_adj <- extract_pvalues(stats_adj)
    df_fluid_adj = stats_adj$data
  
    # Save p-values
    if(high_detectability) {
      write_xlsx(pvals_adj,
                 paste0(output_dir, "/", prefix, "_adjusted_age_sex_pvalues_", fluid, "_high_detectability.xlsx"))
    } else{
      write_xlsx(pvals_adj,
                 paste0(output_dir, "/", prefix, "_adjusted_age_sex_pvalues_", fluid, ".xlsx"))
    }
   
    results[[fluid]] <- list(
      pvals = pvals,
      pvals_adjusted = pvals_adj
    )
    
    # Top N plot (violin)
    p_top_n <- plot_top_proteins_violin(df, stats, td, fluid, top_n = top_n, subtype = subtype)

    if(!high_detectability) {
      filespath = paste0(output_dir, "/boxplots_", fluid, "/")
      if (!file.exists(filespath)) {
        dir.create(filespath, recursive = TRUE)
      }
      width = max(ceiling(top_n / 5) * 0.6 + 6, 15)
      height = max(ceiling(top_n / 5) * 0.3 + 6, 12)
      ggsave(paste0(filespath, "/", prefix, "_Top",top_n,"_LOD_", fluid, ".pdf"),
             p_top_n, width = width, height = height)
    }
    
    # p_top15_adj <- plot_top_proteins_violin(df_fluid_adj, stats_adj, td, fluid,adjusted = TRUE)
    
    # if(!high_detectability) {
    #   filespath = paste0(output_dir, "/boxplots_", fluid, "/")
    #   if (!file.exists(filespath)) {
    #     dir.create(filespath, recursive = TRUE)
    #   }
    #   ggsave(paste0(filespath, "/", prefix, "_Top15_LOD_", fluid, "_adj.pdf"),
    #          p_top15_adj, width = 15, height = 18)
    # }
    
    # All proteins (multi-page PDF)
    targets <- unique(df$Target)
  
    if(!high_detectability) {
      filespath = paste0(output_dir, "/boxplots_", fluid, "/")
      if (!file.exists(filespath)) {
        dir.create(filespath, recursive = TRUE)
      }
      
      pdf(paste0(filespath, "/", prefix, "_ALLproteins_LOD_", fluid, "_others.pdf"),
          width = 6, height = 5.7)
      for (t in targets) {
        p <- plot_single_protein_violin(df, stats, td, fluid, t, subtype = subtype)
        print(p)
      }
      dev.off()
    }
    
    # targets <- unique(df_fluid_adj$Target)
   
    # if(!high_detectability) {
    #   filespath = paste0(output_dir, "/boxplots_", fluid, "/")
    #   if (!file.exists(filespath)) {
    #     dir.create(filespath, recursive = TRUE)
    #   }
    #   pdf(paste0(filespath, "/", prefix, "_ALLproteins_LOD_", fluid, "_others_adj.pdf"),
    #       width = 6, height = 5.7)
    #   for (t in targets) {
    #     p <- plot_single_protein_violin(df_fluid_adj, stats_adj, td, fluid, t,adjusted = TRUE)
    #     print(p)
    #   }
    #   dev.off()
    # }
    
  }
  
  return(results)
}

# Clean datasets (all samples OR PGMC subset)
prepare_dataset <- function(protein_data, sample_map) {
  protein_data %>%
    filter(SampleType == "Sample") %>%
    select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
    left_join(sample_map %>% rename(SampleName = `PatientID`), by = "SampleName") %>%
    filter(!is.na(type)) %>%
    mutate(
      type   = factor(type),
      sex    = factor(sex),
      #center = factor(center),
      age    = as.numeric(age)
    )
}

# Run ANOVA + pairwise tests for one fluid (CSF/SERUM/PLASMA)
format_p <- function(p, digits = 3) {
  ifelse(
    p < 10^(-digits),
    paste0("<", formatC(10^(-digits), format = "f", digits = digits)),
    formatC(p, format = "f", digits = digits)
  )
}

run_stats_for_fluid <- function(df, fluid, value_col = "NPQ", 
                                adjust = FALSE,
                                covariates = c("age", "sex")) {
  
  df_f <- df %>%
    filter(SampleMatrixType == fluid) %>%
    filter(!is.na(.data[[value_col]]))
  
  # formula for raw data or adjusted data 
  fmla <- if(adjust) {
    as.formula(paste(value_col, "~ type +", paste(covariates, collapse = " + ")))
  } else {
    as.formula(paste(value_col, "~ type"))
  }
  
  # ANOVA
  anova_res <- df_f %>%
    group_by(Target) %>%
    rstatix::anova_test(fmla, effect.size = "partial_eta_squared") %>%
    ungroup() %>%
    as_tibble() %>%
    mutate(Fluid = fluid, Adjustment = ifelse(adjust, "ADJ", "UNADJ")) %>% 
    filter(Effect == "type")
  
  
  if(nrow(anova_res) == 0) {
    warning("No ANOVA results for fluid: ", fluid)
    anova_res <- NULL
  }
  
  # Welch t-tests
  if(adjust){
    # Covariate-adjusted t-tests using linear model residuals
    pairwise_res <- df_f %>%
      group_by(Target) %>%
      do({
        sub_df <- .
        # Fit linear model with type + covariates
        fit <- lm(as.formula(paste(value_col, "~ type +", paste(covariates, collapse = " + "))), data = sub_df)
        # Compute residuals + group effects
        emmeans_fit <- try(emmeans::emmeans(fit, specs = "type"), silent = TRUE)
        if(inherits(emmeans_fit, "try-error")) return(tibble())
        pw <- emmeans::contrast(emmeans_fit, method = "pairwise") %>% as_tibble() %>%
          tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
          rename(p = p.value) %>%
          mutate(Target = unique(sub_df$Target))
        pw
      }) %>%
      ungroup() %>%
      group_by(group1, group2) %>%
      adjust_pvalue(method = "BH") %>%
      ungroup() %>%
      add_significance() %>%
      mutate(Fluid = fluid, Adjustment = "ADJ")
  } else {
    # Regular Welch t-tests without adjustment
    pairwise_res <- df_f %>%
      group_by(Target) %>%
      rstatix::t_test(as.formula(paste(value_col, "~ type")), var.equal = FALSE) %>%
      ungroup() %>%
      group_by(group1, group2) %>%
      adjust_pvalue(method = "BH") %>%
      ungroup() %>%
      add_significance() %>%
      mutate(Fluid = fluid, Adjustment = "UNADJ")
  }
  
  if(nrow(pairwise_res) == 0) {
    warning("No pairwise results for fluid: ", fluid)
    pairwise_res <- NULL
  }
  
  if(adjust){
    df_f <- df_f %>%
      group_by(Target) %>%
      mutate(NPQ_adj = local({
        sub_df <- cur_data()
        ok <- complete.cases(sub_df[, c(value_col, covariates)])
        res <- rep(NA_real_, nrow(sub_df))
        if(sum(ok) >= length(covariates) + 1){
          f <- as.formula(paste(value_col, "~", paste(covariates, collapse = " + ")))
          fit <- lm(f, data = sub_df[ok, , drop = FALSE])
          X <- model.matrix(fit)
          beta <- coef(fit)
          drop_cols <- colnames(X) %in% unlist(lapply(covariates, function(v) grep(paste0("^", v), colnames(X), value = TRUE)))
          res[ok] <- sub_df[[value_col]][ok] - as.numeric(X[, drop_cols, drop = FALSE] %*% beta[drop_cols])
        }
        res
      })) %>%
      ungroup()
  }
  
  pairwise_res <- pairwise_res %>%
    mutate(
      p_label = format_p(p, digits = 3)
    )
  
  list(
    anova = anova_res,
    pairwise = pairwise_res,
    data = df_f
  )
}


# Extract p-values for all targets
extract_pvalues <- function(stats_list) {
  
  pw <- stats_list$pairwise %>%
    mutate(
      comparison = paste(group1, group2, sep = "_")
    ) 
  
  pval_wide = pw %>%
    select(Fluid, Target, comparison, p) %>%
    pivot_wider(
      names_from = comparison,
      values_from = p,
      names_glue = "pvalue_{comparison}"
    )
  
  padj_wide = pw %>%
    select(Fluid, Target, comparison, p.adj) %>%
    pivot_wider(
      names_from = comparison,
      values_from = p.adj,
      names_glue = "padj_{comparison}"
    )
  
  signif_wide = pw %>%
    select(Fluid, Target, comparison, p.adj.signif) %>%
    pivot_wider(
      names_from = comparison,
      values_from = p.adj.signif,
      names_glue = "padj_signif_{comparison}"
    )
  
  combined <- pval_wide %>%
    left_join(padj_wide, by = c("Fluid","Target")) %>%
    left_join(signif_wide, by = c("Fluid","Target")) 
  
  comparisons <- unique(pw$comparison)
  
  ordered_cols <- c(
    "Fluid", "Target",
    unlist(lapply(comparisons, function(cmp) {
      c(
        paste0("pvalue_", cmp),
        paste0("padj_", cmp),
        paste0("padj_signif_", cmp)
      )
    }))
  )
  
  final_table <- combined %>% select(any_of(ordered_cols))
  
  pw_final = stats_list$anova %>% select(Target,pvalue_anova = p) %>%
    left_join(final_table)  %>%
    arrange(pvalue_anova)
  
  pw_final = pw_final[,c(1,3,2,4:ncol(pw_final))]
 
  return(pw_final)
}

# Top N significantly changing proteins 
plot_top_proteins_violin <- function(df, stats_list, td, fluid, top_n = 15, subtype = FALSE) {
  
  # Check ANOVA results
  anova_res <- stats_list$anova
  
  if(is.null(anova_res) || nrow(anova_res) == 0) {
    warning("No ANOVA results for ", fluid, "; skipping plot.")
    return(NULL)
  }
  
  # Select top proteins
  top_proteins_df <- anova_res %>%
    arrange(p) %>%
    slice(1:min(top_n, nrow(.)))
  
  top_proteins <- top_proteins_df$Target

  if (subtype == FALSE) {
    df_top <- df %>% filter(SampleMatrixType == fluid, 
                            Target %in% top_proteins,
                            type %in% c("ALS","CTRL"))
  } else {
    df_top <- df %>% filter(SampleMatrixType == fluid, 
                            Target %in% top_proteins,
                            type %in% c("CTRL","alpha","beta"))
  }
  
  if(nrow(df_top) == 0) {
    warning("No data for top proteins in ", fluid, "; skipping plot.")
    return(NULL)
  }
  
  # Pairwise p-values
  pairwise_res <- get_pairwise_sig(stats_list, fluid, top_proteins, df_top)
  
  # # LOD per protein (original LOD-plate)
  # lod_df <- df_top %>%
  #   group_by(Target) %>%
  #   summarise(
  #     LOD_val = get_lod(Target[1], fluid, td),
  #     max_val = max(NPQ, na.rm = TRUE)
  #   ) %>%
  #   mutate(y_position_text = LOD_val * 1.02) %>%
  #   ungroup()
  # 
  
  # LOD per protein (Project-LOD)
  lod_df <- df_top %>%
    left_join(td) %>%
    group_by(Target) %>%
    summarise(
      LOD_val = ProjectLOD %>% unique(),
      max_val = max(NPQ, na.rm = TRUE)
    ) %>%
    mutate(y_position_text = LOD_val * 1.02) %>%
    ungroup()
  
  ALS_CTRL <- c("ALS","CTRL")
  df_top <- df_top %>%
    mutate(type = case_when(
      type %in% ALS_CTRL ~ factor(type, levels = c("CTRL","ALS")),
      TRUE ~ factor(type, levels = c("alpha", "beta"))
    ))
  df_top$type = droplevels(df_top$type)

  # plot
  p <- ggplot(df_top, aes(x = type, y = NPQ, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(data = subset(df_top, type != "others"),
                aes(color = subtype), width = 0.15, alpha = 0.5, size = 2) +
    geom_jitter(data = subset(df_top, type == "others"),
                aes(color = subtype), position = position_dodge(width = 0.8), alpha = 0.5, size = 2) +
    facet_wrap(~Target, scales = "free_y", ncol = min(top_n, 5)) +
    scale_fill_manual(values  = c('CTR' = 'black',  
                                    'ALS' = '#D73027',
                                    'alpha' = '#FC8D62',
                                    'beta' = '#8DA0CB')) +
    scale_color_manual(values  = c('CTR' = 'black',  
                                    'ALS' = '#D73027',
                                    'alpha' = '#FC8D62',
                                    'beta' = '#8DA0CB'))  +
    labs(
      x = "Group",
      y = "NPQ",
      title = paste("Top", top_n, "Protein Expression -", fluid)
      ) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title.x = element_blank(),
      text = element_text(size = 15),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 14),          # Axis tick labels
      strip.text = element_text(size = 16, face = "bold"),  # Facet labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "none"
    ) 

  
  
  # Add LOD lines and labels
  p <- p +
    geom_hline(
      data = lod_df %>% tidyr::unnest(LOD_val),
      aes(yintercept = LOD_val),
      linetype = "dashed",
      color = "gray55",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = lod_df  %>% tidyr::unnest(LOD_val) %>% mutate(label = paste("LOD =", signif(LOD_val, 3))),
      aes(x = 2, y = y_position_text, label = label),
      color = "gray55",
      size = 4,
      hjust = 0,
      inherit.aes = FALSE
    )
    
  # Add pairwise p-values only if available
  if (!is.null(pairwise_res) && nrow(pairwise_res) > 0) {
    p <- p + stat_pvalue_manual(
      pairwise_res,
      xmin = "group1",
      xmax = "group2",
      label = "p_label",
      y.position = "y.position",
      size = 5,
      bracket.nudge.y = 0.02 * max(df_top$NPQ, na.rm = TRUE),
      inherit.aes = FALSE
    )
  }
  
  return(p)
}

# Get pairwise p-values filtered by cutoff for plotting
get_pairwise_sig <- function(stats_list, fluid, proteins, df_top, p_cutoff = 0.1) {
  
  pairwise_df <- stats_list$pairwise %>%
    filter(Target %in% proteins, p <= p_cutoff)
  
  if(nrow(pairwise_df) == 0) return(pairwise_df)
  
  # Max NPQ per protein 
  max_vals <- df_top %>%
    group_by(Target) %>%
    summarise(
      max_y = max(NPQ, na.rm = TRUE),
      y_step = diff(range(NPQ, na.rm = TRUE)) * 0.08,
      .groups = "drop"
    )
  
  # Assign incremental y positions
  pairwise_df <- pairwise_df %>%
    left_join(max_vals, by = "Target") %>%
    group_by(Target) %>%
    mutate(
      comparison_index = row_number(),
      y.position = max_y + comparison_index * y_step
    ) %>%
    ungroup()
  
  return(pairwise_df)
}

# Single proteins
plot_single_protein_violin <- function(df, stats_list, td, fluid, protein, subtype = FALSE) {
  
  df_prot <- df %>% filter(SampleMatrixType == fluid, 
                           Target == protein,
                           type %in% c("ALS","CTRL","alpha","beta"))
  anova_p <- stats_list$anova %>% filter(Target == protein) %>% pull(p)
  pairwise_res <- get_pairwise_sig(stats_list, fluid, protein, df_top = df_prot)
  
  # Determine type order
  if (subtype == FALSE) {
    order_types <- c("CTRL","ALS")
  } else {
    order_types <- c("CTRL","alpha","beta")
  }
  
  df_prot = df_prot %>%
    mutate(type = case_when(
      type %in% order_types ~ factor(type, levels = order_types),
      TRUE ~ factor(type, levels = order_types)
    ))
  
  # # LOD for this protein (original LOD-plate)
  # LOD_val <- get_lod(protein, fluid, td) %>% unique()
  
  # LOD for this protein (Project-LOD)
  LOD_val = get_project_lod(protein, fluid, td) %>% unique()
  lod_df <- df_prot %>%
    summarise(
      LOD_val = get_project_lod(protein, fluid, td),
      #max_val = ifelse(adjusted, max(NPQ_adj, na.rm = TRUE),max(NPQ, na.rm = TRUE))
      max_val = max(NPQ, na.rm = TRUE)
    ) %>%
    mutate(y_position_text = LOD_val * 1.02) %>%
    ungroup()
  
  max_val <- max(df_prot$NPQ, na.rm = TRUE)
  y_text <- max_val * 1.02
  
  p <- ggplot(df_prot, aes(x = type, y = NPQ, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(data = subset(df_prot, type != "others"),
                aes(color = subtype), width = 0.15, alpha = 0.5, size = 2) +
    geom_jitter(data = subset(df_prot, type == "others"),
                aes(color = subtype), position = position_dodge(width = 0.8), alpha = 0.5, size = 2) +
    scale_fill_manual(values  = c('CTR' = 'black',  
                                  'ALS' = '#D73027',
                                  'alpha' = '#FC8D62',
                                  'beta' = '#8DA0CB')) +
    scale_color_manual(values  = c('CTR' = 'black',  
                                   'ALS' = '#D73027',
                                   'alpha' = '#FC8D62',
                                   'beta' = '#8DA0CB')) +
    labs(
      x = "Group",
      y = "NPQ",
      title = paste(protein),
      subtitle = paste("ANOVA p =", signif(anova_p, 3),
                       "; LOD =", signif(LOD_val,3))) +
    theme_test() + 
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      axis.title.x = element_blank(),
      text = element_text(size = 15),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 14),          # Axis tick labels
      strip.text = element_text(size = 16, face = "bold"),  # Facet labels
      plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
      plot.subtitle = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "none"
    )
  
  # p <- p +
  #   annotation_custom(
  #     grob = grid::rectGrob(
  #       gp = grid::gpar(fill = "white", col = "black", lwd = 1)
  #     ),
  #     ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf
  #   ) +
  #   annotation_custom(
  #     grob = grid::textGrob(
  #       label = protein,
  #       gp = grid::gpar(fontsize = 16, fontface = "bold")
  #     ),
  #     ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf
  #   ) +
  #   theme(
  #     plot.margin = margin(t = 30),  # extra room for strip
  #     strip.background = element_blank(),
  #     strip.text = element_blank()
  #   )
    
  if (nrow(pairwise_res) > 0) {
    p <- p + stat_pvalue_manual(
      pairwise_res,
      xmin = "group1",
      xmax = "group2",
      label = "p_label",
      y.position = "y.position",
      bracket.nudge.y = 0.05 * max_val,
      size = 6,
      inherit.aes = FALSE
    )
  }
  
  # Add LOD line and label (once only; use scalar LOD_val to avoid overlap)
  lod_scalar <- if (is.list(lod_df$LOD_val)) unlist(lod_df$LOD_val)[1] else lod_df$LOD_val[1]
  p <- p +
    geom_hline(
      yintercept = lod_scalar,
      linetype = "dashed",
      color = "gray55",
      inherit.aes = FALSE
    ) +
    annotate(
      "text",
      x = 2,
      y = lod_scalar * 1.02,
      label = paste0("LOD = ", signif(lod_scalar, 3)),
      color = "gray55",
      hjust = 0,
      size = 4,
      inherit.aes = FALSE
    )

  p <- p + scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)))
  
  return(p)
}

# Get LOD for a given protein and fluid 
get_lod <- function(protein, fluid, td) {
  lod_val <- td %>%
    filter(Target == protein, SampleMatrixType == fluid) %>%
    pull(TargetLOD_NPQ)
  if(length(lod_val) == 0) return(NA) else return(lod_val)
}

get_project_lod <- function(protein, fluid, td) {
  lod_val <- td %>%
    filter(Target == protein, SampleMatrixType == fluid) %>%
    pull(ProjectLOD)
  if(length(lod_val) == 0) return(NA) else return(lod_val)
}

#=================Parse command line arguments=================
option_list = list(
  make_option(c("-i", "--input"), type="character", default="Results/00_Initialization/npq_counts.xlsx", help="Input NPQ counts file"),
  make_option(c("-o", "--output"), type="character", default="Results/02_Analysis_Groups", help="Output analysis groups path"),
  make_option(c("-m","--metadata"), type="character", default="Results/00_Initialization/all_participants_IDs.xlsx", help="Input all participants IDs file"),
  make_option(c("-t","--target_detectability_extra"), type="character", default="Results/01_Data_Mining/target_detectability_extra.xlsx", help="Input target detectability file with project-LOD"),
  make_option(c("-d","--td"), type="character", default="Results/01_Data_Mining/td.xlsx", help="Input target detectability file"),
  make_option(c("-s","--subtype"), type="logical", default=FALSE, help="If TRUE, analyze subtype (alpha, beta, CTRL) instead of ALS and CTRL"),
  make_option(c("-n","--top_n"), type="integer", default=15, help="Top N proteins to plot")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#=================options setting=================
if (is.null(opt$input)) {
  print("NO INPUT NPQ COUNTS FILE SUPPLIED, EXITING!")
  stop("Please provide the input NPQ counts file path!")
} else {
  input_file <- opt$input
  protein_data <- read_excel(input_file)
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

if (is.null(opt$target_detectability_extra)) {
  print("NO TARGET DETECTABILITY FILE PATH SUPPLIED, EXITING!")
  stop("Please provide the target detectability file path!")
} else {
  target_detectability_extra_file <- opt$target_detectability_extra
  target_detectability_extra <- read_excel(target_detectability_extra_file)
}

if (is.null(opt$td)) {
  print("NO TD FILE PATH SUPPLIED, EXITING!")
  stop("Please provide the td file path!")
} else {
  td_file <- opt$td
  td <- read_excel(td_file)
}

if (is.null(opt$subtype)) {
  print("NO SUBTYPE OPTION SUPPLIED, EXITING!")
  stop("Please provide the subtype option!")
} else {
  subtype <- opt$subtype
}

if (is.null(opt$top_n)) {
  print("NO TOP N OPTION SUPPLIED, EXITING!")
  stop("Please provide the top n option!")
} else {
  top_n <- opt$top_n
}

###############################################################################
# Run pipeline
###############################################################################
## 1. All samples
# If subtype = TRUE, set type to alpha or beta from k2; otherwise keep as is
if (subtype) {
samples_ID_type <- samples_ID_type %>%
  mutate(type = case_when(
    type == "ALS" & k2 == "alpha" ~ "alpha",
    type == "ALS" & k2 == "beta"  ~ "beta",
    TRUE ~ type
  ))
}

results_ALL <- run_full_pipeline(
  protein_data    = protein_data,
  sample_map      = samples_ID_type %>% 
    mutate(subtype = type) %>%
    dplyr::rename(PatientID = Tube_ID),
  td = td %>% left_join(target_detectability_extra %>% select(Target,ProjectLOD) %>% distinct()) %>%
    dplyr::select(SampleMatrixType,Target,ProjectLOD) %>% distinct(),
  prefix          = "ALLsamples",
  output_dir      = output_dir,
  subtype         = subtype
)

saveRDS(results_ALL, paste0(output_dir, "/results_ALL.rds"))