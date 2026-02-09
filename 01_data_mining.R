### Libraries
suppressMessages(library("optparse"))
suppressMessages(library("readxl"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("ggExtra"))
suppressMessages(library("cowplot"))
#=================Functions=================
# Collect sample counts per fluid
count_samples_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid, !is.na(type)) %>%
    distinct(SampleName, type) %>%
    group_by(type) %>%
    summarise(nr_samples = n(), .groups = "drop") %>%
    mutate(biofluid = fluid)
}

make_corr <- function(df1, df2, label1, label2, filename) {
  merged <- inner_join(df1, df2, by = c("Target", "UniProtID"))
  p <- corr_plot(merged,
                 paste0("NPQ_", label1),
                 paste0("NPQ_", label2),
                 paste("Correlation", label1, "vs", label2))
  
  pdf(filename)
  print(p)
  dev.off()
  p
}

get_mean_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid) %>%
    distinct(Target, UniProtID, NPQ) %>%
    group_by(Target, UniProtID) %>%
    summarise(NPQ = mean(NPQ), .groups = "drop") %>%
    rename_with(~ paste0("NPQ_", fluid), "NPQ")
}

# Correlation plot generator
corr_plot <- function(data, x, y, title) {
  
  cor_test <- cor.test(data[[x]], data[[y]], use = "complete.obs")
  r_val <- round(cor_test$estimate, 2)
  p_val <- ifelse(cor_test$p.value < 0.001,
                  formatC(cor_test$p.value, format = "e", digits = 2),
                  round(cor_test$p.value, 3))
  
  cor_text <- paste0("italic(R) == ", r_val,
                     " * ',' ~ italic(p) == ", p_val)
  
  p <- ggscatter(data, x = x, y = y,
                 add = "reg.line",
                 add.params = list(color = "blue",fill = "lightblue"),
                 conf.int = TRUE) +
    annotate("text",
             x = min(data[[x]], na.rm = TRUE),
             y = max(data[[y]], na.rm = TRUE),
             hjust = 0,
             label = cor_text,
             parse = TRUE, size = 5) +
    labs(title = title, x = x, y = y) +
    theme( plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
           axis.title = element_text(size = 16), 
           axis.text = element_text(size = 14), 
           legend.title = element_text(size = 14), 
           legend.text = element_text(size = 14))
  
  ggMarginal(p, type = "histogram", color = "lightblue", fill = "lightblue")
}


# etectability by fluid and plate
plot_detectability <- function(df, fluid, by_plate = FALSE) {
  
  df_f <- df %>% filter(SampleMatrixType == fluid)
  
  if (by_plate) {
    df_f <- df_f %>%
      group_by(PlateID, Target) %>%
      mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
             mean_pct   = mean(TargetDetectability, na.rm = TRUE)) %>%
      ungroup()
  } else {
    df_f <- df_f %>%
      group_by(Target) %>%
      mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
             mean_pct   = mean(TargetDetectability, na.rm = TRUE)) %>%
      ungroup()
  }
  
  df_f <- df_f %>%
    mutate(Target = factor(Target,
                           levels = unique(Target[order(mean_value, decreasing = TRUE)])
    ))
  
  p <- ggplot(df_f,
              aes(x = Target, y = TargetDetectability_value,
                  fill = mean_pct)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_gradient(low = "blue", high = "orange") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Target Detectability -", fluid),
      x = "Target",
      y = "Detectability (NPQ - LOD)",
      fill = "Detectability (%)"
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  if (by_plate) p <- p + facet_wrap(~PlateID, ncol = 1)
  
  p
}

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
#=================Parse command line arguments=================
option_list = list(
  make_option(c("-i", "--input"), type="character", default="Results/00_Initialization/npq_counts.xlsx", help="Input NPQ counts file"),
  make_option(c("-o", "--output"), type="character", default="Results/01_Data_Mining", help="Output data mining path"),
  make_option(c("-m","--metadata"), type="character", default="Results/00_Initialization/all_participants_IDs.xlsx", help="Input all participants IDs file"),
  make_option(c("-t","--target_detectability"), type="character", default="Results/00_Initialization/target_detectability.xlsx", help="Input target detectability file")
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

if (is.null(opt$metadata)) {
  print("NO METADATA FILE PATH SUPPLIED, EXITING!")
  stop("Please provide the metadata file path!")
} else {
  metadata_file <- opt$metadata
  samples_ID_type <- read_excel(metadata_file)
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

if (is.null(opt$target_detectability)) {
  print("NO TARGET DETECTABILITY FILE PATH SUPPLIED, EXITING!")
  stop("Please provide the target detectability file path!")
} else {
  target_detectability_file <- opt$target_detectability
  target_detectability <- read_excel(target_detectability_file)
}
#=================Data Processing=================
# Merge protein data with sample type
protein_data_IDs <- protein_data %>%
  filter(SampleType == "Sample") %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  left_join(samples_ID_type %>% rename(SampleName = `Tube_ID`), by = "SampleName")

# target_detectability

# get count table of APOE
table_APOE = protein_data_IDs %>%
    rename(subtype = genetics)  %>%
    mutate(
    subtype = ifelse(subtype == "CTRL" | is.na(subtype), "No subtype",subtype)) %>%
  filter(Target == "APOE4") %>%
  mutate(APOE_status =  ifelse(NPQ > 10, "Carrier",ifelse(NPQ <= 10, "Non-carrier",NA)))

carrier_table <- table_APOE %>%
  group_by(SampleMatrixType, type, subtype, APOE_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = APOE_status,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total = Carrier + `Non-carrier`
  )

writexl::write_xlsx(carrier_table,paste0(output_dir, "/APOE4_carriers.xlsx"))

###############################################
# Sample Counts Across Fluids
###############################################

fluids <- c("SERUM", "PLASMA", "CSF", "TEARS")

sample_counts <- bind_rows(
  lapply(fluids, function(f) count_samples_per_fluid(protein_data_IDs, f))
)

writexl::write_xlsx(sample_counts,paste0(output_dir, "/samples_biofluid_overview.xlsx"))

###############################################
# Correlation Between Fluids
###############################################
df_SERUM  <- get_mean_per_fluid(protein_data_IDs, "SERUM")
df_PLASMA <- get_mean_per_fluid(protein_data_IDs, "PLASMA")
df_CSF    <- get_mean_per_fluid(protein_data_IDs, "CSF")
df_TEARS  <- get_mean_per_fluid(protein_data_IDs, "TEARS")

p1 = make_corr(df_PLASMA, df_SERUM,  "PLASMA", "SERUM", paste0(output_dir, "/correlation_plasma_serum.pdf"))
p2 = make_corr(df_SERUM,  df_CSF,    "SERUM",  "CSF",   paste0(output_dir, "/correlation_serum_CSF.pdf"))
p3 = make_corr(df_PLASMA, df_CSF,    "PLASMA", "CSF",   paste0(output_dir, "/correlation_plasma_CSF.pdf"))
p4 = make_corr(df_TEARS,  df_SERUM,    "TEARS",  "SERUM",   paste0(output_dir, "/correlation_tears_serum.pdf"))
p5 = make_corr(df_TEARS,  df_PLASMA,    "TEARS",  "PLASMA",   paste0(output_dir, "/correlation_tears_plasma.pdf"))
p6 = make_corr(df_TEARS,  df_CSF,    "TEARS",  "CSF",   paste0(output_dir, "/correlation_tears_csf.pdf"))

combined <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, align = "h") 
height = length(combined)/2 * 2
width = length(combined)/2 * 3
pdf(paste0(output_dir, "/correlation_combined.pdf"), width = width, height = height) 
print(combined) 
dev.off()

###############################################
#Target Detectability per Fluid & Plate
###############################################
td <- protein_data  %>%
  dplyr::select(PlateID, SampleMatrixType, Target, NPQ) %>%
  left_join(target_detectability %>% rename(TargetLOD_NPQ = Target_LOD, TargetDetectability = Target_Detectability), by = c("Target","PlateID")) %>%
  mutate(TargetDetectability_value = NPQ - TargetLOD_NPQ,
         TargetDetectability =
           as.numeric(TargetDetectability) * 100) %>%
  mutate(TargetDetectability_value =
           ifelse(Target %in% c("APOE", "CRP"),
                  abs(TargetDetectability_value),
                  TargetDetectability_value))

writexl::write_xlsx(td,paste0(output_dir, "/td.xlsx"))

fluids <- c("CSF", "PLASMA", "SERUM", "TEARS")

# Save plots
for (fluid in fluids) {
  p1 <- plot_detectability(td, fluid, by_plate = FALSE)
  p2 <- plot_detectability(td, fluid, by_plate = TRUE)
  
  pdf(paste0(output_dir, "/TargetDetectability_", fluid, ".pdf"),
                  width = 28, height = 6)
  print(p1)
  dev.off()
  
  pdf(paste0(output_dir, "/TargetDetectability_byPlate_", fluid, ".pdf"),
                  width = 28, height = 8)
  print(p2)
  dev.off()
}

#################################################################
# NPQ Distributions by Plate original data
#################################################################
protein_long_alltogether <- protein_data %>%
  mutate(
    SampleMatrixType = factor(
      SampleMatrixType,
      levels = c("PLASMA", "SERUM", "CSF", "TEARS", "CONTROL")
    )
  ) %>%
  dplyr::select(
    Target, SampleMatrixType, PlateID, NPQ
  ) %>%
  rename(value = NPQ)

unique_targets <- unique(protein_data$Target)

for (target in unique_targets) {
  
  df_t_alltogether <- protein_long_alltogether %>%
    filter(Target == target)
  
  p_alltogether <- ggboxplot(
    df_t_alltogether,
    x = "SampleMatrixType",
    y = "value",
    color = "PlateID",
    add = "jitter"
  ) +
    scale_color_manual(values = my_palette) +
    labs(
      title = target,
      y = "NPQ"
    ) +
    theme_bw()
  
  if (!file.exists(paste0(output_dir, "/NPQ_fluid_plate"))) {
    dir.create(paste0(output_dir, "/NPQ_fluid_plate"), recursive = TRUE)
  }
  pdf(paste0(output_dir, "/NPQ_fluid_plate/NPQ_", target, ".pdf"),
      width = 10, height = 5)
  print(p_alltogether)
  dev.off()
}

#################################################################
# Project-LOD of original data
#################################################################
# Project-LOD non-normalized data
df_reads <- protein_data %>%
  filter(SampleType == "NC") %>%
  mutate(
    reads = 2^NPQ - 1
  )

lod_linear <- df_reads %>%
  group_by(Target) %>%
  summarise(
    mean_reads = mean(reads, na.rm = TRUE),
    sd_reads   = sd(reads, na.rm = TRUE),
    lod_reads  = mean_reads + 3 * sd_reads,
    .groups = "drop"
  )

lod_project <- lod_linear %>%
  mutate(
    LOD_NPQ = log2(lod_reads + 1)
  ) %>%
  select(Target, LOD_NPQ)

## Attached to original target detectability
target_detectability_extra = target_detectability %>%
  #dplyr::rename(Target = TargetName) %>%
  #left_join(lod_project_norm) %>%
  left_join(lod_project) %>%
  dplyr::rename(original_TargetLOD = Target_LOD,
                #ProjectLOD_norm = LOD_NPQ_norm,
                ProjectLOD = LOD_NPQ) %>%
  arrange(Target)

writexl::write_xlsx(target_detectability_extra,paste0(output_dir, "/target_detectability_extra.xlsx"))
#################################################################
#  Check targets' detectability across samples
#################################################################
protein_data_with_lod <- protein_data %>%
  left_join(target_detectability_extra %>%
              dplyr::select(Target,ProjectLOD) %>% 
              distinct()) %>%
  mutate(below_lod = ifelse(Target %in% c("APOE","CRP"),
                NPQ > ProjectLOD,
               NPQ < ProjectLOD))

detectability_summary <- protein_data_with_lod %>%
  group_by(SampleMatrixType, Target) %>%
  summarise(
    n_samples = n(),
    n_below_lod = sum(below_lod, na.rm = TRUE),
    frac_below_lod = n_below_lod / n_samples,
    .groups = "drop"
  ) %>%
  filter(SampleMatrixType != "CONTROL") %>%
  arrange(Target,frac_below_lod) %>%
  mutate(percent_below_lod = frac_below_lod * 100,
         detectability = ifelse(frac_below_lod>0.5,"low","high"))

writexl::write_xlsx(detectability_summary,paste0(output_dir, "/detectability_summary.xlsx"))

targets_below_LOD = detectability_summary %>%
  select(SampleMatrixType, Target, percent_below_lod) %>%
  tidyr::pivot_wider(
    names_from = SampleMatrixType,
    values_from = percent_below_lod
  )

writexl::write_xlsx(targets_below_LOD,paste0(output_dir, "/targets_below_LOD.xlsx"))