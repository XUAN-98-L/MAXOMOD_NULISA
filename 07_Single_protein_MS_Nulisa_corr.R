# MS (Validation cohort, Supp table 13) vs NULISA (NPQ) per protein — paired by patient ID
# Same individuals: VC MS columns (CSF1…) map via Patient.ID in clinical_data_with_cluster_VC.csv;
# NULISA uses ID from all_participants_IDs (Tube_ID differs from VC Tube.ID but patient matches).

suppressMessages(library(optparse))
library(readxl)
library(dplyr)
library(tidyr)

root <- "/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD"
default_init <- file.path(root, "CNS_immune/Results/CNS_panel/Without_tears/00_Initialization")
default_out <- file.path(root, "CNS_immune/Results/CNS_panel/Without_tears/07_MS_NULISA_correlation")

option_list <- list(
  make_option(
    c("-i", "--init"),
    type = "character",
    default = default_init,
    help = "Folder with npq_counts.xlsx and all_participants_IDs.xlsx [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = default_out,
    help = "Output directory for correlation CSV [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

init_dir <- opt$init
output_dir <- opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

npq_counts <- read_excel(file.path(init_dir, "npq_counts.xlsx"))
npq_counts <- npq_counts %>% filter(SampleMatrixType == "CSF")

all_participants_IDs <- read_excel(file.path(init_dir, "all_participants_IDs.xlsx"))
NULISA_CNS_CSF <- npq_counts %>%
  left_join(all_participants_IDs, by = c("SampleName" = "Tube_ID"))

VC <- read.csv(
  "/Users/xliu2942/Documents/Projects/MAXOMOD/14324_0_source_data_169635_t9mdsb/Supp table 13 - intensity_imputed_log2transf_norm_VC.csv",
  check.names = FALSE
)

inter <- intersect(npq_counts$Target, VC$Protein)
npq_counts <- npq_counts %>% filter(Target %in% inter)
VC <- VC %>% filter(Protein %in% inter)
rownames(VC) <- VC$Protein
VC <- VC %>% select(-Protein)

# Patient.ID → MS matrix column (CSF1, …)
vc_clin <- read.csv(file.path(root, "Data/clinical_data_with_cluster_VC.csv"), check.names = FALSE)
vc_map <- vc_clin %>%
  transmute(
    patient_id = as.character(Patient.ID),
    ms_col = as.character(`CSF ID`)
  ) %>%
  filter(ms_col %in% colnames(VC)) %>%
  distinct(patient_id, .keep_all = TRUE)

nul_ids <- unique(as.character(NULISA_CNS_CSF$ID))
common_ids <- sort(intersect(nul_ids, vc_map$patient_id))

if (length(common_ids) < 3L) {
  stop(
    "Need at least 3 paired subjects. Tube_ID does not match VC MS columns; ",
    "pairing uses Patient.ID. Found ", length(common_ids), " common IDs."
  )
}

ms_for_patient <- function(pid) {
  vc_map$ms_col[match(pid, vc_map$patient_id)]
}

# NULISA: one NPQ per Target × patient (ID)
nul_wide <- NULISA_CNS_CSF %>%
  filter(Target %in% inter, as.character(ID) %in% common_ids) %>%
  group_by(Target, ID) %>%
  summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = ID, values_from = NPQ)

corr_rows <- lapply(inter, function(prot) {
  ms_y <- vapply(common_ids, function(pid) {
    mc <- ms_for_patient(pid)
    if (is.na(mc) || !mc %in% colnames(VC)) {
      return(NA_real_)
    }
    as.numeric(VC[prot, mc])
  }, numeric(1))

  nu_row <- nul_wide %>% filter(Target == prot)
  if (nrow(nu_row) == 0L) {
    nu_y <- rep(NA_real_, length(common_ids))
  } else {
    nu_y <- as.numeric(unlist(nu_row[, common_ids, drop = FALSE]))
  }

  ok <- is.finite(ms_y) & is.finite(nu_y)
  n_pair <- sum(ok)
  if (n_pair < 3L) {
    return(tibble(
      Target = prot,
      n = n_pair,
      spearman_r = NA_real_,
      spearman_p = NA_real_,
      pearson_r = NA_real_,
      pearson_p = NA_real_
    ))
  }
  ms_ok <- ms_y[ok]
  nu_ok <- nu_y[ok]
  ct_s <- suppressWarnings(cor.test(ms_ok, nu_ok, method = "spearman", exact = FALSE))
  ct_p <- suppressWarnings(cor.test(ms_ok, nu_ok, method = "pearson"))
  tibble(
    Target = prot,
    n = n_pair,
    spearman_r = unname(ct_s$estimate),
    spearman_p = ct_s$p.value,
    pearson_r = unname(ct_p$estimate),
    pearson_p = ct_p$p.value
  )
})

corr_tbl <- bind_rows(corr_rows)
out_path <- file.path(output_dir, "MS_NULISA_VC_per_protein_correlation.csv")
write.csv(corr_tbl, out_path, row.names = FALSE)
message("Paired subjects (Patient.ID): ", length(common_ids))
message("Wrote ", normalizePath(out_path, mustWork = FALSE))
