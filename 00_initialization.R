### Libraries
suppressMessages(library("dplyr"))
suppressMessages(library("readxl"))
suppressMessages(library("optparse"))
#=================Functions=================
#=================Parse command line arguments=================
option_list = list(
  make_option(c("-i", "--input"), type="character", default="Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx", help="Input NPQ counts file"),
  make_option(c("-o", "--output"), type="character", default="Results/00_Initialization", help="Output initialization path"),
  make_option(c("-m","--metadata"), type="character", default="Data/NULISA_MAXOMOD_TearALS_final.xlsx", help="Input metadata file"),
  make_option(c("-c","--cluster"), type="character", default="Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv", help="Clustering allocation file from MAXOMOD project Discovery and Validation cohorts")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#=================options setting=================
if (is.null(opt$input)) {
  print("NO INPUT NPQ COUNTS FILE SUPPLIED, EXITING!")
  stop("Please provide the input NPQ counts file path!")
} else {
  input_file <- opt$input
  npq_counts <- read_excel(input_file,sheet = "NPQ Counts")
  # Aβ38, Aβ40, Aβ42 display error
  npq_counts$Target[npq_counts$Target == "AÎ²38"] = "Aβ38"
  npq_counts$Target[npq_counts$Target == "AÎ²40"] = "Aβ40"
  npq_counts$Target[npq_counts$Target == "AÎ²42"] = "Aβ42"
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
}

if (is.null(opt$cluster)) {
  print("NO CLUSTER FILE PATH SUPPLIED,current directory will be used!")
  stop("Please provide the cluster file path!")
} else {
  cluster_file <- opt$cluster
  cluster_file_DC <- strsplit(cluster_file, ",")[[1]][1]
  cluster_file_VC <- strsplit(cluster_file, ",")[[1]][2]
  cluster_info_DC <- read.csv(cluster_file_DC, sep = ",",check.names = FALSE)
  cluster_info_VC <- read.csv(cluster_file_VC, sep = ",",check.names = FALSE)
}

#=================Data Processing=================
# target detectability
target_detectability <- read_excel(input_file, 
           sheet = "Target Detectability") %>% rename(Target_Detectability = "Target Detectability",Target_LOD = "Target LOD (NPQ)")
writexl::write_xlsx(target_detectability, paste0(output_dir, "/target_detectability.xlsx"))

# information on sample IDs
sample_ID_info <- read_excel(input_file, 
          sheet = "Sample Information") %>%
  dplyr::select("Sample Name","PlateID","Well Position","Patient ID (from manifest - for Plate 01)") %>% rename(SampleName = "Sample Name",PlateID = "PlateID",WellPosition = "Well Position", Mapping_ID = "Patient ID (from manifest - for Plate 01)")

# IDS of patients from different fluids
CSF <- read_excel(metadata_file, 
          sheet = "CSF") %>% dplyr::select("Patient","ID","Tube-ID","Material","Rack position Nick","Rack number") %>% rename(Tube_ID = "Tube-ID",Rack_Position_Nick = "Rack position Nick",Rack_number = "Rack number")
Serum <- read_excel(metadata_file, 
          sheet = "Serum") %>% dplyr::select("Patient","ID","Tube-ID","Material","Rack position Nick","Rack number") %>% rename(Tube_ID = "Tube-ID",Rack_Position_Nick = "Rack position Nick",Rack_number = "Rack number")
Plasma <- read_excel(metadata_file, 
          sheet = "Plasma") %>% dplyr::select("Patient","ID","Tube-ID","Material","Rack position Nick","Rack number") %>% rename(Tube_ID = "Tube-ID",Rack_Position_Nick = "Rack position Nick",Rack_number = "Rack number")
Tears <- read_excel(metadata_file, 
          sheet = "Tears") %>% dplyr::select("Patient","ID","Tube-ID","Material","Rack position Nick","Rack number") %>% rename(Tube_ID = "Tube-ID",Rack_Position_Nick = "Rack position Nick",Rack_number = "Rack number")

all_participants_IDs = do.call("rbind",
                               list(CSF,
                                    Serum,
                                    Plasma,
                                    Tears))

# remove interial control from npq_counts
#npq_counts = npq_counts %>% filter(SampleMatrixType != "CONTROL")

#Assign Tear fluid SampleName back into Tube_ID in Tears
Tears = Tears %>% left_join(sample_ID_info , by = c("Patient" = "Mapping_ID")) %>%
    dplyr::filter(!is.na(SampleName))

# Replace npq_counts$SampleName with Tears_sample$Tube_ID where SampleName matches
# Take only the first Tube_ID for each SampleName+PlateID combination to avoid many-to-many relationship (Tears has multiple Tube_ID for the same SampleName+PlateID, TearsL or TearsR. Now it's all mapped to TearsL.)
Tears_sample <- Tears %>% 
    dplyr::select(SampleName, PlateID, Tube_ID) %>%
    group_by(SampleName, PlateID) %>%
    slice_head(n = 1) %>%
    ungroup()

npq_counts <- npq_counts %>% 
    left_join(Tears_sample, by = c("SampleName" = "SampleName", "PlateID" = "PlateID")) %>%
    mutate(SampleName = ifelse(!is.na(Tube_ID), Tube_ID, SampleName)) %>%
    dplyr::select(-Tube_ID)

# leave only the participants IDs that are in the npq_counts
all_participants_IDs <- all_participants_IDs %>% dplyr::filter(Tube_ID %in% npq_counts$SampleName)

# merge clinical info back into all_participants_IDs
# Sample_info from MAXOMOD project Validation cohort
cluster_info_VC <- cluster_info_VC %>% dplyr::select("Patient.ID","disease","sex","age","Nfl","pNFh","genetics","progression_rate","progression_group","slow_vital_capacity","onset","limb","age_at_onset","disease_duration","batchid","ECAS","tube_id","CSF ID","k2","k3") %>% 
    rename(ID = "Patient.ID",CSF_ID = "CSF ID")

# Sample_info from MAXOMOD project Discovery cohort
cluster_info_DC <- cluster_info_DC %>% dplyr::select("Patient ID","Tube ID","Group","Sex","Age at collection","NfL (pg/ml)","pNFh (pg/ml)","Genetic","ALS progresion rate per month (delta ALSFRS-R / days *30)","Progression group","slow vital capacity (%)","Disease onset (location of first wekness)","Limb Onset","Age at onset","Disease duration (until collection in days)","center","ECAS","Tube ID","CSF ID","k2","k3") %>% 
    rename(ID = "Patient ID",CSF_ID = "CSF ID", tube_id = "Tube ID", disease = "Group", sex = "Sex", age = "Age at collection", Nfl = "NfL (pg/ml)", pNFh = "pNFh (pg/ml)", genetics = "Genetic", progression_rate = "ALS progresion rate per month (delta ALSFRS-R / days *30)", progression_group = "Progression group", slow_vital_capacity = "slow vital capacity (%)", onset = "Disease onset (location of first wekness)", limb = "Limb Onset", age_at_onset = "Age at onset", disease_duration = "Disease duration (until collection in days)", batchid = "center", ECAS = "ECAS", k2 = "k2", k3 = "k3") %>%
    dplyr::select(colnames(cluster_info_VC))

# clean-up for na or NA, change into NA values 
cluster_info_VC <- cluster_info_VC %>% mutate(across(where(is.character), ~ ifelse(. == "na" | . == "NA", NA, .)))
cluster_info_DC <- cluster_info_DC %>% mutate(across(where(is.character), ~ ifelse(. == "na" | . == "NA", NA, .)))

cluster_info_VC$cohort = "VC"
cluster_info_DC$cohort = "DC"

cluster_info = do.call("rbind",
                       list(cluster_info_VC,
                            cluster_info_DC))

# for cluster_info$genetics, if it's "Genom neg.", or "Panel neg.", change into "negative"; if it's "not performed" or "not_performed", change into NA
cluster_info$genetics = ifelse(cluster_info$genetics == "Genom neg." | cluster_info$genetics == "Panel neg.", "negative",
                               ifelse(cluster_info$genetics == "not performed" | cluster_info$genetics == "not_performed", NA, cluster_info$genetics))

# convert to character
cluster_info$ID = as.character(cluster_info$ID)
all_participants_IDs <- all_participants_IDs %>% left_join(cluster_info, by = c("ID" = "ID"))
all_participants_IDs$Tube_ID = as.character(all_participants_IDs$Tube_ID)

# add type to all_participants_IDs (to each row, assign a unique number, group by disease)
#all_participants_IDs$type = ave(all_participants_IDs$Tube_ID, all_participants_IDs$disease, FUN = seq_along)
all_participants_IDs$type = toupper(all_participants_IDs$disease)

writexl::write_xlsx(all_participants_IDs, paste0(output_dir, "/all_participants_IDs.xlsx"))

# change SampleMatrixType OTHER into TEARS
npq_counts$SampleMatrixType = ifelse(npq_counts$SampleMatrixType == "OTHER", "TEARS", npq_counts$SampleMatrixType)
writexl::write_xlsx(npq_counts, paste0(output_dir, "/npq_counts.xlsx"))