# Prepare data
Rscript Script/00_initialization.R --input "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx" --output "VC_only/Results/00_Initialization" --metadata "Data/NULISA_MAXOMOD_TearALS_final.xlsx" --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort "VC"

# opt$input = "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx"
# opt$output = "VC_only/Results/00_Initialization"
# opt$metadata = "Data/NULISA_MAXOMOD_TearALS_final.xlsx"
# opt$cluster = "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv"
# opt$cohort = "VC"

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input VC_only/Results/00_Initialization/npq_counts.xlsx --output VC_only/Results/01_Data_Mining --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --target_detectability VC_only/Results/00_Initialization/target_detectability.xlsx

# opt$input = "VC_only/Results/00_Initialization/npq_counts.xlsx"
# opt$output = "VC_only/Results/01_Data_Mining"
# opt$metadata = "VC_only/Results/00_Initialization/all_participants_IDs.xlsx"
# opt$target_detectability = "VC_only/Results/00_Initialization/target_detectability.xlsx"

# Analyze data, calculate ANOVA & pairwise tests for each fluid
Rscript Script/02_analysis_groups.R --input VC_only/Results/00_Initialization/npq_counts.xlsx --output VC_only/Results/02_Analysis_Groups --subtype FALSE --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra VC_only/Results/01_Data_Mining/target_detectability_extra.xlsx --td VC_only/Results/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input VC_only/Results/00_Initialization/npq_counts.xlsx --output VC_only/Results/02_Analysis_Groups_subtypes --subtype TRUE --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra VC_only/Results/01_Data_Mining/target_detectability_extra.xlsx --td VC_only/Results/01_Data_Mining/td.xlsx --top_n 20

###### RUN PCA VISUALIZATION ######
# rm -rf VC_only/Results/03_PCA_v*
Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "plateid"

# adjust for age and sex
Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_adjusted --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "type"

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_adjusted_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "subtype"

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_adjusted_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_adjusted_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input VC_only/Results/01_Data_Mining/protein_data_IDs.xlsx --output VC_only/Results/03_PCA_vis_adjusted_subtypes --metadata VC_only/Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts VC_only/Results/00_Initialization/npq_counts.xlsx --group_by "plateid"

#############################only keep the targets that are DEGs in the VC cohort#############################
# Prepare data
Rscript Script/00_initialization_onlyDGEsInVC.R --input "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx" --output "VC_only_onlyDGEsInVC/Results/00_Initialization" --metadata "Data/NULISA_MAXOMOD_TearALS_final.xlsx" --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort "VC"

# opt$input = "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx"
# opt$output = "VC_only_onlyDGEsInVC/Results/00_Initialization"
# opt$metadata = "Data/NULISA_MAXOMOD_TearALS_final.xlsx"
# opt$cluster = "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv"
# opt$cohort = "VC"

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input VC_only_onlyDGEsInVC/Results/00_Initialization/npq_counts.xlsx --output VC_only_onlyDGEsInVC/Results/01_Data_Mining --metadata VC_only_onlyDGEsInVC/Results/00_Initialization/all_participants_IDs.xlsx --target_detectability VC_only_onlyDGEsInVC/Results/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input VC_only_onlyDGEsInVC/Results/00_Initialization/npq_counts.xlsx --output VC_only_onlyDGEsInVC/Results/02_Analysis_Groups --subtype FALSE --metadata VC_only_onlyDGEsInVC/Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra VC_only_onlyDGEsInVC/Results/01_Data_Mining/target_detectability_extra.xlsx --td VC_only_onlyDGEsInVC/Results/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input VC_only_onlyDGEsInVC/Results/00_Initialization/npq_counts.xlsx --output VC_only_onlyDGEsInVC/Results/02_Analysis_Groups_subtypes --subtype TRUE --metadata VC_only_onlyDGEsInVC/Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra VC_only_onlyDGEsInVC/Results/01_Data_Mining/target_detectability_extra.xlsx --td VC_only_onlyDGEsInVC/Results/01_Data_Mining/td.xlsx --top_n 20