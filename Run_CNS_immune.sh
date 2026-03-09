#################CNSDisease panel#################
# Prepare data
Rscript Script/00_initialization.R --input "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx" --output "CNS_immune/Results/CNS_panel/With_tears/00_Initialization" --metadata "Data/NULISA_MAXOMOD_TearALS_final.xlsx" --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort "VC" --tears TRUE

# opt$input = "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx"
# opt$output = "CNS_immune/Results/CNS_panel/With_tears/00_Initialization_customised_targets"
# opt$metadata = "Data/NULISA_MAXOMOD_TearALS_final.xlsx"
# opt$cluster = "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv"
# opt$cohort = "VC"
# opt$tears = "TRUE"

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Results/CNS_panel/With_tears/00_Initialization/target_detectability.xlsx

# Analyze data, calculate ANOVA & pairwise tests for each fluid
Rscript Script/02_analysis_groups.R --input CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/CNS_panel/With_tears/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/CNS_panel/With_tears/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/td.xlsx --top_n 20

###### RUN PCA VISUALIZATION ######
# rm -rf CNS_immune/Results/CNS_panel/With_tears/03_PCA_v*
Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"

# adjust for age and sex
Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_adjusted --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "type"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "subtype"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/With_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/With_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/With_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/With_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"



#################CNSDisease panel without tears#################
# Prepare data
Rscript Script/00_initialization.R --input "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx" --output "CNS_immune/Results/CNS_panel/Without_tears/00_Initialization" --metadata "Data/NULISA_MAXOMOD_TearALS_final.xlsx" --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort "VC" --tears FALSE

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/target_detectability.xlsx

# Analyze data, calculate ANOVA & pairwise tests for each fluid
Rscript Script/02_analysis_groups.R --input CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/td.xlsx --top_n 20

###### RUN PCA VISUALIZATION ######
# rm -rf CNS_immune/Results/CNS_panel/Without_tears/03_PCA_v*
Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_adjusted --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "type"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "subtype"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/CNS_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/CNS_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/CNS_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"

#################Immune panel#################
# doesn't have tears
Rscript Script/00_initialization.R --input "Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx" --output "CNS_immune/Results/Immune_panel/Without_tears/00_Initialization" --metadata "Data/NULISA_MAXOMOD_TearALS_final.xlsx" --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort "VC"

Rscript Script/01_data_mining.R --input CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/target_detectability.xlsx

# opt$input = "VC_only/Results/immune/00_Initialization/npq_counts.xlsx"
# opt$output = "VC_only/Results/immune/01_Data_Mining"
# opt$metadata = "VC_only/Results/immune/00_Initialization/all_participants_IDs.xlsx"
# opt$target_detectability = "VC_only/Results/immune/00_Initialization/target_detectability.xlsx"

# Analyze data, calculate ANOVA & pairwise tests for each fluid
Rscript Script/02_analysis_groups.R --input CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/td.xlsx --top_n 20

###### RUN PCA VISUALIZATION ######
# rm -rf CNS_immune/Results/Immune_panel/Without_tears/03_PCA_v*
Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx

# PCA part does not work now due to the target “IFNA2” does not appear in Appendix 2: NULISAseq™ Inflammation Panel 250. However, it is present in the file “P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx” for plasma, but not for serum or CSF.

# opt$input = "CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx"
# opt$output = "CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis"
# opt$metadata = "CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx"
# opt$subtype = FALSE
# opt$seed = 123
# opt$label = FALSE
# opt$fluids = "CSF,SERUM,PLASMA"
# opt$adjust = FALSE
# opt$npq_counts = "CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx"
# opt$group_by = "type"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"

# # adjust for age and sex
Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "type"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "subtype"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"


#################################### 00_initialization_customised_targets.R ###################################
# DEGs from DC or VC
# cp /Users/xliu2942/Documents/Projects/MAXOMOD/MAXOMOD_CSF/Validation/03_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv VC_Differential_expression_analysis_for_k2.csv
# cp /Users/xliu2942/Documents/Projects/MAXOMOD/MAXOMOD_CSF/Discovery/03_Differential_expression_analysis_subclusters/Differential_expression_analysis_for_k2.csv DC_Differential_expression_analysis_for_k2.csv
# # alpha or beta DEGs from DC or VC
# cp /Users/xliu2942/Documents/Projects/MAXOMOD/MAXOMOD_CSF/Discovery/04_Vis_Differential_expression_analysis_subclusters/all_important_protein_names_k2_alpha.csv DC_alpha_Differential_expression_analysis_for_k2.csv
# cp /Users/xliu2942/Documents/Projects/MAXOMOD/MAXOMOD_CSF/Discovery/04_Vis_Differential_expression_analysis_subclusters/all_important_protein_names_k2_beta.csv DC_beta_Differential_expression_analysis_for_k2.csv
# cp /Users/xliu2942/Documents/Projects/MAXOMOD/MAXOMOD_CSF/Validation/04_Vis_Differential_expression_analysis_subclusters/all_important_protein_names_k2_alpha.csv VC_alpha_Differential_expression_analysis_for_k2.csv
# cp /Users/xliu2942/Documents/Projects/MAXOMOD/MAXOMOD_CSF/Validation/04_Vis_Differential_expression_analysis_subclusters/all_important_protein_names_k2_beta.csv VC_beta_Differential_expression_analysis_for_k2.csv

# CNS panel, for VC DEGs
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx --output CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs VC_only

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/01_Data_Mining --metadata CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/VC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20


# CNS panel, for DC DEGs
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx --output CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs DC_only

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/01_Data_Mining --metadata CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/DC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20


# CNS panel, for common DEGs in both cohort
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx --output CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs both

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/both_DEGs/01_Data_Mining --metadata CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/both_DEGs/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/both_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/both_DEGs/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/both_DEGs/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/CNS_panel/both_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/both_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/both_DEGs/01_Data_Mining/td.xlsx --top_n 20




# Immune panel, for VC DEGs
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx --output CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs VC_only

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/01_Data_Mining --metadata CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/VC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20


# Immune panel, for DC DEGs
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx --output CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs DC_only

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/01_Data_Mining --metadata CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/DC_only_DEGs/01_Data_Mining/td.xlsx --top_n 20


# Immune panel, for common DEGs in both cohort
# only two proteins are detected in both cohorts
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx --output CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs both

# # Analyze data, calculated correlation between fluids & sample detectability
# Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/both_DEGs/01_Data_Mining --metadata CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/target_detectability.xlsx

# Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/both_DEGs/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/both_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/both_DEGs/01_Data_Mining/td.xlsx --top_n 20

# Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/both_DEGs/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/Immune_panel/both_DEGs/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/both_DEGs/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/both_DEGs/01_Data_Mining/td.xlsx --top_n 20



#########WGCNA module#########
# CNS panel, for turquoise module
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx --output CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs NULL --genelist Data/turquoise.txt

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/turquoise/01_Data_Mining --metadata CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/turquoise/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/turquoise/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/turquoise/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/turquoise/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/CNS_panel/turquoise/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/turquoise/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/turquoise/01_Data_Mining/td.xlsx --top_n 20

# CNS panel, for blue module
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx --output CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs NULL --genelist Data/blue.txt

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/blue/01_Data_Mining --metadata CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/blue/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/blue/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/blue/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/CNS_panel/blue/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/CNS_panel/blue/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/CNS_panel/blue/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/CNS_panel/blue/01_Data_Mining/td.xlsx --top_n 20

# Immune panel, for turquoise module
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx --output CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs NULL --genelist Data/turquoise.txt

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/turquoise/01_Data_Mining --metadata CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/turquoise/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/turquoise/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/turquoise/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/turquoise/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/Immune_panel/turquoise/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/turquoise/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/turquoise/01_Data_Mining/td.xlsx --top_n 20

# CNS panel, for blue module
Rscript Script/00_initialization_customised_targets.R --input Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx --output CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization --metadata Data/NULISA_MAXOMOD_TearALS_final.xlsx --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort VC --tears FALSE --DEGs NULL --genelist Data/blue.txt

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/blue/01_Data_Mining --metadata CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/all_participants_IDs.xlsx --target_detectability CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/blue/02_Analysis_Groups --subtype FALSE --metadata CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/blue/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/blue/01_Data_Mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/npq_counts.xlsx --output CNS_immune/Specific_proteins/Immune_panel/blue/02_Analysis_Groups_subtypes --subtype TRUE --metadata CNS_immune/Specific_proteins/Immune_panel/blue/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra CNS_immune/Specific_proteins/Immune_panel/blue/01_Data_Mining/target_detectability_extra.xlsx --td CNS_immune/Specific_proteins/Immune_panel/blue/01_Data_Mining/td.xlsx --top_n 20

#################################### Find which proteins are commen with MS-based analysis ###################################
# library(readxl)
# CNS = read_excel("/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/CNS_immune/Results/CNS_panel/With_tears/00_Initialization/target_detectability.xlsx")
# Immune = read_excel("/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/target_detectability.xlsx")

# DC = read.csv("/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/Data/DC_Differential_expression_analysis_for_k2.csv")
# VC = read.csv("/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/Data/VC_Differential_expression_analysis_for_k2.csv")

# # Find the common proteins between CNS and Immune
# output_dir = "/Users/xliu2942/Documents/Projects/NULISA_MAXOMOD/Data"
# common_proteins_CNS_DC = intersect(CNS$Target, DC$name)
# common_proteins_CNS_VC = intersect(CNS$Target, VC$name)
# common_proteins_Immune_DC = intersect(Immune$Target, DC$name)
# common_proteins_Immune_VC = intersect(Immune$Target, VC$name)
# #column name is "GeneID"
# library(data.table)
# write.table(common_proteins_CNS_DC, file.path(output_dir, "common_proteins_CNS_DC.txt"), quote = FALSE, row.names = FALSE, col.names = "GeneID")
# write.table(common_proteins_CNS_VC, file.path(output_dir, "common_proteins_CNS_VC.txt"), quote = FALSE, row.names = FALSE, col.names = "GeneID")
# write.table(common_proteins_Immune_DC, file.path(output_dir, "common_proteins_Immune_DC.txt"), quote = FALSE, row.names = FALSE, col.names = "GeneID")
# write.table(common_proteins_Immune_VC, file.path(output_dir, "common_proteins_Immune_VC.txt"), quote = FALSE, row.names = FALSE, col.names = "GeneID")

# # Print the common proteins
# print(paste("Common proteins between CNS and DC:", length(common_proteins_CNS_DC)))
# print(paste("Common proteins between CNS and VC:", length(common_proteins_CNS_VC)))
# print(paste("Common proteins between Immune and DC:", length(common_proteins_Immune_DC)))
# print(paste("Common proteins between Immune and VC:", length(common_proteins_Immune_VC)))


bash Script/Run_commen_between_cohorts.sh