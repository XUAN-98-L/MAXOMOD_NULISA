#################CNSDisease panel#################
# Prepare data
Rscript Script/00_initialization.R --input "Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx" --output "CNS_immune/Results/CNS_panel/With_tears/00_Initialization" --metadata "Data/NULISA_MAXOMOD_TearALS_final.xlsx" --cluster "Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv" --cohort "VC" --tears TRUE

# opt$input = "Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx"
# opt$output = "VC_only/Results/00_Initialization_immune"
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

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "sex"

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust FALSE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"

# # adjust for age and sex
# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "type"

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "subtype"

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "cohort"

# Rscript Script/03_PCA_vis.R --input CNS_immune/Results/Immune_panel/Without_tears/01_Data_Mining/protein_data_IDs.xlsx --output CNS_immune/Results/Immune_panel/Without_tears/03_PCA_vis_adjusted_subtypes --metadata CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA" --adjust TRUE --npq_counts CNS_immune/Results/Immune_panel/Without_tears/00_Initialization/npq_counts.xlsx --group_by "plateid"


