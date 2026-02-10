# Prepare data
Rscript Script/00_initialization.R

# Analyze data, calculated correlation between fluids & sample detectability
Rscript Script/01_data_mining.R --input Results/00_Initialization/npq_counts.xlsx --output Results/01_Data_Mining --metadata Results/00_Initialization/all_participants_IDs.xlsx --target_detectability Results/00_Initialization/target_detectability.xlsx

# Analyze data, calculate ANOVA & pairwise tests for each fluid
Rscript Script/02_analysis_groups.R --input Results/00_Initialization/npq_counts.xlsx --output Results/02_Analysis_Groups --subtype FALSE --metadata Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra Results/01_Data_Mining/target_detectability_extra.xlsx --td Results/01_Data_mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input Results/00_Initialization/npq_counts.xlsx --output Results/02_Analysis_Groups_subtypes --subtype TRUE --metadata Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra Results/01_Data_Mining/target_detectability_extra.xlsx --td Results/01_Data_mining/td.xlsx --top_n 20

###### RUN PCA VISUALIZATION ######
# rm -rf Results/03_PCA_v*
Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts Results/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts Results/00_Initialization/npq_counts.xlsx

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust FALSE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "plateid"

# adjust for age and sex
Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_adjusted --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype FALSE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "type"

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_adjusted_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "subtype"

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_adjusted_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "sex"

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_adjusted_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "cohort"

Rscript Script/03_PCA_vis.R --input Results/01_Data_Mining/protein_data_IDs.xlsx --output Results/03_PCA_vis_adjusted_subtypes --metadata Results/00_Initialization/all_participants_IDs.xlsx --subtype TRUE --seed 123 --label FALSE --fluids "CSF,SERUM,PLASMA,TEARS" --adjust TRUE --npq_counts Results/00_Initialization/npq_counts.xlsx --group_by "plateid"
