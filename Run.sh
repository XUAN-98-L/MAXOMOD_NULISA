Rscript Script/00_initialization.R

Rscript Script/01_data_mining.R --input Results/00_Initialization/npq_counts.xlsx --output Results/01_Data_Mining --metadata Results/00_Initialization/all_participants_IDs.xlsx --target_detectability Results/00_Initialization/target_detectability.xlsx

Rscript Script/02_analysis_groups.R --input Results/00_Initialization/npq_counts.xlsx --output Results/02_Analysis_Groups --subtype FALSE --metadata Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra Results/01_Data_Mining/target_detectability_extra.xlsx --td Results/01_Data_mining/td.xlsx --top_n 20

Rscript Script/02_analysis_groups.R --input Results/00_Initialization/npq_counts.xlsx --output Results/02_Analysis_Groups_subtypes --subtype TRUE --metadata Results/00_Initialization/all_participants_IDs.xlsx --target_detectability_extra Results/01_Data_Mining/target_detectability_extra.xlsx --td Results/01_Data_mining/td.xlsx --top_n 20
