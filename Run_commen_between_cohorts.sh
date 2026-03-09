#!/usr/bin/env bash
# Run pipeline for each common-proteins genelist: CNS_DC, CNS_VC, Immune_DC, Immune_VC.
# CNS genelists -> CNS panel; Immune genelists -> Immune panel.
set -e

DATA_DIR="Data"
SCRIPT_DIR="Script"
CLUSTER_FILES="Data/clinical_data_with_cluster_DC.csv,Data/clinical_data_with_cluster_VC.csv"
METADATA="Data/NULISA_MAXOMOD_TearALS_final.xlsx"
CNS_INPUT="Data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx"
IMMUNE_INPUT="Data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx"

run_one() {
  local genelist="$1"
  local base_name panel input_xlsx
  base_name=$(basename "$genelist" .txt)
  if [[ "$base_name" == *"CNS"* ]]; then
    panel="CNS_panel"
    input_xlsx="$CNS_INPUT"
  elif [[ "$base_name" == *"Immune"* ]]; then
    panel="Immune_panel"
    input_xlsx="$IMMUNE_INPUT"
  else
    echo "Skip: cannot determine panel from $genelist (expected CNS or Immune in name)"
    return 0
  fi
  local base_out="CNS_immune/Specific_proteins/${panel}/${base_name}"
  local init_out="${base_out}/00_Initialization"
  local mining_out="${base_out}/01_Data_Mining"
  local groups_out="${base_out}/02_Analysis_Groups"
  local groups_sub_out="${base_out}/02_Analysis_Groups_subtypes"

  echo "========== $base_name ($panel) =========="
  Rscript "$SCRIPT_DIR/00_initialization_customised_targets.R" \
    --input "$input_xlsx" \
    --output "$init_out" \
    --metadata "$METADATA" \
    --cluster "$CLUSTER_FILES" \
    --cohort VC \
    --tears FALSE \
    --DEGs NULL \
    --genelist "$genelist"

  Rscript "$SCRIPT_DIR/01_data_mining.R" \
    --input "${init_out}/npq_counts.xlsx" \
    --output "$mining_out" \
    --metadata "${init_out}/all_participants_IDs.xlsx" \
    --target_detectability "${init_out}/target_detectability.xlsx"

  Rscript "$SCRIPT_DIR/02_analysis_groups.R" \
    --input "${init_out}/npq_counts.xlsx" \
    --output "$groups_out" \
    --subtype FALSE \
    --metadata "${init_out}/all_participants_IDs.xlsx" \
    --target_detectability_extra "${mining_out}/target_detectability_extra.xlsx" \
    --td "${mining_out}/td.xlsx" \
    --top_n 20

  Rscript "$SCRIPT_DIR/02_analysis_groups.R" \
    --input "${init_out}/npq_counts.xlsx" \
    --output "$groups_sub_out" \
    --subtype TRUE \
    --metadata "${init_out}/all_participants_IDs.xlsx" \
    --target_detectability_extra "${mining_out}/target_detectability_extra.xlsx" \
    --td "${mining_out}/td.xlsx" \
    --top_n 20
}

for f in \
  "${DATA_DIR}/common_proteins_CNS_DC.txt" \
  "${DATA_DIR}/common_proteins_CNS_VC.txt" \
  "${DATA_DIR}/common_proteins_Immune_DC.txt" \
  "${DATA_DIR}/common_proteins_Immune_VC.txt"; do
  if [[ -f "$f" ]]; then
    run_one "$f"
  else
    echo "Skip: file not found $f"
  fi
done

echo "Done."
