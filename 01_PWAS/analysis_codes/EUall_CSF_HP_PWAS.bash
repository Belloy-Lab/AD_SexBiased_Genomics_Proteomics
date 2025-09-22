#!/bin/bash

# Exit the script on any command failure
set -e

# Function to run an Rscript command and track failures
run_rscript() {
    local job_label="$1"
    local r_command="$2"

    echo "============================================================"
    echo "Starting job: $job_label"
    echo "Start time: $(date)"
    echo "Command:"
    echo "$r_command"
    echo "------------------------------------------------------------"

    # Run the R script
    if eval "$r_command"; then
        echo "Job $job_label completed successfully at $(date)"
    else
        echo "ERROR: Job $job_label failed at $(date)" >&2
        exit 1
    fi

    echo "------------------------------------------------------------"
    echo "Finished job: $job_label"
    echo "============================================================"
    echo ""
}

# ----------------------
# Begin job submissions
# ----------------------

run_rscript "Job 1" "Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.no_APOE.shared_var.hg38_sex-strat_NGI_CSF_intersected.sumstats \
    --weights /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt/NGI_CSF_female_cis_HP-weights.pos \
    --weights_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
    --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
    --chr 16 \
    --out /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/Sex/Female/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38.sex-strat_NGI_CSF_weights.HP.dat"

run_rscript "Job 2" "Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.no_APOE.shared_var.hg38_sex-strat_NGI_CSF_intersected.sumstats \
    --weights /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt/NGI_CSF_male_cis_HP-weights.pos \
    --weights_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
    --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
    --chr 16 \
    --out /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/Sex/Male/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38.sex-strat_NGI_CSF_weights.HP.dat"

run_rscript "Job 3" "Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.no_APOE.shared_var.hg38_sex-strat_NGI_CSF_intersected.sumstats \
    --weights /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt/NGI_CSF_male_cis_HP-weights.pos  \
    --weights_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
    --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
    --chr 16 \
    --out /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/Sex/FemaleGWAS_MalePW/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38.sex-strat_NGI_CSF_weights-Males.HP.dat"

run_rscript "Job 4" "Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.no_APOE.shared_var.hg38_sex-strat_NGI_CSF_intersected.sumstats \
    --weights /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt/NGI_CSF_female_cis_HP-weights.pos  \
    --weights_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/pwas_cis_wgt \
    --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files/LD/LD/CSF_proteomics_3107_samples_maf_0032_geno_1_chr_ \
    --chr 16 \
    --out /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/Sex/MaleGWAS_FemalePW/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38.sex-strat_NGI_CSF_weights-Females.HP.dat"

run_rscript "Job 5" "Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.no_APOE.shared_var.hg38_NGI_CSF_intersected.sumstats \
    --weights /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_Weights/NGI_CSF_cis_HP-weights.pos \
    --weights_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_Weights \
    --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_LD_Reference_Files/CSF_proteomics_LONGS_PPMI_Stanford_GARFIELD_samples_geno_0.1_mac_10_3506_samples_chr \
    --chr 16 \
    --out /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/NonSex/Female/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38.NGI_CSF_weights.HP.dat"

run_rscript "Job 6" "Rscript /storage1/fs1/belloy/Active/01_References/02_Commands_and_Tutorials/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.no_APOE.shared_var.hg38_NGI_CSF_intersected.sumstats \
    --weights /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_Weights/NGI_CSF_cis_HP-weights.pos \
    --weights_dir /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_Weights \
    --ref_ld_chr /storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_LD_Reference_Files/CSF_proteomics_LONGS_PPMI_Stanford_GARFIELD_samples_geno_0.1_mac_10_3506_samples_chr \
    --chr 16 \
    --out /storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/CSF/NonSex/Male/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38.NGI_CSF_weights.HP.dat"

echo "All jobs completed successfully!"