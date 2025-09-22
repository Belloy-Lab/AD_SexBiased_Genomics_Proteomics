**Sex Strat AD PWAS**

## Tables

###

```bash
# These variables must be set prior to submitting jobs or opening an interactive session. 
export LSF_DOCKER_VOLUMES="/storage1/fs1/belloy/Active:/storage1/fs1/belloy/Active \
/storage2/fs1/belloy2/Active:/storage2/fs1/belloy2/Active /scratch1/fs1/belloy:/scratch1/fs1/belloy $HOME:$HOME"
export CONDA_ENVS_DIRS="/storage1/fs1/belloy/Active/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/belloy/Active/conda/pkgs/"
export PATH="/opt/conda/bin:$PATH"
export LSF_DOCKER_ENTRYPOINT=/bin/bash
```


```bash
bsub -g /$USER/compute-belloy -J SSgenelist -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.err \
-q subscription -R 'rusage[mem=32GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Figures_Tables/sex-specific_gene_lists.R \
        --brain_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/Brain/Sex/ \
        --csf_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/CSF/Sex/ \
        --out_dir "/storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/
```

###
```bash
bsub -g /$USER/compute-belloy -J BrainCombZ -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.err \
-q subscription -R 'rusage[mem=32GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Figures_Tables/Brain_sex-specific-gene-combZ_AFRad.R \
        --EUall_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/ \
        --AFR_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --f_EUall_stats ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.exclude_APOE_region.shared_var.hg19_1KG_EURintersected.txt \
        --m_EUall_stats ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.exclude_APOE_region.shared_var.hg19_1KG_EURintersected.txt \
        --f_AFR_stats AFRad_Females_case_control.full.hg19_AFR-ad_LDintersected.txt \
        --m_AFR_stats AFRad_Males_case_control.full.hg19_AFR-ad_LDintersected.txt \
        --EUall_f_sex_het ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_fsg_sex-het.txt \
        --EUall_m_sex_het ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_msg_sex-het.txt \
        --AFR_f_sex_het AFRad_Females_AD_cc.full.hg19_AFRad_LD_PWAS_W23_ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var-fsg_sex-het.txt \
        --AFR_m_sex_het AFRad_Males_AD_cc.full.hg19_AFRad_LD_PWAS_W23_ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var-msg_sex-het.txt
```
```bash
bsub -g /$USER/compute-belloy -J CSFcombZ -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.err \
-q subscription -R 'rusage[mem=32GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Figures_Tables/CSF_sex-specific-gene-combZ_AFRad.R \
        --EUall_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/ \
        --AFR_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/AFR \
        --fs_EUall_stats ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.no_APOE.shared_var.hg38_sex-strat_NGI_CSF_intersected.txt \
        --ms_EUall_stats ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.no_APOE.shared_var.hg38_sex-strat_NGI_CSF_intersected.txt \
        --fn_EUall_stats ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.no_APOE.shared_var.hg38_NGI_CSF_intersected.txt \
        --mn_EUall_stats ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.no_APOE.shared_var.hg38_NGI_CSF_intersected.txt \
        --fs_AFR_stats AFRad_Females_case_control.full.hg38_sex-strat_NGI_CSF_intersected.txt \
        --ms_AFR_stats AFRad_Males_case_control.full.hg38_sex-strat_NGI_CSF_intersected.txt \
        --fn_AFR_stats AFRad_Females_case_control.full.hg38_NGI_CSF_intersected.txt \
        --mn_AFR_stats AFRad_Males_case_control.full.hg38_NGI_CSF_intersected.txt \
        --EUall_f_sex_het ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_NGI-CSFc_fsg_wHP_sex-het.txt \
        --EUall_m_sex_het ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_NGI-CSFc_msg_sex-het.txt \
        --AFR_f_sex_het AFRad_Females_AD_cc.full.hg38_PWAS_NGI-CSFc_ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var-fsg_wHP_sex-het.txt \
        --AFR_m_sex_het AFRad_Males_AD_cc.full.hg38_PWAS_NGI-CSFc_ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var-msg_sex-het.txt
```

### EUR AD PWAS sex-specific hits and correspnding AD PWAS results across the EU_all analyses by sex
Sex-specific AD PWAS results for the sex-specific genes identified in the principal European sex-specific AD PWAS including the UKB cohort and corresponding findings from the European sensitivity AD PWAS. Findings were considered consistent in sensitivity analyses if the discovery sex maintained P<0.05 and the opposite sex maintained P>0.05. The most significant findings across primary and secondary discoveries are reported when applicable.

```bash
bsub -g /$USER/compute-belloy -J TS20 -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.err \
-q subscription -R 'rusage[mem=32GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Figures_Tables/TableS20.R \
        --dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/ \
        --female_top_gene ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes.txt \
        --male_top_gene ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes.txt \
        --noUKB_brain_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/Brain/Sex/ \
        --noUKB_CSF_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_noUKB/CSF/Sex/ \
        --noUKB_female_primary_brain_PWASfile ADGC_ADSP_FinnGen_Females_cc.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt \
        --noUKB_male_primary_brain_PWASfile ADGC_ADSP_FinnGen_Males_cc.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_noMHC_ext2Mb_non-strat_sex-strat_W23_weights.txt \
        --noUKB_female_primary_CSF_PWASfile ADGC_ADSP_FinnGen_Females_cc.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt \
        --noUKB_male_primary_CSF_PWASfile ADGC_ADSP_FinnGen_Males_cc.gen090.noAPOE.shared_var.hg38_PWAS_wHP_noMHC_ext2Mb_non-strat_sex-strat_CSFcis_weights.txt \
        --Known_locus_file /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/SupportingFiles/AD_Risk_Loci_consensus_2024.txt
```



### EUR AD PWAS sex-specific hits, and corresponding AFR PWAS of AD findings and cross-ancestry consistency results
Sex-specific AD PWAS results for the sex-specific genes identified in the principal European sex-specific AD PWAS, and corresponding findings in the sex-stratified AFR AD PWAS and the combined sample-size weighted results. Observations of sex heterogeneity concistency are indicated. EUR findings were considered consistent in African ancestry analyses if the P-value in the discovery sex improved after meta-analyzing the EUR and AFR findings, while the opposite sex maintained P>0.05. The most significant findings across the primary and secondary discoveries are reported when applicable.
```bash
bsub -g /$USER/compute-belloy -J TS23 -n 1 -N \
-o /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.out \
-e /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/logs/EU_all_sex-specific-gene-list_making.%J.err \
-q subscription -R 'rusage[mem=32GB] span[hosts=1]' -G compute-belloy-t1 -sla compute-belloy-t1 -a 'docker(dmr07083/fusion-project:4.3.2)' \
Rscript /storage2/fs1/belloy2/Active/04_Code/$USER/PWAS/Figures_Tables/TableS23.R \
        --dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/Figures_Tables/ \
        --female_top_gene ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes.txt \
        --male_top_gene ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes.txt \
        --PD_dir /storage2/fs1/belloy2/Active/05_Projects/$USER/PWAS/EU_all/ \
        --fb_combZ ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_fsg_sex-het_comb-AFRad_full-Z.txt \
        --mb_combZ ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_msg_sex-het_comb-AFRad_full-Z.txt \
        --fc_combZ ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_NGI-CSFc_fsg_sex-het_comb-AFRad_full-Z.txt \
        --mc_combZ ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg38_PWAS_NGI-CSFc_msg_sex-het_comb-AFRad_full-Z.txt
```
