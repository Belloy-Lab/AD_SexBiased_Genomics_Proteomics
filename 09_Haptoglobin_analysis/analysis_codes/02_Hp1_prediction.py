#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 8 2025

@author: talozzil
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Train final model using top 20 most frequently selected SNPs from retrain model
and predict HP1 status in all the cohorts
"""

import pandas as pd
import numpy as np
from sklearn.svm import SVC
import joblib
from collections import Counter
import os
from sklearn.metrics import r2_score

# Load selected SNPs from retrain model
snp_df = pd.read_csv("selected_snps_by_perm.csv")
snp_lists = snp_df["selected_snps"].str.split(";")
flat_snp_list = [snp for sublist in snp_lists for snp in sublist]
snp_counts = Counter(flat_snp_list)
top_20_snps = [snp for snp, _ in snp_counts.most_common(20)]


# Load genotype data
metadata_cols = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
base_path = "/mnt/barry/home/lia/HP/HP_SNPs_ADSP_ADRC"
adsp_raw = pd.read_csv(f"{base_path}/HP_imp_ADSP_WGS.raw", delim_whitespace=True).drop(columns=metadata_cols)
topmed_raw = pd.read_csv(f"{base_path}/HP_imp_TOPMed.raw", delim_whitespace=True).drop(columns=metadata_cols)
eur_raw = pd.read_csv(f"{base_path}/HP_haplotype_variant_PLINK_file_3506_EUR_samples.raw", delim_whitespace=True).drop(columns=metadata_cols)
UKB_raw = pd.read_csv("/mnt/barry/home/lia/HP/UKB_for_HP1_impute/UKB_HP1_imputation.raw", delim_whitespace=True).drop(columns=metadata_cols)
f05_raw = pd.read_csv("/mnt/barry/home/lia/HP/f05_genoEUR/f05_genoEUR_HPprotGWAS.raw", delim_whitespace=True).drop(columns=metadata_cols)

#%%

# Normalize SNP names by removing minor allele suffix
def strip_minor_allele(col):
    return col.rsplit('_', 1)[0] if '_' in col else col

# Extract positions from SNP IDs
bim_topmed = pd.read_csv(f"{base_path}/HP_imp_TOPMed.bim", delim_whitespace=True, header=None,
                         names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
bim_adsp = pd.read_csv(f"{base_path}/HP_imp_ADSP_WGS.bim", delim_whitespace=True, header=None,
                       names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
bim_eur = pd.read_csv(f"{base_path}/HP_haplotype_variant_PLINK_file_3506_EUR_samples.bim", delim_whitespace=True, header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
bim_UKB = pd.read_csv("/mnt/barry/home/lia/HP/UKB_for_HP1_impute/UKB_HP1_imputation.bim", delim_whitespace=True, header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])


bim_f05 = pd.read_csv("/mnt/barry/home/lia/HP/f05_genoEUR/f05_genoEUR_HPprotGWAS.bim", delim_whitespace=True, header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])

# Map positions to SNP IDs
pos_to_rsid_topmed = dict(zip(bim_topmed["POS"], bim_topmed["SNP"]))
pos_to_rsid_adsp = dict(zip(bim_adsp["POS"], bim_adsp["SNP"]))
pos_to_snp_eur = dict(zip(bim_eur["POS"], bim_eur["SNP"]))
pos_to_snp_UKB = dict(zip(bim_UKB["POS"], bim_UKB["SNP"]))
pos_to_snp_f05 = dict(zip(bim_f05["POS"], bim_f05["SNP"]))

# Extract position from SNP ID
def extract_position(snp_id):
    try:
        return int(snp_id.split(":")[1])
    except:
        return None

snp_pos = [extract_position(snp) for snp in top_20_snps if extract_position(snp) is not None]

snp_df = pd.DataFrame({"original_snp": top_20_snps, "position": snp_pos})
snp_df["adsp"] = snp_df["position"].map(pos_to_rsid_adsp)
snp_df["topmed"] = snp_df["position"].map(pos_to_rsid_topmed)
snp_df["eur"] = snp_df["position"].map(pos_to_snp_eur)
snp_df["UKB"] = snp_df["position"].map(pos_to_snp_UKB)
snp_df["f05"] = snp_df["position"].map(pos_to_snp_f05)

snp_df = snp_df.dropna(subset=['f05']).reset_index(drop=True)


top_20_snps_no_nan = snp_df.dropna(subset=['f05']).reset_index(drop=True)
top_20_snps_no_nan=top_20_snps_no_nan.original_snp

# Load ADRC phenotype and matched SNP matrix
LRS = pd.read_excel("SNPs_close_to_SV/LRS_subgroup.xlsx", index_col=0).T
SRS = pd.read_excel("SNPs_close_to_SV/SRS_subgroup_thr_LD01.xlsx", index_col=0)

common_samples = LRS.index.intersection(SRS.columns)
LRS = LRS.loc[common_samples]
SRS = SRS.loc[:, common_samples]

ADRC = SRS.loc[SRS.index.intersection(top_20_snps_no_nan)].T
ADRC = np.nan_to_num(ADRC, nan=np.nanmean(ADRC, axis=0))
y = np.array(LRS['HP1_del'])

# Train final model
model = SVC(kernel='poly', degree=3, C=1.0, probability=True)
model.fit(ADRC, y)
joblib.dump(model, "final_model_top20_f05_retrain_snps.pkl")
print("✔ Final model trained on ADRC using top 20 retrain SNPs.")

ADRC_predict = model.predict(ADRC)

# Compute R² on ADRC prediction
r2_adrc = r2_score(LRS['HP1_del'], ADRC_predict)
print(f"✔ R² score for ADRC prediction: {r2_adrc:.4f}")

summary = []

# ADRC  statistics (training data)
n_alleles = 2 * len(y)
n_alt = np.sum(y) * 2
freq = n_alt / n_alleles
maf = min(freq, 1 - freq)
n_pred_0 = (y == 0).sum()
n_pred_1 = (y == 1).sum()
n_pred_2 = (y == 2).sum() if (y == 2).any() else 0

summary.append({
    "Cohort": "ADRC_training",
    "MAF": maf,
    "n_samples": len(y),
    "n_pred_0": n_pred_0,
    "n_pred_1": n_pred_1,
    "n_pred_2": n_pred_2,
    "missing_individuals": 0
})

y=ADRC_predict

# ADRC  statistics (predicting data)
n_alleles = 2 * len(y)
n_alt = np.sum(y) * 2
freq = n_alt / n_alleles
maf = min(freq, 1 - freq)
n_pred_0 = (y == 0).sum()
n_pred_1 = (y == 1).sum()
n_pred_2 = (y == 2).sum() if (y == 2).any() else 0

summary.append({
    "Cohort": "ADRC_predict",
    "MAF": maf,
    "n_samples": len(y),
    "n_pred_0": n_pred_0,
    "n_pred_1": n_pred_1,
    "n_pred_2": n_pred_2,
    "missing_individuals": 0
})


# Define prediction function and summary collection

snp_maf_records = []


def compute_snp_maf_per_snp(df, label):
    for snp in snp_df["original_snp"]:
        if snp in df.columns:
            geno = df[snp].dropna().astype(float)
            freq = geno.sum() / (2 * len(geno))
            maf = min(freq, 1 - freq)
            snp_maf_records.append({"SNP": snp, "Cohort": label, "MAF": maf})


def predict_and_save(model, raw_df, fam_path, rename_dict, output_path, original_snps, label):
    fam_df = pd.read_csv(fam_path, delim_whitespace=True, header=None,
                         names=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    sample_ids = fam_df["IID"]
    raw_df.index = sample_ids

    raw_df.columns = raw_df.columns.str.replace(r"_[ACGT]+$", "", regex=True)
    df_renamed = raw_df.rename(columns=rename_dict)
    df_renamed = df_renamed.loc[:, original_snps]

    compute_snp_maf_per_snp(df_renamed, label)

    n_missing_rows = df_renamed.isna().any(axis=1).sum()
    print(f"{label}: Found {n_missing_rows} rows with NaN values before prediction.")

    df_renamed.to_csv(f"prediction/filtered_input_{label}.csv")

    df_clean = df_renamed.dropna()
    sample_ids_clean = df_clean.index
    preds = model.predict(df_clean)

    pd.DataFrame({"sample_id": sample_ids_clean, "HP1_pred": preds}).to_csv(output_path, index=False)

    n_alleles = 2 * len(preds)
    n_alt = np.sum(preds) * 2
    freq = n_alt / n_alleles
    maf = min(freq, 1 - freq)
    n_pred_0 = (preds == 0).sum()
    n_pred_1 = (preds == 1).sum()
    n_pred_2 = (preds == 2).sum() if (preds == 2).any() else 0

    summary.append({
        "Cohort": label,
        "MAF": maf,
        "n_samples": len(preds),
        "n_pred_0": n_pred_0,
        "n_pred_1": n_pred_1,
        "n_pred_2": n_pred_2,
        "missing_individuals": n_missing_rows
    })
    print(f"✔ Predictions completed for {label}. MAF: {maf:.4f}")


# Predict for all cohorts
predict_and_save(
    model,
    topmed_raw,
    os.path.join(base_path, "HP_imp_TOPMed.fam"),
    dict(zip(snp_df["topmed"], snp_df["original_snp"])),
    "prediction/HP1_predictions_TOPMed_v2.csv",
    snp_df["original_snp"],
    "TOPMed"
)

predict_and_save(
    model,
    eur_raw,
    os.path.join(base_path, "HP_haplotype_variant_PLINK_file_3506_EUR_samples.fam"),
    dict(zip(snp_df["eur"], snp_df["original_snp"])),
    "prediction/HP1_predictions_EUR_v2.csv",
    snp_df["original_snp"],
    "EUR"
)

predict_and_save(
    model,
    adsp_raw,
    os.path.join(base_path, "HP_imp_ADSP_WGS.fam"),
    dict(zip(snp_df["adsp"], snp_df["original_snp"])),
    "prediction/HP1_predictions_ADSP_v2.csv",
    snp_df["original_snp"],
    "ADSP"
)

predict_and_save(
    model,
    UKB_raw,
    "/mnt/barry/home/lia/HP/UKB_for_HP1_impute/UKB_HP1_imputation.fam",
    dict(zip(snp_df["UKB"], snp_df["original_snp"])),
    "prediction/HP1_predictions_UKB_v2.csv",
    snp_df["original_snp"],
    "UKB"
)

predict_and_save(
    model,
    f05_raw,
    "/mnt/barry/home/lia/HP/f05_genoEUR/f05_genoEUR_HPprotGWAS.fam",
    dict(zip(snp_df["f05"], snp_df["original_snp"])),
    "predictionHP1_predictions_f05_v2.csv",
    snp_df["original_snp"],
    "f05"
)


# SNP-level MAF from training ADRC
adrc_df = SRS.loc[SRS.index.intersection(top_20_snps)].T
compute_snp_maf_per_snp(adrc_df, "ADRC_train")

# Save final summaries
pd.DataFrame(summary).to_csv("prediction/HP1_prediction_summary_MAF_f05.csv", index=False)
pd.DataFrame(snp_maf_records).to_csv("prediction/SNP_MAF_summary_v2.csv", index=False)

