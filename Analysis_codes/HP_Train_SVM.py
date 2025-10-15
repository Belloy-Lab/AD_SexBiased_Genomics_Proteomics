#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 7 2025

@author: talozzil
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Train a model on ADRC using SNPs shared across all cohorts, selecting 20 SNPs per permutation
"""

import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif
import joblib


#%% Load BIM files
bim_topmed = pd.read_csv(
    "/Path/to/HP_imp_TOPMed.bim",
    delim_whitespace=True,
    header=None,
    names=["CHR", "SNP", "CM", "POS", "A1", "A2"]
)
bim_adsp = pd.read_csv(
    "/Path/to/HP_imp_ADSP_WGS.bim",
    delim_whitespace=True,
    header=None,
    names=["CHR", "SNP", "CM", "POS", "A1", "A2"]
)
bim_eur = pd.read_csv(
    "/Path/to/HP_haplotype_variant_PLINK_file_3506_EUR_samples.bim",
    delim_whitespace=True,
    header=None,
    names=["CHR", "SNP", "CM", "POS", "A1", "A2"]
)

bim_UKB = pd.read_csv(
    "/Path/to/UKB_for_HP1_impute/UKB_HP1_imputation.bim",
    delim_whitespace=True,
    header=None,
    names=["CHR", "SNP", "CM", "POS", "A1", "A2"]
)

#%% Map positions to SNP IDs
pos_to_rsid_topmed = dict(zip(bim_topmed["POS"], bim_topmed["SNP"]))
pos_to_rsid_adsp = dict(zip(bim_adsp["POS"], bim_adsp["SNP"]))
pos_to_snp_eur = dict(zip(bim_eur["POS"], bim_eur["SNP"]))
pos_to_snp_UKB = dict(zip(bim_UKB["POS"], bim_UKB["SNP"]))

#%% Extract position from SNP ID

def extract_position(snp_id):
    try:
        return int(snp_id.split(":")[1])
    except:
        return None

#%% Load ADRC data
LRS = pd.read_excel("/Path/to/LRS_subgroup.xlsx", index_col=0).T
SRS = pd.read_excel("/Path/to/SRS_subgroup_thr_LD01.xlsx", index_col=0)

common_samples = LRS.index.intersection(SRS.columns)
LRS = LRS.loc[common_samples]
SRS = SRS.loc[:, common_samples]

#%% Define shared SNPs across cohorts
snp_positions = [extract_position(snp) for snp in SRS.index if extract_position(snp) is not None]
snp_df = pd.DataFrame({"original_snp": SRS.index, "position": snp_positions})
snp_df["adsp"] = snp_df["position"].map(pos_to_rsid_adsp)
snp_df["topmed"] = snp_df["position"].map(pos_to_rsid_topmed)
snp_df["eur"] = snp_df["position"].map(pos_to_snp_eur)
snp_df["UKB"] = snp_df["position"].map(pos_to_snp_UKB)
snp_df = snp_df.dropna()
shared_snps = snp_df["original_snp"].tolist()

SRS = SRS.loc[shared_snps]
X_full = np.array(SRS)
X_full = np.nan_to_num(X_full, nan=np.nanmean(X_full, axis=1, keepdims=True)).T

snp_names = SRS.index.to_numpy()
y = np.array(LRS['HP1_del'])

#%% Permutation training with data-driven SNP selection
n_permutations = 100
k_best_snps = 20
results = []
snp_tracking = []

for i in range(n_permutations):
    print(f"\n--- Permutation {i+1}/{n_permutations} ---")

    X_train, X_test, y_train, y_test = train_test_split(
        X_full, y, test_size=0.3, stratify=y, random_state=i
    )

    # Step 1: Variance Threshold
    vt = VarianceThreshold(threshold=0.01)
    X_train_vt = vt.fit_transform(X_train)
    X_test_vt = vt.transform(X_test)
    vt_mask = vt.get_support()
    snps_after_vt = snp_names[vt_mask]

    # Step 2: SelectKBest
    selector = SelectKBest(score_func=f_classif, k=min(k_best_snps, X_train_vt.shape[1]))
    X_train_sel = selector.fit_transform(X_train_vt, y_train)
    X_test_sel = selector.transform(X_test_vt)

    final_snps = snps_after_vt[selector.get_support()]
    snp_tracking.append({
        "perm": i+1,
        "selected_snps": list(final_snps)
    })

    model = SVC(kernel='poly', degree=3, C=1.0)
    model.fit(X_train_sel, y_train)

    joblib.dump(model, f"train/model_perm_{i+1}.pkl")

    yhat_train = model.predict(X_train_sel)
    yhat_test = model.predict(X_test_sel)

    train_acc = accuracy_score(y_train, yhat_train)
    test_acc = accuracy_score(y_test, yhat_test)

    print(f"Train Accuracy: {train_acc:.3f} | Test Accuracy: {test_acc:.3f}")

    results.append({
        "perm": i+1,
        "train_acc": train_acc,
        "test_acc": test_acc,
        "n_snps": len(final_snps)
    })

#%% Save results
results_df = pd.DataFrame(results)
results_df.to_csv("permutation_results.csv", index=False)

snp_df = pd.DataFrame([
    {"perm": item["perm"], "selected_snps": ";".join(item["selected_snps"])}
    for item in snp_tracking
])
snp_df.to_csv("selected_snps_by_perm.csv", index=False)

print("\n✔ ADRC permutation-based model training complete.")
