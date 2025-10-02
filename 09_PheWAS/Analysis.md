**AD Sex-Biased Genomics & Proteomics**

### Scripts
- `01_Train_SVM.py` — Train SVM classifier.
- `02_Hp1_prediction.py` — Predict HP1.
- `03_CSF_phewas_code.R` — Run CSF protein PheWAS.
- `04_CSF_phewas_append_results.R` — Collate/append CSF per-protein results into a single table.
- `03_plasma_phewas_code.R` — Run plasma protein PheWAS.
- `04_plasma_phewas_append_results.R` — Collate/append plasma per‑protein results into a single table.


### How to run
Edit path variables at the top of each script, then:

```bash
python analysis_codes/01_Train_SVM.py

python analysis_codes/02_Hp1_prediction.py

Rscript analysis_codes/03_CSF_phewas_code.R
Rscript analysis_codes/04_CSF_phewas_append_results.R

Rscript analysis_codes/03_plasma_phewas_code.R
Rscript analysis_codes/04_plasma_phewas_append_results.R
```

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
