**AD Sex-Biased Genomics & Proteomics**

# xQTL Colocalization Pipeline
The aim of this pipeline is to utilize a wide range of QTL datasets to run QTL colocalization analyses with a target dataset. In the example codes provided, the target dataset is a non-stratified GWAS of Alzheimer's disease.

## Requirements
* R Packages: [coloc](https://cran.r-project.org/web/packages/coloc/index.html).
 * Target dataset (GWAS summary stats for specific phenotype)
 * QTL datasets - see QTL_datasets.md for list of possible QTL datasets and tissues

 
## Workflow
* Create input reference file for the bash scipt with the following variables:
  1. study - Study from which QTL dataset comes from
  2. tissue - specific tissue/cell/brain region the QTL dataset investigates
  3. filepath - path to QTL dataset on local or remote server
  4. qtl_type - what QTL type you are investigating (eQTL, pQTL, sQTL, mQTL, haQTL, caQTL)
     
* Create gene reference file to loop through in R script
  1. chrom - chromosome number the locus of interest falls on
  2. bp_38_start - starting basepair position for larger window approach (GRch38)
  3. bp_38_end - ending basepair position for larger window approach (GRch38)
  4. bp_38_med - middle base pair position for smaller window approach (GRch38)
  5. locus_index - locus of interest for QTL colocalization analysis
  6. stratum - indicate which stratum oof the target dataset that will be used
  7. discovery - (optional) used to label output data and locus compare plots

## Datasets 
### eQTLs
* eQTL Catalogue eQTLs: 69 tissue/cell types across 10 studies (Gtex_v8, CommonMind, Braineac2, Young 2019, Kasela 2017, BLUEPRINT, CEDAR, Fairfax 2014, Quach 2016, Nedelec 2016)
  * Access at: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv
  * quant_method: 'ge'
  * Build: GRCh38
 
* Fujita 2022 eQTLS: 7 cell types
  * Access at: https://www.synapse.org/Synapse:syn52363777
  * Build: GRCh38

* Wingo 2023 eQTLs: 1 tissue (DLPFC) across 3 stratum (Non-stratifed, Female, Male)
  * Access at: https://www.synapse.org/Synapse:syn51150434/wiki/621280
  * Section: 'Brain cis-eQTLS'
  * Build: GRCh37 --> need to liftover to Build GRCh38

* MetaBrain eQTLs: 5 brain regions (Basal Ganglia, Cerebellum, Cortex, Hippocampus, Spinal Cord)
  * Access at: https://download.metabrain.nl/files.html
  * Build: GRCh38
 
* BrainMeta eQTLs: Multiple brain regions
  * Access at: [https://download.metabrain.nl/files.html](https://yanglab.westlake.edu.cn/software/smr/#DataResource)
  * Build: GRCh37


### sQTLs
* eQTL Catalogue sQTLs: 62 tissue/cell types across 6 studies (Gtex_v8, CommonMind, Braineac2, BLUEPRINT, Quach 2016, Nedelec 2016)
  * Access at: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv
  * quant_method: 'leafcutter'
  * Build: GRCh38
 

### pQTLS 
* Wingo 2023 pQTLs: 1 tissue (DLPFC) across 3 stratum (Non-stratifed, Female, Male)
  * Access at: https://www.synapse.org/Synapse:syn51150434/wiki/621280
  * Section: 'Brain cis-pQTLS'
  * Build: GRCh37 --> need to liftover to Build GRCh38

* ARIC pQTL: 1 tissue (plasma)
  * Access at: http://nilanjanchatterjeelab.org/pwas/

## Data formats and Organization
### Genes or SNPs input format
This CSV file defines the genomic regions or loci of interest to be used as input for the xQTL mapping pipeline. Each row corresponds to a specific locus identified in prior analyses (e.g., GWAS, PWAS) and includes chromosomal coordinates, the associated gene name or locus label, and stratification metadata.

Column descriptions:
- chrom: Chromosome number (e.g., 1, 2, ..., 22)
- bp_38_sta: Start position of the locus (hg38 genome build)
- bp_38_end: End position of the locus (hg38 genome build)
- bp_38_me: Midpoint or center position of the locus (hg38)
- locus_ind: Locus index or gene symbol associated with the region
- stratum: Stratification variable (e.g., Female, Male)
- discovery: Tissue or data type in which the locus was discovered (e.g., Brain, CSF)

Table : Example input gene list.

| chrom | bp_38_sta | bp_38_end | bp_38_me   | locus_ind | stratum | discovery |
|-------|-----------|-----------|------------|-----------|---------|-----------|
| 1     | 32278265  | 33313442  | 32795853.5 | YARS1     | Female  | Brain     |
| 1     | 150500305 | 151504836 | 151002570.5| MINDY1    | Female  | Brain     |
| 1     | 151277945 | 152287266 | 151782605.5| TDRKH     | Female  | Brain     |


Once the list is finalized, please move the CSV file into the 04_Code directory on the server.
Let me know when you're ready to document the next input or step!

### QTL Dataset Format and Organization
This section outlines the structure and requirements for incorporating QTL datasets into the analysis pipeline.

🧾 1. Master QTL Metadata Table (CSV format)
The pipeline uses a centralized metadata file to register and manage all QTL datasets. Each row corresponds to a single QTL dataset and includes study information, tissue context, file location, and QTL type.

| Column   | Description                                                              |
|----------|--------------------------------------------------------------------------|
| study    | Name of the dataset or study (e.g., BrainMeta, eQTLgen)                  |
| tissue   | Tissue type (e.g., DLPFC, blood, microglia, adipose)                     |
| filepath | Full path to the gzipped or plain-text QTL file                          |
| qtl_type | Type of QTL (e.g., eQTL, pQTL, sQTL, mQTL, caQTL, etc.)                  |

```csv
study,tissue,filepath,qtl_type
BrainMeta,DLPFC,/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/BrainMeta/eQTL/DLPFC/DLPFC_eQTL.txt,eQTL
eQTLgen,blood,/storage1/fs1/belloy/Active/02_Data/02_Processed/QTL/eQTLgen/eQTL/blood/blood_eQTL.csv,eQTL
```

📄 2. QTL Data File Format (gzipped or plain TSV/CSV)
Each QTL dataset file should follow a standardized tabular format to ensure compatibility with the pipeline. The expected format includes SNP-level statistics and gene mappings.

Each QTL dataset file should follow this column format:
| Column                     | Description                                                             |
|----------------------------|-------------------------------------------------------------------------|
| rsID                       | SNP rsID identifier                                                     |
| CHR                        | Chromosome number                                                       |
| BP                         | Base-pair position (hg38)                                               |
| A1, A2                     | Effect and reference alleles                                            |
| EAF                        | Effect allele frequency                                                 |
| beta, se, pvalue           | Effect size, standard error, and p-value for association                |
| molecular_trait_object_id | Identifier of the molecular trait (e.g., gene ENSG ID)                  |
| gene_id                    | Gene identifier                                                         |
| posID, id12, id21          | Alternative SNP identifiers                                             |
| N                          | Sample size                                                             |

 ```tsv
 rsID	CHR	BP	A1	A2	EAF	beta	se	pvalue	molecular_trait_object_id	gene_id	posID	id12	id21	N
rs201327123	1	14677	G	A	0.0342	-0.0782	0.3951	0.8437	ENSG00000177757	ENSG00000177757	1:14677	1:14677:G:A	1:14677:A:G	73
rs201327123	1	14677	G	A	0.0342	0.1309	0.2411	0.5892	ENSG00000187583	ENSG00000187583	1:14677	1:14677:G:A	1:14677:A:G	73
```
Both .tsv and .tsv.gz formats are supported.

✅ Integration Instructions
Format each QTL dataset to match the required structure above.
Save the dataset as a .tsv or .tsv.gz file.
Add a new entry for the dataset in the master CSV metadata table, ensuring:
Filepath is absolute
qtl_type is correctly specified (e.g., eQTL, pQTL, etc.)
Confirm that the metadata file is complete and accessible to the pipeline to ensure a successful run.

---
**Citation:** see [main repository README](../README.md) 
**License:** see [main repository README](../README.md)
