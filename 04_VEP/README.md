**Sex Strat AD PWAS**

## Variant Effect Predictor
This repo contains analysis methods and R scripts for performing Variant Effect Predictor (VEP) analysis. The workflow prepares VCF files for the Ensembl web VEP tool, processes the VEP output, and generates final summary tables and figures.

### Workflow

#### Preprocessing (R script)
- Takes variant data.  
- Formats input into VCF files compatible with the Ensembl web VEP tool.  

#### VEP Analysis (Ensembl Web Tool)
- Upload the VCF files generated in preprocessing.  
- Download the annotation results from the Ensembl VEP web interface.  

#### Post-processing (R script)
- Reads the VEP output (`.txt` format).  
- Cleans and summarizes annotations.  
- Produces the final tables and figures.  