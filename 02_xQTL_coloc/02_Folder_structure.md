**Sex Strat AD xQTL**

## Project Directory Structure for Analysis
This project involves GWAS summary statistics from EU_all dataset (sex startified) and 169 xQTL (includin X eQTL, X pQTL, X sQTL, X mQTL, X haQTL and X caQTL) datasets. To streamline analysis and keep outputs organized, we define a structured folder hierarchy for each analysis scenario.

### Set Working Directory
Navigate to your project directory (choose the appropriate path based on your system). Please replace '$USER' with your actual username or preferred folder name. Also, make sure to replace '$USER' with your actual username or preferred folder name in all R and Bash scripts located in the analysis_codes directory.

```bash
cd /storage1/fs1/belloy/Active/05_Projects/$USER/
mkdir xQTL
cd xQTL

# OR

cd /storage2/fs1/belloy2/Active/05_Projects/$USER/
mkdir xQTL
cd xQTL
```

### Create Analysis Folder Structure
This directory structure organizes the results of various xQTL analysis.

Run the following commands to set up the full directory layout:
```bash
mkdir -p logs
mkdir -p ABF
mkdir -p ABF_LCplots
mkdir -p SUSIE
mkdir -p SUSIE_LCplots
```

### Set Up Code Directory
Navigate to 04_Code directory from storage1 or storage2 and set up project-specific code folders:
```bash
cd /storage1/fs1/belloy/Active/04_Code/$USER/
mkdir xQTL
cd xQTL

# OR

cd /storage2/fs1/belloy2/Active/04_Code/$USER/
mkdir xQTL
cd xQTL
```

Please download the analysis_codes folder from this repository and copy its entire contents into the 04_Code directory on your server.

If any folder paths differ from the default setup, make sure to update them accordingly in all R and Bash scripts.

### You're All Set
After running the above commands, you’ll have a clean, organized folder layout ready for:  
- Intermediate outputs  
- Final results  
- Project-specific scripts  
This structure ensures reproducibility, clarity, and ease of collaboration across datasets and analysis stages.

__END__