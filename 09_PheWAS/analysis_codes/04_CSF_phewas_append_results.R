## PheWAS Results append for CSF

library(data.table)
library(tidyverse)
#library(ggpubr)



############################################
# Function to extract the analyte label
extract_analyte_label <- function(file_path) {
  # Use regular expression to match the pattern and extract the desired part
  pattern <- "region\\.(.*?)\\.glm"
  match <- regmatches(file_path, regexpr(pattern, file_path))
  # Extract the specific part
  analyte_label <- sub("region\\.", "", sub("\\.glm", "", match))
  return(analyte_label)
}

#################################################################################
# Append non-conditioned results
#nc_dir = "/Volumes/belloy/Active/05_Projects/noahc/opt_QTL/post_hoc/HP/phewas/output_data/non_conditioned/"
nc_dir = "Output_data/new_non_conditioned/"
nc_files <- list.files(nc_dir, pattern = ".linear$", full.names = TRUE)
print(length(nc_files))
#################################################################################

nc_results = data.frame()

for (i in nc_files) {
  
  # Index file
  file = i
  print(file)
  
  # Extract analyte label
  analyte_label <- extract_analyte_label(file)
  print(paste0("Reading in file: ", analyte_label))
  
  # Read in indexed file - create analyte variable and P-value
  file1 = fread(file, header = TRUE) 
  n_variants = n_distinct(file1$ID)
  
  file1 = file1 %>% 
    mutate(analyte = analyte_label,
           P = 10^(-LOG10_P),
           num_variant= n_variants)
  
  # Append file results to final results
  nc_results = rbind(nc_results, file1)
  
}

print(paste0("Completed data append for non-conditioned directory"))



# # Get top analyte for every variant
# top_nc = nc_results %>%
#   group_by(ID) %>%
#   filter(P == min(P)) %>%
#   ungroup()
# 
# rm(nc_results)


#################################################################
# Read in HP1 Results
hp = fread("/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/post_hoc/HP/phewas/HP1_genotype_phewas_CSF_pQTL.txt")

# Create reference df for gene name
ref = hp %>% 
  mutate(analyte = SomaID) %>% 
  dplyr::select(analyte, EntrezGeneSymbol) 


# Merge gene name reference to top_nc results
top_nc2 = merge(nc_results, ref, by = "analyte") %>%
  arrange(POS)


# Write out file
fwrite(top_nc2, "/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/post_hoc/HP/phewas/output_data/HP_CSF_pheWAS_non-conditioned_appended.txt", sep = "\t")
print(paste0("Wrote out data for non-conditioned"))

rm(nc_results)
rm(top_nc2)

##################################################################################################################################################################
##################################################################################################################################################################
# Function to extract the analyte label
extract_analyte_label2 <- function(file_path) {
  # Use regular expression to match the pattern and extract the desired part
  pattern <- "region_cond_HP1\\.(.*?)\\.glm"
  match <- regmatches(file_path, regexpr(pattern, file_path))
  # Extract the specific part
  analyte_label <- sub("region_cond_HP1\\.", "", sub("\\.glm", "", match))
  return(analyte_label)
}

#################################################################################
# Append non-conditioned results
# con_dir = "/Volumes/belloy/Active/05_Projects/noahc/opt_QTL/post_hoc/HP/phewas/output_data/conditioned/"
con_dir = "/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/post_hoc/HP/phewas/output_data/new_conditioned/"
con_files <- list.files(con_dir, pattern = ".linear$", full.names = TRUE)
print(length(con_files))

# Set empty df to start
con_results = data.frame()

for (i in con_files) {
  
  # Index file
  file = i
  
  # Extract analyte label
  analyte_label <- extract_analyte_label2(file)
  print(paste0("Reading in file: ", analyte_label))
  
  # Read in indexed file - create analyte variable and P-value
  file1 = fread(file, header = TRUE) 
  n_variants = n_distinct(file1$ID)
  
  file1 = file1 %>% 
    mutate(analyte = analyte_label,
           P = 10^(-LOG10_P),
           num_variant = n_variants)
  
  con_results = rbind(con_results, file1)
  
}

print(paste0("Completed data append for conditioned"))


# Get top analyte for each variant
# top_con = con_results %>% 
#   group_by(ID) %>%
#   filter(P == min(P)) %>% 
#   ungroup()
# 
# rm(con_results)



# Merge gene name reference to top_nc results
top_con2 = merge(con_results, ref, by = "analyte") %>% 
  rename(Gene = EntrezGeneSymbol) %>%
  arrange(POS)


# Write out file
fwrite(top_con2, "/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/post_hoc/HP/phewas/output_data/HP_CSF_pheWAS_conditioned_appended.txt", sep = "\t")
print(paste0("Wrote out data for conditioned"))

