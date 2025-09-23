### Preprocessing/Analysis code for Autosomal QTL analysis (eQTL Catalogue eQTL ABF)

print(paste("Starting preprocessing for tissue:", study, tissue, qtl_type, sep = " "))


##########################################
## pQTL file list
source_dir = as.character(filepath)
pQTL_files = list.files(source_dir, pattern = "$", full.names = TRUE)


###################################
## read in gene list
all_suggestive_snps = fread(as.character(gene_list))
print(all_suggestive_snps)


##########################################
final_results = data.frame()
colocalization_results = data.frame()


for (i in 1:nrow(all_suggestive_snps)) {
  
  # Define variables in gene_list
  chrom = as.integer(all_suggestive_snps[i, chrom])
  bp_38_start = as.double(all_suggestive_snps[i, bp_38_start])
  bp_38_end = as.double(all_suggestive_snps[i, bp_38_end])
  bp_38_med = as.double(all_suggestive_snps[i, bp_38_med])
  locus_index = all_suggestive_snps[i, locus_index]
  stratum = all_suggestive_snps[i, stratum]
  discovery = all_suggestive_snps[i, discovery]
  range = 1e6
  
  # Print input parameters
  print(paste("Starting Data Merge of:", locus_index, "in", study, tissue, qtl_type, sep = " "))
  
  
  # Assign the corresponding data frame to ad_df based on the value in the sample column
  if (stratum == "Female") {
    print(paste("AD stratum is Female"))
    ad_df <- Female_rb
  } else if (stratum == "Male") {
    print(paste("AD stratum is Male"))
    ad_df <- Male_rb
  } else if (stratum == "SexHet") {
    print(paste("AD stratum is Sex-Het"))
    ad_df <- SexHet_rb
  } else if (stratum == "NonStrat") {
    print(paste("AD stratum is Non-Stratified"))
    ad_df <- NS_rb
  } else if (stratum == "apoe44") {
    print(paste("AD stratum is APOE4/4"))
    ad_df <- apoe44
  } else if (stratum == "apoe4_neg") {
    print(paste("AD stratum is APOE4-"))
    ad_df <- apoe4_neg
  } else if (stratum == "apoe4_pos") {
    print(paste("AD stratum is APOE4+"))
    ad_df <- apoe4_pos
  }
  
  # Filter AD df to correct chromosome
  ad_df2 = ad_df %>% 
    filter(CHR == chrom)
  
  
  # Filter ref_seq to chromosome
  ref_seq_2 = ref_seq %>% 
    filter(type == "gene", seqnames == paste0("chr", chrom)) %>% 
    group_by(gene_name) %>% 
    filter(width == max(width)) %>% 
    ungroup() %>% 
    arrange(start) %>% 
    dplyr::select(start, end, width, gene_name, gene_id, transcript_id, gene_type) %>% 
    mutate(gene_id = str_extract(gene_id, "^[^.]+"))
  
  setDT(ref_seq_2)
  
  
  # Filter ref_seq to corresponding gene of interest - get +/- 2 MB range around gene - and extract all genes
  ref_seq_3 = ref_seq_2 %>% 
    filter(gene_name == locus_index)
  
  if (nrow(ref_seq_3) == 0) {
    message1 = paste("Gene:", locus_index, "is not present in ref_seq", sep = " ")
    print(message1)
    next
  }
  
  
  gene_name <- ref_seq_3$gene_name
  lower_bound <- ref_seq_3$start-2e6
  upper_bound <- ref_seq_3$end+2e6
  
  
  # Filter ref_seq to contain only rows within +/- 2Mb range around ref_seq
  ref_seq_4 = ref_seq_2 %>% 
    filter(start >= lower_bound & end <= upper_bound)
  
  final_data = data.frame()
  
  
  # Loop through each row in ref_seq_4 --> extract gene_name, start, and end --> see if gene name is in ref_file --> extract analyte file corresponding to gene
  for (j in 1:nrow(ref_seq_4)) {
    gene_n <- ref_seq_4$gene_name[j]
    
    # Build the pattern to search for within the filepaths
    pattern <- paste0("/", gene_n, "_")
    
    # Find the filepath that contains the gene name followed by an underscore
    matched_files <- grep(pattern, pQTL_files, value = TRUE)
    
    if(length(matched_files) > 0) {
      matched_files = matched_files[1]
      # If files were matched, print them or perform other actions
      file_x <- basename(matched_files)
      print(paste("The following directories were found for gene", gene_n, ":"))
      print(matched_files)
    } else {
      # If no files were matched, print a message or perform other actions
      print(paste("No directories found for gene", gene_n))
      next
    }
    
    filepath_ref = paste0(matched_files, "/discovery_chr", chrom, "_", file_x, ".gz")
    
    # Read the file for specific analyte
    df <- fread(filepath_ref)
    df2 = df %>% 
      dplyr::rename(CHR = CHROM, BP = GENPOS) %>% 
      filter(CHR == chrom,
             BP >= lower_bound & BP <= upper_bound) %>% 
      mutate(gene_name = gene_n,
             locus = locus_index)
    
    final_data <- rbind(final_data, df2)
    print(paste("Appended:", gene_n, "to final data frame", sep = " "))
    
  } 
  
  
  message(paste("Appended dataset for:", locus_index, "complete.", sep = " "))
  
  if(nrow(final_data) == 0) {
    message2 = paste("No files for", locus_index, ":")
    print(message2)
    next
  } 
  
  
  ##########################################  
  ## Make appropriate variable for data merge with AD data 
  df_3 = final_data %>% 
    dplyr::rename(A1_FREQ = A1FREQ) %>% 
    mutate(posID = paste(CHR, BP, sep = ":"),
           id12 = paste(CHR, BP, ALLELE0, ALLELE1, sep = ":"),
           id21 = paste(CHR, BP, ALLELE1, ALLELE0, sep = ":"),
           pvalue = 10^(-LOG10P))
  
  
  ##########################################  
  ### Merge AD data with UKB data  
  merged_12 = merge(df_3, ad_df2, by.x = "id12", by.y = "posID") 
  merged_21 = merge(df_3, ad_df2, by.x = "id21", by.y = "posID") %>% 
    mutate(BETA.x = BETA.x*-1,
           A1_FREQ = 1 - A1_FREQ)  
  
  
  # Combine the two datasets 
  merged_updated = rbind(merged_12, merged_21) %>% 
    arrange(CHR.x, BP.x) 
  print(paste("There are", n_distinct(merged_updated$SNP), "in first merge", sep = " "))
  
  
  
  # Create appropriate variable names for efficient QTL colocalization    
  merged_updated <- merged_updated %>%
    filter(abs(A1FREQ - A1_FREQ) <= 0.1) %>%
    dplyr::rename(CHR = CHR.x, BP = BP.x, AD_P = P, AD_EAF = A1FREQ, AD_BETA = BETA.y, AD_SE = SE.y, AD_N = N_incl,
                  QTL_P = pvalue, QTL_EAF = A1_FREQ, QTL_BETA = BETA.x, QTL_SE = SE.x, QTL_N = N, ALLELE0 = ALLELE0.x, ALLELE1 = ALLELE1.x) %>%
    dplyr::select(CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, Z, AD_P, AD_N,
                  QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, locus, gene_name, s)
  print(paste("There are", n_distinct(merged_updated$SNP), "after filter", sep = " "))
  
  
  
  ##########################################
  # Assign variables for caQTL colocalization
  merged <- merged_updated
  name_dset <- paste(study, tissue, "chr", chrom, qtl_type, "RB", stratum, sep = "_")
  print(LC_dir)
  print(LC_threshold)
  
  # Remove unneeded objects from environment
  rm(merged_12)
  rm(merged_21)
  rm(ad_df)
  rm(ad_df2)
  rm(ref_seq_2)
  gc()
  
  print(system("free -h"))
  
  
  ## Source caQTL function
  source("/storage2/fs1/belloy2/Active/04_Code/$USER/xQTL/functions/autosome_pQTL_function_abf.R")  
  
  # Run caQTL colocalization
  print(paste("Starting analysis for locus", locus_index, "in tissue:", tissue, sep = " "))
  
  result <- pQTL_colocalization(merged = merged_updated, all_suggestive_snps, name_dset, LC_dir, LC_threshold)
  
  
  # Check to see if results came out of function
  if (is.null(result)) {
    message = paste("No results for:", locus_index, "in tissue:", tissue, qtl_type, sep = " ")
    print(message)
    next
  }
  
  
  # Accumulate results
  colocalization_results <- rbind(colocalization_results, result)
  
  # Append the current SNP details results to the final results, if non-empty
  if (nrow(colocalization_results) > 0) {
    final_results <- rbind(final_results, colocalization_results)
  }
  
  
}



## Arrange Final Results
final_results <- final_results %>% 
  arrange(desc(PP4)) %>% 
  group_by(locus, gene_name, PP3, PP4) %>% 
  distinct(PP4, .keep_all = T)
final_results = as.data.frame(final_results)


# Write out final results to output file
fwrite(final_results, paste(out_dir, study, tissue, qtl_type, "abf_hg38.csv", sep = "_"))


print(paste("Completed QTL colocalization ABF analysis for:", study, tissue, sep = " "))

