
###################################
## Read in QTL data with consistent format
print("Reading in QTL dataset")
df = fread(as.character(filepath))


# Read in gene list
all_suggestive_snps = fread(as.character(gene_list)) 
print(all_suggestive_snps)


##########################################
# Set output results to empty dataframe
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
  print(paste("Filtering AD GWAS to chromosome:", chrom, sep = " "))
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
  
  
  
  ###################################
  # Filter QTL file to referenced chromosome
  df_2 = df %>% 
    filter(CHR == chrom) %>% 
    arrange(BP)
  
  
  ###################################
  # Merge eQTL data with ref_seq data
  df_3 = merge(df_2, ref_seq_2, all.x = T, by = "gene_id") %>% 
    arrange(BP)
  
  if (study == "MetaBrain" | study == "eQTLgen") {
    
    df_3 = df_3 %>% 
      mutate(gene_name = GeneSymbol)
    
  }
  
  
  ######################
  print("Merging AD GWAS with QTL dataset")
  if (study == "kosoy" | study == "eQTLgen" | study == "SingleBrain"){
    
    ###################################
    ## Merge QTL file and AD File by ID12 and ID 21
    merged_12 = merge(df_3, ad_df2, by.x = "id12", by.y = "posID") 
    merged_21 = merge(df_3, ad_df2, by.x = "id21", by.y = "posID") %>% 
      mutate(beta = beta*-1)
    
    
    # Combine the two datasets 
    merged_updated = rbind(merged_12, merged_21) %>% 
      arrange(CHR.x, BP.x) 
    print(paste("There are", n_distinct(merged_updated$SNP), "in first merge", sep = " "))
    
    
    # Create appropriate variable names for efficient QTL colocalization
    merged_updated <- merged_updated %>%
      dplyr::rename(CHR = CHR.x, BP = BP.x, AD_P = P, AD_EAF = A1FREQ, AD_BETA = BETA, AD_SE = SE, AD_N = N_incl,
                    QTL_P = pvalue, QTL_BETA = beta, QTL_SE = se, QTL_N = N) %>%
      mutate(QTL_EAF = AD_EAF) %>% 
      dplyr::select(CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, Z, AD_P, AD_N,
                    QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, gene_name, molecular_trait_object_id, s)
    print(paste("There are", n_distinct(merged_updated$SNP), "after filter", sep = " "))   
    
  } else {
    
    ###################################
    ## Merge QTL file and AD File by ID12 and ID 21
    merged_12 = merge(df_3, ad_df2, by.x = "id12", by.y = "posID") 
    merged_21 = merge(df_3, ad_df2, by.x = "id21", by.y = "posID") %>% 
      mutate(beta = beta*-1,
             EAF = 1 - EAF)
    
    
    # Combine the two datasets 
    merged_updated = rbind(merged_12, merged_21) %>% 
      arrange(CHR.x, BP.x) 
    print(paste("There are", n_distinct(merged_updated$SNP), "in first merge", sep = " "))
    
    
    # Create appropriate variable names for efficient QTL colocalization
    merged_updated <- merged_updated %>%
      filter(abs(A1FREQ - EAF) <= 0.1) %>%
      dplyr::rename(CHR = CHR.x, BP = BP.x, AD_P = P, AD_EAF = A1FREQ, AD_BETA = BETA, AD_SE = SE, AD_N = N_incl,
                    QTL_P = pvalue, QTL_EAF = EAF, QTL_BETA = beta, QTL_SE = se, QTL_N = N) %>%
      dplyr::select(CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, Z, AD_P, AD_N,
                    QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, gene_name, molecular_trait_object_id, s)
    print(paste("There are", n_distinct(merged_updated$SNP), "after filter", sep = " "))    
    
    
  }
  
  merged = merged_updated
  
  ##########################################
  # Assign variables for eQTL colocalization
  name_dset <- paste(study, tissue, "chr", chrom, qtl_type, "RB", stratum, sep = "_") # Naming convention for QTL results output table
  print(LC_dir) # Locus Compare output directory
  print(LC_threshold) # PP4 threshold for LC creation
  
  # Remove unneeded objects from environment
  rm(merged_12)
  rm(merged_21)
  rm(ad_df)
  rm(ad_df2)
  rm(df_3)
  rm(ref_seq_2)
  gc()
  
  print(system("free -h"))
  
  
  ### Detrmine which colocalization to run:
  
  if (qtl_type == "eQTL") {
    ## Source eQTL SuSiE function - update function with appropriate input paramters
    source("/storage2/fs1/belloy2/Active/04_Code/sivas/xQTL/functions/autosome_eQTL_function_susie.R") 
    # Run eQTL colocalization
    print(paste("Starting analysis for locus", locus_index, "in tissue:", tissue, sep = " "))
    result <- eQTL_colocalization(merged = merged_updated, all_suggestive_snps, name_dset, LC_dir, LC_threshold)
    
  } else if (qtl_type == "sQTL") {
    ## Source sQTL SuSiE function - update function with appropriate input paramters
    source("/storage2/fs1/belloy2/Active/04_Code/sivas/xQTL/functions/autosome_sQTL_function_susie.R") 
    # Run sQTL colocalization
    print(paste("Starting analysis for locus", locus_index, "in tissue:", tissue, sep = " "))
    result <- sQTL_colocalization(merged = merged_updated, all_suggestive_snps, name_dset, LC_dir, LC_threshold)
    
  } else if (qtl_type == "pQTL") {
    ## Source pQTL SuSiE function - update function with appropriate input paramters
    source("/storage2/fs1/belloy2/Active/04_Code/sivas/xQTL/functions/autosome_pQTL_function_susie.R") 
    # Run pQTL colocalization
    print(paste("Starting analysis for locus", locus_index, "in tissue:", tissue, sep = " "))
    result <- pQTL_colocalization(merged = merged_updated, all_suggestive_snps, name_dset, LC_dir, LC_threshold)
    
  }
  
  
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
  
  
  # Write out final results to output file
  fwrite(colocalization_results, paste(out_dir, study, tissue, qtl_type, i, "susie_hg38.csv", sep = "_"))
  
  
  if (i > 1) {
    # Delete old result file
    old_i = i-1
    file.remove(paste(out_dir, study, tissue, qtl_type, old_i, "susie_hg38.csv", sep = "_"))
    
  }
  
  
}


## Arrange Final Results
final_results <- final_results %>% 
  arrange(desc(PP4)) %>% 
  group_by(locus, gene_name, PP3, PP4) %>% 
  distinct(PP4, .keep_all = T)
final_results = as.data.frame(final_results)


# Write out final results to output file
fwrite(final_results, paste(out_dir, study, tissue, qtl_type, "susie_hg38.csv", sep = "_"))

file.remove(paste(out_dir, study, tissue, qtl_type, i, "susie_hg38.csv", sep = "_"))

# Completed work flow
print(paste("Completed QTL colocalization SuSiE analysis for:", study, tissue, sep = " "))


