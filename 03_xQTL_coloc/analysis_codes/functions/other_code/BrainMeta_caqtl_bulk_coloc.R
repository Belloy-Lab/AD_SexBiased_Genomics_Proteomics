### Preprocessing/Analysis code for Autosomal QTL analysis (BrainMeta caQTL ABF)


print(paste("Starting preprocessing for tissue:", study, tissue, qtl_type, sep = " "))


###################################
## read in tissue qtl file
df = fread(as.character(filepath))
print(head(df))


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
  
  
  
  ###################################
  # Filter QTL file to referenced chromosome - Make ID12 and ID21
  df_2 = df %>% 
    filter(CHR == chrom) %>% 
    mutate(CHR = chrom,
           posID = paste(CHR, BP, sep = ":"),
           id12 = paste(CHR, BP, ALLELE0, ALLELE1, sep = ":"),
           id21 = paste(CHR, BP, ALLELE1, ALLELE0, sep = ":")) %>% 
    arrange(BP)
  
  
  
  ###################################  
  ## Merge QTL file and AD File by ID12 and ID 21
  merged_12 = merge(df_2, ad_df2, by.x = "id12", by.y = "posID") 
  merged_21 = merge(df_2, ad_df2, by.x = "id21", by.y = "posID") %>% 
    mutate(BETA.x = BETA.x*-1)
  
  # Combine the two datasets 
  merged_updated = rbind(merged_12, merged_21) %>% 
    arrange(CHR.x, BP.x) 
  print(paste("There are", n_distinct(merged_updated$SNP), "in first merge", sep = " "))
  
  
  # Create appropriate variable names for efficient QTL colocalization
  merged_updated <- merged_updated %>%
    dplyr::rename(CHR = CHR.x, BP = BP.x, ALLELE1 = ALLELE1.x, ALLELE0 = ALLELE0.x, AD_P = P, AD_EAF = A1FREQ, 
                  AD_BETA = BETA.y, AD_SE = SE.y, AD_N = N_incl,
                  QTL_P = p, QTL_BETA = BETA.x, QTL_SE = SE.x, QTL_N = N, Peak = Probe, PeakPos = Probe_bp) %>% 
    mutate(QTL_EAF = AD_EAF) %>% 
    dplyr::select(CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, Z, AD_P, AD_N,
                  QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, Peak, PeakPos, s)
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
  source("/storage2/fs1/belloy2/Active/04_Code/$USER/xQTL/functions/autosome_caQTL_function_abf.R")  
  
  # Run caQTL colocalization
  print(paste("Starting analysis for locus", locus_index, "in tissue:", tissue, sep = " "))
  
  result <- caQTL_colocalization(merged, all_suggestive_snps, name_dset, out_dir)
  
  
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

