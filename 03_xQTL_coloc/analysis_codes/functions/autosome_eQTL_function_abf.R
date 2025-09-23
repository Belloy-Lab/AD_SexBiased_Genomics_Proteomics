eQTL_colocalization <- function(merged, all_suggestive_snps, name_dset, LC_dir, LC_threshold) {
  
  print(paste("Running colocalization on", name_dset, sep = " "))
  
  bp = as.integer(bp_38_med)
  bp_s = as.integer(bp_38_start)
  bp_e = as.integer(bp_38_end)
  chrom = as.integer(chrom)
  range = 1e6
  
  print(chrom)
  print(bp_s)
  print(bp_e)
  print(bp)
  print(range)
  print(locus_index)
  print(study)
  print(tissue)
  print(discovery)
  print(LC_dir)
  print(LC_threshold)
  
  
  if (discovery != "pleiotropy") {
    
    colocalizing_genes <- merged %>% 
      filter(BP >= bp - range & BP <= bp + range) %>% 
      distinct(gene_name, .keep_all = T) %>% 
      arrange(gene_name) %>% 
      filter(!is.na(gene_name))
    
  } else if (discovery == "pleiotropy"){
    
    colocalizing_genes <- merged %>% 
      filter(BP >= bp_s - range & BP <= bp_e + range) %>% 
      distinct(gene_name, .keep_all = T) %>% 
      arrange(gene_name) %>% 
      filter(!is.na(gene_name))
    
  }
  
  # Check if colocalizing_genes is empty
  if (nrow(colocalizing_genes) == 0) {
    print("There are 0 genes for the basepair position in the merged file")
    return(NULL)
  }
  
  colocalization_results <- data.frame()
  
  # Loop through each colocalizing gene
  if (nrow(colocalizing_genes) > 0) {
    print(paste("Running colocalization on", locus_index, sep = " "))
    
    for (gene_row in 1:nrow(colocalizing_genes)) {
      
      gene_n = colocalizing_genes$gene_name[gene_row]
      print(paste("Running colocalization on", gene_n, sep = " "))
      
      
      if (discovery != "pleiotropy") {
        
        
        merged_c = merged %>% 
          filter(BP >= (bp - range) & BP <= (bp + range)) %>%
          filter(AD_BETA != 0,
                 QTL_BETA != 0,
                 gene_name == gene_n, AD_P > 0, QTL_P > 0) %>% 
          mutate(avarbeta_d = AD_SE^2,
                 varbeta_t2 = QTL_SE^2) %>% 
          distinct(SNP, .keep_all = T)
        
        # ## Prepare datasets - AD GWAS
        # ad_sub = merged %>% 
        #   mutate(varbeta_ad = AD_SE^2) %>% 
        #   filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        #   filter(varbeta_ad != 0, !is.na(SNP)) %>% 
        #   filter(gene_name == gene_n, AD_P > 0) %>% 
        #   distinct(SNP, .keep_all = T)
        # 
        # ## Prepare datasets - eQTL
        # qtl_sub <- merged %>%
        #   mutate(varbeta_qtl = QTL_SE^2) %>% 
        #   filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        #   filter(varbeta_qtl != 0, !is.na(SNP)) %>% 
        #   filter(gene_name == gene_n, QTL_P > 0) %>% 
        #   distinct(SNP, .keep_all = T)
        
        
      } else if (discovery == "pleiotropy"){
        
        merged_c <- merged %>% 
          filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
          filter(AD_BETA != 0,
                 QTL_BETA != 0,
                 gene_name == gene_n, AD_P > 0, QTL_P > 0) %>% 
          mutate(varbeta_ad = AD_SE^2,
                 varbeta_t2 = QTL_SE^2) %>% 
          distinct(SNP, .keep_all = T)
        
      }
      
      
      # Skip the iteration if either ad_sub or qtl_sub is empty
      if (nrow(merged_c) == 0) {
        next
      }
      
      
      D1 <- list(
        pvalues = merged_c$AD_P,
        type = "cc",
        snp = merged_c$SNP,
        position = merged_c$BP,
        N = merged_c$AD_N,
        MAF = merged_c$AD_EAF,
        s = merged_c$s[1]
      )
      
      check_dataset(D1)
      
      D2 <- list(
        pvalues = merged_c$QTL_P,
        type = "quant",
        snp = merged_c$SNP,
        position = merged_c$BP,
        N = merged_c$QTL_N,
        MAF = merged_c$QTL_EAF
      )
      
      check_dataset(D2)
      
      coloc_results = coloc.abf(dataset1 = D1, dataset2 = D2)
      
      
      pp_values <- coloc_results$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")]
      
      # Replace NaN values with 0
      pp_values[is.nan(pp_values)] <- 0
      
      # Get colocalization summary and add to the results data frame
      colocalization_row <- data.frame(
        CHROM = chrom,
        BP = bp,
        dset = name_dset,
        locus = locus_index,
        gene_name = gene_n,
        molecular_trait_id = "",
        PP0 = pp_values["PP.H0.abf"],
        PP1 = pp_values["PP.H1.abf"],
        PP2 = pp_values["PP.H2.abf"],
        PP3 = pp_values["PP.H3.abf"],
        PP4 = pp_values["PP.H4.abf"],
        lc_path = "",
        qtl_type = "eQTL",
        method = "abf",
        hit1 = 0,
        hit2 = 0,
        n_snps = n_distinct(merged_c$SNP),
        discovery = discovery,
        stringsAsFactors = FALSE
      )
      
      colocalization_results <- rbind(colocalization_results, colocalization_row)
      
      # Check if PP4 >= 0.70 and create locus compare plot
      if (colocalization_row$PP4 >= LC_threshold) {
        
        # Get variables for outcome variable 1 (Alzheimer's disease) - extract rsID and p-value
        ad_comp <- merged_c %>%
          dplyr::select(SNP, AD_P) %>%
          dplyr::rename(rsid = SNP, pval = AD_P) %>% 
          arrange(pval)
        
        # Get variables for outcome variable 2 (QTL) - extract rsID and p-value
        qtl_comp <- merged_c %>%
          dplyr::select(SNP, QTL_P) %>%
          dplyr::rename(rsid = SNP, pval = QTL_P) %>% 
          arrange(pval)
        
        # Define lead SNP (Top SNP in AD dataset)
        rsID1 = ad_comp$rsid[1]
        rsID2 = qtl_comp$rsid[1]
        
        # Check if snp is in the intersection of ad_comp and qtl_comp
        if (!(rsID1 %in% intersect(ad_comp$rsid, qtl_comp$rsid))) {
          next
        }
        
        # Split name_dset into components
        t1 <- paste("AD", stratum, "GWAS", sep = " ")
        t2 <- paste(study, tissue, gene_n, qtl_type, sep = " ")
        
        
        # Continue to next iteration in case of error
        tryCatch({
          lc <- locuscompare(in_fn1 = ad_comp, in_fn2 = qtl_comp, snp = rsID1, title = t1, title2 = t2, genome = "hg38", lz_ylab_linebreak = T)
          file_name <- paste(LC_dir, colocalization_row$dset, "_COLOC_", colocalization_row$locus, "_", colocalization_row$gene_name, "_", rsID1, "_", discovery, "_topAD.pdf", sep = "")
          
          # Save file
          pdf(file_name, height = 6, width = 12)
          print(lc)
          dev.off()
        }, error = function(e) {
          # Print error message if needed (optional)
          message("Skipping iteration due to error: ", e$message)
        })
        
        
        # Check if snp is in the intersection of ad_comp and qtl_comp
        if (!(rsID2 %in% intersect(ad_comp$rsid, qtl_comp$rsid))) {
          next
        }
        
        # Split name_dset into components
        t1 <- paste("AD", stratum, "GWAS", sep = " ")
        t2 <- paste(study, tissue, gene_n, qtl_type, sep = " ")
        
        
        # Continue to next iteration in case of error
        tryCatch({
          lc <- locuscompare(in_fn1 = ad_comp, in_fn2 = qtl_comp, snp = rsID2, title = t1, title2 = t2, genome = "hg38", lz_ylab_linebreak = T)
          file_name <- paste(LC_dir, colocalization_row$dset, "_COLOC_", colocalization_row$locus, "_", colocalization_row$gene_name, "_", rsID2, "_", discovery, "_topQTL.pdf", sep = "")
          
          # Save file
          pdf(file_name, height = 6, width = 12)
          print(lc)
          dev.off()
        }, error = function(e) {
          # Print error message if needed (optional)
          message("Skipping iteration due to error: ", e$message)
        })
        
      }
      
    }
    
  }
  return(colocalization_results)   
  
}


