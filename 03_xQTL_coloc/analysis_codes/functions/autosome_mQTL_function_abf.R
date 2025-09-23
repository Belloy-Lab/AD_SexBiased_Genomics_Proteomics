mQTL_colocalization <- function(merged, all_suggestive_snps, name_dset, out_dir) {
  
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
  
  if (discovery != "pleiotropy") {
  
  colocalizing_genes <- merged %>% 
    filter(BP >= bp - range & BP <= bp + range) %>% 
    distinct(CpG, .keep_all = T) 
  
  } else if (discovery == "pleiotropy"){
    
    colocalizing_genes <- merged %>% 
      filter(BP >= bp_s - range & BP <= bp_e + range) %>% 
      distinct(CpG, .keep_all = T) 
    
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
      
      gene_n = colocalizing_genes$CpG[gene_row]
      gene_pos = colocalizing_genes$CpGpos[gene_row]
      print(paste("Running colocalization on", gene_n, sep = " "))
      
      
      if (discovery != "pleiotropy") {
        
      ## Prepare datasets - AD GWAS
      ad_sub = merged %>% 
        mutate(varbeta_ad = AD_SE^2) %>% 
        filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        filter(varbeta_ad != 0, !is.na(SNP)) %>% 
        filter(CpG == gene_n,
               AD_P > 0) %>% 
        distinct(SNP, .keep_all = T)
      
      ## Prepare datasets - mQTL
      qtl_sub <- merged %>%
        mutate(varbeta_qtl = QTL_SE^2) %>% 
        filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        filter(varbeta_qtl != 0, !is.na(SNP)) %>% 
        filter(CpG == gene_n,
               QTL_P > 0) %>% 
        distinct(SNP, .keep_all = T)
      
      } else if (discovery == "pleiotropy"){
        
        ## Prepare datasets - AD GWAS
        ad_sub = merged %>% 
          mutate(varbeta_ad = AD_SE^2) %>% 
          filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
          filter(varbeta_ad != 0, !is.na(SNP)) %>% 
          filter(CpG == gene_n, AD_P > 0) %>% 
          distinct(SNP, .keep_all = T)
        
        ## Prepare datasets - mQTL
        qtl_sub <- merged %>%
          mutate(varbeta_qtl = QTL_SE^2) %>% 
          filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
          filter(varbeta_qtl != 0, !is.na(SNP)) %>% 
          filter(CpG == gene_n,
                 QTL_P > 0) %>% 
          distinct(SNP, .keep_all = T)
        
        
      }
      
      # Skip the iteration if either ad_sub or qtl_sub is empty
      if (nrow(ad_sub) == 0 || nrow(qtl_sub) == 0) {
        next
      }
      
      
      D1 <- list(
        pvalues = ad_sub$AD_P,
        type = "cc",
        snp = ad_sub$SNP,
        position = ad_sub$BP,
        N = ad_sub$AD_N,
        MAF = ad_sub$AD_EAF,
        s = ad_sub$s[1]
      )
      
      check_dataset(D1)
      
      D2 <- list(
        pvalues = qtl_sub$QTL_P,
        type = "quant",
        snp = qtl_sub$SNP,
        position = qtl_sub$BP,
        N = qtl_sub$QTL_N,
        MAF = qtl_sub$QTL_EAF
      )
      
      check_dataset(D2)
      
      coloc_results = coloc.abf(dataset1 = D1, dataset2 = D2)
      
      # Get colocalization summary and add to the results data frame
      colocalization_row <- data.frame(
        CHROM = chrom,
        BP = bp,
        dset = name_dset,
        locus = locus_index,
        gene_name = gene_n,
        molecular_trait_id = gene_pos,
        PP0 = coloc_results$summary["PP.H0.abf"],
        PP1 = coloc_results$summary["PP.H1.abf"],
        PP2 = coloc_results$summary["PP.H2.abf"],
        PP3 = coloc_results$summary["PP.H3.abf"],
        PP4 = coloc_results$summary["PP.H4.abf"],
        lc_path = "",
        qtl_type = "mQTL",
        method = "abf",
        hit1 = 0,
        hit2 = 0, 
        n_snps = n_distinct(ad_sub$SNP),
        discovery = discovery,
        stringsAsFactors = FALSE
      )
      
      colocalization_results <- rbind(colocalization_results, colocalization_row)
      
      # Check if PP4 >= 0.70 and create locus compare plot
      if (colocalization_row$PP4 >= LC_threshold) {
        
        # Get variables for outcome variable 1 (Alzheimer's disease) - extract rsID and p-value
        ad_comp <- ad_sub %>%
          dplyr::select(SNP, AD_P) %>%
          dplyr::rename(rsid = SNP, pval = AD_P) %>% 
          arrange(pval)
        
        # Get variables for outcome variable 2 (QTL) - extract rsID and p-value
        qtl_comp <- qtl_sub %>%
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


