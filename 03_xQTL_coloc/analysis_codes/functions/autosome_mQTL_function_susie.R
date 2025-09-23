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
        
        ## Prepare datasets - AD GWAS & QTL
        merged_c = merged %>% 
          mutate(varbeta_ad = AD_SE^2,
                 varbeta_t2 = QTL_SE^2) %>% 
          filter(BP >= (bp - range) & BP <= (bp + range)) %>%
          filter(varbeta_ad != 0, varbeta_t2 != 0, !is.na(SNP)) %>% 
          filter(CpG == gene_n,
                 QTL_P > 0) %>% 
          distinct(SNP, .keep_all = T)
        
        
        
      } else if (discovery == "pleiotropy"){
        
        ## Prepare datasets - AD GWAS
        merged_c = merged %>% 
          mutate(varbeta_ad = AD_SE^2,
                 varbeta_t2 = QTL_SE^2) %>% 
          filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
          filter(varbeta_ad != 0, varbeta_t2 != 0, !is.na(SNP)) %>% 
          filter(CpG == gene_n,
                 QTL_P > 0) %>% 
          distinct(SNP, .keep_all = T)
        
        
      }
      
      
      snps_for_ld <- list(merged_c$SNP)
      
      # Check if colocalizing_genes is empty
      if (nrow(merged_c) == 0) {
        print("There are 0 SNPs for this gene in the merged file")
        next
      }
      
      
      pl <- "plink1.9"
      in_fold <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_reference_genotype_data/EU_TOPMed/"
      in_file <- paste0("TOPMed_ch", chrom)
      
      fold_w <- "/storage1/fs1/belloy/Active/05_Projects/noahc/opt_QTL/Xchr_LD/mQTL/"
      snp_file <- paste0(fold_w, study, "_", tissue, "_mQTL_", discovery, "_", locus_index, ".txt")
      fwrite(snps_for_ld, snp_file, sep = "\t", col.names = FALSE)
      
      
      file_ld <- paste("LD_matrix_", study, tissue, "mQTL", discovery, locus_index, sep = "_")
      
      command <- paste(pl, 
                       "--bfile", paste(in_fold, in_file, sep = ""),
                       "--allow-no-sex",
                       "--extract", snp_file,
                       "--keep-allele-order",
                       "--make-bed --out", paste(fold_w, file_ld, ".data", sep = ""),
                       sep = " ")
      system(command)
      
      command <- paste(pl, 
                       "--bfile", paste(fold_w, file_ld, ".data", sep = ""),
                       "--allow-no-sex",
                       "--keep-allele-order",
                       "--r --matrix",
                       "--out", paste(fold_w, file_ld, sep = ""),
                       sep = " ")
      system(command)
      
      bim <- fread(paste(fold_w, file_ld, ".data.bim", sep = ""))
      snps <- unlist(bim$V2)
      
      
      LD <- fread(paste(fold_w, file_ld, ".ld", sep = ""))
      LD <- as.matrix(LD)
      colnames(LD) <- snps
      rownames(LD) <- snps
      
      print(system("free -h"))
      # Remove LD Matrix
      file.remove(paste0(fold_w, paste0(file_ld, ".data.bed")))
      file.remove(paste0(fold_w, paste0(file_ld, ".data.bim")))
      file.remove(paste0(fold_w, paste0(file_ld, ".data.fam")))
      file.remove(paste0(fold_w, paste0(file_ld, ".ld")))
      file.remove(paste0(fold_w, paste0(file_ld, ".log")))
      file.remove(paste0(fold_w, paste0(file_ld, ".nosex")))
      file.remove(paste0(fold_w, paste0(file_ld, ".data.log")))
      file.remove(paste0(fold_w, paste0(file_ld, ".data.nosex")))
      file.remove(snp_file)
      
      
      # AD trait coloc.susie preperation
      ADl <- as.list(merged_c[,c("AD_BETA","varbeta_ad","SNP","BP","AD_EAF")])
      names(ADl)[1:5] <- c("beta","varbeta","snp","position","MAF")
      names(ADl$beta) <- ADl$snp
      names(ADl$varbeta) <- ADl$snp
      ADl$N <- median(merged$AD_N, na.rm = TRUE)
      
      
      matching_snps_ADl <- ADl$snp %in% snps
      ADl$beta <- ADl$beta[matching_snps_ADl]
      ADl$varbeta <- ADl$varbeta[matching_snps_ADl]
      ADl$snp <- ADl$snp[matching_snps_ADl]
      ADl$position <- ADl$position[matching_snps_ADl]
      ADl$MAF <- ADl$MAF[matching_snps_ADl]
      ADl$LD <- LD
      ADl$type <- "cc"
      cc <- ADl
      
      common_snps <- intersect(cc$snp, colnames(cc$LD))
      cc_index <- match(common_snps, cc$snp)
      ld_index <- match(common_snps, colnames(cc$LD))
      
      cc$beta <- cc$beta[cc_index]
      cc$varbeta <- cc$varbeta[cc_index]
      cc$snp <- cc$snp[cc_index]
      cc$position <- cc$position[cc_index]
      cc$MAF <- cc$MAF[cc_index]
      cc$LD <- cc$LD[ld_index, ld_index]
      
      
      # Trait 2 coloc.susie preperation
      traitl <- as.list(merged_c[,c("QTL_BETA","varbeta_t2","SNP","BP","QTL_EAF")])
      names(traitl)[1:5] <- c("beta","varbeta","snp","position","MAF")
      names(traitl$beta) <- traitl$snp
      names(traitl$varbeta) <- traitl$snp
      traitl$N <- round(median(merged$QTL_N, na.rm = TRUE))
      
      matching_snps_traitl <- traitl$snp %in% snps
      traitl$beta <- traitl$beta[matching_snps_traitl]
      traitl$varbeta <- traitl$varbeta[matching_snps_traitl]
      traitl$snp <- traitl$snp[matching_snps_traitl]
      traitl$position <- traitl$position[matching_snps_traitl]
      traitl$MAF <- traitl$MAF[matching_snps_traitl]
      traitl$LD <- LD
      traitl$type <- "quant"
      eq <- traitl
      
      common_snps <- intersect(eq$snp, colnames(eq$LD))
      eq_index <- match(common_snps, eq$snp)
      ld_index <- match(common_snps, colnames(eq$LD))
      
      eq$beta <- eq$beta[eq_index]
      eq$varbeta <- eq$varbeta[eq_index]
      eq$snp <- eq$snp[eq_index]
      eq$position <- eq$position[eq_index]
      eq$MAF <- eq$MAF[eq_index]
      eq$LD <- eq$LD[ld_index, ld_index]
      
      
      if (length(eq$LD) == 1 | length(cc$LD) == 1){
        next
      }
      
      print(paste("Starting to run coloc susie on:", gene_n, sep = " "))
      print(system("free -h"))
      ## Run coloc.susie with sequenced coverages
      coverage_values <- seq(0.95, 0.05, by = -0.1)
      cs_length_S1 <- 0
      cs_length_S2 <- 0
      for (coverage in coverage_values) {
        tryCatch({
          S1 <- runsusie(cc, coverage = coverage, max_iter = 10000)
          cs_length_S1 <- length(S1$sets$cs)
          # print(system("free -h"))
          
          S2 <- runsusie(eq, coverage = coverage, max_iter = 10000)
          cs_length_S2 <- length(S2$sets$cs)
          # print(system("free -h"))
          
          if (cs_length_S1 > 0 & cs_length_S2 > 0) {
            break
          }
        }, error = function(e) {
          # Print the error message if you want to
          message("Error encountered: ", e)
          # Continue to the next iteration
        })
      }
      
      gc()
      
      
      if (!exists("S1") || !exists("S2")) {
        print(paste("S1 or S2 does not exist for: ", gene_n, sep = " "))
        next
      }
      
      if (class(S1) != "susie" || class(S2) != "susie") {
        print(paste("S1 or S2 failed to converge for: ", gene_n, sep = " "))
        next
      }
      
      
      # Run coloc.susie and return the maximum posterior probability
      out <- coloc.susie(S1, S2)
      
      if (is.null(out$summary$PP.H4.abf)) {
        print("No coloc results")
        coloc_PP0 = 0
        coloc_PP1 = 0
        coloc_PP2 = 0
        coloc_PP3 = 0
        coloc_PP4 = 0
        coloc_hit1 = 0
        coloc_hit2 = 0
      }else{
        print("Coloc results found")
        print(out$summary)
        index_row <- which.max(out$summary$PP.H4.abf)
        coloc_PP0 = out$summary$PP.H0.abf[index_row]
        coloc_PP1 = out$summary$PP.H1.abf[index_row]
        coloc_PP2 = out$summary$PP.H2.abf[index_row]
        coloc_PP3 = out$summary$PP.H3.abf[index_row]
        coloc_PP4 = out$summary$PP.H4.abf[index_row]
        
        if (is.na(coloc_PP4)) {
          coloc_PP4 = 0
        }
        
        coloc_hit1 = out$summary$hit1[index_row]
        coloc_hit2 = out$summary$hit2[index_row]
      }
      
      
      colocalization_row <- data.frame(
        CHROM = chrom,
        BP = bp,
        dset = name_dset,
        locus = locus_index,
        gene_name = gene_n,
        molecular_trait_id = "",
        PP0 = coloc_PP0,
        PP1 = coloc_PP1,
        PP2 = coloc_PP2,
        PP3 = coloc_PP3,
        PP4 = coloc_PP4,
        lc_path = "",
        qtl_type = "mQTL",
        method = "susie",
        hit1 = coloc_hit1,
        hit2 = coloc_hit2,
        n_snps = nrow(LD),
        discovery = discovery,
        # tissue = tissue,
        # analyte = analyte,
        stringsAsFactors = FALSE
      )
      
      colocalization_results <- rbind(colocalization_results, colocalization_row)
      
      
      # Check if PP4 >= 0.70 and create locus compare plot
      if (coloc_PP4 >= LC_threshold) {
        
        ad_comp <- merged_c %>%
          dplyr::select(SNP, AD_P) %>%
          rename(rsid = SNP, pval = AD_P) %>% 
          arrange(pval)
        
        qtl_comp <- merged_c %>%
          dplyr::select(SNP, QTL_P) %>%
          rename(rsid = SNP, pval = QTL_P) %>% 
          arrange(pval)
        
        rsID = ad_comp$rsid[1]
        
        # Check if snp from colocalization_row$rsID is in the intersection of ad_comp and qtl_comp
        if (!(rsID %in% intersect(ad_comp$rsid, qtl_comp$rsid))) {
          next
        }
        
        # Split name_dset into components
        t1 <- paste("AD", stratum, "GWAS", sep = " ")
        t2 <- paste(study, tissue, gene_n, qtl_type, sep = " ")
        
        
        # Continue to next iteration in case of error
        tryCatch({
          lc <- locuscompare(in_fn1 = ad_comp, in_fn2 = qtl_comp, snp = rsID, title = t1, title2 = t2, genome = "hg38", lz_ylab_linebreak = TRUE)
          file_name <- paste(LC_dir, colocalization_row$dset, "_COLOC_", colocalization_row$locus, "_", colocalization_row$gene_name, "_", rsID, "_", discovery, "_SuSiE.pdf", sep = "")
          
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
      
      
