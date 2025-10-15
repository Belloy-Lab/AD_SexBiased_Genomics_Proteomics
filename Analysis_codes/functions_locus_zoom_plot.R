find_best_snp <- function(data_list,signal,chr,bp_start,bp_end){
    data <- data_list[[signal]]
    AD_temp <- data %>%
        dplyr::filter(CHR == chr) %>%
        dplyr::filter(BP >= bp_start & BP <= bp_end)

    AD_temp$P <- as.numeric(AD_temp$P)
    AD_temp$BP <- as.numeric(AD_temp$BP)
    best_gwas_snp <- AD_temp %>%
        dplyr::filter(!is.na(P)) %>%
        dplyr::arrange(P, abs(BP - mean(BP))) %>%
        dplyr::slice(1)
    chr = best_gwas_snp$CHR
    bp = best_gwas_snp$BP
    snp = best_gwas_snp$SNP
    print(paste0("Best SNP at CHR: ", chr, " BP: ", bp))
    return(c(chr,bp,snp))
}

subset_pQTL <- function(pQTL,Gene_id){
    pQTL_temp <- pQTL %>%
        dplyr::filter(gene_id == Gene_id)
    print(paste0("Number of pQTL is ",nrow(pQTL_temp)))
    return(pQTL_temp)
}

merge_dataset <- function(AD_female,AD_male,pQTL_female,pQTL_male,pQTL_combined,chr,bp,range){
    range <- as.numeric(range)
    chr <- as.numeric(chr)
    bp <- as.numeric(bp)
    AD_female_data_sub <- AD_female %>% dplyr::filter(CHR == chr & BP >= bp - range & BP <= bp + range)
    AD_male_data_sub <- AD_male %>% dplyr::filter(CHR == chr & BP >= bp - range & BP <= bp + range)
    pQTL_female_data_sub <- pQTL_female %>% dplyr::filter(CHR == chr & BP >= bp - range & BP <= bp + range)
    pQTL_male_data_sub <- pQTL_male %>% dplyr::filter(CHR == chr & BP >= bp - range & BP <= bp + range)
    pQTL_combined_data_sub <- pQTL_combined %>% dplyr::filter(CHR == chr & BP >= bp - range & BP <= bp + range)

    AD_female_data_sub_12 <- merge(AD_female_data_sub,pQTL_combined_data_sub[,c("id12","gene_id")],by.x = "id12",by.y = "id12") %>% dplyr::select(-gene_id)
    AD_female_data_sub_21 <- merge(AD_female_data_sub,pQTL_combined_data_sub[,c("id21","gene_id")],by.x = "id12",by.y = "id21") %>% dplyr::select(-gene_id)
    AD_female_data_sub_merged <- rbind(AD_female_data_sub_12,AD_female_data_sub_21) %>% dplyr::select(-id12, -id21)

    print(paste0("Number of AD GWAS Female SNPs: ", nrow(AD_female_data_sub_merged)))

    AD_male_data_sub_12 <- merge(AD_male_data_sub,pQTL_combined_data_sub[,c("id12","gene_id")],by.x = "id12",by.y = "id12") %>% dplyr::select(-gene_id)
    AD_male_data_sub_21 <- merge(AD_male_data_sub,pQTL_combined_data_sub[,c("id21","gene_id")],by.x = "id12",by.y = "id21") %>% dplyr::select(-gene_id)
    AD_male_data_sub_merged <- rbind(AD_male_data_sub_12,AD_male_data_sub_21) %>% dplyr::select(-id12, -id21)

    print(paste0("Number of AD GWAS Male SNPs: ", nrow(AD_male_data_sub_merged)))

    AD_female_data_sub_merged <- AD_female_data_sub_merged %>% 
        dplyr::mutate(
            id12 = paste(CHR, BP, A1, A2, sep = ":"),
            id21 = paste(CHR, BP, A2, A1, sep = ":"))
    
    AD_male_data_sub_merged <- AD_male_data_sub_merged %>% 
        dplyr::mutate(
            id12 = paste(CHR, BP, A1, A2, sep = ":"),
            id21 = paste(CHR, BP, A2, A1, sep = ":"))

    pQTL_female_data_sub_12 <- merge(AD_female_data_sub_merged[,c("id12","SNP")],pQTL_female_data_sub,by.x = "id12",by.y = "id12") %>% dplyr::select(-id12, -id21)
    pQTL_female_data_sub_21 <- merge(AD_female_data_sub_merged[,c("id12","SNP")],pQTL_female_data_sub,by.x = "id12",by.y = "id21") %>% dplyr::select(-id12.y, -id12)
    pQTL_female_data_sub_merged <- rbind(pQTL_female_data_sub_12,pQTL_female_data_sub_21)

    print(paste0("Number of pQTL Female SNPs: ", nrow(pQTL_female_data_sub_merged)))

    pQTL_male_data_sub_12 <- merge(AD_male_data_sub_merged[,c("id12","SNP")],pQTL_male_data_sub,by.x = "id12",by.y = "id12") %>% dplyr::select(-id12, -id21)
    pQTL_male_data_sub_21 <- merge(AD_male_data_sub_merged[,c("id12","SNP")],pQTL_male_data_sub,by.x = "id12",by.y = "id21") %>% dplyr::select(-id12.y, -id12)
    pQTL_male_data_sub_merged <- rbind(pQTL_male_data_sub_12,pQTL_male_data_sub_21)

    print(paste0("Number of pQTL Male SNPs: ", nrow(pQTL_male_data_sub_merged)))

    pQTL_combined_data_sub_12 <- merge(AD_female_data_sub_merged[,c("id12","SNP")],pQTL_combined_data_sub,by.x = "id12",by.y = "id12") %>% dplyr::select(-id12, -id21)
    pQTL_combined_data_sub_21 <- merge(AD_female_data_sub_merged[,c("id12","SNP")],pQTL_combined_data_sub,by.x = "id12",by.y = "id21") %>% dplyr::select(-id12.y, -id12)
    pQTL_combined_data_sub_merged <- rbind(pQTL_combined_data_sub_12,pQTL_combined_data_sub_21)

    print(paste0("Number of pQTL Combined SNPs: ", nrow(pQTL_combined_data_sub_merged)))

    merged_data_list <- list(
        AD_female    = AD_female_data_sub_merged,
        AD_male      = AD_male_data_sub_merged,
        pQTL_female  = pQTL_female_data_sub_merged,
        pQTL_male    = pQTL_male_data_sub_merged,
        pQTL_combined= pQTL_combined_data_sub_merged
    )
    common_id <- Reduce(intersect, list(
        AD_female_data_sub_merged$SNP,
        AD_male_data_sub_merged$SNP,
        pQTL_female_data_sub_merged$SNP,
        pQTL_male_data_sub_merged$SNP,
        pQTL_combined_data_sub_merged$SNP
    ))
    return(list(
        merged_data_list = merged_data_list,
        common_id       = common_id
    ))
}

calculate_LD <- function(common_id,snp,gene){
    in_file <- "TOPMed_ch"

    fold_w <- "/Path/to/data_temp_save/LD_matrix_for_plots/"
    file_snp <- "common_id"
    pl <- "plink1.9"
    in_fold <- "/Path/to/EU_TOPMed/data/in/PLINK/format/separeted/byCHR/"
    in_file <- "TOPMed_ch"
    keep_file <- "/Path/to/EU_nHISP_TOPMed_list_Unrel.txt"
    print(paste0("common_id length: ", length(common_id)))
    print(paste0("common_id saved to ", paste0(fold_w, file_snp,"_",snp,"_",gene,".txt")))

    fwrite(data.frame(common_id), paste0(fold_w, file_snp,"_",snp,"_",gene,".txt"), sep = "\t", col.names = FALSE)

    print(paste0("Best SNP: ", snp))
    print(paste0("common_id saved to ", paste0(fold_w, file_snp,"_",snp,".txt")))
    in_file_chr <- paste0(in_file,chr)

    command <- paste(pl, 
                    "--bfile", paste0(in_fold, in_file_chr),
                    "--allow-no-sex",
                    "--extract", paste0(fold_w, file_snp,"_",snp,"_",gene,".txt"),
                    "--keep-allele-order",
                    "--make-bed --out", paste0(fold_w, paste0(snp,"_",gene, ".data")),
                    sep = " ")
    system(command)

        command <- paste(pl, 
                    "--bfile", paste0(fold_w, paste0(snp,"_",gene, ".data")),
                    "--allow-no-sex",
                    "--keep-allele-order",
                    "--keep", keep_file,
                    "--r2 yes-really --matrix",
                    "--out", paste0(fold_w, snp,"_",gene),
                    sep = " ")
    system(command)

    return(c(paste0(fold_w, paste0(snp,"_",gene,".data.bim")),paste0(fold_w, paste0(snp,"_",gene, ".ld"))))
}

read_LD <- function(merged_data_list,snp,bim_path,ld_path){
    AD_female_merged <- merged_data_list[["AD_female"]]
    AD_male_merged <- merged_data_list[["AD_male"]]
    pQTL_female_merged <- merged_data_list[["pQTL_female"]]
    pQTL_male_merged <- merged_data_list[["pQTL_male"]]
    pQTL_combined_merged <- merged_data_list[["pQTL_combined"]]
    # Load the BIM file and corresponding LD matrix
    bim_file <- bim_path
    if (file.exists(bim_file)) {
    bim <- fread(bim_file)
    snps_name <- unlist(bim$V2)  # Extract the SNP names

    LD_file <- ld_path
    LD <- fread(LD_file)
    LD <- as.matrix(LD)  # Convert LD data to a matrix

    # Check if dimensions of SNP names and LD matrix match
    if (!is.null(LD) && length(snps_name) == ncol(LD) && length(snps_name) == nrow(LD)) {
        colnames(LD) <- snps_name
        rownames(LD) <- snps_name

        LD_df <- as.data.frame(LD)  # Convert LD matrix to data frame
        LD_df$V1 <- snps_name
    } else {
        warning("Mismatch between SNPs and LD matrix dimensions")
    }
    } else {
    warning("BIM file not found, skipping iteration")
    }

    print(paste0("Number of SNPs in LD matrix: ", nrow(LD_df)))

    ld_values <- LD_df[,c("V1",snp)] %>% dplyr::rename(SNP = V1, r2 = snp) %>% arrange(desc(SNP))

    AD_female_data_sub_merged <- AD_female_merged %>% dplyr::left_join(ld_values, by = "SNP")
    AD_male_data_sub_merged <- AD_male_merged %>% dplyr::left_join(ld_values, by = "SNP")
    pQTL_female_data_sub_merged <- pQTL_female_merged %>% dplyr::left_join(ld_values, by = "SNP")
    pQTL_male_data_sub_merged <- pQTL_male_merged %>% dplyr::left_join(ld_values, by = "SNP")
    pQTL_combined_data_sub_merged <- pQTL_combined_merged %>% dplyr::left_join(ld_values, by = "SNP")

    merged_data_list_with_LD <- list()
    merged_data_list_with_LD[["AD_female"]] <- AD_female_data_sub_merged
    merged_data_list_with_LD[["AD_male"]] <- AD_male_data_sub_merged
    merged_data_list_with_LD[["pQTL_female"]] <- pQTL_female_data_sub_merged
    merged_data_list_with_LD[["pQTL_male"]] <- pQTL_male_data_sub_merged
    merged_data_list_with_LD[["pQTL_combined"]] <- pQTL_combined_data_sub_merged

    return(merged_data_list_with_LD)
}

find_AD_ylim <- function(AD_female_data_sub_merged,AD_male_data_sub_merged,snp,bp){
    bp <- as.numeric(bp)
    bp_start <- min(AD_female_data_sub_merged$BP)-1
    bp_end <- max(AD_female_data_sub_merged$BP)+1
    loc1 <- locus(data = AD_female_data_sub_merged, chrom = "CHR", pos = "BP", p = "P", index_snp = snp, LD = "r2", flank = c(bp-bp_start,bp_end-bp), ens_db = "EnsDb.Hsapiens.v86")
    loc2 <- locus(data = AD_male_data_sub_merged, chrom = "CHR", pos = "BP", p = "P", index_snp = snp, LD = "r2", flank = c(bp-bp_start,bp_end-bp), ens_db = "EnsDb.Hsapiens.v86")
    y_max_AD <- min(c(min(loc1$data$P),min(loc2$data$P)))
    y_max_AD <- ceiling(-log10(y_max_AD))
    print(paste0("y_max_AD: ", y_max_AD))

    return(y_max_AD)
}

find_pQTL_ylim <- function(pQTL_female_data_sub_merged,pQTL_male_data_sub_merged,pQTL_combined_data_sub_merged,snp,bp){
    bp <- as.numeric(bp)
    bp_start <- min(pQTL_female_data_sub_merged$BP)-1
    bp_end <- max(pQTL_female_data_sub_merged$BP)+1
    loc3 <- locus(data = pQTL_female_data_sub_merged, chrom = "CHR", pos = "BP", p = "P", index_snp = snp, LD = "r2", flank = c(bp-bp_start,bp_end-bp), ens_db = "EnsDb.Hsapiens.v86")
    loc4 <- locus(data = pQTL_male_data_sub_merged, chrom = "CHR", pos = "BP", p = "P", index_snp = snp, LD = "r2", flank = c(bp-bp_start,bp_end-bp), ens_db = "EnsDb.Hsapiens.v86")
    loc5 <- locus(data = pQTL_combined_data_sub_merged, chrom = "CHR", pos = "BP", p = "P", index_snp = snp, LD = "r2", flank = c(bp-bp_start,bp_end-bp), ens_db = "EnsDb.Hsapiens.v86")
    y_max_pQTL <- min(c(min(loc3$data$P),min(loc4$data$P),min(loc5$data$P)))
    y_max_pQTL <- ceiling(-log10(y_max_pQTL))
    print(paste0("y_max_pQTL: ", y_max_pQTL))

    return(y_max_pQTL)
}

locus_sub_plot_38 <- function(data, snp, bp, bp_start, bp_end, y_max, y_label, plot_title = NULL) {
    bp       <- as.numeric(bp)
    bp_start <- as.numeric(bp_start)
    bp_end   <- as.numeric(bp_end)
    y_max    <- as.numeric(y_max)

    y_label <- bquote(bold(.(y_label)) ~ -log[10](P))
    loc0 <- locus(data = data, chrom = "CHR", pos = "BP", p = "P", 
                index_snp = snp, LD = "r2", 
                flank = c(bp - bp_start, bp_end - bp), 
                ens_db = "EnsDb.Hsapiens.v86")
    loc0 <- link_recomb(loc0, genome = "hg38")

    p <- gg_scatter(loc0, pcutoff = 5e-08, align = TRUE, xticks = FALSE, 
                    cex.lab = 1.25, cex.axis = 1.1, font.axis = 2, font.lab = 2,
                    ylim = c(0, y_max), legend_pos = NULL, ylab = y_label) +
        theme(
            axis.title.y = element_text(size = 13, face = "bold"),     
            axis.title.y.right = element_text(size = 10, face = "plain"),
            plot.title = element_text(
                size = 14, 
                face = "bold",
                hjust = 0.5,  
                margin = margin(b = 10)
            )
        )

    if (!is.null(plot_title)) {
        p <- p + labs(title = plot_title)
    }

    p <- p + 
         geom_point(data = subset(loc0$data, SNP == loc0$index_snp), 
                    mapping = aes(x = BP / 1e6, y = logP), 
                    shape = 23,       
                    fill = "purple",  
                    color = "black",  
                    size = 3,        
                    inherit.aes = FALSE)  

    return(p)
}

genetracks_sub_plot <- function(loc0, Gene_name1, Gene_name2,maxrows_num) {
    maxrows_num <- as.numeric(maxrows_num)
    print(paste0("Gene_name1: ", Gene_name1))
    print(paste0("Gene_name2: ", Gene_name2))
    print(paste0("maxrows_num: ", maxrows_num))
  g <- gg_genetracks(loc0, 
             highlight = c(Gene_name1, Gene_name2), 
             filter_gene_biotype = "protein_coding", 
             gene_col = "skyblue", 
             showExons = TRUE, 
             xticks = TRUE,
             maxrows = maxrows_num)

  return(g)
}

cleanup_LD <- function(snp,gene){
    fold_w <- "/Path/to/data_temp_save/LD_matrix_for_plots/"
    file_snp <- "common_id"
    file.remove(paste0(fold_w, file_snp,"_",snp,"_",gene,".txt"))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".data.bed")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".data.bim")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".data.fam")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".ld")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".log")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".nosex")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".data.log")))
    file.remove(paste0(fold_w, paste0(snp,"_",gene, ".data.nosex")))
    print("Files removed.")
    return(0)
}

locus_compare_plot_38 <- function(AD_data, pQTL_data, chr, index_snp, ad_title, pqtl_title,position) {
  AD_data_P <- AD_data %>% 
    dplyr::select(all_of(c("SNP", "P", "CHR", "BP"))) %>% 
    dplyr::mutate(P = as.numeric(P)) %>% 
    dplyr::rename(rsid = SNP, pval1 = P)
  
  pQTL_data_P <- pQTL_data %>% 
    dplyr::select(all_of(c("SNP", "P", "CHR", "BP"))) %>% 
    dplyr::mutate(P = as.numeric(P)) %>% 
    dplyr::rename(rsid = SNP, pval2 = P)
  
  merged_data <- merge(AD_data_P, pQTL_data_P, by = "rsid")
  
  merged_data <- merged_data %>% 
    dplyr::mutate(
      logp1 = -log10(pval1),
      logp2 = -log10(pval2),
      chr = CHR.x,
      pos = BP.x
    ) %>% 
    dplyr::mutate(label = if_else(rsid == index_snp, rsid, ""))

  merged_data <- merged_data %>% 
    dplyr::mutate(label = NA_character_)

  ld_snp <- retrieve_LD(chr, index_snp, 'EUR')
  color_df <- locuscomparer::assign_color(merged_data$rsid, snp = index_snp, ld = ld_snp)
  shape_vec <- ifelse(merged_data$rsid == index_snp, 23, 21) %>% stats::setNames(merged_data$rsid)
  size_vec <- ifelse(merged_data$rsid == index_snp, 3, 2) %>% stats::setNames(merged_data$rsid)
  
  sc <- locuscomparer::make_scatterplot(
    merged = merged_data,
    title1 = ad_title,
    title2 = pqtl_title,
    color = color_df,
    shape = shape_vec,
    size = size_vec,
    legend = TRUE,
    legend_position = position
  )

  sc$layers <- sc$layers[!sapply(sc$layers, function(x) {
        inherits(x$geom, "GeomTextRepel") | inherits(x$geom, "GeomText")
    })]

  if (any(grepl("geom_point", sapply(sc$layers, function(x) class(x$geom)[1])))) {
    sc$layers <- sc$layers[-which(sapply(sc$layers, function(x) 
      inherits(x$geom, "GeomPoint") & 
      !is.null(x$data$label) & 
      any(x$data$label != "")))]
  }

  sc <- sc +
    labs(x = bquote(bold(.(ad_title)) ~ -log[10](P)),
         y = bquote(bold(.(pqtl_title)) ~ -log[10](P))) +
    theme(
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13)
    )
  
  return(sc)
}

S14_plots <- function(pQTL_female, pQTL_male, pQTL_combined, 
                      AD_GWAS_female, AD_GWAS_male, data_list,
                      signal1, signal2, bp_start, bp_end, chr,
                      Gene_id, Gene_name1, Gene_name2, range,index,level,maxrows_num,position,ld_table,out_dir) {

    pQTL_female_filtered1   <- subset_pQTL(pQTL_female, Gene_id)
    pQTL_male_filtered1     <- subset_pQTL(pQTL_male, Gene_id)
    pQTL_combined_filtered1 <- subset_pQTL(pQTL_combined, Gene_id)

    best_snp1 <- find_best_snp(data_list, signal1, chr, bp_start, bp_end)
    mr1 <- merge_dataset(
        AD_GWAS_female, AD_GWAS_male,
        pQTL_female_filtered1, pQTL_male_filtered1,
        pQTL_combined_filtered1,
        best_snp1[1], best_snp1[2], range
    )
    data_list1  <- mr1$merged_data_list
    common_id1  <- mr1$common_id
    best_snp_actual1 <- find_best_snp(data_list1, signal1, chr, bp_start, bp_end)
    LD_results1 <- calculate_LD(common_id1, best_snp_actual1[3], Gene_name2)
    merged_data_list_with_LD1 <- read_LD(
        data_list1,
        best_snp_actual1[3],
        LD_results1[1],
        LD_results1[2]
    )
    y_lim_AD1   <- find_AD_ylim(merged_data_list_with_LD1[["AD_female"]],
                                merged_data_list_with_LD1[["AD_male"]],
                                best_snp_actual1[3], best_snp_actual1[2])
    y_lim_pQTL1 <- find_pQTL_ylim(merged_data_list_with_LD1[["pQTL_female"]],
                                merged_data_list_with_LD1[["pQTL_male"]],
                                merged_data_list_with_LD1[["pQTL_combined"]],
                                best_snp_actual1[3], best_snp_actual1[2])

    pQTL_female_filtered2   <- subset_pQTL(pQTL_female, Gene_id)
    pQTL_male_filtered2     <- subset_pQTL(pQTL_male, Gene_id)
    pQTL_combined_filtered2 <- subset_pQTL(pQTL_combined, Gene_id)

    best_snp2 <- find_best_snp(data_list, signal2, chr, bp_start, bp_end)

    mr2         <- merge_dataset(
                    AD_GWAS_female, AD_GWAS_male,
                    pQTL_female_filtered2, pQTL_male_filtered2,
                    pQTL_combined_filtered2,
                    best_snp2[1], best_snp2[2], range
                )
    data_list2  <- mr2$merged_data_list
    common_id2  <- mr2$common_id

    best_snp_actual2 <- find_best_snp(data_list2, signal2, chr, bp_start, bp_end)

    LD_results2 <- calculate_LD(common_id2, best_snp_actual2[3], Gene_name2)
    merged_data_list_with_LD2 <- read_LD(
    data_list2,
    best_snp_actual2[3],
    LD_results2[1],
    LD_results2[2]
    )
    y_lim_AD2   <- find_AD_ylim(merged_data_list_with_LD2[["AD_female"]],
                                merged_data_list_with_LD2[["AD_male"]],
                                best_snp_actual2[3], best_snp_actual2[2])
    y_lim_pQTL2 <- find_pQTL_ylim(merged_data_list_with_LD2[["pQTL_female"]],
                                merged_data_list_with_LD2[["pQTL_male"]],
                                merged_data_list_with_LD2[["pQTL_combined"]],
                                best_snp_actual2[3], best_snp_actual2[2])
    if (level == "Primary") {
        if (signal2 == "pQTL_female") {
            pqtl_subtitle <- "Top pQTL female Signal"
        } else if (signal2 == "pQTL_male") {
            pqtl_subtitle <- "Top pQTL male Signal"
        }
    } else {
        pqtl_subtitle <- "Top pQTL combined Signal"
    }

    if(signal1 == "AD_female"){
        ad_subtitle <- "Top AD female Signal"
    }else if(signal1 == "AD_male"){
        ad_subtitle <- "Top AD male Signal"
    }
                 
    p_adf1 <- locus_sub_plot_38(merged_data_list_with_LD1[["AD_female"]],
                                best_snp_actual1[3], best_snp_actual1[2],
                                bp_start, bp_end, y_lim_AD1, "AD female ",plot_title = ad_subtitle)
    p_adf2 <- locus_sub_plot_38(merged_data_list_with_LD2[["AD_female"]],
                                best_snp_actual2[3], best_snp_actual2[2],
                                bp_start, bp_end, y_lim_AD2, "AD female ",plot_title = pqtl_subtitle)
    # AD GWAS male
    p_adm1 <- locus_sub_plot_38(merged_data_list_with_LD1[["AD_male"]],
                                best_snp_actual1[3], best_snp_actual1[2],
                                bp_start, bp_end, y_lim_AD1, "AD male ")
    p_adm2 <- locus_sub_plot_38(merged_data_list_with_LD2[["AD_male"]],
                                best_snp_actual2[3], best_snp_actual2[2],
                                bp_start, bp_end, y_lim_AD2, "AD male ")
    # pQTL female
    p_pqf1 <- locus_sub_plot_38(merged_data_list_with_LD1[["pQTL_female"]],
                                best_snp_actual1[3], best_snp_actual1[2],
                                bp_start, bp_end, y_lim_pQTL1, "pQTL female ")
    p_pqf2 <- locus_sub_plot_38(merged_data_list_with_LD2[["pQTL_female"]],
                                best_snp_actual2[3], best_snp_actual2[2],
                                bp_start, bp_end, y_lim_pQTL2, "pQTL female ")
    # pQTL male
    p_pqm1 <- locus_sub_plot_38(merged_data_list_with_LD1[["pQTL_male"]],
                                best_snp_actual1[3], best_snp_actual1[2],
                                bp_start, bp_end, y_lim_pQTL1, "pQTL male ")
    p_pqm2 <- locus_sub_plot_38(merged_data_list_with_LD2[["pQTL_male"]],
                                best_snp_actual2[3], best_snp_actual2[2],
                                bp_start, bp_end, y_lim_pQTL2, "pQTL male ")
    # pQTL combined
    p_pqc1 <- locus_sub_plot_38(merged_data_list_with_LD1[["pQTL_combined"]],
                                best_snp_actual1[3], best_snp_actual1[2],
                                bp_start, bp_end, y_lim_pQTL1, "pQTL combined ")
    p_pqc2 <- locus_sub_plot_38(merged_data_list_with_LD2[["pQTL_combined"]],
                                best_snp_actual2[3], best_snp_actual2[2],
                                bp_start, bp_end, y_lim_pQTL2, "pQTL combined ")

    p_gene1 <- genetracks_sub_plot(
    loc = {
        tmp <- locus(data = merged_data_list_with_LD1[["AD_female"]],
                    chrom = "CHR", pos = "BP", p = "P",
                    index_snp = best_snp_actual1[3], LD = "r2",
                    flank = c(as.numeric(best_snp_actual1[2]) - bp_start,
                                bp_end - as.numeric(best_snp_actual1[2])),
                    ens_db = "EnsDb.Hsapiens.v86")
        link_recomb(tmp, genome = "hg38")
    },
    Gene_name1, Gene_name2,maxrows_num
    )

    p_gene2 <- genetracks_sub_plot(
    loc = {
        tmp <- locus(data = merged_data_list_with_LD2[["AD_female"]],
                    chrom = "CHR", pos = "BP", p = "P",
                    index_snp = best_snp_actual2[3], LD = "r2",
                    flank = c(as.numeric(best_snp_actual2[2]) - bp_start,
                                bp_end - as.numeric(best_snp_actual2[2])),
                    ens_db = "EnsDb.Hsapiens.v86")
        link_recomb(tmp, genome = "hg38")
    },
    Gene_name1, Gene_name2,maxrows_num
    )
    p_adf1 <- p_adf1 + theme(plot.margin = margin(t = 5, r = 8.5, b = 5, l = 5))
    p_adf2 <- p_adf2 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 8.5))
    p_adm1 <- p_adm1 + theme(plot.margin = margin(t = 5, r = 8.5, b = 5, l = 5))
    p_adm2 <- p_adm2 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 8.5))
    p_gene1 <- p_gene1 + theme(plot.margin = margin(t = 5, r = 40, b = 5, l = 5))
    p_gene2 <- p_gene2 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 40))
    level <- as.character(level)
    if (level == "Primary") {

        if (signal2 == "pQTL_female") {
            pQTL_index <- "pQTL_female"
            pQTL_label <- "pQTL female "
        } else if (signal2 == "pQTL_male") {
            pQTL_index <- "pQTL_male"
            pQTL_label <- "pQTL male "
        }
    } else {
        pQTL_index <- "pQTL_combined"
        pQTL_label <- "pQTL combined "
    }
        
    if (signal1 == "AD_female") {
        AD_index_sc <- "AD_female"
        AD_label <- "AD female "
    } else if (signal1 == "AD_male") {
        AD_index_sc <- "AD_male"
        AD_label <- "AD male "
    }
    

    p_sc1 <- locus_compare_plot_38(merged_data_list_with_LD1[[AD_index_sc]],
                                    merged_data_list_with_LD1[[pQTL_index]],
                                    chr, best_snp_actual1[3],
                                    AD_label, pQTL_label,position)
    p_sc2 <- locus_compare_plot_38(merged_data_list_with_LD2[[AD_index_sc]],
                                    merged_data_list_with_LD2[[pQTL_index]],
                                    chr, best_snp_actual2[3],
                                    AD_label, pQTL_label,position)
    p_sc1 <- p_sc1 + theme(plot.margin = margin(t = 5, r = 30, b = 5, l = 5))
    p_sc2 <- p_sc2 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30))                        

    final_plot <- 
        (p_adf1 | p_adf2) /
        (p_adm1 | p_adm2) /
        (p_pqf1 | p_pqf2) /
        (p_pqm1 | p_pqm2) /
        (p_pqc1 | p_pqc2) /
        (p_gene1 | p_gene2) /

        (p_sc1 | p_sc2) +

        plot_layout(
            nrow = 7,
            heights = c(1, 1, 1, 1, 1, 0.5, 2)
        )

    final_plot_with_margin <- ggdraw(final_plot) +
    theme(plot.margin = margin(0, 0.5, 1, 0.5, unit = "inches"))
    index <- as.character(index)
    title_S14 <- ggdraw() + 
    draw_label(paste0("S14.", index), 
                size = 24, 
                fontface = "bold",
                x = 0.02, 
                hjust = 0, 
                vjust = 1
    ) +
    theme(plot.margin = margin(1, 0, 0, 0.5, unit = "inches"))  
    if(signal1 == "AD_female"){
        title_index <- " Female "
    }else if(signal1 == "AD_male"){
        title_index <- " Male "
    }
    title_global <- ggdraw() + 
    draw_label(
        paste0(Gene_name2,title_index,level," Discovery"),
        size = 18,  
        fontface = 'bold', 
        x = 0.5, 
        hjust = 0.5,
        vjust = 0,
        lineheight = 0.8  
    ) +
    theme(plot.margin = margin(0.5, 0, 0, 0, unit = "inches")) 


    final_plot_with_title <- plot_grid(
    title_S14,
    title_global,
    final_plot_with_margin,
    ncol = 1,
    rel_heights = c(0.04, 0.03, 0.93)  
    )

    ggsave(
    filename = paste0(out_dir,index,"_",Gene_name2, "_s14_supplement_plot.pdf"),
    plot = final_plot_with_title,
    width = 8.27 * 2, 
    height = 11.69 * 2
    )
    print(paste0("Saved plot to ", paste0(out_dir,index,"_",Gene_name2, "_s14_supplement_plot.pdf")))

    cleanup_LD(best_snp_actual1[3], Gene_name2)
    cleanup_LD(best_snp_actual2[3], Gene_name2)

    return(final_plot_with_title)
}


read_in_csf <- function(Gene_id){
    csf_female = fread(paste0("/Path/to/Female/CSF/pQTL/files/separated/by/Gene_id.glm.linear.gz"))

    csf_male = fread(paste0("/Path/to/Male/CSF/pQTL/files/separated/by/Gene_id.glm.linear.gz"))

    csf_combined = fread(paste0("/Path/to/Combined/CSF/pQTL/files/separated/by/Gene_id.glm.linear.gz"))

    csf_female <- csf_female %>%
    dplyr::rename(CHR = `#CHROM`, BP = `POS`, A1 = `A1`, A2 = `REF`) %>%
    dplyr::mutate(
        id12 = paste(CHR, BP, A1, A2, sep = ":"),
        id21 = paste(CHR, BP, A2, A1, sep = ":")) %>%
    dplyr::mutate(
        P = 10^(-LOG10_P)
    ) %>%
    dplyr::mutate(
        gene_id = Gene_id
    )

    csf_male <- csf_male %>%
    dplyr::rename(CHR = `#CHROM`, BP = `POS`, A1 = `A1`, A2 = `REF`) %>%
    dplyr::mutate(
        id12 = paste(CHR, BP, A1, A2, sep = ":"),
        id21 = paste(CHR, BP, A2, A1, sep = ":")) %>%
    dplyr::mutate(
        P = 10^(-LOG10_P)
    ) %>%
    dplyr::mutate(
        gene_id = Gene_id
    )

    csf_combined <- csf_combined %>%
    dplyr::rename(CHR = `#CHROM`, BP = `POS`, A1 = `A1`, A2 = `REF`) %>%
    dplyr::mutate(
        id12 = paste(CHR, BP, A1, A2, sep = ":"),
        id21 = paste(CHR, BP, A2, A1, sep = ":")) %>%
    dplyr::mutate(
        P = 10^(-LOG10_P)
    ) %>%
    dplyr::mutate(
        gene_id = Gene_id
    )

    data_list <- list()
    data_list[["csf_female"]] <- csf_female
    data_list[["csf_male"]] <- csf_male
    data_list[["csf_combined"]] <- csf_combined
    return(data_list)
}
