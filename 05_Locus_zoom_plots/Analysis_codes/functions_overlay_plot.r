find_max_ylim_AD <- function(data, snps) {
  loc <- locus(data = data, chrom = "CHR", pos = "BP", p = "P", index_snp = snps, fix_window = 1e6, ens_db = "EnsDb.Hsapiens.v86")
  return(ceiling(max(-log10(as.numeric(loc$data$P)), na.rm = TRUE)))
}

merge_ld_values <- function(df, ld1, ld2, ld3) {
  merged_ld1 <- inner_join(df, ld1, by = "SNP") %>%
    mutate(r2_1 = r2) %>%
    as.data.frame()

  merged_ld2 <- inner_join(df, ld2, by = "SNP") %>%
    mutate(r2_2 = r2) %>%
    as.data.frame()

  merged_ld3 <- inner_join(df, ld3, by = "SNP") %>%
    mutate(r2_3 = r2) %>%
    as.data.frame()

    merged_withld_123 <- merge(
    merged_ld1,
    merged_ld2[c("SNP", "r2_2")],
    by = "SNP"
    )

    merged_withld_123 <- merge(
    merged_withld_123,
    merged_ld3[c("SNP", "r2_3")],
    by = "SNP"
    )

  return(merged_withld_123)
}

select_data_merge_LD <- function(rsid,rsid_2,rsid_3,chr,bp,df_f,df_m,df_1,df_2,df_3,df_4,df_5,fold_w,file_snp,in_file,pl,in_fold,keep_file,window_range) {
    print("start")
    df_f_1mb <- df_f %>% 
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range)  %>% 
        dplyr::mutate(id12 = paste(CHR, BP, ALLELE1, ALLELE0, sep = ":")) %>%
        dplyr::mutate(id21 = paste(CHR, BP, ALLELE0, ALLELE1, sep = ":"))

    df_m_1mb <- df_m %>% 
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range) %>%
        dplyr::mutate(id12 = paste(CHR, BP, ALLELE1, ALLELE0, sep = ":")) %>%
        dplyr::mutate(id21 = paste(CHR, BP, ALLELE0, ALLELE1, sep = ":"))

    df_1_1mb <- df_1 %>% 
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range) %>%
        dplyr::mutate(id = paste(CHR, BP, A1, A2, sep = ":"))

    df_2_1mb <- df_2 %>% 
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range) %>%
        dplyr::mutate(id = paste(CHR, BP, A1, A2, sep = ":"))

    df_3_1mb <- df_3 %>%
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range) %>%
        dplyr::mutate(id = paste(CHR, BP, A1, A2, sep = ":"))

    df_4_1mb <- df_4 %>%
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range) %>%
        dplyr::mutate(id = paste(CHR, BP, A1, A2, sep = ":"))

    df_5_1mb <- df_5 %>%
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window_range & BP <= bp + window_range) %>%
        dplyr::mutate(id = paste(CHR, BP, A1, A2, sep = ":"))

    print("df cut done")

    merged_df_1_1mb_12 <- merge(df_1_1mb,df_f_1mb[,c("id12","SNP")],by.x = "id",by.y = "id12")
    merged_df_1_1mb_21 <- merge(df_1_1mb,df_m_1mb[,c("id21","SNP")],by.x = "id",by.y = "id21")
    merged_df_1_1mb <- rbind(merged_df_1_1mb_12,merged_df_1_1mb_21)

    merged_df_2_1mb_12 <- merge(df_2_1mb, df_f_1mb[, c("id12", "SNP")], by.x = "id", by.y = "id12")
    merged_df_2_1mb_21 <- merge(df_2_1mb, df_m_1mb[, c("id21", "SNP")], by.x = "id", by.y = "id21")
    merged_df_2_1mb <- rbind(merged_df_2_1mb_12, merged_df_2_1mb_21)

    merged_df_3_1mb_12 <- merge(df_3_1mb, df_f_1mb[, c("id12", "SNP")], by.x = "id", by.y = "id12")
    merged_df_3_1mb_21 <- merge(df_3_1mb, df_m_1mb[, c("id21", "SNP")], by.x = "id", by.y = "id21")
    merged_df_3_1mb <- rbind(merged_df_3_1mb_12, merged_df_3_1mb_21)

    merged_df_4_1mb_12 <- merge(df_4_1mb, df_f_1mb[, c("id12", "SNP")], by.x = "id", by.y = "id12")
    merged_df_4_1mb_21 <- merge(df_4_1mb, df_m_1mb[, c("id21", "SNP")], by.x = "id", by.y = "id21")
    merged_df_4_1mb <- rbind(merged_df_4_1mb_12, merged_df_4_1mb_21)

    merged_df_5_1mb_12 <- merge(df_5_1mb, df_f_1mb[, c("id12", "SNP")], by.x = "id", by.y = "id12")
    merged_df_5_1mb_21 <- merge(df_5_1mb, df_m_1mb[, c("id21", "SNP")], by.x = "id", by.y = "id21")
    merged_df_5_1mb <- rbind(merged_df_5_1mb_12, merged_df_5_1mb_21)

    common_id <- intersect(df_f_1mb$SNP,df_m_1mb$SNP)
    merged_df_f_1mb <- df_f_1mb %>% dplyr::filter(SNP %in% common_id)
    merged_df_m_1mb <- df_m_1mb %>% dplyr::filter(SNP %in% common_id)
    print("df merge done")

    fwrite(data.frame(common_id), paste0(fold_w, rsid, file_snp), sep = "\t", col.names = FALSE)

    print("write snps list done")
    in_file_chr <- paste0(in_file,chr)
    
    command <- paste(pl, 
                 "--bfile", paste0(in_fold, in_file_chr),
                 "--allow-no-sex",
                 "--extract", paste0(fold_w, rsid, file_snp),
                 "--keep-allele-order",
                 "--make-bed --out", paste0(fold_w, paste0(rsid, ".data")),
                 sep = " ")
    system(command)

    command <- paste(pl, 
                 "--bfile", paste0(fold_w, paste0(rsid, ".data")),
                 "--allow-no-sex",
                 "--keep-allele-order",
                 "--keep", keep_file,
                 "--r2 yes-really --matrix",
                 "--out", paste0(fold_w, rsid),
                 sep = " ")
    system(command) 

    # Load the BIM file and corresponding LD matrix
    bim_file <- paste0(fold_w, paste0(rsid, ".data.bim"))
    if (file.exists(bim_file)) {
    bim <- fread(bim_file)
    snps_name <- unlist(bim$V2)  # Extract the SNP names
    
    LD_file <- paste0(fold_w, rsid, ".ld")
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

    print("LD done")

    # Ensure correct selection syntax
    ld_values <- LD_df %>%
    dplyr::select(V1, all_of(rsid)) %>%
    dplyr::rename(SNP = V1, r2 = rsid) %>%
    arrange(desc(SNP))

    ld_values_2 <- LD_df %>%
    dplyr::select(V1, all_of(rsid_2)) %>%
    dplyr::rename(SNP = V1, r2 = rsid_2) %>%
    arrange(desc(SNP))

    ld_values_3 <- LD_df %>%
    dplyr::select(V1, all_of(rsid_3)) %>%
    dplyr::rename(SNP = V1, r2 = rsid_3) %>%
    arrange(desc(SNP))

    print("ld values select done")

    merged_female_withld <- merge_ld_values(merged_df_f_1mb,ld_values,ld_values_2,ld_values_3)
    print("ld values merge done for ad female")
    merged_male_withld <- merge_ld_values(merged_df_m_1mb,ld_values,ld_values_2,ld_values_3)
    print("ld values merge done for ad male")
    merged_1_withld <- merge_ld_values(merged_df_1_1mb,ld_values,ld_values_2,ld_values_3)
    merged_2_withld <- merge_ld_values(merged_df_2_1mb,ld_values,ld_values_2,ld_values_3)
    merged_3_withld <- merge_ld_values(merged_df_3_1mb,ld_values,ld_values_2,ld_values_3)
    merged_4_withld <- merge_ld_values(merged_df_4_1mb,ld_values,ld_values_2,ld_values_3)
    merged_5_withld <- merge_ld_values(merged_df_5_1mb,ld_values,ld_values_2,ld_values_3)

    list <- list()
    list[["female"]] <- merged_female_withld
    list[["male"]] <- merged_male_withld
    list[["1"]] <- merged_1_withld
    list[["2"]] <- merged_2_withld
    list[["3"]] <- merged_3_withld
    list[["4"]] <- merged_4_withld
    list[["5"]] <- merged_5_withld

    return(list)
}

get_diff_groups <- function(df, col1, col2) {
  df <- df %>% dplyr::mutate(diff = .data[[col1]] - .data[[col2]])
  list(
    "diff_0.2"   = df %>% dplyr::filter(abs(diff) <= 0.2),
    "diff_0.4_1" = df %>% dplyr::filter(diff > 0.2 & diff <= 0.4),
    "diff_0.4_2" = df %>% dplyr::filter(diff < -0.2 & diff >= -0.4),
    "diff_0.6_1" = df %>% dplyr::filter(diff > 0.4 & diff <= 0.6),
    "diff_0.6_2" = df %>% dplyr::filter(diff < -0.4 & diff >= -0.6),
    "diff_0.8_1" = df %>% dplyr::filter(diff > 0.6 & diff <= 0.8),
    "diff_0.8_2" = df %>% dplyr::filter(diff < -0.6 & diff >= -0.8)
  )
}

add_diff_layers <- function(p, diff_groups, color_mapping) {
  p + purrr::imap(diff_groups, ~ {
    if (nrow(.x) > 0) {
      geom_point(
        data = .x,
        aes(x = BP / 1e6, y = logP, size = max_r2 * 3),
        fill = color_mapping[.y],
        shape = 21,
        color = "white",
        alpha = 1,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    }
  }) +
    theme(legend.position = "none")
}

get_ggplot <- function(max_ylim_ad,bp,merged_female_withld,merged_df,rsid,rsid_2,rsid_3,y_lable,margin_y){
    index_color_list <- c("#ff00fa", "#ffc900", "#00ffdb")
    LD_color_list <- c(
    "#d5b14e","#bfb560","#aaba73","#95be85","#80c397","#6ac7a9","#55ccbc",
        "#d1ae55","#b9b369","#a2b87d","#8bbd91","#74c2a6","#5dc7b9",

        "#d5a063","#bfa475","#aaa988","#95ad9a","#80b2ac","#6ab6be",
            "#d19c6c","#b9a180","#a2a694","#8baba8","#74b0bc",

            "#d58f78","#bf948a","#aa989c","#959daf","#80a1c1",
                "#d18a83","#b98f97","#a294ab","#8b98be",

                "#d57e8d","#bf839f","#aa87b1","#958cc3",
                    "#d17799","#b97cad","#a281c1",

                    "#d56ea2","#bf72b4","#aa77c6",
                        "#c479a6","#b86ac4",

                        "#d55db6","#bf61c9",
                            "#d153c7",
                            
                            "#d54ccd"
    )

    bp_start <- min(merged_female_withld$BP)
    bp_end <- max(merged_female_withld$BP)

    max_ylim <- ceiling(max(-log10(as.numeric(merged_df$P)), na.rm = TRUE))
    min_ylim <- floor(min(-log10(as.numeric(merged_df$P)), na.rm = TRUE))

    # Generate locus plot
    loc1 <- locus(data = merged_df, chrom = "CHR", pos = "BP", p = "P", index_snp = rsid, LD = "r2_1", flank = c(bp - bp_start, bp_end - bp), ens_db = "EnsDb.Hsapiens.v86")
    loc1 <- link_recomb(loc1, genome = "hg38")

    if (y_lable == "AD Female GWAS" | y_lable == "AD Male GWAS"){
        max_ylim <- max_ylim_ad
        p1b <- gg_scatter(loc1, pcutoff = 5e-08, ylim = c(min_ylim, max_ylim), labels = NULL, xticks = FALSE, cex.lab = 1, ylab = y_lable,legend_pos=NULL,recomb_col = "#2d2d2d",
                  LD_scheme = c("white", "white", "white", "white", "white", "white", "white")) 
    }else{
        p1b <- gg_scatter(loc1, pcutoff = 5e-1000, ylim = c(min_ylim, max_ylim), labels = NULL, xticks = FALSE, cex.lab = 1, ylab = "", legend_pos = NULL, recomb_col = "#2d2d2d",
                LD_scheme = c("white", "white", "white", "white", "white", "white", "white")) + 
            theme(
                axis.text.y = element_text(color = "white"),
                axis.ticks.y = element_line(color = "black"),
                axis.ticks.y.right = element_line(color = "black"),
                axis.text.y.right = element_text(color = "black")
            )
    }

    p1b <- p1b + 
        theme(
            axis.title.y = element_text(
            size = 7, 
            face = "bold", 
            margin = margin(r = margin_y)
            ),
            axis.title.y.right = element_text(size = 7, face = "plain"),
            plot.title = element_text(
            size = 5, 
            face = "bold",
            hjust = 0.5, 
            margin = margin(b = 10)
        )
    )

    p1b_fixed <- p1b
    p1b_fixed$layers <- p1b_fixed$layers[!sapply(p1b_fixed$layers, function(layer) inherits(layer$geom, "GeomPoint"))]

    loc1_0_subset <- loc1$data[loc1$data$ld < 0.2,]
    loc2_0_subset <- loc1$data[loc1$data$r2_2 < 0.2,]
    loc3_0_subset <- loc1$data[loc1$data$r2_3 < 0.2,]
    snp_loc1 <- loc1_0_subset$SNP
    snp_loc2 <- loc2_0_subset$SNP
    snp_loc3 <- loc3_0_subset$SNP

    common_snps <- Reduce(intersect, list(snp_loc1, snp_loc2, snp_loc3))
    loc_0_subset <- loc1_0_subset[loc1_0_subset$SNP %in% common_snps,]

    p1b_fixed <- p1b_fixed +
        geom_point(
            data = subset(loc1$data, SNP %in% loc_0_subset$SNP),
            aes(x = BP / 1e6, y = logP), 
            shape = 21, 
            fill = "white", 
            color = "#3c3c3c", 
            size = 1, 
            alpha = 0.25,
            inherit.aes = FALSE
        )
    
    loc1$data$max_r2 <- pmax(
        ifelse(loc1$data$ld   == 1, NA, loc1$data$ld),
        ifelse(loc1$data$r2_2 == 1, NA, loc1$data$r2_2),
        ifelse(loc1$data$r2_3 == 1, NA, loc1$data$r2_3),
        na.rm = TRUE
        )
    loc1_filtered <- loc1$data %>%
    dplyr::filter(
        ld >= 0.2 & ld <= 1,
        r2_2 >= 0.2 & r2_2 <= 1,
        r2_3 < 0.2
    )
    loc1_diff_groups <- get_diff_groups(loc1_filtered, "ld", "r2_2")
    color_mapping1 <- c(
    "diff_0.2"   = "#958cc3",
    "diff_0.4_1" = "#80a1c1",
    "diff_0.4_2" = "#aa77c6",
    "diff_0.6_1" = "#6ab6be",
    "diff_0.6_2" = "#bf61c9",
    "diff_0.8_1" = "#55ccbc",
    "diff_0.8_2" = "#d54ccb"
    )
    p1b_fixed <- add_diff_layers(p1b_fixed, loc1_diff_groups, color_mapping1)

    loc2_filtered <- loc1$data %>%
    dplyr::filter(
        ld < 0.2,
        r2_2 >= 0.2 & r2_2 <= 1,
        r2_3 >= 0.2 & r2_3 <= 1
    )
    loc2_diff_groups <- get_diff_groups(loc2_filtered, "r2_2", "r2_3")
    color_mapping2 <- c(
    "diff_0.2"   = "#d57e8d",
    "diff_0.4_1" = "#d56ea2",
    "diff_0.4_2" = "#d58f78",
    "diff_0.6_1" = "#d55db6",
    "diff_0.6_2" = "#d5a063",
    "diff_0.8_1" = "#d54ccd",
    "diff_0.8_2" = "#d5b14e"
    )
    p1b_fixed <- add_diff_layers(p1b_fixed, loc2_diff_groups, color_mapping2)

    loc3_filtered <- loc1$data %>%
    dplyr::filter(
        ld >= 0.2 & ld <= 1,
        r2_2 < 0.2,
        r2_3 >= 0.2 & r2_3 <= 1
    )
    loc3_diff_groups <- get_diff_groups(loc3_filtered, "ld", "r2_3")
    color_mapping3 <- c(
    "diff_0.2"   = "#95be85",
    "diff_0.4_1" = "#80c397",
    "diff_0.4_2" = "#aaba73",
    "diff_0.6_1" = "#6ac7a9",
    "diff_0.6_2" = "#bfb560",
    "diff_0.8_1" = "#55ccbc",
    "diff_0.8_2" = "#d5b14e"
    )
    p1b_fixed <- add_diff_layers(p1b_fixed, loc3_diff_groups, color_mapping3)

    loc_all_filtered <- loc1$data %>%
    dplyr::filter(
        ld >= 0.2 & ld < 1,
        r2_2 >= 0.2 & r2_2 <1,
        r2_3 >= 0.2 & r2_3 < 1
    ) %>% distinct(SNP, .keep_all = TRUE)

    # print(loc_all_filtered)

    # stop("stop here")

    snp_color_mapping <- c(
        "rs10444736" = "#b98f97",
        "rs112213845" = "#b98f97",
        "rs12147398" = "#b98f97",
        "rs17127085" = "#b98f97",
        "rs17182391" = "#b98f97",
        "rs17182412" = "#95ad9a",
        "rs56175115" = "#95ad9a",
        "rs56917431" = "#b98f97",
        "rs59923225" = "#b98f97",
        "rs61985808" = "#b98f97",
        "rs61985848" = "#b98f97",
        "rs61987082" = "#b98f97",
        "rs7145805" = "#b98f97",
        "rs7157889" = "#95ad9a",
        "rs7158675" = "#95ad9a",
        "rs8008784" = "#95ad9a",
        "rs8009020" = "#b98f97",
        "rs8016448" = "#b98f97",
        "rs8018050" = "#b98f97"
        )

    for (i in 1:nrow(loc_all_filtered)) {
        temp_LD3 <- loc_all_filtered[i, ]
        color_LD3 <- snp_color_mapping[temp_LD3$SNP]
        p1b_fixed <- p1b_fixed +
            geom_point(
                data = subset(loc1$data, SNP %in% temp_LD3$SNP),
                aes(x = BP / 1e6, y = logP, size = max_r2 * 3),
                shape = 21,
                fill = color_LD3,
                color = "white",
                alpha = 1,
                inherit.aes = FALSE,
                show.legend = FALSE
            )

    }

    loc1_subset <- loc1$data[loc1$data$ld >= 0.2 & loc1$data$ld < 1 & loc1$data$r2_2 < 0.2 & loc1$data$r2_3 < 0.2, ]
    loc2_subset <- loc1$data[loc1$data$r2_2 >= 0.2 & loc1$data$r2_2 < 1 & loc1$data$ld < 0.2 & loc1$data$r2_3 < 0.2, ]
    loc3_subset <- loc1$data[loc1$data$r2_3 >= 0.2 & loc1$data$r2_3 < 1 & loc1$data$ld < 0.2 & loc1$data$r2_2 < 0.2, ]

    p1b_fixed <- p1b_fixed +
    geom_point(
        data = subset(loc1$data, SNP %in% loc1_subset$SNP),
        aes(x = BP / 1e6, y = logP, size = r2 * 3),
        shape = 21,
        fill = "#00ffdb",
        color = "white",
        alpha = 1,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    geom_point(
        data = subset(loc1$data, SNP %in% loc2_subset$SNP),
        aes(x = BP / 1e6, y = logP, size = r2_2 * 3),
        shape = 21,
        fill = "#ff00fa",
        color = "white",
        alpha = 1,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    geom_point(
        data = subset(loc1$data, SNP %in% loc3_subset$SNP),
        aes(x = BP / 1e6, y = logP, size = r2_3 * 3),
        shape = 21,
        fill = "#ffc900",
        color = "white",
        alpha = 1,
        inherit.aes = FALSE,
        show.legend = FALSE
    ) 

    p1b_fixed <- p1b_fixed  + geom_point(data = subset(loc1$data, SNP == rsid), 
                                        mapping = aes(x = BP / 1e6, y = logP), 
                                        shape = 23,       
                                        fill = "#00ffdb",  
                                        color = "#000000",  
                                        alpha = 0.75,
                                        size = 4, 
                                        inherit.aes = FALSE) +
    geom_point(data = subset(loc1$data, SNP == rsid_2), 
                mapping = aes(x = BP / 1e6, y = logP), 
                shape = 23,       
                fill = "#ff00fa",  
                color = "#000000",  
                alpha = 0.75,
                size = 4,      
                inherit.aes = FALSE)+  
    geom_point(data = subset(loc1$data, SNP == rsid_3), 
                mapping = aes(x = BP / 1e6, y = logP), 
                shape = 23,       
                fill = "#ffc900",  
                color = "#000000",  
                alpha = 0.75,
                size = 4,      
                inherit.aes = FALSE) 

    if (rsid %in% loc1_filtered$SNP) {
        loc1_filtered_rsid <- loc1_filtered %>% dplyr::filter(SNP == rsid)
        loc1_diff_groups_rsid <- get_diff_groups(loc1_filtered_rsid, "ld", "r2_2")
        p1b_fixed <- add_diff_layers(p1b_fixed, loc1_diff_groups_rsid, color_mapping1)
    }else if (rsid %in% loc3_filtered$SNP) {
        loc3_filtered_rsid <- loc3_filtered %>% dplyr::filter(SNP == rsid)
        loc3_diff_groups_rsid <- get_diff_groups(loc3_filtered_rsid, "ld", "r2_3")
        p1b_fixed <- add_diff_layers(p1b_fixed, loc3_diff_groups_rsid, color_mapping3)
    }

    if (rsid_2 %in% loc2_filtered$SNP) {
        loc2_filtered_rsid_2 <- loc2_filtered %>% dplyr::filter(SNP == rsid_2)
        loc2_diff_groups_rsid_2 <- get_diff_groups(loc2_filtered_rsid_2, "r2_2", "r2_3")
        p1b_fixed <- add_diff_layers(p1b_fixed, loc2_diff_groups_rsid_2, color_mapping2)
    }else if (rsid_2 %in% loc1_filtered$SNP) {
        loc1_filtered_rsid_2 <- loc1_filtered %>% dplyr::filter(SNP == rsid_2)
        loc1_diff_groups_rsid_2 <- get_diff_groups(loc1_filtered_rsid_2, "ld", "r2_2")
        p1b_fixed <- add_diff_layers(p1b_fixed, loc1_diff_groups_rsid_2, color_mapping1)
    }

    if (rsid_3 %in% loc2_filtered$SNP) {
        loc2_filtered_rsid_3 <- loc2_filtered %>% dplyr::filter(SNP == rsid_3)
        loc2_diff_groups_rsid_3 <- get_diff_groups(loc2_filtered_rsid_3, "r2_2", "r2_3")
        p1b_fixed <- add_diff_layers(p1b_fixed, loc2_diff_groups_rsid_3, color_mapping2)
    }else if (rsid_3 %in% loc3_filtered$SNP) {
        loc3_filtered_rsid_3 <- loc3_filtered %>% dplyr::filter(SNP == rsid_3)
        loc3_diff_groups_rsid_3 <- get_diff_groups(loc3_filtered_rsid_3, "ld", "r2_3")
        p1b_fixed <- add_diff_layers(p1b_fixed, loc3_diff_groups_rsid_3, color_mapping3)
    }
    p1b_fixed <- p1b_fixed + scale_size_continuous(range = c(1, 3))

    p1b_front <- p1b
    p1b_front$layers <- p1b_front$layers[!sapply(p1b_front$layers, function(layer) inherits(layer$geom, "GeomPoint"))]

    p1b_front <- p1b_front +
    theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
    )
    
    combined_plot <- ggdraw() +
    draw_plot(p1b_fixed, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(p1b_front, x = 0, y = 0, width = 1, height = 1)  

    ggsave(paste0("/storage1/fs1/belloy/Active/05_Projects/chenyu.y/overlay_plot/fig/sub1_",rsid,"_",y_lable,".png"), combined_plot, width = 6.6*1.5, height = 4*1.5, units = "cm", dpi = 600)
    print(paste0("save plot to ","/storage1/fs1/belloy/Active/05_Projects/chenyu.y/overlay_plot/fig/sub1_",rsid,"_",y_lable,".png"))
    return(combined_plot)
}

get_ggplot_compare <- function(merged_female_withld, merged_df, signal, index_snp_choose,ylab,final,margin_y,axis_numbers,axis_num_text) {
    
    if (signal == "left") {
        point_color <- "#ff00fa"
        merged_female_withld <- merged_female_withld[, !names(merged_female_withld) %in% c("r2","r2_1", "r2_3")]
        merged_female_withld$r2 <- merged_female_withld$r2_2
        merged_df <- merged_df[, !names(merged_df) %in% c("r2","r2_1", "r2_3")]
        merged_df$r2 <- merged_df$r2_2
    } else if (signal == "right") {
        point_color <- "#ffc900"
        merged_female_withld <- merged_female_withld[, !names(merged_female_withld) %in% c("r2","r2_1", "r2_2")]
        merged_female_withld$r2 <- merged_female_withld$r2_3
        merged_df <- merged_df[, !names(merged_df) %in% c("r2","r2_1", "r2_2")]
        merged_df$r2 <- merged_df$r2_3
    } else if (signal == "middle") {
        point_color <- "#00ffdb"
        merged_female_withld <- merged_female_withld[, !names(merged_female_withld) %in% c("r2","r2_2", "r2_3")]
        merged_female_withld$r2 <- merged_female_withld$r2_1
        merged_df <- merged_df[, !names(merged_df) %in% c("r2","r2_2", "r2_3")]
        merged_df$r2 <- merged_df$r2_1
    }
    
    merged_data <- dplyr::inner_join(merged_female_withld, merged_df, by = "SNP", 
                                     suffix = c("_female", "_other"))
    
    xlim <- ceiling(max(-log10(merged_data$P_female)))
    xlim_min <- floor(min(-log10(merged_data$P_female)))
    ylim <- ceiling(max(-log10(merged_data$P_other)))
    ylim_min <- min(-log10(merged_data$P_other))
    
    out_data <- merged_data %>% dplyr::filter(r2_female < 0.2)
    pplot <- ggplot(out_data, aes(x = -log10(P_female), y = -log10(P_other), size = 0.5)) +
        geom_point(alpha = 0.25, shape = 21, fill = "#ffffff", color = "#3c3c3c", show.legend = FALSE) +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA),
            axis.title.x = element_text(size = 7, face = "bold"),
            axis.title.y = element_text(size = 7, face = "bold"),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.text.x = element_text(size = 8, face = "bold"),
            axis.text.y = element_text(
                size = 8, 
                face = "bold", 
                margin = margin(r = margin_y, unit = "pt")
                )
        ) +
        labs(x = "", y = "")
    
    r2_data <- merged_data %>% dplyr::filter(r2_female >= 0.2)
    pplot <- pplot + 
        geom_point(data = r2_data, 
                   aes(x = -log10(P_female), y = -log10(P_other), size = r2_female * 3),
                   shape = 21, fill = point_color, color = "white", alpha = 1, show.legend = FALSE)
    
    index_data <- merged_data %>% dplyr::filter(SNP == index_snp_choose)
    pplot <- pplot + 
            geom_point(data = index_data, 
                    aes(x = -log10(P_female), y = -log10(P_other)),
                    size = 3.5, shape = 23, fill = point_color, color = "black", alpha = 1, show.legend = FALSE) +
        scale_size_continuous(range = c(0.7, 2.7))

    if(ylim < 5){step_axis_y <- 1}else if(ylim < 10){step_axis_y <- 2}else{step_axis_y <- 5}
    if(xlim < 5){step_axis_x <- 1}else if(xlim < 10){step_axis_x <- 2}else{step_axis_x <- 5}

    pplot <- pplot +
        coord_cartesian(xlim = c(xlim_min, xlim), ylim = c(ylim_min, ylim)) +
        scale_x_continuous(breaks = seq(xlim_min, xlim, by = step_axis_x)) +
        scale_y_continuous(breaks = axis_numbers, labels = axis_num_text)

    pplot <- pplot +
        theme(
            axis.ticks.x = element_line(color = "black"),
            axis.ticks.y = element_line(color = "black"),
        )

    ggsave(paste0("/storage1/fs1/belloy/Active/05_Projects/chenyu.y/overlay_plot/fig/sub1_compare_",rsid,"_", ylab, ".png"), 
        pplot, width = 4*1.5, height = 4*1.5, units = "cm", dpi = 600)
    print("save compare plot")
    return(pplot)
}



get_ggplot_genetrack <-function(maxrows_num,bp,merged_female_withld,merged_df,rsid,rsid_2,rsid_3,...){
    highlight_genes <- c(...)
    bp_start <- min(merged_female_withld$BP)
    bp_end <- max(merged_female_withld$BP)
    # Generate locus plot
    loc1 <- locus(data = merged_df, chrom = "CHR", pos = "BP", p = "P", index_snp = rsid, LD = "r2_1", flank = c(bp - bp_start, bp_end - bp), ens_db = "EnsDb.Hsapiens.v86")
    loc1 <- link_recomb(loc1, genome = "hg38")

    p1b <- gg_genetracks(loc1, 
             highlight = highlight_genes, 
             filter_gene_biotype = "protein_coding", 
             gene_col = "black", 
             exon_col = "black",
             exon_border = "black",
             showExons = TRUE, 
             cex.text = 0.5,
             xticks = TRUE,
             maxrows = maxrows_num)
    p1b <- p1b + theme(text = element_text(face = "bold"))
    ggsave(paste0("/storage1/fs1/belloy/Active/05_Projects/chenyu.y/overlay_plot/fig/",rsid,"_sub1_genetrack.png"), 
            p1b, width = 6.7*1.75, height = 4.2*1.75,units = "cm", dpi = 600)
    return(p1b)
}


read_in_csf <- function(Gene_id){
    csf_female = fread(paste0("/storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files",
                            "/other/pQTL-female-results/pqtl.",Gene_id,".glm.linear.gz"))

    csf_male = fread(paste0("/storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_sex-stratified_pQTL_files",
                            "/other/pQTL-male-results/pqtl.",Gene_id,".glm.linear.gz"))

    csf_combined = fread(paste0("/storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files",
                            "/Other_Files/02_no_na_no_palindromics/",Gene_id,".glm.linear.gz"))

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
