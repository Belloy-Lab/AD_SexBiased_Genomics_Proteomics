process_AD_data <- function(df,chr,bp,window){
    df_1mb <- df %>% 
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window & BP <= bp + window) %>%
        dplyr::mutate(id12 = paste(CHR, BP, ALLELE1, ALLELE0, sep = ":")) %>%
        dplyr::mutate(id21 = paste(CHR, BP, ALLELE0, ALLELE1, sep = ":"))

    return(df_1mb)
}


process_qtl_data <- function(df,chr,bp,window){
    df_1mb <- df %>% 
        dplyr::filter(CHR == chr) %>% 
        dplyr::filter(BP >= bp - window & BP <= bp + window) %>%
        dplyr::mutate(id = paste(CHR, BP, A1, A2, sep = ":"))  %>% 
        dplyr::mutate(P = pvalue)

    return(df_1mb)
}

merge_qtl <- function(merged_female_1mb,merged_male_1mb,qtl_df_1mb,rsid){
    merged_qtl_df_1mb_12 <- merge(qtl_df_1mb,merged_female_1mb[,c("id12","SNP")],by.x = "id",by.y = "id12")
    merged_qtl_df_1mb_21 <- merge(qtl_df_1mb,merged_male_1mb[,c("id21","SNP")],by.x = "id",by.y = "id21")
    merged_qtl_df_1mb <- rbind(merged_qtl_df_1mb_12,merged_qtl_df_1mb_21)


    common_id <- intersect(merged_female_1mb$SNP, merged_male_1mb$SNP)

    merged_qtl_df_1mb <- merged_qtl_df_1mb %>% dplyr::filter(SNP %in% common_id)

    fwrite(data.frame(common_id), paste0(fold_w, rsid, file_snp), sep = "\t", col.names = FALSE)

    return(merged_qtl_df_1mb)
}

calculate_ld <- function(fold_w, file_snp,rsid, pl,in_fold,in_file,chr){
    in_file_chr <- paste0(in_file,chr)
    keep_file <- "/storage1/fs1/belloy/Active/01_References/01_Files/EU_nHISP_TOPMed_list_Unrel.txt"
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
    file.remove(paste0(fold_w, paste0(rsid, ".ld")))
    file.remove(paste0(fold_w, paste0(rsid, ".log")))
    file.remove(paste0(fold_w, paste0(rsid, ".data.bed")))
    file.remove(paste0(fold_w, paste0(rsid, ".data.bim")))
    file.remove(paste0(fold_w, paste0(rsid, ".data.fam")))
    file.remove(paste0(fold_w, paste0(rsid, ".data.log")))
    file.remove(paste0(fold_w, paste0(rsid, ".nosex")))
    file.remove(paste0(fold_w, paste0(rsid, ".data.nosex")))
    file.remove(paste0(fold_w, paste0(rsid, "_common_id.txt")))

    return(LD_df)
}

get_snp_ld <- function(LD_df, rsid, df) {
  ld_values <- LD_df %>%
    dplyr::select(V1, all_of(rsid)) %>%       # safe selection of external var
    dplyr::rename(SNP = V1) %>%                # rename V1 → SNP
    dplyr::rename_with(~ "r2", all_of(rsid)) %>%  # rename the selected rsid column → r2
    dplyr::arrange(desc(SNP))
  
  df <- as.data.frame(df)
  df_withld <- dplyr::inner_join(df, ld_values, by = "SNP")
  return(df_withld)
}

subplot_scatter_top <- function(female_90_1mb_withld,bp,max_ylim,df,rsid,ylabel,top_title_suffix){
    bp_start <- min(female_90_1mb_withld$BP)
    bp_end <- max(female_90_1mb_withld$BP)

    loc1 <- locus(data = df, chrom = "CHR", pos = "BP", p = "P", index_snp = rsid, LD = "r2", 
                    flank = c(bp - bp_start, bp_end - bp), ens_db = "EnsDb.Hsapiens.v86")
    # loc1 <- link_recomb(loc1, genome = "hg38")
    bw_url <- "http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recomb1000GAvg.bw"
    tmp <- tempfile(fileext = ".bw")
    download.file(bw_url, tmp, mode = "wb")
    recomb_gr <- import.bw(tmp)
    loc1 <- link_recomb(loc1, recomb = recomb_gr)

    p <- gg_scatter(loc1, pcutoff = 5e-8, ylim = c(0, max_ylim), labels = NULL, xticks = FALSE, cex.lab = 1, ylab = ylabel,legend_pos=NULL) + 
        theme(
            axis.title.y = element_text(size = 8, face = "bold"),     
            axis.title.y.right = element_text(size = 8, face = "plain"),
            plot.title = element_text(
                size = 14, 
                face = "bold",
                hjust = 0.5, 
                margin = margin(b = 10)
            )
        ) + 
        geom_point(data = subset(loc1$data, SNP == loc1$index_snp), 
                mapping = aes(x = BP / 1e6, y = logP), 
                shape = 23,       
                fill = "purple",  
                color = "black",  
                size = 4,        
                inherit.aes = FALSE)+
        labs(title = paste0("Index snp: ", rsid,top_title_suffix))
    p <- p + 
        theme(
            legend.background = element_rect(fill = "transparent")
        )

    return(p)
}

subplot_scatter_AD <- function(female_90_1mb_withld,bp,max_ylim,df,rsid,ylabel){
    bp_start <- min(female_90_1mb_withld$BP)
    bp_end <- max(female_90_1mb_withld$BP)

    loc1 <- locus(data = df, chrom = "CHR", pos = "BP", p = "P", index_snp = rsid, LD = "r2", 
                    flank = c(bp - bp_start, bp_end - bp), ens_db = "EnsDb.Hsapiens.v86")
    # loc1 <- link_recomb(loc1, genome = "hg38")
    bw_url <- "http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recomb1000GAvg.bw"
    tmp <- tempfile(fileext = ".bw")
    download.file(bw_url, tmp, mode = "wb")
    recomb_gr <- import.bw(tmp)
    loc1 <- link_recomb(loc1, recomb = recomb_gr)

    p <- gg_scatter(loc1, pcutoff = 5e-8, ylim = c(0, max_ylim), labels = NULL, xticks = FALSE, cex.lab = 1, ylab = ylabel,legend_pos=NULL) +
            theme(
                axis.title.y = element_text(size = 8, face = "bold"),     
                axis.title.y.right = element_text(size = 8, face = "plain"),
                plot.title = element_text(
                    size = 14, 
                    face = "bold",
                    hjust = 0.5, 
                    margin = margin(b = 10)
                ) 
            ) +
            geom_point(data = subset(loc1$data, SNP == loc1$index_snp), 
                    mapping = aes(x = BP / 1e6, y = logP), 
                    shape = 23,       
                    fill = "purple",  
                    color = "black",  
                    size = 4,        
                    inherit.aes = FALSE)

    return(p)
}

subplot_scatter_qtl <- function(female_90_1mb_withld,bp,df,rsid,ylabel,inorout){
    bp_start <- min(female_90_1mb_withld$BP)
    bp_end <- max(female_90_1mb_withld$BP)
    p_clean   <- as.numeric(df$P); p_clean <- p_clean[is.finite(p_clean) & p_clean > 0]
    min_ylim  <- min(-log10(p_clean), na.rm=TRUE); max_ylim <- ceiling(max(-log10(p_clean), na.rm=TRUE))

    loc0 <- locus(data = df, index_snp = rsid, flank = c(bp - bp_start, bp_end - bp),
                 chrom = "CHR", pos = "BP", p = "P", LD = "r2", labs = "SNP", ens_db = "EnsDb.Hsapiens.v86")
    bw_url <- "http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recomb1000GAvg.bw"
    tmp <- tempfile(fileext = ".bw")
    download.file(bw_url, tmp, mode = "wb")
    recomb_gr <- import.bw(tmp)
    loc0 <- link_recomb(loc0, recomb = recomb_gr)

    p <- gg_scatter(loc0, pcutoff = 9e-800, labels = NULL, xticks = FALSE, cex.lab = 1, ylim = c(min_ylim, max_ylim),ylab = ylabel,legend_pos=NULL) +
            theme(
                axis.title.y = element_text(size = 8, face = "bold"),     
                axis.title.y.right = element_text(size = 8, face = "plain"),
                plot.title = element_text(
                    size = 14, 
                    face = "bold",
                    hjust = 0.5, 
                    margin = margin(b = 10)
                ) 
            ) +
            geom_point(data = subset(loc0$data, SNP == loc0$index_snp), 
                        mapping = aes(x = BP / 1e6, y = logP), 
                        shape = 23,       
                        fill = "purple",  
                        color = "black",  
                        size = 4,        
                        inherit.aes = FALSE)

    if(inorout =="in"){
        return(p)
    }else if(inorout =="out"){
        p$layers <- Filter(function(l){
        !(inherits(l$geom, "GeomPoint") &&
            !is.null(l$aes_params$shape) &&
            l$aes_params$shape == 23)
        }, p$layers)
        return(p)
    }

}


subplot_genetrack <- function(female_90_1mb_withld,bp,rsid,size,...){
    bp_start <- min(female_90_1mb_withld$BP)
    bp_end <- max(female_90_1mb_withld$BP)
    loc0 <- locus(data = female_90_1mb_withld, index_snp = rsid, flank = c(bp - bp_start, bp_end - bp),
                chrom = "CHR", pos = "BP", p = "P", LD = "r2", labs = "SNP", ens_db = "EnsDb.Hsapiens.v86")
    g1a <- gg_genetracks(loc0,
            highlight = c(...),
            cex.text = size,
            filter_gene_biotype = "protein_coding", 
            gene_col = "skyblue", 
            showExons = TRUE, 
            xticks = TRUE)
}



subplot_compare <- function(AD_data, pQTL_data, chr, index_snp, ad_title, pqtl_title,position) {
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

    length <- nrow(AD_data)
    ld_snp <- data.frame(SNP_A = rep(index_snp, length),SNP_B = AD_data$SNP, R2 = AD_data$r2)
    color_df <- locuscomparer::assign_color(merged_data$rsid, snp = index_snp, ld = ld_snp)
    shape_vec <- ifelse(merged_data$rsid == index_snp, 23, 21) %>% stats::setNames(merged_data$rsid)
    size_vec <- ifelse(merged_data$rsid == index_snp, 3, 2) %>% stats::setNames(merged_data$rsid)
    
    sc <- make_scatterplot_ycy(  
        merged = merged_data,
        title1 = ad_title,
        title2 = pqtl_title,
        color = color_df,
        shape = shape_vec,
        size = size_vec,
        legend = FALSE,
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
        labs(x = bquote(bold(.(ad_title)) ~ log[10](P)),
            y = bquote(bold(.(pqtl_title)) ~ log[10](P))) +
        theme(
        axis.title.x = element_text(size = 8, face = "bold"),
        axis.title.y = element_text(size = 8, face = "bold")
        )

    
    return(sc)
}

subplot_blank <- function(text){
    sc <- ggplot() +
    # set coordinate limits so that (0.5,0.5) is the center
    xlim(0, 1) + 
    ylim(0, 1) +
    # remove all axes, gridlines, ticks, etc.
    theme_void() +
    # force the panel background to solid white
    theme(panel.background = element_rect(fill = "white", colour = "white")) +
    # add your message at the center
    annotate(
        "text",
        x = 0.5, y = 0.5,
        label = text,
        size = 4,             # adjust text size as needed
        hjust = 0.5,          # horizontal centering
        vjust = 0.5           # vertical centering
    )
    return(sc)
}






make_scatterplot_ycy = function (merged, title1, title2, color, shape, size, legend = TRUE, legend_position = c('bottomright','topright','topleft')) {

    p = ggplot(merged, aes(logp1, logp2)) +
        geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
        geom_point(data = merged[merged$label != "",],
                   aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)) +
        xlab(bquote(.(title1) ~ log[10] * '(P)')) +
        ylab(bquote(.(title2) ~ log[10] * '(P)')) +
        scale_fill_manual(values = color, guide = "none") +
        scale_shape_manual(values = shape, guide = "none") +
        scale_size_manual(values = size, guide = "none") +
        ggrepel::geom_text_repel(aes(label = label))+
        theme_classic()

    if (legend == TRUE) {
        legend_position = match.arg(legend_position)
        if (legend_position == 'bottomright'){
            legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
        } else if (legend_position == 'topright'){
            legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
        } else {
            legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
        }

        p = ggdraw(p) +
            geom_rect(data = legend_box,
                      aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                      color = "black",
                      fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
            draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
            draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
            draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
            draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
    }

    return(p)
}