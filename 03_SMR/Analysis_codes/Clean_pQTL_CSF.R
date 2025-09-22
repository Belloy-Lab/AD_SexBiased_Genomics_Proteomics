#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(vroom)
  library(dplyr)
})

# Define command-line options
option_list <- list(
  make_option("--pqtl_female_dir", type = "character", help = "Directory path to female pQTL files"),
  make_option("--pqtl_male_dir", type = "character", help = "Directory path to male pQTL files"),
  make_option("--pqtl_both_dir", type = "character", help = "Directory path to both-sex pQTL files"),
  make_option("--gwasdir", type = "character", help = "Directory path to GWAS files"),
  make_option("--female_gwas", type = "character", help = "Female GWAS summary statistics filename"),
  make_option("--male_gwas", type = "character", help = "Male GWAS summary statistics filename"),
  make_option("--gene", type = "character", help = "Full path to candidate gene list txt file"),
  make_option("--output_dir", type = "character", default = ".", help = "Output directory [default: current dir]")
)


# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Check for missing required options
required_opts <- c("pqtl_female_dir", "pqtl_male_dir", "pqtl_both_dir", "gwasdir", 
                   "female_gwas", "male_gwas", "gene", "output_dir")
missing_opts <- required_opts[!sapply(required_opts, function(x) !is.null(opt[[x]]))]
if (length(missing_opts) > 0) {
  stop("Missing required arguments: ", paste(missing_opts, collapse = ", "))
}

# Step 1 Fileter SNPs under significant threshold in the cis-region
candidate_gene_list <- fread(file.path(opt$gene))
candidate_gene_list$pos <- (candidate_gene_list$Start_hg38+candidate_gene_list$End_hg38)/2

# selecting IV from NGI CSF Female/Male/Both GWAS pQTL datasets
candidate_gene_list_NGI <- filter(candidate_gene_list,candidate_gene_list$Tissue=='CSF')

# [1] NGI Female
dir_path <- file.path(opt$pqtl_female_dir, paste0("pqtl.",candidate_gene_list_NGI$Analyte,".glm.linear.gz"))
NGI_pQTL_female_final_all <- vroom(dir_path[1])
NGI_pQTL_female_final_all <- NGI_pQTL_female_final_all[0,]

for(i in 1:nrow(candidate_gene_list_NGI)) {
  NGI_pQTL_female <- vroom(dir_path[i])
  NGI_pQTL_female <- as.data.frame(NGI_pQTL_female)
  NGI_pQTL_female$P <- 10^(-NGI_pQTL_female$LOG10_P)
  
  # add gene name
  NGI_pQTL_female$updated_Gene <- candidate_gene_list_NGI$updated_Gene[i]
  NGI_pQTL_female$Gene_chr <- candidate_gene_list_NGI$CHR[i]
  NGI_pQTL_female$Gene_pos <- candidate_gene_list_NGI$pos[i]
  
  # # selecting threshold P < 5e-8
  NGI_pQTL_female_filtered <- NGI_pQTL_female # not select
  
  # ensure cis region (±2MB)
  NGI_pQTL_female_final <- filter(NGI_pQTL_female_filtered,
                                  NGI_pQTL_female_filtered$`#CHROM`==NGI_pQTL_female_filtered$Gene_chr &
                                  NGI_pQTL_female_filtered$POS >= NGI_pQTL_female_filtered$Gene_pos-1000000 &
                                  NGI_pQTL_female_filtered$POS <= NGI_pQTL_female_filtered$Gene_pos+1000000)
  
  NGI_pQTL_female_final_all <- rbind(NGI_pQTL_female_final_all,NGI_pQTL_female_final)
}

##########################################################################################
# [2] NGI Male
##########################################################################################
dir_path <- file.path(opt$pqtl_male_dir, paste0("pqtl.",candidate_gene_list_NGI$Analyte,".glm.linear.gz"))
NGI_pQTL_male_final_all <- vroom(dir_path[1])
NGI_pQTL_male_final_all <- NGI_pQTL_male_final_all[0,]

for (i in 1:nrow(candidate_gene_list_NGI)) {
  NGI_pQTL_male <- vroom(dir_path[i])
  NGI_pQTL_male <- as.data.frame(NGI_pQTL_male)
  NGI_pQTL_male$P <- 10^(-NGI_pQTL_male$LOG10_P)
  
  # add gene name
  NGI_pQTL_male$updated_Gene <- candidate_gene_list_NGI$updated_Gene[i]
  NGI_pQTL_male$Gene_chr <- candidate_gene_list_NGI$CHR[i]
  NGI_pQTL_male$Gene_pos <- candidate_gene_list_NGI$pos[i]
  
  # # selecting threshold P < 5e-8
  NGI_pQTL_male_filtered <- NGI_pQTL_male # not select
  
  
  # ensure cis region (±2MB)
  NGI_pQTL_male_final <- filter(NGI_pQTL_male_filtered,
                                NGI_pQTL_male_filtered$`#CHROM`==NGI_pQTL_male_filtered$Gene_chr &
                                NGI_pQTL_male_filtered$POS >= NGI_pQTL_male_filtered$Gene_pos-1000000 &
                                NGI_pQTL_male_filtered$POS <= NGI_pQTL_male_filtered$Gene_pos+1000000)
  
  NGI_pQTL_male_final_all <- rbind(NGI_pQTL_male_final_all,NGI_pQTL_male_final)
}

##########################################################################################
# [3] NGI Both
##########################################################################################
dir_path <- file.path(opt$pqtl_both_dir, paste0(candidate_gene_list_NGI$Analyte,".glm.linear.gz"))
NGI_pQTL_both_final_all <- vroom(dir_path[1])
NGI_pQTL_both_final_all <- NGI_pQTL_both_final_all[0,]

for (i in 1:nrow(candidate_gene_list_NGI)) {
  NGI_pQTL_both <- vroom(dir_path[i])
  NGI_pQTL_both <- as.data.frame(NGI_pQTL_both)
  NGI_pQTL_both$P <- 10^(-NGI_pQTL_both$LOG10_P)
  
  # add gene name
  NGI_pQTL_both$updated_Gene <- candidate_gene_list_NGI$updated_Gene[i]
  NGI_pQTL_both$Gene_chr <- candidate_gene_list_NGI$CHR[i]
  NGI_pQTL_both$Gene_pos <- candidate_gene_list_NGI$pos[i]
  
  # # selecting threshold P < 5e-8
  NGI_pQTL_both_filtered <- NGI_pQTL_both # not select
  
  # ensure cis region (±2MB)
  NGI_pQTL_both_final <- filter(NGI_pQTL_both_filtered,
                                NGI_pQTL_both_filtered$`#CHROM`==NGI_pQTL_both_filtered$Gene_chr &
                                NGI_pQTL_both_filtered$POS >= NGI_pQTL_both_filtered$Gene_pos-1000000 &
                               NGI_pQTL_both_filtered$POS <= NGI_pQTL_both_filtered$Gene_pos+1000000)
  
  NGI_pQTL_both_final_all <- rbind(NGI_pQTL_both_final_all,NGI_pQTL_both_final)
}

###### Step 2 exclude SNPs not exist in the AD GWAS and update rsid in the pQTL ##########
NGI_pQTL_female_final_all$gene_id <- paste0(NGI_pQTL_female_final_all$updated_Gene,"_female")
NGI_pQTL_male_final_all$gene_id <- paste0(NGI_pQTL_male_final_all$updated_Gene,"_male")
NGI_pQTL_both_final_all$gene_id <- paste0(NGI_pQTL_both_final_all$updated_Gene,"_both")

# keep colnames identical
NGI_pQTL_female_final_all <- select(NGI_pQTL_female_final_all, "#CHROM",POS,ID,ALT,REF,BETA,SE,P,updated_Gene,Gene_chr,Gene_pos,gene_id)
NGI_pQTL_male_final_all <- select(NGI_pQTL_male_final_all,colnames(NGI_pQTL_female_final_all))
NGI_pQTL_both_final_all <- select(NGI_pQTL_both_final_all, "#CHROM",POS,ID,ALT,REF,BETA,SE,P,updated_Gene,Gene_chr,Gene_pos,gene_id)

NGI_pQTL_combine <- rbind(NGI_pQTL_female_final_all,NGI_pQTL_male_final_all,NGI_pQTL_both_final_all)

# Load GWAS summary stats
trait_name <- c("Female", "Male")
trait_dir <- c(file.path(opt$gwasdir, opt$female_gwas),
               file.path(opt$gwasdir, opt$male_gwas))

trait_all <- vroom(trait_dir[1]) |> mutate(trait = trait_name[1])
trait_all <- select(trait_all, CHR, BP, SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N_incl, trait)

for (i in 2:length(trait_name)) {
  trait <- vroom(trait_dir[i]) |> mutate(trait = trait_name[i])
  trait <- select(trait, all_of(colnames(trait_all)))
  trait_all <- bind_rows(trait_all, trait)
}

trait_all$chr_pos_a12 <- paste0(trait_all$CHR,":", trait_all$BP,":", trait_all$ALLELE1,":", trait_all$ALLELE0)
NGI_pQTL_combine$chr_pos_a12 <- paste0(NGI_pQTL_combine$`#CHROM`,":", NGI_pQTL_combine$POS,":", NGI_pQTL_combine$ALT,":", NGI_pQTL_combine$REF)
NGI_pQTL_combine$chr_pos_a21 <- paste0(NGI_pQTL_combine$`#CHROM`,":", NGI_pQTL_combine$POS,":", NGI_pQTL_combine$REF,":", NGI_pQTL_combine$ALT)

NGI_pQTL_combine <- NGI_pQTL_combine %>%
  filter(chr_pos_a12 %in% trait_all$chr_pos_a12 | chr_pos_a21 %in% trait_all$chr_pos_a12)

# [2] update SNP
# add rsid
trait_all_distinct <- dplyr::select(trait_all,SNP,chr_pos_a12)
trait_all_distinct <- distinct(trait_all_distinct,chr_pos_a12,.keep_all = T)

NGI_pQTL_combine <- left_join(NGI_pQTL_combine,trait_all_distinct,by = c("chr_pos_a12"="chr_pos_a12"))
NGI_pQTL_combine <- left_join(NGI_pQTL_combine,trait_all_distinct,by = c("chr_pos_a21"="chr_pos_a12"))
NGI_pQTL_combine$SNP <- NA
NGI_pQTL_combine <- mutate(NGI_pQTL_combine,SNP = coalesce(SNP,SNP.x,SNP.y))

NGI_pQTL_combine$chr_pos_a12 <- NULL
NGI_pQTL_combine$chr_pos_a21 <- NULL
NGI_pQTL_combine$SNP.x <- NULL
NGI_pQTL_combine$SNP.y <- NULL
NGI_pQTL_combine$ID <- NULL

# Step 3 Make besd-format files for SMR analysis
my_flist <- distinct(NGI_pQTL_combine,NGI_pQTL_combine$gene_id,.keep_all=T)
my_flist <- my_flist[,c("#CHROM","gene_id","Gene_pos","updated_Gene")]

my_flist$GeneticDistance <- 0 #
my_flist$Orientation <- NA #

esd_dir <- file.path(opt$output_dir, "Combined_esd_files")
if (!dir.exists(esd_dir)) {
  dir.create(esd_dir, recursive = TRUE)
}

my_flist$PathOfEsd <- file.path(esd_dir, paste0(my_flist$gene_id, ".esd"))

my_flist <- select(my_flist, "#CHROM", gene_id, GeneticDistance, Gene_pos, updated_Gene, Orientation, PathOfEsd)
colnames(my_flist) <- c("Chr", "ProbeID", "GeneticDistance", "ProbeBp", "Gene", "Orientation", "PathOfEsd")

# Save as Linux-readable format
write.table(my_flist, file.path(opt$output_dir, "CSF_my.flist"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8", eol = "\n")

for (probe_ID in unique(NGI_pQTL_combine$gene_id)) {
  esd <- filter(NGI_pQTL_combine,gene_id==probe_ID)
  esd$MAF <- NA # NGI don't have MAF, add NA
  esd_save_name <- file.path(esd_dir, paste0(esd$gene_id[1],".esd")) # save path should be identical to the my_flist # nolint
  esd <- select(esd,
                "#CHROM",SNP,POS,ALT,REF,MAF,BETA,SE,P)
  colnames(esd) <- c("Chr","SNP","Bp","A1","A2","Freq","Beta","se","p")
  write.table(esd,esd_save_name,sep = "\t", quote = FALSE, row.names = FALSE)
}
