#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(vroom)
  library(dplyr)
})

# Define command-line options
option_list <- list(
  make_option("--pqtl_dir", type = "character", help = "Directory path to brain pQTL files"),
  make_option("--pqtl_female", type = "character", help = "Female brain pQTL file"),
  make_option("--pqtl_male", type = "character", help = "Male brain pQTL file"),
  make_option("--pqtl_both", type = "character", help = "non-sex stratified brain pQTL file"),
  make_option("--gwasdir", type = "character", help = "Directory path to GWAS files"),
  make_option("--female_gwas", type = "character", help = "Female GWAS summary statistics filename"),
  make_option("--male_gwas", type = "character", help = "Male GWAS summary statistics filename"),
  make_option("--gene", type = "character", help = "Full path to candidate gene list txt file"),
  make_option("--output_dir", type = "character", default = ".", help = "Output directory [default: current dir]")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

candidate_gene_list <- fread(file.path(opt$gene))
candidate_gene_list$pos <- (candidate_gene_list$Start_hg38+candidate_gene_list$End_hg38)/2
candidate_gene_list_wingo <- filter(candidate_gene_list,candidate_gene_list$Tissue=='Brain')
colnames(candidate_gene_list_wingo)[2] <- "gene_id"

# [1] wingo Female
wingo_pQTL_female <- vroom(file.path(opt$pqtl_dir, opt$pqtl_female)) ## load female wingo dataset
wingo_pQTL_female <- as.data.frame(wingo_pQTL_female)
wingo_pQTL_female <- semi_join(wingo_pQTL_female,candidate_gene_list_wingo,by='gene_id')

# add gene name
candidate_gene_list_wingo_add_info <- candidate_gene_list_wingo[,c("gene_id","updated_Gene","CHR","pos")]
colnames(candidate_gene_list_wingo_add_info)[3]="Gene_chr"
colnames(candidate_gene_list_wingo_add_info)[4]="Gene_pos"
wingo_pQTL_female <- left_join(wingo_pQTL_female,candidate_gene_list_wingo_add_info,by="gene_id")

wingo_pQTL_female_filtered <- wingo_pQTL_female

# ensure cis region (±2MB)
wingo_pQTL_female_final <- filter(wingo_pQTL_female_filtered,
                                  wingo_pQTL_female_filtered$CHR==wingo_pQTL_female_filtered$Gene_chr &
                                  wingo_pQTL_female_filtered$BP >= wingo_pQTL_female_filtered$Gene_pos-1000000 &
                                  wingo_pQTL_female_filtered$BP <= wingo_pQTL_female_filtered$Gene_pos+1000000)

# [2] wingo male
wingo_pQTL_male <- vroom(file.path(opt$pqtl_dir, opt$pqtl_male)) ## load male wingo dataset
wingo_pQTL_male <- as.data.frame(wingo_pQTL_male)
wingo_pQTL_male <- semi_join(wingo_pQTL_male,candidate_gene_list_wingo,by='gene_id')

# add gene name
candidate_gene_list_wingo_add_info <- candidate_gene_list_wingo[,c("gene_id","updated_Gene","CHR","pos")]
colnames(candidate_gene_list_wingo_add_info)[3]="Gene_chr"
colnames(candidate_gene_list_wingo_add_info)[4]="Gene_pos"
wingo_pQTL_male <- left_join(wingo_pQTL_male,candidate_gene_list_wingo_add_info,by="gene_id")

wingo_pQTL_male_filtered <- wingo_pQTL_male # not select

# ensure cis region (±2MB)
wingo_pQTL_male_final <- filter(wingo_pQTL_male_filtered,
                                  wingo_pQTL_male_filtered$CHR==wingo_pQTL_male_filtered$Gene_chr &
                                    wingo_pQTL_male_filtered$BP >= wingo_pQTL_male_filtered$Gene_pos-1000000 &
                                    wingo_pQTL_male_filtered$BP <= wingo_pQTL_male_filtered$Gene_pos+1000000)

# [3] wingo both
wingo_pQTL_both <- vroom(file.path(opt$pqtl_dir, opt$pqtl_both)) ## load both wingo dataset
wingo_pQTL_both <- as.data.frame(wingo_pQTL_both)
wingo_pQTL_both <- semi_join(wingo_pQTL_both,candidate_gene_list_wingo,by='gene_id')

# add gene name
candidate_gene_list_wingo_add_info <- candidate_gene_list_wingo[,c("gene_id","updated_Gene","CHR","pos")]
colnames(candidate_gene_list_wingo_add_info)[3]="Gene_chr"
colnames(candidate_gene_list_wingo_add_info)[4]="Gene_pos"
wingo_pQTL_both <- left_join(wingo_pQTL_both,candidate_gene_list_wingo_add_info,by="gene_id")

wingo_pQTL_both_filtered <- wingo_pQTL_both # not select

# ensure cis region (±2MB)
wingo_pQTL_both_final <- filter(wingo_pQTL_both_filtered,
                                  wingo_pQTL_both_filtered$CHR==wingo_pQTL_both_filtered$Gene_chr &
                                    wingo_pQTL_both_filtered$BP >= wingo_pQTL_both_filtered$Gene_pos-1000000 &
                                    wingo_pQTL_both_filtered$BP <= wingo_pQTL_both_filtered$Gene_pos+1000000)

# Step 2 exclude SNPs not exist in the AD GWAS and update rsid in the pQTL ##########
wingo_pQTL_female_final$gene_id <- paste0(wingo_pQTL_female_final$gene_id,"_female")
wingo_pQTL_male_final$gene_id <- paste0(wingo_pQTL_male_final$gene_id,"_male")
wingo_pQTL_both_final$gene_id <- paste0(wingo_pQTL_both_final$gene_id,"_both")

# [1] exclude SNP
wingo_pQTL_combine <- rbind(wingo_pQTL_female_final,wingo_pQTL_male_final,wingo_pQTL_both_final)
dim(wingo_pQTL_combine)
colnames(wingo_pQTL_combine)

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
wingo_pQTL_combine$chr_pos_a12 <- paste0(wingo_pQTL_combine$CHR,":", wingo_pQTL_combine$BP,":", wingo_pQTL_combine$A1,":", wingo_pQTL_combine$A2)
wingo_pQTL_combine$chr_pos_a21 <- paste0(wingo_pQTL_combine$CHR,":", wingo_pQTL_combine$BP,":", wingo_pQTL_combine$A2,":", wingo_pQTL_combine$A1)

wingo_pQTL_combine <- wingo_pQTL_combine %>%
  filter(chr_pos_a12 %in% trait_all$chr_pos_a12 | chr_pos_a21 %in% trait_all$chr_pos_a12)

# [2] update SNP, add rsid
trait_all_distinct <- dplyr::select(trait_all,SNP,chr_pos_a12)
trait_all_distinct <- distinct(trait_all_distinct,chr_pos_a12,.keep_all = T)

colnames(wingo_pQTL_combine)[1] <- "rsid"
wingo_pQTL_combine <- left_join(wingo_pQTL_combine,trait_all_distinct,by = c("chr_pos_a12"="chr_pos_a12"))
wingo_pQTL_combine <- left_join(wingo_pQTL_combine,trait_all_distinct,by = c("chr_pos_a21"="chr_pos_a12"))
wingo_pQTL_combine$SNP <- NA
wingo_pQTL_combine <- mutate(wingo_pQTL_combine,SNP = coalesce(SNP,SNP.x,SNP.y))

wingo_pQTL_combine$chr_pos_a12 <- NULL
wingo_pQTL_combine$chr_pos_a21 <- NULL
wingo_pQTL_combine$SNP.x <- NULL
wingo_pQTL_combine$SNP.y <- NULL
wingo_pQTL_combine$rsid <- NULL

# Step 3 Make besd-format files for SMR analysis
my_flist <- distinct(wingo_pQTL_combine,wingo_pQTL_combine$gene_id,.keep_all=T)
my_flist <- my_flist[,c("CHR","gene_id","Gene_pos","updated_Gene")]

my_flist$GeneticDistance <- 0 #
my_flist$Orientation <- NA #

esd_dir <- file.path(opt$output_dir, "Combined_esd_files")
if (!dir.exists(esd_dir)) {
  dir.create(esd_dir, recursive = TRUE)
}

my_flist$PathOfEsd <- file.path(esd_dir, paste0(my_flist$gene_id,".esd"))

my_flist <- select(my_flist, CHR, gene_id, GeneticDistance, Gene_pos, updated_Gene, Orientation, PathOfEsd)
colnames(my_flist) <- c("Chr", "ProbeID", "GeneticDistance", "ProbeBp", "Gene", "Orientation", "PathOfEsd")

# Save as Linux-readable format
write.table(my_flist, file.path(opt$output_dir, "Brain_my.flist"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8", eol = "\n")

for (probe_ID in unique(wingo_pQTL_combine$gene_id)) {
  esd <- filter(wingo_pQTL_combine,gene_id==probe_ID)
  esd_save_name <- file.path(esd_dir, paste0(esd$gene_id[1],".esd"))  # save path should be identical to the my_flist
  esd <- select(esd,
                CHR,SNP,BP,A1,A2,MAF,BETA,SE,P)
  colnames(esd) <- c("Chr","SNP","Bp","A1","A2","Freq","Beta","se","p")
  # the folder "3_esd_files" should be added manually
  write.table(esd,esd_save_name,sep = "\t", quote = FALSE, row.names = FALSE)
}