suppressPackageStartupMessages({
    library(optparse)
    library(vroom)
    library(dplyr)
    library(data.table)
    library(stringr)
})

# Define command-line options
option_list <- list(
  make_option("--gene", type = "character", help = "Full path to candidate gene list txt file"),
  make_option("--smr_res_female", type = "character", help = "Female SMR output file"),
  make_option("--smr_res_male", type = "character", help = "Male SMR output file"),
  make_option("--work_dir", type = "character", default = ".", help = "Working directory [default: current dir]"),
  make_option("--output", type = "character", default = "SMR_res.txt", help = "Output filename [default: SMR_res.txt]", metavar = "filename")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# candidate gene list
gene_list <- fread(opt$gene)
gene_list$a <- paste0(gene_list$Gender,"_",gene_list$Discovery)

gene_list$a <- gsub("Female_Primary","F_Female",gene_list$a)
gene_list$a <- gsub("Female_Secondary","Both_Female",gene_list$a)
gene_list$a <- gsub("Female_Both","F_Female",gene_list$a)
gene_list$a <- gsub("Male_Secondary","Both_Male",gene_list$a)

gene_list$a <- paste0(gene_list$updated_Gene,"_",gene_list$a)


# combine smr results
res_female <- vroom(file.path(opt$work_dir, opt$smr_res_female))
res_female$`AD Gender` <- "Female"

res_male <- vroom(file.path(opt$work_dir, opt$smr_res_male))
res_male$`AD Gender` <- "Male"

SMR_res <- rbind(res_female,res_male)
SMR_res <- SMR_res %>%
  mutate(
    sex_pQTL = str_extract(probeID, "(?<=_)[^_]+"),  # Extract the part after the first underscore
    sex_pQTL = case_when(
      sex_pQTL == "female" ~ "F",
      sex_pQTL == "male" ~ "M",
      sex_pQTL == "both" ~ "Both",
      TRUE ~ sex_pQTL  # Keep as-is if not matched (optional)
    )
  )
SMR_res$a <- paste0(SMR_res$Gene,"_",SMR_res$sex_pQTL,"_",SMR_res$`AD Gender`)

#
SMR_res <- semi_join(SMR_res,gene_list,by='a')
SMR_res$a <- NULL
SMR_res$fdr_SMR <- p.adjust(SMR_res$p_SMR,method='fdr')

write.table(SMR_res, file.path(opt$work_dir, opt$output), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
