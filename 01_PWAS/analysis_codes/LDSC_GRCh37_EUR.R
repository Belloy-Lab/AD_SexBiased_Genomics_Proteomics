suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(data.table) 
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(methods)
})

# ---- options ----
option_list <- list(
  make_option(c("-t", "--tDir"), type = "character", help = "Path to the phased 1000 Genomes reference panel (PSAM and PGEN files)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tDir)) {
  stop("Please provide --tDir (top-level input directory).", call. = FALSE)
}

# Input
Fold_1000G <- opt$tDir
Fil_1000G <- "1000g_EUR"

# 1KG output, please make sure to include your path to store outputs from your run of making EUR LD panel.
tDir_1000G <- paste0(tDir, "1KG_EUR/")

# genetic map file
hg19map_fold <- paste0(tDir, "genetic_map_b37/") # recombination rate genetic map files per chr (gr37)

################################################################################ 
# identify EUR subjects in 1KG

# update fam file with ancestry
dem_1000G <- read.table(paste(Fold_1000G,"phase3_corrected.psam",sep=""))
dem_1000G <- dem_1000G[,c(2,1,5)]
colnames(dem_1000G) <- c("SID","IID","sPOP")
dem_1000G$sPOP <- as.character(dem_1000G$sPOP)
dem_1000G <- dem_1000G[which(dem_1000G$sPOP=="EUR"),1:2]
names(dem_1000G)[2] <- "#IID"
dem_1000G <- dplyr::select(dem_1000G, -SID)
write.table(dem_1000G,paste(tDir_1000G,"1K_EUR_IDs.txt",sep=""),
            row.names = F, col.names = F, quote = F
)

################################################################################ 
## Extract EUR subject from 1KG

# get subjects
plink2 <- "/usr/bin/plink2"
arg <- paste(plink2, 
                 "--pgen", paste(Fold_1000G, "all_phase3.pgen", sep=""),
                 "--pvar", paste(Fold_1000G, "all_phase3.pvar", sep=""),
                 "--psam", paste(Fold_1000G,"phase3_corrected.psam",sep=""),
                 "--keep",paste(tDir_1000G,"1K_EUR_IDs.txt",sep=""),
                 "--allow-no-sex",
                 "--allow-extra-chr --chr 1-22",
                 "--maf 0.005",
                 "--max-alleles 2",
                 "--make-bed --out", paste(tDir_1000G, Fil_1000G, sep=""),
                 sep=" ")
    cat("Executing command:", arg, "\n")
system(arg) 

# update centimorgan map in bim file
plink19 <- "/usr/bin/plink1.9"
arg <- paste(plink19, 
                 "--bfile", paste(tDir_1000G, Fil_1000G, sep=""),
                 "--allow-no-sex",
                 "--cm-map",paste(hg19map_fold,"genetic_map_chr@_combined_b37.txt",sep=""),
                 "--make-bed --out", paste(tDir_1000G, Fil_1000G, "_cm", sep=""),
                 sep=" ")
    cat("Executing command:", arg, "\n")
system(arg) 

#split the plink files by chr
for (ch in 1:22) { 
arg <- paste(plink2, 
                 "--bfile", paste(tDir_1000G, Fil_1000G, "_cm", sep=""),
                 "--chr", ch,
                 "--make-bed --out", paste(tDir_1000G, Fil_1000G, "_cm_ch", ch, sep=""),
                 sep=" ")
    cat("Executing command:", arg, "\n")
system(arg) 
}

#create a tab delimited file of the updated plink files with cm that are not separated by chr
bim <- fread(paste(tDir_1000G, Fil_1000G, "_cm.bim", sep=""))
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "ALT", "REF")
bim_tab <- subset(bim, select = c(CHR, SNP, ALT, REF, BP))
bim_tab$posID <- paste(bim_tab$CHR, bim_tab$BP, sep=':')
fwrite(bim_tab, paste(tDir_1000G, Fil_1000G, "_cm.tab", sep=""),
        row.names = F, col.names = T, quote = F, sep = '\t'
)

#create snplist
snplist <- subset(bim, select = c(SNP, ALT, REF))
names(snplist)[2] <- "A1"
names(snplist)[3] <- "A2"
fwrite(snplist, paste(tDir_1000G, Fil_1000G, "_cm.snplist", sep=""),
        row.names = F, col.names = T, quote = F, sep= '\t'
)