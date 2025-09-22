################################################################# Danielle Reid
##################################################################### 05/28/24
###create an LD ref panel for EUR ancestry from 1000G build 37

################################################################################ 
## libraries

library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(methods)
library(optparse)

################################################################################
option_list = list(
  make_option("--PATH_plink1.9", action="store", default="plink1.9", type='character',
              help="Path to plink1.9 executable [%default]"),
  make_option("--PATH_plink2", action="store", default="plink2", type='character',
              help="Path to plink2 executable [%default]")
)

opt = parse_args(OptionParser(option_list=option_list))

if ( system( paste(opt$PATH_plink1.9,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink 1.9 could not be executed, set with --PATH_plink1.9\n" , sep='', file=stderr() )
	q()
}

if ( system( paste(opt$PATH_plink2,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink2 could not be executed, set with --PATH_plink2\n" , sep='', file=stderr() )
	q()
}

################################################################################ 
## Input

# general inputs
tDir <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/"

# 1KG inputs 
Fold_1000G <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/"
Fil_1000G <- "1000g_EUR"

# 1KG output, please make sure to include your path to store outputs from your run of making EUR LD panel.
tDir_1000G <- paste(tDir,"1KG_EUR/",sep="")

# genetic map file
hg19map_fold <- "/storage1/fs1/belloy/Active/02_Data/01_Incoming/1000G_Grch37/genetic_map_b37/"

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
arg <- paste(opt$PATH_plink2, 
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
arg <- paste(opt$PATH_plink1.9, 
                 "--bfile", paste(tDir_1000G, Fil_1000G, sep=""),
                 "--allow-no-sex",
                 "--cm-map",paste(hg19map_fold,"genetic_map_chr@_combined_b37.txt",sep=""),
                 "--make-bed --out", paste(tDir_1000G, Fil_1000G, "_cm", sep=""),
                 sep=" ")
    cat("Executing command:", arg, "\n")
system(arg) 

#split the plink files by chr
for (ch in 1:22) { 
arg <- paste(opt$PATH_plink2, 
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