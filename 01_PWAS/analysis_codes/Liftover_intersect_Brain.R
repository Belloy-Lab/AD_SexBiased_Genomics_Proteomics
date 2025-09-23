################################################################# Danielle Reid
##################################################################### 12/03/24
## Lift over to hg19 and intersect with 1KG EUR LD panel #######################
################################################################################

## Enable command line arguments
args = commandArgs(TRUE)

################################################################################
## Libraries

library(data.table)
library(plyr)
library(dplyr)
library(openxlsx)
library(tidyr)
library(methods)
library(optparse)

#############################################################
option_list = list(
  make_option("--ld_fold", action="store", default=NA, type='character',
              help="Path to LD files [required]"),
  make_option("--snplist", action="store", default=NA, type='character',
              help="LD SNPlist file name [required]"),        
  make_option("--ld_tab", action="store", default=NA, type='character',
              help="LD tab file name [required]"),
  make_option("--ref_fold", action="store", default=NA, type='character',
              help="Path to chain file directory [required]"),
  make_option("--chain", action="store", default=NA, type='character',
              help="Liftover chain file name [required]"),        
  make_option("--dir", action="store", default=NA, type='character',
              help="Path to original summary stats directory [required]"),
  make_option("--out_dir", action="store", default=NA, type='character',
              help="Output directory for processed summary stats [required]"),
  make_option("--stats", action="store", default=NA, type='character',
              help="Summary statistics filename [required]"),
  make_option("--rsID", action="store", default=NA, type='integer',
              help="Original rsID column number [required]"),
  make_option("--EA", action="store", default=NA, type='character',
              help="Original REF/effect allele or A1 column name [required]"),
  make_option("--NEA", action="store", default=NA, type='character',
              help="Original ALT/non-effect allele or A2 column name [required]"),
  make_option("--file_name", action="store", default=NA, type='character',
              help="Output base filename [required]")            
)

opt = parse_args(OptionParser(option_list=option_list))

# Validate the inputs
if (is.na(opt$rsID) || is.na(opt$EA) || is.na(opt$NEA)) {
  stop("All options --rsID, --EA, --NEA are required and must be valid.")
}

#path to directory where the LD files are located
ld_fold = opt$ld_fold
ld_snp_list = opt$snplist; ld_pos_list = opt$ld_tab
ld_snp = paste(ld_fold,ld_snp_list,sep=""); ld_pos = paste(ld_fold,ld_pos_list,sep="")

ref_fold = opt$ref_fold
chain = opt$chain
dir = opt$dir
out_dir = opt$out_dir
ldfile = fread(ld_pos,header=T);

#create new column name posID in LD SNP positions file by concatenating the CHR and BP column separated by a :
ldfile$posID = paste(ldfile$CHR,ldfile$BP,sep=":")

sumstats = fread(paste(dir, opt$stats, sep = ""), header = T)      #read in file
sumstats = na.omit(sumstats)       #remove NAs if any
EA = opt$EA; NEA = opt$NEA

# Ensure that the column indices and names are valid
if (opt$rsID <= 0 || opt$rsID > ncol(sumstats)) {
  stop("rsID column index is out of bounds.")
}

if (!opt$EA %in% colnames(sumstats)) {
  stop("EA column name is not found in the data frame.")
}

if (!opt$NEA %in% colnames(sumstats)) {
  stop("NEA column name is not found in the data frame.")
}

rsID_col = opt$rsID
names(sumstats)[rsID_col] <- "rsIDhg38"    #rename hg38 SNP column
sumstats$posID = paste(sumstats$CHR,sumstats$BP,sep=":")       #create posID column

sumstats$SNP = paste(sumstats$CHR,sumstats$BP,sumstats[[EA]],sumstats[[NEA]],sep=":")   #CHR:BP:ALT:REF (CHR:BP:EA:NEA)
sumstats = sumstats[order(CHR, BP)]       # sort chr/bp

## create bed file
bed = sumstats[,c("CHR","BP","BP","SNP")]
setnames(bed, 1:4, c("CHR","START","END","NAME"));
bed = bed[, CHR:=as.character(CHR)]; bed = bed[, CHR:=paste0("chr", CHR)]
bed[,"END"] = bed[,"END"]+1
bed = bed[, START:=as.integer(START)]; bed = bed[, END:=as.integer(END)]
bed_file = paste(out_dir, opt$file_name, "to_lift_hg38.txt" ,sep="")
fwrite(bed, bed_file, col.names=F,row.names=F,quote=F, sep='\t')

## Execute liftOver
liftOver = paste(ref_fold, 'liftOver', sep="")
lifted = paste(out_dir,opt$file_name, "lifted_hg19.txt",sep="")
not_lifted = paste(out_dir, opt$file_name, "not_lifted_hg38tohg19.txt",sep="")
chain_file = paste(ref_fold, chain, sep="")
cmd = paste(liftOver, bed_file, chain_file, lifted, not_lifted, sep= " ")
system(cmd)

## Read in lifted variants
bedl = fread(lifted, sep="\t", header=F)
setnames(bedl, 1:4, c("CHRl","START","END","SNP"));

#take the original summary statistics and merge it with the lifted coordinates dataframe only including SNP, START, and CHRl columns by the SNP column (adds START and CHRl columns)
sumstats = inner_join(x = sumstats, y = bedl[,c("SNP","START","CHRl")], by = "SNP")
sumstats$CHRl = gsub("chr","",sumstats$CHRl); pos = which(sumstats$CHR != sumstats$CHRl)
if (any(pos)) {sumstats = sumstats[-pos,]}
sumstats[,"BP"] = sumstats[,"START"]; sumstats = sumstats[,1:(dim(sumstats)[2]-2)]

sumstats$posID = paste(sumstats$CHR,sumstats$BP,sep=":")
sumstats = as.data.frame(sumstats); sumstats = sumstats[,c(1:(dim(sumstats)[2]-1)),]

fwrite(sumstats,paste(out_dir, opt$file_name, "hg19.txt",sep=""),
col.names=T, sep="\t"
)

sumstats = sumstats[which(sumstats$posID %in% ldfile$posID),]
sumstats = merge(sumstats,ldfile[,c("posID","SNP")],by="posID")
fwrite(sumstats,paste(out_dir, opt$file_name, "hg19_intersected.txt",sep=""),
       col.names=T, sep="\t"
)

