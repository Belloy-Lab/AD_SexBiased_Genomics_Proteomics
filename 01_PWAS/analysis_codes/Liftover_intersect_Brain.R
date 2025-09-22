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

################################################################################
## argument checks

# Validate the inputs
if (is.na(opt$rsID) || is.na(opt$EA) || is.na(opt$NEA)) {
  stop("All options --rsID, --EA, --NEA are required and must be valid.")
}

################################################################################
## General inputs/outputs

#path to directory where the LD files are located
ld_fold = opt$ld_fold

#variable for LD SNP list and LD position list files
ld_snp_list = opt$snplist; ld_pos_list = opt$ld_tab

#makes variable for path to LD SNP list and LD position list files
ld_snp = paste(ld_fold,ld_snp_list,sep=""); ld_pos = paste(ld_fold,ld_pos_list,sep="")

#path to directory where chain file is located
ref_fold = opt$ref_fold

#chain file to liftover from hg38 to hg19 obtained from ucsc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
chain = opt$chain

# directory path to summary stats
dir = opt$dir

# directory path for outputs
out_dir = opt$out_dir

################################################################################
## Get LD reference variant list

#variable for reading the LD SNP positions file and indicates there is a header, while
ldfile = fread(ld_pos,header=T);

#create new column name posID in LD SNP positions file by concatenating the CHR and BP column separated by a :
ldfile$posID = paste(ldfile$CHR,ldfile$BP,sep=":")

################################################################################
## grab summary stats files

## process summary stats files (reformat and liftover)

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
#create bed variable using CHR, BP, BP, and SNP columns from sumstats dataframe
bed = sumstats[,c("CHR","BP","BP","SNP")]
#set column 1-4 names in the bed file to CHR, START, END, and NAME, while
setnames(bed, 1:4, c("CHR","START","END","NAME"));
#convert values in bed dataframe in CHR column to character type and assign it back to the CHR column and then prefix each value in CHR column with "chr"
bed = bed[, CHR:=as.character(CHR)]; bed = bed[, CHR:=paste0("chr", CHR)]
#add 1 to each value in the END column for the bed dataframe
bed[,"END"] = bed[,"END"]+1
#convert the START and END values to integers
bed = bed[, START:=as.integer(START)]; bed = bed[, END:=as.integer(END)]
# save bed file
bed_file = paste(out_dir, opt$file_name, "to_lift_hg38.txt" ,sep="")
fwrite(bed, bed_file, col.names=F,row.names=F,quote=F, sep='\t')

## Execute liftOver
#create variable for activating lifOver using the path for liftOver
liftOver = paste(ref_fold, 'liftOver', sep="")
#create a variable including output file path and file name for the lifted coordinates to hg19
lifted = paste(out_dir,opt$file_name, "lifted_hg19.txt",sep="")
#create a variable including the output file path and file name for the unmapped coordinates to hg19
not_lifted = paste(out_dir, opt$file_name, "not_lifted_hg38tohg19.txt",sep="")
#create a variable that prints the path to the liftover directory and calls the chain file
chain_file = paste(ref_fold, chain, sep="")
#create a variable that activates liftOver and points to the created bed file and the chain file to generate the output lifted and unmapped files
cmd = paste(liftOver, bed_file, chain_file, lifted, not_lifted, sep= " ")
#tells the system to execute the cmd variable that performs the liftover
system(cmd)

## Read in lifted variants
bedl = fread(lifted, sep="\t", header=F)
#create a header with 4 columns named CHRl, START, END, SNP
setnames(bedl, 1:4, c("CHRl","START","END","SNP"));

## Merge with original sumstats df and update coordinates
#take the original summary statistics and merge it with the lifted coordinates dataframe only including SNP, START, and CHRl columns by the SNP column (adds START and CHRl columns)
sumstats = inner_join(x = sumstats, y = bedl[,c("SNP","START","CHRl")], by = "SNP")
#use gsub function to remove the prefix "chr" from the values in CHRl column of the sumstats dataframe and create a vector "pos" containing indices of rows where values in CHR of sumstats are not equal to values in CHRl and if they are not equal remove those rows
sumstats$CHRl = gsub("chr","",sumstats$CHRl); pos = which(sumstats$CHR != sumstats$CHRl)
if (any(pos)) {sumstats = sumstats[-pos,]}
#create a new column named "BP" in sumstats dataframe with values from START column and remove last two columns from dataframe (essentially updates BP to lifted and removes columns START and CHRl)
sumstats[,"BP"] = sumstats[,"START"]; sumstats = sumstats[,1:(dim(sumstats)[2]-2)]

## Update posID to replace old CHR and BP with lifted CHR and BP
#add posID column and input values from CHR and BP separating them with a ":"
sumstats$posID = paste(sumstats$CHR,sumstats$BP,sep=":")
# set sumstats df to data frame, drop "SNP" column
sumstats = as.data.frame(sumstats); sumstats = sumstats[,c(1:(dim(sumstats)[2]-1)),]

#save the lifted and cleaned summary statistics file
fwrite(sumstats,paste(out_dir, opt$file_name, "hg19.txt",sep=""),
col.names=T, sep="\t"
)

################################################################################
## intersect with LD variants
#filters sumstats dataframe to only include rows where posIDs match between the sumstats and the LD SNP positions
sumstats = sumstats[which(sumstats$posID %in% ldfile$posID),]

#merge the sumstats dataframe with the LD positions dataframe subsetted to only the posID and SNP column by posID, this gives you lifted rsIDs
sumstats = merge(sumstats,ldfile[,c("posID","SNP")],by="posID")

# save the newly lifted, cleaned, and intersected tab delimited summary statistics file with column names to the output directory
fwrite(sumstats,paste(out_dir, opt$file_name, "hg19_intersected.txt",sep=""),
       col.names=T, sep="\t"
)

