##prep hg38 sex-strat ADGC_ADSP_UKB_FinnGen_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var summary statistics txt files for mungestats by intersecting with Dan's LD ref panel
#####################Danielle Reid
#################### 12/04/24

#load libraries
library(data.table)
library(plyr)
library(dplyr)
library(openxlsx)
library(readr)

#create path for inputs
dir = "/storage1/fs1/belloy/Active/02_Data/01_Incoming/Belloy_2024_GWAS_AD_november_update/"

#create a dataframe of the hg38 summary stats files
F_stats = fread(paste(dir, "ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var", sep = ""), header = T, sep = "\t")
M_stats = fread(paste(dir, "ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var", sep = ""), header = T, sep = "\t")
F_stats = as.data.frame(F_stats); M_stats = as.data.frame(M_stats)
print("loaded Female hg38 meta sumstats"); print(nrow(F_stats))
print("loaded Male hg38 meta sumstats"); print(nrow(M_stats))

#remove NAs
F_stats = na.omit(F_stats)
print("nrow of Female meta sumstats after removing NAs"); nrow(F_stats)
M_stats = na.omit(M_stats)
print("nrow of Male meta sumstats after removing NAs"); nrow(M_stats)

#change SNP column names
names(F_stats)[2] = "ID"
names(M_stats)[2] = "ID"

#add tmpID column as CHR:BP:A2(ALT):A1(REF) to match IDs with Dan's data (column 5: column 6) (MB sumstats ALLELE0 = ALT, ALLELE1= REF)
F_stats$tmpID = paste(F_stats$CHR,F_stats$BP,F_stats$ALLELE0,F_stats$ALLELE1,sep=":")
M_stats$tmpID = paste(M_stats$CHR,M_stats$BP,M_stats$ALLELE0,M_stats$ALLELE1,sep=":")

##intersect the summary stats files with Dan Western's LD ref panel
#create variable for Dan's folder and load merged bim file into a dataframe and add a header for mungestats
Dan_dir = "/storage1/fs1/belloy/Active/02_Data/01_Incoming/NGI_CSF_pQTL_files/PWAS_LD_Reference_Files/"
bim_dat = fread(paste(Dan_dir, "chr1-22.txt", sep = ""), header = F, sep ='\t')
header_names = c("CHR", "ID", "POS", "BP", "A2", "A1")
colnames(bim_dat) = header_names

#add tmpID column in Dan's LD ref panel as CHR:BP:A1(ref):A2(alt)- this is because some alleles are flipped in his ID names, so this will serve as a temporary ID for matching
bim_dat$tmpID = paste(bim_dat$CHR,bim_dat$BP,bim_dat$A1,bim_dat$A2,sep=":")
bim_dat$tmpID2 = paste(bim_dat$CHR,bim_dat$BP,bim_dat$A2,bim_dat$A1,sep=":")

#intersect sumstats with Dan's LD ref panel where CHR:BP:REF:ALT match
F_stats1 = F_stats[which(F_stats$tmpID %in% bim_dat$tmpID),]
print("nrow of Female meta sumstats after intersection 1"); nrow(F_stats1)

M_stats1 = M_stats[which(M_stats$tmpID %in% bim_dat$tmpID),]
print("nrow of Male meta sumstats after intersection 1"); nrow(M_stats1)

# intersect sumstats with Dan's LD ref panel by tmpID2
F_stats2 = F_stats[which(F_stats$tmpID %in% bim_dat$tmpID2),]
print("nrow of Female meta sumstats after intersection 2"); nrow(F_stats2)

M_stats2 = M_stats[which(M_stats$tmpID %in% bim_dat$tmpID2),]
print("nrow of Male meta sumstats after intersection 2"); nrow(M_stats2)

#concatenate the results into one dataframe each
F = rbind(F_stats1, F_stats2)
print("total nrow of Female meta sumstats after intersection"); nrow(F)

M = rbind(M_stats1, M_stats2)
print("total nrow of Male meta sumstats after intersection"); nrow(M)

#add SNP column with Dan's ID names where the values match
F$SNP <- ifelse(F$tmpID %in% bim_dat$tmpID, 
                          bim_dat$ID[match(F$tmpID, bim_dat$tmpID)], 
                          bim_dat$ID[match(F$tmpID, bim_dat$tmpID2)]
)
M$SNP = ifelse(M$tmpID %in% bim_dat$tmpID, 
                          bim_dat$ID[match(M$tmpID, bim_dat$tmpID)], 
                          bim_dat$ID[match(M$tmpID, bim_dat$tmpID2)]
)

out_dir = "/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/"

#save the intersected summary stats files
fwrite(F,paste(out_dir,"ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.no_APOE.shared_var.hg38_NGI_CSF_intersected.txt",sep=""),
       col.names=T, sep="\t"
)
fwrite(M,paste(out_dir,"ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.no_APOE.shared_var.hg38_NGI_CSF_intersected.txt",sep=""),
       col.names=T, sep="\t"
)
