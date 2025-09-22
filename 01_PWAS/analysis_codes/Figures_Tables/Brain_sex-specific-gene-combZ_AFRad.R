############################### Combine W23 PWAS Z-scores for AFR full and ADGC_ADSP_UKB_FinnGen_case_control_reduced_bias.gwama.clean.gen090.no_APOE.shared_var
###################################################### Danielle Reid - 02/04/2025
#############################################################################################

#load the packages
library(data.table)
library(plyr)
library(dplyr)

#############################################################
option_list = list(
  make_option("--EUall_dir", action="store", default=NA, type='character',
              help="EU_all PWAS analysis directory path [required]"),
  make_option("--AFR_dir", action="store", default=NA, type='character',
              help="AFR admixed PWAS priamry analysis directory path [required]"),
  make_option("--f_EUall_stats", action="store", default=NA, type='character',
              help="EU_all 1KGP intersected female GWAS sumstats file [required]"),
  make_option("--m_EUall_stats", action="store", default=NA, type='character',
              help="EU_all 1KGP intersected male GWAS sumstats file [required]"),
  make_option("--f_AFR_stats", action="store", default=NA, type='character',
              help="AFR LD intersected female GWAS sumstats file [required]"),
  make_option("--m_AFR_stats", action="store", default=NA, type='character',
              help="AFR LD intersected male GWAS sumstats file [required]"),
  make_option("--EUall_f_sex_het", action="store", default=NA, type='character',
              help="EU_all female sex-specific heterogeneity file [required]"),
  make_option("--EUall_m_sex_het", action="store", default=NA, type='character',
              help="EU_all male sex-specific heterogeneity file [required]"),
  make_option("--AFR_f_sex_het", action="store", default=NA, type='character',
              help="AFR female sex-specific heterogeneity file [required]"),
  make_option("--AFR_m_sex_het", action="store", default=NA, type='character',
              help="AFR male sex-specific heterogeneity file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#############################################################
f_dir = "Brain/Sex/Female/"
m_dir = "Brain/Sex/Male/"

# load sex-specific genes with sex-het stats for meta and AFR
f_meta = fread(paste0(opt$EUall_dir, f_dir, opt$EUall_f_sex_het))
m_meta = fread(paste0(opt$EUall_dir, m_dir, opt$EUall_m_sex_het))
f_AFR = fread(paste0(opt$AFR_dir, f_dir, opt$AFR_f_sex_het))
m_AFR = fread(paste0(opt$AFR_dir, m_dir, opt$AFR_m_sex_het))
# f_meta = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Brain/Sex/Female/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_fsg_sex-het.txt")
# m_meta = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/Brain/Sex/Male/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_msg_sex-het.txt")
# f_AFR = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/AFR/Brain/Sex/Female/AFRad_Females_AD_cc.full.hg19_AFRad_LD_PWAS_W23_ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var-fsg_sex-het.txt")
# m_AFR = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/AFR/Brain/Sex/Male/AFRad_Males_AD_cc.full.hg19_AFRad_LD_PWAS_W23_ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var-msg_sex-het.txt")
f_meta = as.data.frame(f_meta); f_AFR = as.data.frame(f_AFR); m_meta = as.data.frame(m_meta); m_AFR = as.data.frame(m_AFR)

# load meta and AFR full GWAS summary stats
f_meta_stats = fread(paste0(opt$EUall_dir, opt$f_EUall_stats))
m_meta_stats = fread(paste0(opt$EUall_dir, opt$m_EUall_stats))
# f_meta_stats = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.exclude_APOE_region.shared_var.hg19_1KG_EURintersected.txt")
# m_meta_stats = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/EU_all/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.exclude_APOE_region.shared_var.hg19_1KG_EURintersected.txt")
f_meta_stats = as.data.frame(f_meta_stats); m_meta_stats = as.data.frame(m_meta_stats)

f_AFR_stats = fread(paste0(opt$AFR_dir, opt$f_AFR_stats))
m_AFR_stats = fread(paste0(opt$AFR_dir, opt$m_AFR_stats))
# f_AFR_stats = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/AFR/AFRad_Females_case_control.full.hg19_AFR-ad_LDintersected.txt")
# m_AFR_stats = fread("/storage2/fs1/belloy2/Active/05_Projects/username/PWAS/AFR/AFRad_Males_case_control.full.hg19_AFR-ad_LDintersected.txt")
f_AFR_stats = as.data.frame(f_AFR_stats); m_AFR_stats = as.data.frame(m_AFR_stats)

#### FEMALES

# obtain sample N for females meta using the max N in the GWAS
f_meta$N.meta <- NA
f_meta$N.meta <- max(f_meta_stats$N_incl)

# obtain sample N for AFR females using the max N in the GWAS
f_meta$N.AFR <- NA
f_meta$N.AFR <- max(f_AFR_stats$N_incl)

# obtain PD.PWAS.Z for AFR females for matching gene
f_meta$AFR.Female.PD.PWAS.Z <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Female.PD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# obtain PD.PWAS.P for AFR females for matching gene
f_meta$AFR.Female.PD.PWAS.P <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Female.PD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# obtain SD.PWAS.Z for AFR females for matching gene
f_meta$AFR.Female.SD.PWAS.Z <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Female.SD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# obtain SD.PWAS.P for AFR females for matching gene
f_meta$AFR.Female.SD.PWAS.P <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Female.SD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# calculate combined PD PWAS Z-scores and p-values
f_meta$comb.Female.PD.PWAS.Z = (sqrt(f_meta$N.meta)*f_meta$Female.PD.PWAS.Z + sqrt(f_meta$N.AFR)*f_meta$AFR.Female.PD.PWAS.Z) / sqrt(f_meta$N.meta + f_meta$N.AFR)
f_meta$comb.Female.PD.PWAS.P = 2*pnorm(abs(f_meta$comb.Female.PD.PWAS.Z), lower.tail = FALSE)

# calculate combined SD PWAS Z-scores and p-values
f_meta$comb.Female.SD.PWAS.Z = (sqrt(f_meta$N.meta)*f_meta$Female.SD.PWAS.Z + sqrt(f_meta$N.AFR)*f_meta$AFR.Female.SD.PWAS.Z) / sqrt(f_meta$N.meta + f_meta$N.AFR)
f_meta$comb.Female.SD.PWAS.P = 2*pnorm(abs(f_meta$comb.Female.SD.PWAS.Z), lower.tail = FALSE)

## add male AFR info
# obtain PD.PWAS.Z for AFR Males for matching gene
f_meta$AFR.Male.PD.PWAS.Z <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Male.PD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# obtain PD.PWAS.P for AFR Males for matching gene
f_meta$AFR.Male.PD.PWAS.P <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Male.PD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# obtain SD.PWAS.Z for AFR Males for matching gene
f_meta$AFR.Male.SD.PWAS.Z <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Male.SD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene)

# obtain SD.PWAS.P for AFR Males for matching gene
f_meta$AFR.Male.SD.PWAS.P <- mapply(function(gene) {
    idx <- which(f_AFR$Gene == gene)
    if (length(idx) >0) {
        return(f_AFR$Male.SD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, f_meta$Gene) 

# obtain sample N for males meta using the max N in the GWAS
f_meta$Male.N.meta <- NA
f_meta$Male.N.meta <- max(m_meta_stats$N_incl)

# obtain sample N for AFR males using the max N in the GWAS
f_meta$Male.N.AFR <- NA
f_meta$Male.N.AFR <- max(m_AFR_stats$N_incl)

# calculate combined PD PWAS Z scores and p-values for males
f_meta$comb.Male.PD.PWAS.Z = (sqrt(f_meta$Male.N.meta)*f_meta$Male.PD.PWAS.Z + sqrt(f_meta$Male.N.AFR)*f_meta$AFR.Male.PD.PWAS.Z) / sqrt(f_meta$Male.N.meta + f_meta$Male.N.AFR)
f_meta$comb.Male.PD.PWAS.P = 2*pnorm(abs(f_meta$comb.Male.PD.PWAS.Z), lower.tail = FALSE)

#calculate combined SD PWAS Z scores and p-values for males
f_meta$comb.Male.SD.PWAS.Z = (sqrt(f_meta$Male.N.meta)*f_meta$Male.SD.PWAS.Z + sqrt(f_meta$Male.N.AFR)*f_meta$AFR.Male.SD.PWAS.Z) / sqrt(f_meta$Male.N.meta + f_meta$Male.N.AFR)
f_meta$comb.Male.SD.PWAS.P = 2*pnorm(abs(f_meta$comb.Male.SD.PWAS.Z), lower.tail = FALSE)

# save
fwrite(f_meta, paste0(opt$EUall_dir, f_dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_fsg_sex-het_comb-AFRad_full-Z.txt"), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, na = NA)

## MALES

# obtain sample N for males meta using the max N in the GWAS
m_meta$N.meta <- NA
m_meta$N.meta <- max(m_meta_stats$N_incl)

# obtain sample N for AFR males using the max N in the GWAS
m_meta$N.AFR <- NA
m_meta$N.AFR <- max(m_AFR_stats$N_incl)

# obtain PD.PWAS.Z for AFR males for matching gene
m_meta$AFR.Male.PD.PWAS.Z <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Male.PD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

#obtain PD.PWAS.P for AFR males for matching gene
m_meta$AFR.Male.PD.PWAS.P <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Male.PD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

# obtain SD.PWAS.Z for AFR males for matching gene
m_meta$AFR.Male.SD.PWAS.Z <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Male.SD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

# obtain SD.PWAS.P for AFR males for matching gene
m_meta$AFR.Male.SD.PWAS.P <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Male.SD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

# calculate combined PD PWAS Z-scores and p-values
m_meta$comb.Male.PD.PWAS.Z = (sqrt(m_meta$N.meta)*m_meta$Male.PD.PWAS.Z + sqrt(m_meta$N.AFR)*m_meta$AFR.Male.PD.PWAS.Z) / sqrt(m_meta$N.meta + m_meta$N.AFR)
m_meta$comb.Male.PD.PWAS.P = 2*pnorm(abs(m_meta$comb.Male.PD.PWAS.Z), lower.tail = FALSE)

# calculate combined SD PWAS Z-scores and p-values
m_meta$comb.Male.SD.PWAS.Z = (sqrt(m_meta$N.meta)*m_meta$Male.SD.PWAS.Z + sqrt(m_meta$N.AFR)*m_meta$AFR.Male.SD.PWAS.Z) / sqrt(m_meta$N.meta + m_meta$N.AFR)
m_meta$comb.Male.SD.PWAS.P = 2*pnorm(abs(m_meta$comb.Male.SD.PWAS.Z), lower.tail = FALSE)

## add female AFR info
# obtain PD.PWAS.Z for AFR Females for matching gene
m_meta$AFR.Female.PD.PWAS.Z <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Female.PD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

#obtain PD.PWAS.P for AFR Females for matching gene
m_meta$AFR.Female.PD.PWAS.P <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Female.PD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

# obtain SD.PWAS.Z for AFR Females for matching gene
m_meta$AFR.Female.SD.PWAS.Z <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Female.SD.PWAS.Z[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

# obtain SD.PWAS.P for AFR Females for matching gene
m_meta$AFR.Female.SD.PWAS.P <- mapply(function(gene) {
    idx <- which(m_AFR$Gene == gene)
    if (length(idx) >0) {
        return(m_AFR$Female.SD.PWAS.P[idx])
    } else {
       return(NA)
    }
}, m_meta$Gene)

# obtain sample N for females meta using the max N in the GWAS
m_meta$Female.N.meta <- NA
m_meta$Female.N.meta <- max(f_meta_stats$N_incl)

#obtain sample N for AFR females using the max N in the GWAS
m_meta$Female.N.AFR <- NA
m_meta$Female.N.AFR <- max(f_AFR_stats$N_incl)

# calculate combined PD PWAS z scores and p-values for females
m_meta$comb.Female.PD.PWAS.Z = (sqrt(m_meta$Female.N.meta)*m_meta$Female.PD.PWAS.Z + sqrt(m_meta$Female.N.AFR)*m_meta$AFR.Female.PD.PWAS.Z) / sqrt(m_meta$Female.N.meta + m_meta$Female.N.AFR)
m_meta$comb.Female.PD.PWAS.P = 2*pnorm(abs(m_meta$comb.Female.PD.PWAS.Z), lower.tail = FALSE)

# calculate combined SD PWAS z scores and p-values for females
m_meta$comb.Female.SD.PWAS.Z = (sqrt(m_meta$Female.N.meta)*m_meta$Female.SD.PWAS.Z + sqrt(m_meta$Female.N.AFR)*m_meta$AFR.Female.SD.PWAS.Z) / sqrt(m_meta$Female.N.meta + m_meta$Female.N.AFR)
m_meta$comb.Female.SD.PWAS.P = 2*pnorm(abs(m_meta$comb.Female.SD.PWAS.Z), lower.tail = FALSE)
	
# save
fwrite(m_meta, paste0(opt$EUall_dir, m_dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_eurLD_PWAS_W23_msg_sex-het_comb-AFRad_full-Z.txt"), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, na = NA)