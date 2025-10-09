## Create Table2 and TableS18 for PWAS paper
## Take Sex-specific results from EUR and make comparisons with AFR output

library(data.table)
library(tidyverse)


### Read in EUR Results - sent from Sathesh
## Female
# Brain Primary
f1b = fread("EU_all/Brain/Female/Primary/EUall_Female_Brain_Primary_sex_specific_gene_list.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Primary",
         discovery_sex = "Female")

# Brain Secondary
f2b = fread("EU_all/Brain/Female/Secondary/EUall_Female_Brain_Secondary_sex_specific_gene_list.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Secondary",
         discovery_sex = "Female")

# CSF Primary
f1c = fread("EU_all/CSF/Female/Primary/EUall_Female_CSF_Primary_sex_specific_gene_list.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Primary",
         discovery_sex = "Female") 

# CSF Secondary
f2c = fread("EU_all/CSF/Female/Secondary/EUall_Female_CSF_Secondary_sex_specific_gene_list.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Secondary",
         discovery_sex = "Female")


## Male
# Brain Primary - NA
# m1b = fread("EU_all/Brain/Male/Primary/EUall_Male_Brain_Primary_sex_specific_gene_list.txt")

# Brain Secondary
m2b = fread("EU_all/Brain/Male/Secondary/EUall_Male_Brain_Secondary_sex_specific_gene_list.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Secondary",
         discovery_sex = "Male")

# CSF Primary - NA
# m1c = fread("EU_all/CSF/Male/Primary/EUall_Male_CSF_Primary_sex_specific_gene_list.txt")

# CSF Secondary
m2c = fread("EU_all/CSF/Male/Secondary/EUall_Male_CSF_Secondary_sex_specific_gene_list.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Secondary",
         discovery_sex = "Male")


### Combine EUR results
eur = rbind(f1b, f2b, f1c, f2c, m2b, m2c) %>% 
  group_by(GENE, ID) %>% 
  mutate(Discovery = ifelse(n() > 1, "Prim+Sec", discovery)) %>% 
  filter(FDR == min(FDR)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  filter(GENE != "PACSIN1", GENE != "CEBPZOS") %>% 
  arrange(CHR, P0)


## Need male Secondary discovery results for YARS1, SCAPER
m2b_ref = fread("EU_all/Brain/Male/Secondary/EUall_Male_Brain_Secondary_merged_with_FDR.txt") %>% 
  dplyr::select(ID, PWAS.Z, PWAS.P, FDR) %>% 
  rename(ref_PWAS.Z = PWAS.Z,
         ref_PWAS.P = PWAS.P,
         ref_PWAS.FDR = FDR)

## Need male mismatach results for INPP5D
# m3b_ref = fread('EU_all/Brain/Male/Opposite/ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.hg19_EUR.Wingo2023WEIGHTS-Females.2.dat') %>% 
m3b_ref = fread('EU_all/Brain/Male/Opposite/EUall.cleaned.brain.male.BrainFemalePW.2.dat') %>% 
  mutate(FDR = p.adjust(TWAS.P, method = "BH")) %>% 
  dplyr::select(ID, TWAS.Z, TWAS.P, FDR) %>% 
  rename(ref_PWAS.Z = TWAS.Z,
         ref_PWAS.P = TWAS.P,
         ref_PWAS.FDR = FDR) %>% 
  filter(ID == "ENSG00000168918.INPP5D")

# Combine m2b_ref and m3b_ref
mb_other = rbind(m2b_ref, m3b_ref)


# Merge significant findings with other reuslts for opposite sex.
eur_full = left_join(eur, mb_other, by = "ID") %>%
  mutate(
    ref_PWAS.Z = if_else(!is.na(ref_PWAS.Z.x), ref_PWAS.Z.x, ref_PWAS.Z.y),
    ref_PWAS.P = if_else(!is.na(ref_PWAS.P.x), ref_PWAS.P.x, ref_PWAS.P.y),
    ref_PWAS.FDR = if_else(!is.na(ref_PWAS.FDR.x), ref_PWAS.FDR.x, ref_PWAS.FDR.y)
  ) %>%
  select(-signs, -ref_PWAS.Z.x, -ref_PWAS.Z.y, -ref_PWAS.P.x, -ref_PWAS.P.y, -ref_PWAS.FDR.x, -ref_PWAS.FDR.y)

####################################################################################################################################################################
####################################################################################################################################################################

### Read in AFR samples
## Female
# Brain Primary
afr_f1b = fread("AFR/Brain/Female/Primary/AFR_Female_Brain_Primary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Primary") %>% 
  rename(AFR_female_Z = PWAS.Z,
         AFR_female_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_female_Z, AFR_female_P)

# Brain Secondary
afr_f2b = fread("AFR/Brain/Female/Secondary/AFR_Female_Brain_Secondary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Secondary") %>% 
  rename(AFR_female_Z = PWAS.Z,
         AFR_female_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_female_Z, AFR_female_P)

# CSF Primary
afr_f1c = fread("AFR/CSF/Female/Primary/AFR_Female_CSF_Primary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Primary") %>% 
  rename(AFR_female_Z = PWAS.Z,
         AFR_female_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_female_Z, AFR_female_P)

# CSF Secondary
afr_f2c = fread("AFR/CSF/Female/Secondary/AFR_Female_CSF_Secondary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Secondary") %>% 
  rename(AFR_female_Z = PWAS.Z,
         AFR_female_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_female_Z, AFR_female_P)

# rbind female samples
afr_female = rbind(afr_f1b, afr_f2b, afr_f1c, afr_f2c)


## Male
# Brain Primary
afr_m1b = fread("AFR/Brain/Male/Primary/AFR_Male_Brain_Primary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Primary") %>% 
  rename(AFR_male_Z = PWAS.Z,
         AFR_male_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_male_Z, AFR_male_P)

# Brain Secondary
afr_m2b = fread("AFR/Brain/Male/Secondary/AFR_Male_Brain_Secondary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "Brain",
         discovery = "Secondary") %>% 
  rename(AFR_male_Z = PWAS.Z,
         AFR_male_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_male_Z, AFR_male_P)

# CSF Primary
afr_m1c = fread("AFR/CSF/Male/Primary/AFR_Male_CSF_Primary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Primary") %>% 
  rename(AFR_male_Z = PWAS.Z,
         AFR_male_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_male_Z, AFR_male_P)

# CSF Secondary
afr_m2c = fread("AFR/CSF/Male/Secondary/AFR_Male_CSF_Secondary_merged_with_FDR.txt") %>% 
  mutate(Tissue = "CSF",
         discovery = "Secondary") %>% 
  rename(AFR_male_Z = PWAS.Z,
         AFR_male_P = PWAS.P) %>% 
  dplyr::select(ID, Tissue, discovery, AFR_male_Z, AFR_male_P)

# rbind male samples
afr_male = rbind(afr_m1b, afr_m2b, afr_m1c, afr_m2c)

####################################################################################################################################################################
####################################################################################################################################################################

## Merge EUR_full with afr female by ID, tissue, and discovery
eur_m1 = left_join(eur_full, afr_female, by = c("ID", "Tissue", "discovery"))

## Merge EUR_m1 with afr male by ID, tissue, and discovery
eur_m2 = left_join(eur_m1, afr_male, by = c("ID", "Tissue", "discovery"))


## Need male Secondary discovery results for YARS1, SCAPER
afr_m2b_ref = fread("AFR/Brain/Male/Secondary/AFR_Male_Brain_Secondary_merged_with_FDR.txt") %>% 
  dplyr::select(ID, PWAS.Z, PWAS.P) %>% 
  rename(AFR_male_Z = PWAS.Z,
         AFR_male_P = PWAS.P)

## Need male mismatach results for INPP5D
#afr_m3b_ref = fread('AFR/Brain/Male/Opposite/AFRad_Males_case_control.full.hg19_AFR-ad_LDintersected.W23.female-WEIGHTS.2.dat') %>% 
afr_m3b_ref = fread('AFR/Brain/Male/Opposite/AFR.cleaned.brain.male.BrainFemalePW.2.dat') %>% 
  dplyr::select(ID, TWAS.Z, TWAS.P) %>% 
  rename(AFR_male_Z = TWAS.Z,
         AFR_male_P = TWAS.P)%>% 
  filter(ID == "ENSG00000168918.INPP5D")

# Combine m2b_ref and m3b_ref
afr_mb_other = rbind(afr_m2b_ref, afr_m3b_ref)


# Merge significant findings with other reuslts for opposite sex.
eur_m3 = left_join(eur_m2, afr_mb_other, by = "ID") %>%
  mutate(
    AFR_male_Z = if_else(!is.na(AFR_male_Z.x), AFR_male_Z.x, AFR_male_Z.y),
    AFR_male_P = if_else(!is.na(AFR_male_P.x), AFR_male_P.x, AFR_male_P.y)) %>%
  select(-AFR_male_Z.x, -AFR_male_Z.y, -AFR_male_P.x, -AFR_male_P.y)


####################################################################################################################################################################
####################################################################################################################################################################

# Read in AD risk loci consensus file for novelty status
nov = fread("SupportingFiles/AD_Risk_Loci_consensus_2024.txt")


# Calculate Update columns for table2 - get novelty status
all1 = eur_m3 %>% 
  mutate(EUR_female_Z = ifelse(discovery_sex == "Female", PWAS.Z, ref_PWAS.Z),
         EUR_female_P = ifelse(discovery_sex == "Female", PWAS.P, ref_PWAS.P),
         EUR_female_P.FDR = ifelse(discovery_sex == "Female", FDR, ref_PWAS.FDR),
         EUR_male_Z = ifelse(discovery_sex == "Male", PWAS.Z, ref_PWAS.Z),
         EUR_male_P = ifelse(discovery_sex == "Male", PWAS.P, ref_PWAS.P),
         EUR_male_P.FDR = ifelse(discovery_sex == "Male", FDR, ref_PWAS.FDR)) %>% 
  rowwise() %>%
  mutate(Novelty = case_when(GENE %in% nov$Locus_consensus ~ "Known locus",
                             any(nov$CHR == CHR & nov$GRCh38_POS >= (P0 - 1e6) & nov$GRCh38_POS <= (P1 + 1e6)) ~ "Known locus",
                             T ~ "Novel locus")) %>% 
  ungroup() %>%
  as.data.frame() %>% 
  rename(Gene = GENE, Start = P0, End = P1) %>% 
  dplyr::select(Gene, CHR, Start, End, Novelty, Discovery, Tissue, discovery_sex, EUR_female_Z, EUR_female_P, EUR_female_P.FDR, EUR_male_Z, EUR_male_P, EUR_male_P.FDR,
                AFR_female_Z, AFR_female_P, AFR_male_Z, AFR_male_P)


# Calculate weighted z-score and p-value for EUR-AFR Meta
all2 = all1 %>% 
  mutate(EUR_AFR_female_Z = (sqrt(636155)*EUR_female_Z + sqrt(5152)*AFR_female_Z) / sqrt(636155+5152),
         EUR_AFR_female_P = 2 * (1 - pnorm(abs(EUR_AFR_female_Z))),
         EUR_AFR_male_Z = (sqrt(507678)*EUR_male_Z + sqrt(2118)*AFR_male_Z) / sqrt(507678+2118),
         EUR_AFR_male_P = 2 * (1 - pnorm(abs(EUR_AFR_male_Z)))) %>% 
  mutate(AFR_sex_het = case_when(discovery_sex == "Female" & EUR_AFR_female_P < EUR_female_P & EUR_AFR_male_P > 0.05 ~ "Yes",
                                 discovery_sex == "Male" & EUR_AFR_male_P < EUR_male_P & EUR_AFR_female_P > 0.05 ~ "Yes",
                                 T ~ "No")) %>% 
  arrange(discovery_sex, CHR, Start) %>% 
  dplyr::select(-discovery_sex)

# Write out final table S18    
fwrite(all2, "~/Desktop/TableS18_new.csv")


# Make table2
all3 = all2 %>% 
  rename(Protein = Gene, CH = CHR, `Gene start` = Start, `Gene end` = End, `Novelty status` = Novelty,
         Female_Z = EUR_female_Z, Female_P = EUR_female_P, Female_FDR_P = EUR_female_P.FDR,
         Male_Z = EUR_male_Z, Male_P = EUR_male_P, Male_FDR_P = EUR_male_P.FDR, `Sex-het. consistent in AFR?` = AFR_sex_het) %>% 
  dplyr::select(Protein, CH, `Gene start`, `Gene end`, `Novelty status`, Discovery, Tissue, Female_Z, Female_P, Female_FDR_P, Male_Z, Male_P, Male_FDR_P, `Sex-het. consistent in AFR?`)
  
# Write out final table 2  
fwrite(all3, "~/Desktop/Table2_new.csv")



## End of script