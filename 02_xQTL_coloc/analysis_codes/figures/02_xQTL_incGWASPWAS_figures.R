## xQTL figure for Sex-stratified GWAS and PWAS Paper
## Combine GWAS and PWAS xQTL results in one figure (Manuscript Figure 3)
# August 5, 2025

library(data.table)
library(tidyverse)
library(patchwork)
library(grid)
library(reshape2)
library(gridExtra)
library(scales)
library(optparse)
option_list <- list(
  make_option("--work_dir", type = "character", help = "Working directory"),
  make_option("--GP_qtl_in", type = "character", help = "Annotated combined PWAS/GWAS CSV"),
  make_option("--mqtl_in", type = "character", help = "Top mQTL CSV"),
  make_option("--haqtl_in", type = "character", help = "Top haQTL CSV"),
  make_option("--caqtl_in", type = "character", help = "Top caQTL CSV"),
  make_option("--gnames_in", type = "character", help = "Updated gene names CSV"),
  make_option("--index_in", type = "character", help = "Index genes CSV"),
  make_option("--plot_out", type = "character", help = "Output figure PNG")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$work_dir)) stop("--work_dir is required")
dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)


##########################################################################################################################################

## Read in appended PWAS data
annotated = fread(file.path(opt$work_dir, "results", opt$annotated_in)) 

## Filter out m/ca/ha data as this will be read back in later with the correct gene names
all = annotated %>% 
  dplyr::filter(qtl_type != "mQTL",
                qtl_type != "haQTL",
                qtl_type != "caQTL") %>% 
  arrange(locus)


## Read in m/ca/ha QTL data with correct names
# Manually updated gene_name column in m/ca/ha QTL files to match last case of gene at locus to align with top of the locus in figure.
mqtl = fread(file.path(opt$work_dir, "results", opt$mqtl_in))
haqtl = fread(file.path(opt$work_dir, "results", opt$haqtl_in))
caqtl = fread(file.path(opt$work_dir, "results", opt$caqtl_in))

tmp = rbind(mqtl, haqtl, caqtl); tmp$molecular_trait_id <- tmp$gene_name
gnames <- fread(file.path(opt$work_dir, "results", opt$gnames_in))

gnames_map <- gnames %>%
  filter(!is.na(gene_name), !is.na(Updated_gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  select(gene_name, Updated_gene_name)

tmp_updated <- tmp %>%
  left_join(gnames_map, by = c("molecular_trait_id" = "gene_name")) %>%
  mutate(gene_name = coalesce(Updated_gene_name, gene_name)) %>%
  select(-Updated_gene_name)


# Combine all qtl data
all1 = rbind(all, tmp_updated) %>% 
  arrange(CHR, BP)

unique(all1$locus)
n_distinct(all1$locus) # N = 56


##########################################################################################################################################
# Create Locus Labels for overview plot - rename genes with aliases - filter out loci not used and genes that did not pass visual inspection
all2 = all1 %>%  
  dplyr::filter(locus != "CEBPZOS",
                locus != "CLIC1",
                locus != "PACSIN1",
                locus != "HLA-DQA2") %>% 
  mutate(stratum_label = case_when(ad_sample == "Female" ~ "F",
                                   ad_sample == "Male" ~ "M"),
         locus_name = case_when(locus == "YARS1" ~ "YARS",
                                locus == "TDRKH" | locus == "MINDY1" ~ "TDRKH/MINDY1",
                                locus == "GBA" ~ "KRCTAP2/GBA1",
                                locus == "CD84" ~ "CD84 (ADAMTS4)",
                                locus == "WDFY1" ~ "WDFY1-AP1S3",
                                locus == "USP4" | locus == "USP19" ~ "USP4/USP19",
                                locus == "CLNK" ~ "CLNK/HS3ST1",
                                locus == "CLU" ~ "CLU/PTK2B",
                                locus == "ANXA11" | locus == "TSPAN14" ~ "TSPAN14/ANXA11",
                                locus == "SLC24A4" ~ "SLC24A4/RIN3",
                                locus == "ETFA" | locus == "SCAPER" | locus == "RCN2" ~ "ETFA/SCAPER/RCN2",
                                locus == "HP" ~ "PMFBP1/HP",
                                locus == "ENO3" | locus == "RABEP1" ~ "ENO3 (SCIMP/RABEP1/CHRNE)",
                                locus == "TOM1L2" ~ "TOM1L2 (MYO15A)",
                                locus == "NSF" ~ "NSF/MAPT/KANSL1",
                                locus == "PHB1" ~ "PHB1 (ABI3)",
                                locus == "SLC44A2" | locus == "CARM1" ~ "SLC44A2/CARM1",
                                locus == "ADAMTS1" ~ "APP/ADAMTS1",
                                T ~ locus),
         gene_name = case_when(gene_name == "GBA" ~ "GBA1",
                               gene_name == "YARS1" ~ "YARS",
                               T ~ gene_name),
         locus_label = case_when(locus_name == "KRCTAP2/GBA1" ~ "KRCTAP2/GBA1 (M, GWAS & CSF)",
                                 locus_name == "INPP5D" ~ "INPP5D (F, GWAS & Brain)",
                                 locus_name == "TSPAN14/ANXA11" ~ "TSPAN14/ANXA11 (F, GWAS & Brain)",
                                 locus_name == "PMFBP1/HP" ~ "PMFBP1/HP (F, GWAS & CSF)",
                                 locus_name == "ENO3 (SCIMP/RABEP1/CHRNE)" ~ "ENO3 (SCIMP/RABEP1/CHRNE) (F, CSF & Brain)",
                                 T ~ paste0(locus_name, " (", stratum_label, ", ", discovery, ")"))) %>% 
  filter(!(locus == "HP" & gene_name == "HP" & dset == "BrainMeta_DLPFC_chr_16_eQTL_RB_Female"),
         !(locus == "HP" & gene_name == "DHODH" & dset == "eQTL_catalogue_commonmind_dlpfc_chr_16_eQTL_RB_Female"),
         !(locus == "HP" & gene_name == "HPR" & dset == "eQTL_catalogue_tibial_nerve_chr_16_eQTL_RB_Female"),
         !(locus == "HP" & gene_name == "HPR" & dset == "BrainMeta_DLPFC_chr_16_eQTL_RB_Female"),
         !(locus == "NCK2" & gene_name == "NCK2" & dset == "eQTLgen_blood_chr_2_eQTL_RB_Female"),
         !(locus == "GBA" & gene_name == "MSTO1"),
         !(locus == "GBA" & gene_name == "FDPS"),
         !(locus == "GBA" & gene_name == "RIT1"),
         !(locus == "CLU" & gene_name == "SCARA3")) %>% 
  arrange(CHR, BP) 


n_distinct(all2$locus_name) # N = 45
n_distinct(all2$locus_label) # N = 45


######################################################################################################################################################################################################################
## Get Best PP4 for locus with GWAS and PWAS support - Create PP4 value which is rounded down for the figure - create locus gene name used to intersect index df
all3 = all2 %>% 
  group_by(locus_label, gene_name, dset) %>% 
  filter(best_PP4 == max(best_PP4)) %>% 
  ungroup() %>% 
  dplyr::mutate(PP4 = round(best_PP4, 2),
                locus_gene = paste(locus_label, gene_name, sep = "-:-")) 

n_distinct(all3$locus_name) # N = 45
n_distinct(all3$locus_label) # N = 45


################################################################################################
## Insert gene and locus lists - filter
index = fread(file.path(opt$work_dir, "results", opt$index_in))

index = index %>% 
  mutate(locus_gene = paste(locus_index, genes, sep = "-:-"))
n_distinct(index$locus_index)


## Combine index file with appended file to only include prioritized genes
all4 = inner_join(all3, index, by = "locus_gene") %>% 
  dplyr::mutate(PP4 = round(best_PP4, 2),
                tissue_type = case_when(tissue_type == "Kosoy 2020 Microglia" ~ "Kosoy 2022 Microglia",
                                        T ~ tissue_type)) %>% 
  arrange(tissue_type)

n_distinct(all4$locus_name) # N = 30



######################################################################################################################################################################################################################
## Create cell_type variable for color in overview plot
unique(all4$tissue_type)
n_distinct(all4$tissue_type)

all4_v1 = all4 %>% 
  mutate(cell_type = case_when(
    tissue_type == unique(all4$tissue_type)[1] ~ "Blood/Plasma",
    tissue_type == unique(all4$tissue_type)[2] ~ "Single Cell",
    tissue_type %in% unique(all4$tissue_type)[3:4] ~ "Brain",
    tissue_type == unique(all4$tissue_type)[5] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[6] ~ "Brain",
    tissue_type == unique(all4$tissue_type)[7] ~ "Single Cell",
    tissue_type %in% unique(all4$tissue_type)[8:11] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[12] ~ "Single Cell",
    tissue_type %in% unique(all4$tissue_type)[13:14] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[15] ~ "Blood/Plasma",
    tissue_type == unique(all4$tissue_type)[16] ~ "Brain",
    tissue_type == unique(all4$tissue_type)[17] ~ "Other Tissue",
    tissue_type == unique(all4$tissue_type)[18] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[19] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[20] ~ "Brain",
    tissue_type %in% unique(all4$tissue_type)[21:23] ~ "CSF",
    tissue_type %in% unique(all4$tissue_type)[24:25] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[26] ~ "Blood/Plasma",
    tissue_type %in% unique(all4$tissue_type)[27:29] ~ "Brain",
    tissue_type == unique(all4$tissue_type)[30] ~ "Single Cell",
    tissue_type == unique(all4$tissue_type)[31] ~ "Blood/Plasma",
    tissue_type == unique(all4$tissue_type)[32] ~ "Brain")) %>% 
  mutate(tissue_qtl = paste(tissue_type, qtl_type, sep = " "),
         qtl_type = case_when(qtl_type == "mQTL" ~ "m/ca/haQTL",
                              qtl_type == "haQTL" ~ "m/ca/haQTL",
                              qtl_type == "caQTL" ~ "m/ca/haQTL",
                              T ~ qtl_type)) %>% 
  group_by(locus_label, gene_name, tissue_qtl) %>% 
  dplyr::filter(PP4 == max(PP4)) %>% 
  arrange(tissue_qtl)


######################################################################################################################################################################################################################
## Create collapsed tissue_qtl variable for x-axis
unique(all4_v1$tissue_qtl)
n_distinct(all4_v1$tissue_qtl)

all4_v2 = all4_v1 %>% 
  mutate(tissue_qtl_c = case_when(
    tissue_qtl == unique(all4_v1$tissue_qtl)[1] ~ "Blood/Plasma pQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[2] ~ "Monocytes/Macrophages eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[3] ~ "Monocytes/Macrophages sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[4] ~ "Brain m/ha/caQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[5] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[6] ~ "Brain m/ha/caQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[7] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[8] ~ "Brain sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[9] ~ "Monocytes/Macrophages eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[10] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[11] ~ "Brain sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[12] ~ "Monocytes/Macrophages eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[13] ~ "Astrocytes eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[14] ~ "Endothelial Cells eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[15] ~ "Excitatory Neurons eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[16] ~ "Inhibitory Neurons eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[17] ~ "Microglia eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[18] ~ "OPCs eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[19] ~ "Oligodendrocytes eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[20] ~ "Blood/Plasma eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[21] ~ "Blood/Plasma sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[22] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[23] ~ "Brain sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[24] ~ "Other Tissue eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[25] ~ "Other Tissue sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[26] ~ "T-Cells eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[27] ~ "Microglia caQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[28] ~ "Microglia eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[29] ~ "Brain eQTL",
    tissue_qtl %in% unique(all4_v1$tissue_qtl)[30:32] ~ "PWAS CSF pQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[33] ~ "Monocytes/Macrophages eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[34] ~ "Monocytes/Macrophages sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[35] ~ "Monocytes/Macrophages eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[36] ~ "Monocytes/Macrophages sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[37] ~ "Blood/Plasma pQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[38] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[39] ~ "PWAS Brain pQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[40] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[41] ~ "PWAS Brain pQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[42] ~ "Brain eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[43] ~ "PWAS Brain pQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[44] ~ "Microglia eQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[45] ~ "Microglia sQTL",
    tissue_qtl == unique(all4_v1$tissue_qtl)[46] ~ "Blood/Plasma eQTL",
    tissue_qtl %in% unique(all4_v1$tissue_qtl)[47:48] ~ "Brain m/ha/caQTL")) %>% 
  group_by(locus_label, gene_name, tissue_qtl_c) %>% 
  dplyr::filter(PP4 == max(PP4)) %>% 
  arrange(locus_label)



######################################################################################################################################################################################################################
# Filter out irrelevant QTL groups
all4_v2 = all4_v2 %>%
  filter(tissue_qtl_c != "Endothelial Cells eQTL" & tissue_qtl_c != "Astrocytes eQTL" & tissue_qtl_c != "Microglia sQTL")

################################################################################################
# Cutomize order of x-axis (QTL dataset)
unique(all4_v2$tissue_qtl_c)
custom_order_tissue_qtl = c("PWAS Brain pQTL", "Brain eQTL", "Brain sQTL", "Brain m/ha/caQTL", "PWAS CSF pQTL", "Blood/Plasma pQTL", "Blood/Plasma eQTL", "Blood/Plasma sQTL",
                            "Other Tissue eQTL", "Other Tissue sQTL", "T-Cells eQTL", "Monocytes/Macrophages eQTL", "Monocytes/Macrophages sQTL",
                            "Microglia eQTL",   "Microglia caQTL",  "Excitatory Neurons eQTL", "Inhibitory Neurons eQTL", 
                            "Oligodendrocytes eQTL", "OPCs eQTL")

# Removed: "Microglia sQTL", "Astrocytes eQTL", "Endothelial Cells eQTL",
missing_items1 <- setdiff(all4_v2$tissue_qtl_c, custom_order_tissue_qtl)
missing_items2 <- setdiff(custom_order_tissue_qtl, all4_v2$tissue_qtl_c)

all4_v2$tissue_qtl_c <- factor(all4_v2$tissue_qtl_c, levels = custom_order_tissue_qtl)


################################################################################################
# Customize order of locus_labels
unique(all4_v2$locus_label)

custom_order_locus = c("ALPL (F, GWAS)", "KRCTAP2/GBA1 (M, GWAS & CSF)", "CD84 (ADAMTS4) (F, CSF)", "NCK2 (F, GWAS)", "MFSD6 (F, Brain)", "WDFY1-AP1S3 (M, GWAS)",
                       "INPP5D (F, GWAS & Brain)", "RFTN1 (M, Brain)", "USP4/USP19 (F, Brain)", "MANF (F, CSF)", "NCK1 (F, Brain)", "METTL14 (F, Brain)",
                       "RAPGEF2 (F, Brain)", "ACSL6 (F, Brain)", "CARMIL1 (M, Brain)", "CDK14 (F, Brain)", "RASA4B (F, Brain)", "SHC3 (F, GWAS)", "SFXN4 (F, Brain)",
                       "PSEN1 (F, Brain)", "ETFA/SCAPER/RCN2 (F, Brain)", "PMFBP1/HP (F, GWAS & CSF)", "SLC7A5 (F, Brain)", "ENO3 (SCIMP/RABEP1/CHRNE) (F, CSF & Brain)", 
                       "TOM1L2 (MYO15A) (F, Brain)", "PHB1 (ABI3) (M, Brain)", "SLC44A2/CARM1 (F, Brain)", "RBCK1 (F, GWAS)", "APP/ADAMTS1 (M, GWAS)", "ARSA (F, Brain)")


missing_items1 <- setdiff(all4_v2$locus_label, custom_order_locus)
missing_items2 <- setdiff(custom_order_locus, all4_v2$locus_label)

all4_v2$locus_label <- factor(all4_v2$locus_label, levels = custom_order_locus)
n_distinct(all4_v2$locus_label) # N = 30


################################################################################################
# Customize order of qtl_types
unique(all4_v2$qtl_type)

custom_order_qtl_type = c("eQTL", "pQTL", "sQTL", "m/ca/haQTL")


missing_items1 <- setdiff(all4_v2$qtl_type, custom_order_qtl_type)
missing_items2 <- setdiff(custom_order_qtl_type, all4_v2$qtl_type)

all4_v2$qtl_type <- factor(all4_v2$qtl_type, levels = custom_order_qtl_type)


######################################################################################################################################################################################################################
# Define the color mapping for each cell_type
color_map <- c(
  "Blood/Plasma" = "#FF69B4",
  "CSF" = "yellow1",
  "Other Tissue" = "#DDA0DD",
  "Single Cell" = "#00CED1",
  "Brain" = "#FFCC66"
)


# Precompute the fill colors based on PP4 and cell_type
all4_v2$fill_color <- sapply(1:nrow(all4_v2), function(i) {
  cell_type <- as.character(all4_v2$cell_type[i])
  pp4 <- all4_v2$PP4[i]
  # Get the base color for this cell_type, default to Brain if not found
  base_color <- color_map[ifelse(cell_type %in% names(color_map), cell_type, "Brain")]
  # Create a gradient from white to the base color
  ramp <- colour_ramp(c("#FFFFFF", base_color))
  # Apply the gradient based on PP4 value (0 to 1)
  ramp(pp4)
})

# Create Base Plot - Uses shapes - Not used in Manuscript Figure 3
p1 = ggplot(data = all4_v2, aes(y = gene_name, x = tissue_qtl_c, label = PP4)) +
  geom_point(aes(size = PP4, colour = cell_type, alpha = PP4, shape = qtl_type), position = position_nudge(x = -0.2)) + 
  geom_text(nudge_x = 0.2, size = 4) + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 0, face = "bold", colour = "black", hjust = 0), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
        axis.text.y = element_text(face = "italic", size = 12, colour = "black"),
        axis.title = element_text(size = 16),
        title = element_text(size = 18, face = "bold"),
        strip.background = element_blank(),
        legend.box.margin = margin(c(0,0,0,0)), 
        legend.margin = margin(c(0,0,0,0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        strip.placement = "outside",
        
        panel.spacing.y = unit(x = 0,units = "points"), 
        panel.border = element_rect(fill = NA, size = 0.4), 
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black"))+
  scale_size_continuous(limits = c(0,1)) +
  scale_alpha_continuous(limits = c(0,1)) +
  scale_x_discrete(position = "bottom") +
  labs(x = "", y = "", colour = "Tissue Type", shape = 'QTL Type', title = "") +
  facet_grid(locus_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 13)) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_colour_manual(values = color_map)

print(p1)


######################################################################################################################################################################################################################
## Version 2 - creating individual cells:

# Precompute fill_color2 for p2: full color if meets threshold, else white
all4_v2$fill_color2 <- sapply(1:nrow(all4_v2), function(i) {
  pp4 <- all4_v2$PP4[i]
  cell_type <- as.character(all4_v2$cell_type[i])
  tissue_qtl_c <- as.character(all4_v2$tissue_qtl_c[i])
  
  # Determine if it's a pQTL
  is_pqtl <- tissue_qtl_c %in% c("PWAS Brain pQTL", "PWAS CSF pQTL")
  
  # Set threshold based on QTL type
  threshold <- ifelse(is_pqtl, 0.4, 0.7)
  
  # Check condition: >= 0.4 for pQTL, > 0.7 for xQTL
  meets_threshold <- ifelse(is_pqtl, pp4 >= 0.4, pp4 >= 0.7)
  
  # Get the base color for this cell_type, default to Brain if not found
  base_color <- color_map[ifelse(cell_type %in% names(color_map), cell_type, "Brain")]
  
  # Apply color logic: full color if meets threshold, else white
  if (meets_threshold) {
    return(base_color)
  } else {
    return("white")
  }
})


# Define the base plot function to avoid code duplication
create_plot <- function(fill_col) {
  ggplot(data = all4_v2, aes(y = gene_name, x = tissue_qtl_c)) +
    geom_tile(aes(fill = .data[[fill_col]]), color = "black", width = 0.8, height = 0.8) +
    geom_text(aes(label = ifelse(fill_color2 != "white", sprintf("%.2f", PP4), "")), fontface = "bold", size = 3.5, color = "black") + # Use for text in cells with PP4 > 0.7
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold", colour = "black", hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = "black", face = "bold"),
      axis.text.y = element_text(face = "italic", size = 12, colour = "black"),
      axis.title = element_text(size = 16),
      title = element_text(size = 18, face = "bold"),
      strip.background = element_blank(),
      legend.box.margin = margin(c(0,0,0,0)),
      legend.margin = margin(c(0,0,0,0)),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      strip.placement = "outside",
      panel.spacing.y = unit(x = 0, units = "points"),
      panel.border = element_rect(fill = NA, size = 0.4),
      panel.grid = element_blank(),
      axis.ticks = element_line(colour = "black")
    ) +
    scale_x_discrete(position = "bottom") +
    labs(x = "", y = "", fill = "Tissue Type") +
    facet_grid(locus_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 13)) +
    scale_fill_identity()
}


# Create p2
p2 <- create_plot("fill_color2")

# Print the plots
print(p2)



################################################################################################################################################################
## Include priority score to the left of the plot

score_plot <- ggplot(data = all4_v2, aes(x = 1, y = gene_name, fill = factor(score), label = score)) +
  geom_tile(color = "white") +
  geom_text(size = 5, face = "bold") +
  scale_fill_manual(values = c("3" = "#d0f0c0", "2" = "#90ee90", "1" = "#32cd32")) +  # Custom shades of green
  labs(x = "Priority Score") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "italic", colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        axis.title.x = element_blank(),
        # axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, angle = 45, face = "bold"),
        axis.title.y = element_blank(),
        plot.margin = margin(l = 10, unit = "pt"),
        panel.spacing = unit(0, "points"),
        legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.y.right = element_text(angle = 0, face = "bold", size = 13, colour = "black")) +
  facet_grid(locus_label ~ ., scales = "free_y", space = "free_y")

print(score_plot)



################################################################################################################################################################
## Include DAL results for each gene

all4_v3 = all4_v2 %>% 
  dplyr::select(locus_label, gene_name, sex.Z, sex.sig, Z.women, women.sig, Z.men, men.sig, AD.sex_het.Zscore, sex_het.sig) %>% 
  group_by(gene_name) %>% 
  distinct(gene_name, .keep_all = TRUE) %>% 
  ungroup() %>% 
  dplyr::select(-tissue_qtl_c)
all4_v3 = as.data.frame(all4_v3)


################################################################################################
# Customize order of locus_labels - use exact same orde as for cell-plot
unique(all4_v3$locus_label)

custom_order_locus = c("ALPL (F, GWAS)", "KRCTAP2/GBA1 (M, GWAS & CSF)", "CD84 (ADAMTS4) (F, CSF)", "NCK2 (F, GWAS)", "MFSD6 (F, Brain)", "WDFY1-AP1S3 (M, GWAS)",
                       "INPP5D (F, GWAS & Brain)", "RFTN1 (M, Brain)", "USP4/USP19 (F, Brain)", "MANF (F, CSF)", "NCK1 (F, Brain)", "METTL14 (F, Brain)",
                       "RAPGEF2 (F, Brain)", "ACSL6 (F, Brain)", "CARMIL1 (M, Brain)", "CDK14 (F, Brain)", "RASA4B (F, Brain)", "SHC3 (F, GWAS)", "SFXN4 (F, Brain)",
                       "PSEN1 (F, Brain)", "ETFA/SCAPER/RCN2 (F, Brain)", "PMFBP1/HP (F, GWAS & CSF)", "SLC7A5 (F, Brain)", "ENO3 (SCIMP/RABEP1/CHRNE) (F, CSF & Brain)", 
                       "TOM1L2 (MYO15A) (F, Brain)", "PHB1 (ABI3) (M, Brain)", "SLC44A2/CARM1 (F, Brain)", "RBCK1 (F, GWAS)", "APP/ADAMTS1 (M, GWAS)", "ARSA (F, Brain)")


missing_items1 <- setdiff(all4_v3$locus_label, custom_order_locus)
missing_items2 <- setdiff(custom_order_locus, all4_v3$locus_label)

all4_v3$locus_label <- factor(all4_v3$locus_label, levels = custom_order_locus)


################################################################################################
## Reshape data to long format
df_long <- all4_v3 %>%
  select(locus_label, gene_name, 
         sex.Z, sex.sig, 
         Z.women, women.sig, 
         Z.men, men.sig, 
         AD.sex_het.Zscore, sex_het.sig) %>%
  pivot_longer(
    cols = c(sex.Z, Z.women, Z.men, AD.sex_het.Zscore),
    names_to = "variable",
    values_to = "Z_score"
  ) %>%
  pivot_longer(
    cols = c(sex.sig, women.sig, men.sig, sex_het.sig),
    names_to = "sig_variable",
    values_to = "sig"
  ) %>%
  # Filter to match Z-score and significance pairs
  filter(
    (variable == "sex.Z" & sig_variable == "sex.sig") |
      (variable == "Z.women" & sig_variable == "women.sig") |
      (variable == "Z.men" & sig_variable == "men.sig") |
      (variable == "AD.sex_het.Zscore" & sig_variable == "sex_het.sig")
  ) %>%
  # Replace NA significance with empty string
  mutate(sig = ifelse(is.na(sig), "", sig)) %>%
  # Order variables for plotting
  mutate(variable = factor(variable, levels = c("sex.Z", "Z.women", "Z.men", "AD.sex_het.Zscore")))



################################################################################################################################################################################################
# Create the plot
d1 = ggplot(df_long, aes(x = variable, y = gene_name, fill = Z_score, label = sig)) +
  geom_tile(color = "gray", size = 0.5, width = 0.8, height = 0.8) +  # Set gray border for each cell
  geom_text(color = "black", size = 4, hjust = 0.5, vjust = 0.75) +  # Center significance symbols
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",  # Blue for negative, red for positive
    midpoint = 0, na.value = "white"           # White for Z-score = 0, grey for NA
  ) +
  scale_x_discrete(
    labels = c(
      "sex.Z" = "Females vs. Males",
      "Z.women" = "AD vs. CO - Females",
      "Z.men" = "AD vs. CO - Males",
      "AD.sex_het.Zscore" = "AD vs. CO - Sex-Het")) +
  coord_cartesian(clip = "off") +  # Disable clipping to extend grid lines
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, color = "black", hjust = 1),  # Rotate x-axis labels
    axis.title = element_blank(),                       # Remove axis titles
    panel.grid.major.x = element_blank(),               # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),               # Remove minor vertical grid lines
    panel.grid.major.y = element_blank(),  # Add horizontal grid lines between facets
    panel.grid.minor.y = element_blank(),               # Remove minor horizontal grid lines
    axis.ticks.y = element_blank(),                     # Remove y-axis tick marks
    axis.line.y = element_blank(),                      # Remove y-axis line
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.y.right = element_text(angle = 0, face = "bold", size = 13, colour = "black"), # Adjust facet labels
    panel.spacing.y = unit(0, "lines")                  # Remove space between facets
  ) +
  facet_grid(locus_label ~ ., scales = "free_y", space = "free_y") +  # Group by locus_label
  labs(fill = "Z-score")  # Legend title

print(d1)


################################################################################################################################################################################################
# Combine the plots using patchwork
combined_plot <- p2 + score_plot + d1 + plot_layout(ncol = 3, widths = c(5, 0.5, 1))

# Display the combined plot
print(combined_plot)


file_name <- opt$plot_out
out_dir <- file.path(opt$work_dir, "results")

file_path <- file.path(out_dir, file_name)
ggsave(filename = file_path, plot = combined_plot, device = "png", width = 30, height = 16, units = "in", dpi = 320)

