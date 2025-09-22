#!/usr/bin/env Rscript

## make Venn Diagrams for manuscript because files need to be vectors
###################### Danielle Reid - 05/23/25
#######################################################

#load the packages
suppressPackageStartupMessages({
	library(plyr)
	library(dplyr)
	library(readr)
	library(data.table)
	library(tidyr)
	library(ggVennDiagram)
	library(ggplot2)
	library(optparse)
	library(data.table)
})

# Define command line options
option_list = list(
  make_option("--dir", type = "character", default = NULL,
              help = "Directory containing input files", metavar = "character"),
  make_option("--QCfile", type = "character", default = NULL,
              help = "Path to first input file (Overlap_W23-norm_NGI-CSF-QC.txt)", metavar = "character"),
  make_option("--weightfile", type = "character", default = NULL,
              help = "Path to second input file (Overlap_cis_NGI-CSF_W23_weights.txt)", metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# dir <- "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/DMR_plots/"
# data <- fread(paste0(dir, "Overlap_W23-norm_NGI-CSF-QC.txt"), sep = '\t'); dat = as.data.frame(dat)
# data2 <- fread(paste0(dir, "Overlap_cis_NGI-CSF_W23_weights.txt"), sep = '\t'); dat = as.data.frame(dat2)

# Assign arguments to variables
dir <- opt$dir
dat <- fread(paste0(dir, opt$QCfile), sep = '\t'); dat = as.data.frame(dat)
setwd(dir)
names(dat)[3] = "W23"
# Create set1: ENSG values that have a value of 1 for Somascan
set1 <- dat %>%
  filter(Somascan == 1) %>%
  pull(ENSG
)
# Create set2: ENSG values that have a value of 1 for W23
set2 <- dat %>%
  filter(W23 == 1) %>%
  pull(ENSG
)


venn <- venn.diagram(
  x = list(set1, set2), category.names = c("CSF" , "Brain"), filename = NULL, imagetype = 'tiff',
  # Numbers
  cex = 1.7, fontface = "bold",
  # Circles
  lwd = 1, scaled = FALSE, fill = c("navy", "gray80"), alpha = c("0.5", "0.5"),
  # Set names
  cat.cex = 2, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(30, -30), cat.dist = c(0.03, 0.04))
pdf(paste0(dir, "Overlap_W23-norm_NGI-CSF-QC.pdf"), height = 6, width = 6)
grid::grid.draw(venn)
dev.off()

## load overlap file for NGI CSF cis proteins vs W23 protein weights
dat2 = fread(paste0(dir, opt$weightfile), sep = '\t'); dat2 = as.data.frame(dat2)

pdf(paste0(dir, "Overlap_cis_NGI-CSF_W23_weights.v1.pdf"))
upset(dat2, sets = c("NGI_CSF_Non-sex-strat_Weights", "NGI_CSF_Female_Weights", "NGI_CSF_Male_Weights", "W23_Non-sex-strat_Weights", "W23_Female_Weights", "W23_Male_Weights"),
        keep.order = T, order.by = "freq", mainbar.y.label = "CSF and Brain Gene/Protein Intersections in Weights",
        sets.x.label = "Proteins/Genes"
)
dev.off()

colnames(dat2) = c("ENSG", "CSF_Non-sex-strat", "CSF_Female", "CSF_Male", "Brain_Non-sex-strat", "Brain_Female", "Brain_Male")

pdf(paste0(dir, "Overlap_cis_NGI-CSF_W23_weights.v2.pdf"))
upset(dat2, sets = c("CSF_Non-sex-strat", "CSF_Female", "CSF_Male", "Brain_Non-sex-strat", "Brain_Female", "Brain_Male"),
        keep.order = T, order.by = "freq", mainbar.y.label = "CSF and Brain Gene/Protein Intersections in Weights",
        sets.x.label = "Proteins/Genes"
)
dev.off()

## make venn diagram of overlap in weights between CSF and brain
#create set : ENSG values that have 1 for anything CSF
c1 <- dat2 %>%
    filter(`CSF_Non-sex-strat` == 1) %>%
    pull(ENSG
)

c2 <- dat2 %>%
    filter(CSF_Female == 1) %>%
    pull(ENSG
)

c3 <- dat2 %>%
    filter(CSF_Male == 1) %>%
    pull(ENSG
)

c = union(c1, c(c2, c3))

set1 <- unique(c)

# create set2: ENSG values that have 1 for anything brain
b1 <- dat2 %>%
    filter(`Brain_Non-sex-strat` == 1) %>%
    pull(ENSG
)

b2 <- dat2 %>%
    filter(Brain_Female == 1) %>%
    pull(ENSG
)

b3 <- dat2 %>%
    filter(Brain_Male == 1) %>%
    pull(ENSG
)

b = union(b1, c(b2, b3))

set2 <- unique(b)

venn2 <- venn.diagram(
  x = list(set1, set2), category.names = c("CSF" , "Brain"), filename = NULL, imagetype = 'tiff',
  # Numbers
  cex = 1.7, fontface = "bold",
  # Circles
  lwd = 1, scaled = FALSE, fill = c("navy", "gray80"), alpha = c("0.5", "0.5"),
  # Set names
  cat.cex = 2, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(30, -30), cat.dist = c(0.03, 0.04))
  
pdf(paste0(dir, "Overlap_cis_NGI-CSF_W23_weights.v3.pdf"), height = 6, width = 6)
grid::grid.draw(venn2)
dev.off()
