#!/usr/bin/env Rscript


############################################################
#create circos plots instead of using the miami plot
## Michael's preferred option 3: 1 circos plot with sex-stratified PWAS results per tissue
## (outer most - inner circle: female brain, male brain, female CSF, male CSF)
# figure 2 dotplot matrix for PWAS manuscript main figure
# written by Danielle M. Read 
############################################################

#load library
#load the packages
suppressPackageStartupMessages({
  library(data.table)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(haven)
  library(gridExtra)
  library(ggpattern)
  library(UpSetR)
  library(stringr)
  library(biomaRt)
  library(CMplot)
  library(readxl)
  library(ggpubr)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
  library(forcats)
})

# Define command line options
option_list = list(
  make_option("--dir", type = "character", default = NULL,
              help = "Directory containing input files", metavar = "character"),
  make_option("--f_ss_genes", type = "character", default = NULL,
              help = "Female top sex-specific genes from both brain/CSF PWAS", metavar = "character"),
  make_option("--m_ss_genes", type = "character", default = NULL,
              help = "Male top sex-specific genes from both brain/CSF PWAS", metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
dir <- opt$dir
setwd(dir)

# dir = "/storage2/fs1/belloy2/Active/05_Projects/sivas/dreidRep/sumstats/Final_plots/"
# f = fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Females_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-female-specific-genes_ExDMR.txt"))
# m = fread(paste0(dir, "ADGC_ADSP_UKB_FinnGen_Males_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_PWAS_noMHC_ext2Mb_top-male-specific-genes_ExDMR.txt"))
f = fread(paste0(dir, opt$f_ss_genes)); m = as.data.frame(m)
m = fread(paste0(dir, opt$f_ss_genes)); f = as.data.frame(f)

f$P0 = as.numeric(f$P0); f$TWAS.P = as.numeric(f$TWAS.P); f$CHR = as.numeric(f$CHR)
m$P0 = as.numeric(m$P0); m$TWAS.P = as.numeric(m$TWAS.P); m$CHR = as.numeric(m$CHR)

# set female and male sig thresholds per tissue
#females
bp = 0.000772
bs = 0.00083
cp = 0.000799
cs = 0.000788

fb.sig = max(bp, bs)
fc.sig = max(cp, cs)

#males
bp = 0.000332
bs = 0.000233
cp = 6.91e-05
cs = 0.000554

mb.sig = max(bp, bs)
mc.sig = max(cp, cs)

# remember to remove CEBPZOS
f = subset(f, Gene != "CEBPZOS") # why are we removing this gene?

# for known AD gene label - PWAS gene(AD gene) (e.g., CD84(ADAMTS4))
f$Gene <- gsub("CD84", "CD84(ADAMTS4)", f$Gene)
f$Gene <- gsub("ENO3", "ENO3(SCIMP)", f$Gene)
f$Gene <- gsub("TOM1L2", "TOM1L2(MYO15A)", f$Gene)
m$Gene <- gsub("PHB1", "PHB1(ABI3)", m$Gene)

# for multiple genes in a loci, only label most sig gene at locus (previously TDRKH/MINDY1 -> TDRKH)
f$topSNP[f$Gene == "MINDY1"] <- 0
f$topSNP[f$Gene == "USP19"] <- 0
f$topSNP[f$Gene == "ANXA11"] <- 0
f$topSNP[f$Gene == "SCAPER"] <- 0
f$topSNP[f$Gene == "RCN2"] <- 0
f$topSNP[f$Gene == "RABEP1"] <- 0
f$topSNP[f$Gene == "CARM1"] <- 0

# to get CMplot to correctly pick a gene for highlighting, create unique names with TWAS.P
f$Unique_ID <- paste(f$Gene, f$TWAS.P, sep = ":")
m$Unique_ID <- paste(m$Gene, m$TWAS.P, sep = ":")

# subset the female and male dataframes per tissue
fb = subset(f, Tissue == "B", select = c(Unique_ID, CHR, P0, TWAS.P, topSNP, label_color, Gene))
fc = subset(f, Tissue == "C", select = c(Unique_ID, CHR, P0, TWAS.P, topSNP, label_color, Gene))
mb = subset(m, Tissue == "B", select = c(Unique_ID, CHR, P0, TWAS.P, topSNP, label_color, Gene))
mc = subset(m, Tissue == "C", select = c(Unique_ID, CHR, P0, TWAS.P, topSNP, label_color, Gene))

# change colnames
colnames(fb) <- c("Unique_ID", "CHR", "POS", "TWAS.P", "topSNP", "label_color", "Gene")
colnames(fc) <- c("Unique_ID", "CHR", "POS", "TWAS.P", "topSNP", "label_color", "Gene")
colnames(mb) <- c("Unique_ID", "CHR", "POS", "TWAS.P", "topSNP", "label_color", "Gene")
colnames(mc) <- c("Unique_ID", "CHR", "POS", "TWAS.P", "topSNP", "label_color", "Gene")

## update blue color to what Michael is now using in Fig. 1, since it is slightly off from the blue i used
# update label_color per sex for top hits
fb$label_color <- ifelse(!is.na(fb$label_color), "#700000", NA)
fc$label_color <- ifelse(!is.na(fc$label_color), "#700000", NA)
mb$label_color <- ifelse(!is.na(mb$label_color), "#3c7899", NA)
mc$label_color <- ifelse(!is.na(mc$label_color), "#3c7899", NA)

# rbind the female and male results
comb_dat = rbind(mc, fc, mb, fb)

# select genes to highlight
x = filter(comb_dat, topSNP > 0)

# Recombine all Unique_ID/CHR/POS info in the same order as comb_dat
base_dat <- rbind(mc[, c("Unique_ID", "CHR", "POS")],
                  fc[, c("Unique_ID", "CHR", "POS")],
                  mb[, c("Unique_ID", "CHR", "POS")],
                  fb[, c("Unique_ID", "CHR", "POS")])

# Extract the trait values
trait1 <- mc$TWAS.P
trait2 <- fc$TWAS.P
trait3 <- mb$TWAS.P
trait4 <- fb$TWAS.P

# Combine into one dataframe
dat <- data.frame(
  Unique_ID = base_dat$Unique_ID,
  CHR = base_dat$CHR,
  POS = base_dat$POS,
  trait1 = c(trait1, rep(NA, length.out = nrow(fc) + nrow(mb) + nrow(fb))),
  trait2 = c(rep(NA, length.out = nrow(mc)), trait2, rep(NA, length.out = nrow(mb) + nrow(fb))),
  trait3 = c(rep(NA, length.out = nrow(mc) + nrow(fc)), trait3, rep(NA, length.out = nrow(fb))),
  trait4 = c(rep(NA, length.out = nrow(mc) + nrow(fc) + nrow(mb)), trait4)
)

highlight_vec <- setNames(x$Unique_ID, x$Gene)

## grey non-sex-specific genes
CMplot(dat,
        type = "p",
        plot.type = "c", 
        LOG10 = TRUE, 
        col = matrix(c("grey60", "grey60", "grey60", "grey60"), nrow = 4, byrow = TRUE),
        pch = c(17, 17, 19, 19),
        band = 0.5,
        ylim = c(0,15),
        r = 0.1,
        cir.axis = F,
        outward = F,
        cir.axis.col = "black",
        cir.chr.h=1,
        chr.den.col = "grey60",
        threshold = list(mc.sig, fc.sig, mb.sig, fb.sig),
        threshold.lwd = c(1.5), 
        threshold.col = rep("black", 4),
        highlight = highlight_vec,
        highlight.col = x$label_color,
        highlight.cex = 1.3, 
        highlight.text.cex = 1.7,
        highlight.text.col = x$label_color,
        chr.labels = as.character(1:22),
        amplify = F,
        signal.line = 1,
        file = "pdf", file.name = "ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var.W23nCSF.hg38_top-sex-specific-genesCIRCOS4",
        file.output = TRUE, verbose = TRUE, width = 6, height = 6
)


## update blue color to what Michael is now using in Fig. 1, since it is slightly off from the blue i used
# update label_color per sex for top hits
fb$label_color <- ifelse(!is.na(fb$label_color), "#700000", NA)
fc$label_color <- ifelse(!is.na(fc$label_color), "#700000", NA)
mb$label_color <- ifelse(!is.na(mb$label_color), "#3c7899", NA)
mc$label_color <- ifelse(!is.na(mc$label_color), "#3c7899", NA)

comb_dat1 = rbind(mc, fc)
comb_dat2 = rbind(mb, fb)

x1 = filter(comb_dat1, topSNP > 0)
x2 = filter(comb_dat2, topSNP > 0)

base_dat1 <- rbind(mc[, c("Unique_ID", "CHR", "POS")],
                  fc[, c("Unique_ID", "CHR", "POS")]
)

base_dat2 <- rbind(mb[, c("Unique_ID", "CHR", "POS")],
                  fb[, c("Unique_ID", "CHR", "POS")]
)

# Extract the trait values
trait1 <- mc$TWAS.P
trait2 <- fc$TWAS.P
trait3 <- mb$TWAS.P
trait4 <- fb$TWAS.P

dat1 <- data.frame(
  Unique_ID = base_dat1$Unique_ID,
  CHR = base_dat1$CHR,
  POS = base_dat1$POS,
  trait1 = c(trait1, rep(NA, length.out = nrow(fc))),
  trait2 = c(rep(NA, length.out = nrow(mc)), trait2)
)

dat2 <- data.frame(
  Unique_ID = base_dat2$Unique_ID,
  CHR = base_dat2$CHR,
  POS = base_dat2$POS,
  trait3 = c(trait3, rep(NA, length.out = nrow(fb))),
  trait4 = c(rep(NA, length.out = nrow(mb)), trait4)
)

highlight1 <- setNames(x1$Unique_ID, x1$Gene)
highlight2 <- setNames(x2$Unique_ID, x2$Gene)

# grey non-sex-specigic genes
CMplot(dat1,
        type = "p",
        plot.type = "c", 
        LOG10 = TRUE, 
        col = matrix(c("grey60", "grey60"), nrow = 2, byrow = TRUE),
        pch = 17,
        band = 0.5,
        ylim = c(0,15),
        r = 0.1,
        cir.axis = F,
        outward = F,
        cir.axis.col = "black",
        cir.chr.h=1,
        chr.den.col = "grey60",
        threshold = list(mc.sig, fc.sig),
        threshold.lwd = c(1.5), 
        threshold.col = rep("black", 2),
        highlight = highlight1,
        highlight.col = x1$label_color,
        highlight.cex = 1.3, 
        chr.labels = as.character(1:22),
        amplify = F,
        signal.line = 1,
        file = "pdf", file.name = "ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var.CSF.hg38_top-sex-specific-genesCIRCOS4.1",
        file.output = TRUE, verbose = TRUE, width = 6, height = 6
)

# grey non-sex-specific genes
CMplot(dat2,
        type = "p",
        plot.type = "c", 
        LOG10 = TRUE, 
        col = matrix(c("grey60", "grey60"), nrow = 2, byrow = TRUE),
        pch = 19,
        band = 0.5,
        ylim = c(0,15),
        r = 0.1,
        cir.axis = F,
        outward = F,
        cir.axis.col = "black",
        cir.chr.h=1,
        chr.den.col = "grey60",
        threshold = list(mb.sig, fb.sig),
        threshold.lwd = c(1.5), 
        threshold.col = rep("black", 2),
        highlight = highlight2,
        highlight.col = x2$label_color,
        highlight.cex = 1.3, 
        chr.labels = as.character(1:22),
        amplify = F,
        signal.line = 1,
        file = "pdf", file.name = "ADGC_ADSP_UKB_FinnGen_cc_rb.gen090.noAPOE.shared_var.W23.hg38_top-sex-specific-genesCIRCOS4.2",
        file.output = TRUE, verbose = TRUE, width = 6, height = 6
)