#### process FUSION association results from integrating sex-stratified GWAS summary stats
#### with Wingo et al. 2023 non-sex-stratified & sex-stratified brain protein weights
#### to determine significant genes/proteins and then perform sex-specific filters to identify sex-specific genes/proteins
#### written by Danielle M. Reid (dmr07083)

#########################
#load R libraries
library(data.table); library(plyr); library(dplyr); library(readr)
library(ggplot2); library(ggtext); library(haven); library(ggrepel)
library(methods); library(optparse); library(tidyr); library(stringr)
library(hrbrthemes); library(naniar); library(ggnewscale)

################################################################################
option_list = list(
  make_option("--non_strat_dir", action="store", default=NA, type='character',
              help="Path to non-sex-stratified Wingo et al. 2023 brain directory containing sex-specific subdirectories with protein imputation files [required].\n
                      (e.g., Female and Male)"),
  make_option("--sex_strat_dir", action="store", default=NA, type='character',
              help="Path to sex-stratified Wingo et al. 2023 brain directory containing sex-specific subdirectories with protein imputation files [required].\n
                      (e.g., Female, Male, FemaleGWAS_MalePW, and MaleGWAS_FemalePW)"),
  make_option("--female_file", action="store", default=NA, type='character',
              help="Female file name [required]"),
  make_option("--male_file", action="store", default=NA, type='character',
              help="Male file name [required]"),
  make_option("--no_mhc", action="store", default=FALSE, type='character',
              help="Exclude files ending with '.dat.MHC' from processing: [default: %default].\n
                     Set this flag to TRUE to exclude processing the Major Histocompatibility Complex (MHC) region [optional]."),
  make_option("--ext_mhc", action="store", default=2, type='numeric',
              help="Extend exclusion area surrounding MHC region by: +-[default: %default]Mb [optional]."),                   
  make_option("--gene_file", action="store", default=NA, type='character',
              help="Gene name replacement file including path [optional]"),
  make_option("--ad_ref", action="store", default=NA, type='character',
              help="AD risk loci reference file for characterizing loci novelty [optional]"),
  make_option("--lift", action="store", default=NA, type='character',
              help="File to liftover genome coordinates from hg19 to hg38 [optional]")                                 
)

opt = parse_args(OptionParser(option_list=option_list))

#########################

# input options
non_dir = opt$non_strat_dir; print(paste("non-sex-stratified W23 directory:", non_dir)); print(paste("constructed non-sex directory: ", file.path(non_dir)))
sex_dir = opt$sex_strat_dir; print(paste("sex-stratified W23 directory:", sex_dir)); print(paste("constructed sex directory: ", file.path(sex_dir)))
f_dir = "Female"; print(paste("constructed F-N directory: ", file.path(non_dir, f_dir)))
print(paste("constructed F-F directory: ", file.path(sex_dir, f_dir)))
m_dir = "Male"; print(paste("constructed M-N directory: ", file.path(non_dir, m_dir)))
print(paste("constructed M-M directory: ", file.path(sex_dir, m_dir)))
fm_dir = "FemaleGWAS_MalePW"; print(paste("constructed F-M directory: ", file.path(sex_dir, fm_dir)))
mf_dir = "MaleGWAS_FemalePW"; print(paste("constructed M-F directory: ", file.path(sex_dir, mf_dir)))
f_f = opt$female_file; print(paste("female file name:", f_f))
m_f = opt$male_file; print(paste("male file name:", m_f))

###################### PWAS results for non-sex-strat #########################
# Females results
print("processing females results for Wingo et al. 2023 non-sex-stratified brain protein weights")

## Make a list of files in directory to merge ".dat" files
file_list <- list.files(path = file.path(non_dir, f_dir), pattern = ".dat", full.names = TRUE)

# Exclude .dat.MHC files if --mhc is TRUE
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_list <- file_list[!grepl(".dat.MHC$", file_list)]
}

#make dataframe; remove PANEL and FILE columns; clean up df before saving
pwas_dat = ldply(file_list, read_delim, show_col_types = FALSE); print("loaded .dat files")
pwas_dat = select(pwas_dat, -PANEL, -FILE); print("removed PANEL and FILE column")
names(pwas_dat)[1] = "ENSG_ID"
pwas_dat$ID = gsub("ENSG00000\\d\\d\\d\\d\\d\\d.", "", as.character(pwas_dat$ENSG_ID))
if (!is.na(opt$gene_file)) {
  pwas_dat$Gene = pwas_dat$ID
  col_to_move = c("ID", "Gene")
  pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))

  # function to replace ID (gene) names that are not in standard Symbol nomenclature
  replace_genes <- function(df, column, replacements) {
    # Update replacements to use word boundaries for regex
    replacements_with_boundaries <- setNames(
      replacements,
      paste0("\\b", names(replacements), "\\b")
    )

    df %>%
      mutate(!!sym(column) := str_replace_all(.data[[column]], replacements_with_boundaries))
  }
  # updated gene name file
  source(opt$gene_file)

  # define the replacement rules as a named vector
  df <- list(pwas_dat = pwas_dat)
  col <- "Gene"

  # apply the function to the dataframe and column if exists
  for (df_name in names(df)) {
    for (col in col) {
      if (col %in% colnames(df[[df_name]])) {
        df[[df_name]] <- replace_genes(df[[df_name]], col, replacements)
      }
    }
  }

  # unpack modified data frame back to individual variable
  list2env(df, .GlobalEnv)
} else {
   col_to_move = "ID"
   pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))
   names(pwas_dat)[2] = "Gene"
}

# if --lift option is not NA, read in the file and add hg38 coordinates for weights
if (!is.na(opt$lift)) {
  coord_dat <- fread(opt$lift, sep = '\t')

  # merge pwas_dat with the coordinate file by ENSG_ID
  pwas_dat <- pwas_dat %>%
    left_join(coord_dat %>% dplyr::select(ENSG_ID, Start_hg38, End_hg38), by = "ENSG_ID")
}

# save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  pwas_file = paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights-raw.txt")
} else {
  pwas_file = paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights-raw.txt")
}
fwrite(pwas_dat, pwas_file, row.names=F,quote=F, sep='\t'); print("saved pwas_dat file: "); print(pwas_file)

## remove NAs, convert p and P0 & P1 to numeric, calculate fdr_p & sort df on fdr_p
new_pwas_dat = na.omit(pwas_dat); print("nrow of new_pwas_dat after removing NAs"); print(nrow(new_pwas_dat))
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

if (!is.na(opt$lift)) {
  new_pwas_dat$Start_hg38 = as.numeric(new_pwas_dat$Start_hg38)
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

# if used opt$ext_mhc remove the rows with a P0 and P1 (gene start and end) within the extended region
if (!is.na(opt$ext_mhc) && as.numeric(opt$ext_mhc) > 0) {
  mhc_start <- 28477797
  mhc_end <- 33448354
  ext_mhc <- as.numeric(opt$ext_mhc) * 1e6
  ext_mhc_start <- mhc_start - ext_mhc
  ext_mhc_end <- mhc_end + ext_mhc
  print(paste("Extended MHC region:", ext_mhc_start, "to", ext_mhc_end))
  new_pwas_dat <- new_pwas_dat[!(new_pwas_dat$CHR == 6 &
                         (new_pwas_dat$P0 >= ext_mhc_start & new_pwas_dat$P0 <= ext_mhc_end)), ]
}
fdr_p = p.adjust(new_pwas_dat$TWAS.P, method = "BH")
new_pwas_dat$fdr_p = fdr_p
new_pwas_dat = new_pwas_dat[order(new_pwas_dat$fdr_p),]

#create a variable for the path for where to save the summarized PWAS results as a txt file and save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  new_pwas_file = paste0(file.path(non_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights.txt")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  new_pwas_file = paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights.txt")
} else {
   new_pwas_file = paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights.txt")
}
fwrite(new_pwas_dat, new_pwas_file, row.names=F,quote=F, sep='\t'); print("saved new_pwas_dat file: "); print(new_pwas_file)

# print sum of genes with fdr_p below 0.05
new_pwas_dat %>% filter(fdr_p <0.05)
print("Number of genes that pass fdr_p < 0.05"); sum(new_pwas_dat$fdr_p <0.05)

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Set y-axis limit above highest SNP -log10(p) by transforming the largest exponent to positive and add 1
ylim = abs(floor(log10(min(new_pwas_dat$TWAS.P)))) +1; print("ylim = "); print(ylim)

## determine sig threshold
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- new_pwas_dat %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # Using `match` to ensure we get the first occurrence if both fdr_p and max TWAS.P are repeated
  index1 <- match(TRUE, new_pwas_dat$fdr_p == result$fdr_p & new_pwas_dat$TWAS.P == result$TWAS.P)
  sig <- new_pwas_dat$TWAS.P[index1]
  print(paste("index row = ", index1))
  print(paste("p signficance threshold = ", sig))
  print(paste("-log10 sig threshold = ", -log10(sig)))
} else {
  index1 <- NA
  sig <- NA
  print("No valid significance threshold found")
}

# Make manhattan plot of PWAS results for Females
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#009ACD", size = 0.5) +
  geom_text_repel(
    data = subset(new_pwas_dat, TWAS.P <= sig), 
    aes(label = Gene, y = -log10(TWAS.P)),
    color = "#009ACD",
    size = 3,
    nudge_y = 0.3,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_manhplot.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights_manhplot.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights_hg38_manhplot.pdf")
} else {
   file_name <- paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights_manhplot.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# Make manhattan plot of PWAS results for Females - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#009ACD", size = 0.5) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_manhplot2.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights_manhplot2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights_hg38_manhplot2.pdf")
} else {
   file_name <- paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights_manhplot2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()


print("Females results completed for non-sex-strat W23 weigthts")

# remove variables
if (!is.na(opt$gene_file)) {
  rm(df, col, df_name)
}

rm(axis_set, col_to_move, data_cum, fdr_p, file_list, index1, manhplot, new_pwas_dat, 
new_pwas_file, pwas_dat, pwas_file, result, sig, ylim
); print("removed female variables")

#######################################################

#Males results
print("processing males results for Wingo et al. 2023 non-sex-stratified brain protein weights")

#Make a list of files in directory to merge ".dat" files
file_list <- list.files(path = file.path(non_dir, m_dir), pattern = ".dat", full.names = TRUE)

#Exclude .dat.MHC files if --mhc is TRUE
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_list <- file_list[!grepl(".dat.MHC$", file_list)]
}

#make dataframe; remove PANEL and FILE columns; clean up df before saving
pwas_dat = ldply(file_list, read_delim, show_col_types = FALSE); print("loaded .dat files")
pwas_dat = select(pwas_dat, -PANEL, -FILE); print("removed PANEL and FILE column")
names(pwas_dat)[1] = "ENSG_ID"
pwas_dat$ID = gsub("ENSG00000\\d\\d\\d\\d\\d\\d.", "", as.character(pwas_dat$ENSG_ID))
if (!is.na(opt$gene_file)) {
  pwas_dat$Gene = pwas_dat$ID
  col_to_move = c("ID", "Gene")
  pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))

  # define the replacement rules as a named vector
  df <- list(pwas_dat = pwas_dat)
  col <- "Gene"

  # apply the function to the dataframe and column if exists
  for (df_name in names(df)) {
    for (col in col) {
      if (col %in% colnames(df[[df_name]])) {
        df[[df_name]] <- replace_genes(df[[df_name]], col, replacements)
      }
    }
  }

  # unpack modified data frame back to individual variable
  list2env(df, .GlobalEnv)
} else {
   col_to_move = "ID"
   pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))
   names(pwas_dat)[2] = "Gene"
}

# if --lift option is not NA, read in the file and add hg38 coordinates for weights
if (!is.na(opt$lift)) {
  coord_dat <- fread(opt$lift, sep = '\t')

  # merge pwas_dat with the coordinate file by ENSG_ID
  pwas_dat <- pwas_dat %>%
    left_join(coord_dat %>% dplyr::select(ENSG_ID, Start_hg38, End_hg38), by = "ENSG_ID")
}

# save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  pwas_file = paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights-raw.txt")
} else {
  pwas_file = paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights-raw.txt")
}
fwrite(pwas_dat, pwas_file, row.names=F,quote=F, sep='\t'); print("saved pwas_dat file: "); print(pwas_file)

## remove NAs, convert p and P0 & P1 to numeric, calculate fdr_p & sort df on fdr_p
new_pwas_dat = na.omit(pwas_dat); print("nrow of new_pwas_dat after removing NAs"); print(nrow(new_pwas_dat))
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

if (!is.na(opt$lift)) {
  new_pwas_dat$Start_hg38 = as.numeric(new_pwas_dat$Start_hg38)
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

# if used opt$ext_mhc remove the rows with a P0 and P1 (gene start and end) within the extended region
if (!is.na(opt$ext_mhc) && as.numeric(opt$ext_mhc) > 0) {
  print(paste("Extended MHC region:", ext_mhc_start, "to", ext_mhc_end))
  new_pwas_dat <- new_pwas_dat[!(new_pwas_dat$CHR == 6 &
                         (new_pwas_dat$P0 >= ext_mhc_start & new_pwas_dat$P0 <= ext_mhc_end)), ]
}
fdr_p = p.adjust(new_pwas_dat$TWAS.P, method = "BH")
new_pwas_dat$fdr_p = fdr_p
new_pwas_dat = new_pwas_dat[order(new_pwas_dat$fdr_p),]

#create a variable for the path for where to save the summarized PWAS results as a txt file and save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  new_pwas_file = paste0(file.path(non_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights.txt")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  new_pwas_file = paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights.txt")
} else {
   new_pwas_file = paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights.txt")
}
fwrite(new_pwas_dat, new_pwas_file, row.names=F,quote=F, sep='\t'); print("saved new_pwas_dat file: "); print(new_pwas_file)

# print sum of genes with fdr_p below 0.05
new_pwas_dat %>% filter(fdr_p <0.05)
print("Number of genes that pass fdr_p < 0.05"); sum(new_pwas_dat$fdr_p <0.05)

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Set y-axis limit above highest SNP -log10(p) by transforming the largest exponent to positive and add 1
ylim = abs(floor(log10(min(new_pwas_dat$TWAS.P)))) +1; print("ylim = "); print(ylim)

## determine sig threshold
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- new_pwas_dat %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # Using `match` to ensure we get the first occurrence if both fdr_p and max TWAS.P are repeated
  index1 <- match(TRUE, new_pwas_dat$fdr_p == result$fdr_p & new_pwas_dat$TWAS.P == result$TWAS.P)
  sig <- new_pwas_dat$TWAS.P[index1]
  print(paste("index row = ", index1))
  print(paste("p signficance threshold = ", sig))
  print(paste("-log10 sig threshold = ", -log10(sig)))
} else {
  index1 <- NA
  sig <- NA
  print("No valid significance threshold found")
}

# Make manhattan plot of PWAS results for Males
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#2E8B57", size = 0.5) +
  geom_text_repel(
    data = subset(new_pwas_dat, TWAS.P <= sig), 
    aes(label = Gene, y = -log10(TWAS.P)),
    color = "#2E8B57",
    size = 3,
    nudge_y = 0.3,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_manhplot.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights_manhplot.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights_hg38_manhplot.pdf")
} else {
   file_name <- paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights_manhplot.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# Make manhattan plot of PWAS results for Males - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#2E8B57", size = 0.5) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights_manhplot2.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights_manhplot2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights_hg38_manhplot2.pdf")
} else {
   file_name <- paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights_manhplot2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

print("Males results completed for non-sex-strat W23 weigthts")

if (!is.na(opt$gene_file)) {
  rm(df, col, df_name)
}

rm(axis_set, col_to_move, data_cum, fdr_p, file_list, index1, manhplot, new_pwas_dat, 
new_pwas_file, pwas_dat, pwas_file, result, sig, ylim
); print("removed male variables")

###################### PWAS results for sex-matched #########################
print("processing females results for Wingo et al. 2023 sex-stratified brain protein weights")

#Make a list of files in directory to merge ".dat" files
file_list <- list.files(path = file.path(sex_dir, f_dir), pattern = ".dat", full.names = TRUE)

# Exclude .dat.MHC files if --mhc is TRUE
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_list <- file_list[!grepl(".dat.MHC$", file_list)]
}

#make dataframe; remove PANEL and FILE columns; clean up df before saving
pwas_dat = ldply(file_list, read_delim, show_col_types = FALSE); print("loaded .dat files")
pwas_dat = select(pwas_dat, -PANEL, -FILE); print("removed PANEL and FILE column")
names(pwas_dat)[1] = "ENSG_ID"
pwas_dat$ID = gsub("ENSG00000\\d\\d\\d\\d\\d\\d.", "", as.character(pwas_dat$ENSG_ID))
if (!is.na(opt$gene_file)) {
  pwas_dat$Gene = pwas_dat$ID
  col_to_move = c("ID", "Gene")
  pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))

  # define the replacement rules as a named vector
  df <- list(pwas_dat = pwas_dat)
  col <- "Gene"

  # apply the function to the dataframe and column if exists
  for (df_name in names(df)) {
    for (col in col) {
      if (col %in% colnames(df[[df_name]])) {
        df[[df_name]] <- replace_genes(df[[df_name]], col, replacements)
      }
    }
  }

  # unpack modified data frame back to individual variable
  list2env(df, .GlobalEnv)
} else {
   col_to_move = "ID"
   pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))
   names(pwas_dat)[2] = "Gene"
}

# if --lift option is not NA, read in the file and add hg38 coordinates for weights
if (!is.na(opt$lift)) {
  coord_dat <- fread(opt$lift, sep = '\t')

  # merge pwas_dat with the coordinate file by ENSG_ID
  pwas_dat <- pwas_dat %>%
    left_join(coord_dat %>% dplyr::select(ENSG_ID, Start_hg38, End_hg38), by = "ENSG_ID")
}

# save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  pwas_file = paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights-raw.txt")
} else {
  pwas_file = paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights-raw.txt")
}
fwrite(pwas_dat, pwas_file, row.names=F,quote=F, sep='\t'); print("saved pwas_dat file: "); print(pwas_file)

#remove NAs, convert p and P0 & P1 to numeric, calculate fdr_p & sort df on fdr_p
new_pwas_dat = na.omit(pwas_dat); print("nrow of new_pwas_dat after removing NAs"); print(nrow(new_pwas_dat))
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

if (!is.na(opt$lift)) {
  new_pwas_dat$Start_hg38 = as.numeric(new_pwas_dat$Start_hg38)
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

# if used opt$ext_mhc remove the rows with a P0 and P1 (gene start and end) within the extended region
if (!is.na(opt$ext_mhc) && as.numeric(opt$ext_mhc) > 0) {
  print(paste("Extended MHC region:", ext_mhc_start, "to", ext_mhc_end))
  new_pwas_dat <- new_pwas_dat[!(new_pwas_dat$CHR == 6 &
                         (new_pwas_dat$P0 >= ext_mhc_start & new_pwas_dat$P0 <= ext_mhc_end)), ]
}
fdr_p = p.adjust(new_pwas_dat$TWAS.P, method = "BH")
new_pwas_dat$fdr_p = fdr_p
new_pwas_dat = new_pwas_dat[order(new_pwas_dat$fdr_p),]

#create a variable for the path for where to save the summarized PWAS results as a txt file and save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  new_pwas_file = paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights.txt")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  new_pwas_file = paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights.txt")
} else {
   new_pwas_file = paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights.txt")
}
fwrite(new_pwas_dat, new_pwas_file, row.names=F,quote=F, sep='\t'); print("saved new_pwas_dat file: "); print(new_pwas_file)

# print sum of genes with fdr_p below 0.05
new_pwas_dat %>% filter(fdr_p <0.05)
print("Number of genes that pass fdr_p < 0.05"); sum(new_pwas_dat$fdr_p <0.05)

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Set y-axis limit above highest SNP -log10(p) by transforming the largest exponent to positive and add 1
ylim = abs(floor(log10(min(new_pwas_dat$TWAS.P)))) +1; print("ylim = "); print(ylim)

## determine sig threshold
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- new_pwas_dat %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # Using `match` to ensure we get the first occurrence if both fdr_p and max TWAS.P are repeated
  index1 <- match(TRUE, new_pwas_dat$fdr_p == result$fdr_p & new_pwas_dat$TWAS.P == result$TWAS.P)
  sig <- new_pwas_dat$TWAS.P[index1]
  print(paste("index row = ", index1))
  print(paste("p signficance threshold = ", sig))
  print(paste("-log10 sig threshold = ", -log10(sig)))
} else {
  index1 <- NA
  sig <- NA
  print("No valid significance threshold found")
}

# Make manhattan plot of PWAS results for Females
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#009ACD", size = 0.5) +
  geom_text_repel(
    data = subset(new_pwas_dat, TWAS.P <= sig), 
    aes(label = Gene, y = -log10(TWAS.P)),
    color = "#009ACD",
    size = 3,
    nudge_y = 0.3,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_manhplot.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights_manhplot.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights_hg38_manhplot.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights_manhplot.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# Make manhattan plot of PWAS results for Females - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#009ACD", size = 0.5) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_manhplot2.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights_manhplot2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights_hg38_manhplot2.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights_manhplot2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

print("Females results completed for sex-strat W23 weights")

# remove variables
if (!is.na(opt$gene_file)) {
  rm(df, col, df_name)
}

rm(axis_set, col_to_move, data_cum, fdr_p, file_list, index1, manhplot, new_pwas_dat, 
new_pwas_file, pwas_dat, pwas_file, result, sig, ylim
); print("removed female variables")

#######################################################

#Males results
print("processing males results for Wingo et al. 2023 sex-stratified brain protein weights")

#Make a list of files in directory to merge ".dat" files 
file_list <- list.files(path = file.path(sex_dir, m_dir), pattern = ".dat", full.names = TRUE)

#Exclude .dat.MHC files if --mhc is TRUE
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_list <- file_list[!grepl(".dat.MHC$", file_list)]
}

#make dataframe; remove PANEL and FILE columns; clean up df before saving
pwas_dat = ldply(file_list, read_delim, show_col_types = FALSE); print("loaded .dat files")
pwas_dat = select(pwas_dat, -PANEL, -FILE); print("removed PANEL and FILE column")
names(pwas_dat)[1] = "ENSG_ID"
pwas_dat$ID = gsub("ENSG00000\\d\\d\\d\\d\\d\\d.", "", as.character(pwas_dat$ENSG_ID))
if (!is.na(opt$gene_file)) {
  pwas_dat$Gene = pwas_dat$ID
  col_to_move = c("ID", "Gene")
  pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))

  # define the replacement rules as a named vector
  df <- list(pwas_dat = pwas_dat)
  col <- "Gene"

  # apply the function to the dataframe and column if exists
  for (df_name in names(df)) {
    for (col in col) {
      if (col %in% colnames(df[[df_name]])) {
        df[[df_name]] <- replace_genes(df[[df_name]], col, replacements)
      }
    }
  }

  # unpack modified data frame back to individual variable
  list2env(df, .GlobalEnv)
} else {
   col_to_move = "ID"
   pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))
   names(pwas_dat)[2] = "Gene"
}

# if --lift option is not NA, read in the file and add hg38 coordinates for weights
if (!is.na(opt$lift)) {
  coord_dat <- fread(opt$lift, sep = '\t')

  # merge pwas_dat with the coordinate file by ENSG_ID
  pwas_dat <- pwas_dat %>%
    left_join(coord_dat %>% dplyr::select(ENSG_ID, Start_hg38, End_hg38), by = "ENSG_ID")
}

# save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  pwas_file = paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights-raw.txt")
} else {
  pwas_file = paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights-raw.txt")
}
fwrite(pwas_dat, pwas_file, row.names=F,quote=F, sep='\t'); print("saved pwas_dat file: "); print(pwas_file)

#remove NAs, convert p and P0 & P1 to numeric, calculate fdr_p & sort df on fdr_p
new_pwas_dat = na.omit(pwas_dat); print("nrow of new_pwas_dat after removing NAs"); print(nrow(new_pwas_dat))
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

if (!is.na(opt$lift)) {
  new_pwas_dat$Start_hg38 = as.numeric(new_pwas_dat$Start_hg38)
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

# if used opt$ext_mhc remove the rows with a P0 and P1 (gene start and end) within the extended region
if (!is.na(opt$ext_mhc) && as.numeric(opt$ext_mhc) > 0) {
  print(paste("Extended MHC region:", ext_mhc_start, "to", ext_mhc_end))
  new_pwas_dat <- new_pwas_dat[!(new_pwas_dat$CHR == 6 &
                         (new_pwas_dat$P0 >= ext_mhc_start & new_pwas_dat$P0 <= ext_mhc_end)), ]
}
fdr_p = p.adjust(new_pwas_dat$TWAS.P, method = "BH")
new_pwas_dat$fdr_p = fdr_p
new_pwas_dat = new_pwas_dat[order(new_pwas_dat$fdr_p),]

#create a variable for the path for where to save the summarized PWAS results as a txt file and save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  new_pwas_file = paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights.txt")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  new_pwas_file = paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights.txt")
} else {
   new_pwas_file = paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights.txt")
}
fwrite(new_pwas_dat, new_pwas_file, row.names=F,quote=F, sep='\t'); print("saved new_pwas_dat file: "); print(new_pwas_file)

# print sum of genes with fdr_p below 0.05
new_pwas_dat %>% filter(fdr_p <0.05)
print("Number of genes that pass fdr_p < 0.05"); sum(new_pwas_dat$fdr_p <0.05)

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}
# Set y-axis limit above highest SNP -log10(p) by transforming the largest exponent to positive and add 1
ylim = abs(floor(log10(min(new_pwas_dat$TWAS.P)))) +1; print("ylim = "); print(ylim)

## determine sig threshold
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- new_pwas_dat %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # Using `match` to ensure we get the first occurrence if both fdr_p and max TWAS.P are repeated
  index1 <- match(TRUE, new_pwas_dat$fdr_p == result$fdr_p & new_pwas_dat$TWAS.P == result$TWAS.P)
  sig <- new_pwas_dat$TWAS.P[index1]
  print(paste("index row = ", index1))
  print(paste("p signficance threshold = ", sig))
  print(paste("-log10 sig threshold = ", -log10(sig)))
} else {
  index1 <- NA
  sig <- NA
  print("No valid significance threshold found")
}

# Make manhattan plot of PWAS results for Males
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#2E8B57", size = 0.5) +
  geom_text_repel(
    data = subset(new_pwas_dat, TWAS.P <= sig), 
    aes(label = Gene, y = -log10(TWAS.P)),
    color = "#2E8B57",
    size = 3,
    nudge_y = 0.3,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_manhplot.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights_manhplot.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights_hg38_manhplot.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights_manhplot.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# Make manhattan plot of PWAS results for Males - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#2E8B57", size = 0.5) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights_manhplot2.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights_manhplot2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights_hg38_manhplot2.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights_manhplot2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

print("Males results completed for sex-strat W23 weights")

# remove variables
if (!is.na(opt$gene_file)) {
  rm(df, col, df_name)
}

rm(axis_set, col_to_move, data_cum, fdr_p, file_list, index1, manhplot, new_pwas_dat, 
new_pwas_file, pwas_dat, pwas_file, result, sig, ylim
); print("removed male variables")

###################### PWAS results for sex-mismatched #########################

# Females GWAS with MALE weights results
print("processing females GWAS with Wingo et al. 2023 male brain protein weights results")

#Make a list of files in directory to merge ".dat" files
file_list <- list.files(path = file.path(sex_dir, fm_dir), pattern = ".dat", full.names = TRUE)

# Exclude .dat.MHC files if --mhc is TRUE
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_list <- file_list[!grepl(".dat.MHC$", file_list)]
}

#make dataframe; remove PANEL and FILE columns; clean up df before saving
pwas_dat = ldply(file_list, read_delim, show_col_types = FALSE); print("loaded .dat files")
pwas_dat = select(pwas_dat, -PANEL, -FILE); print("removed PANEL and FILE column")
names(pwas_dat)[1] = "ENSG_ID"
pwas_dat$ID = gsub("ENSG00000\\d\\d\\d\\d\\d\\d.", "", as.character(pwas_dat$ENSG_ID))
if (!is.na(opt$gene_file)) {
  pwas_dat$Gene = pwas_dat$ID
  col_to_move = c("ID", "Gene")
  pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))

  # define the replacement rules as a named vector
  df <- list(pwas_dat = pwas_dat)
  col <- "Gene"

  # apply the function to the dataframe and column if exists
  for (df_name in names(df)) {
    for (col in col) {
      if (col %in% colnames(df[[df_name]])) {
        df[[df_name]] <- replace_genes(df[[df_name]], col, replacements)
      }
    }
  }

  # unpack modified data frame back to individual variable
  list2env(df, .GlobalEnv)
} else {
   col_to_move = "ID"
   pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))
   names(pwas_dat)[2] = "Gene"
}

# if --lift option is not NA, read in the file and add hg38 coordinates for weights
if (!is.na(opt$lift)) {
  coord_dat <- fread(opt$lift, sep = '\t')

  # merge pwas_dat with the coordinate file by ENSG_ID
  pwas_dat <- pwas_dat %>%
    left_join(coord_dat %>% dplyr::select(ENSG_ID, Start_hg38, End_hg38), by = "ENSG_ID")
}

# save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  pwas_file = paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights-raw.txt")
} else {
  pwas_file = paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights-raw.txt")
}
fwrite(pwas_dat, pwas_file, row.names=F,quote=F, sep='\t'); print("saved pwas_dat file: "); print(pwas_file)

#remove NAs, convert p and P0 & P1 to numeric, calculate fdr_p & sort df on fdr_p
new_pwas_dat = na.omit(pwas_dat); print("nrow of new_pwas_dat after removing NAs"); print(nrow(new_pwas_dat))
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

if (!is.na(opt$lift)) {
  new_pwas_dat$Start_hg38 = as.numeric(new_pwas_dat$Start_hg38)
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

# if used opt$ext_mhc remove the rows with a P0 and P1 (gene start and end) within the extended region
if (!is.na(opt$ext_mhc) && as.numeric(opt$ext_mhc) > 0) {
  print(paste("Extended MHC region:", ext_mhc_start, "to", ext_mhc_end))
  new_pwas_dat <- new_pwas_dat[!(new_pwas_dat$CHR == 6 &
                         (new_pwas_dat$P0 >= ext_mhc_start & new_pwas_dat$P0 <= ext_mhc_end)), ]
}
fdr_p = p.adjust(new_pwas_dat$TWAS.P, method = "BH")
new_pwas_dat$fdr_p = fdr_p
new_pwas_dat = new_pwas_dat[order(new_pwas_dat$fdr_p),]

#create a variable for the path for where to save the summarized PWAS results as a txt file and save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  new_pwas_file = paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_male-weights.txt")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  new_pwas_file = paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights.txt")
} else {
   new_pwas_file = paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights.txt")
}
fwrite(new_pwas_dat, new_pwas_file, row.names=F,quote=F, sep='\t'); print("saved new_pwas_dat file: "); print(new_pwas_file)

# print sum of genes with fdr_p below 0.05
new_pwas_dat %>% filter(fdr_p <0.05)
print("Number of genes that pass fdr_p < 0.05"); sum(new_pwas_dat$fdr_p <0.05)

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Set y-axis limit above highest SNP -log10(p) by transforming the largest exponent to positive and add 1
ylim = abs(floor(log10(min(new_pwas_dat$TWAS.P)))) +1; print("ylim = "); print(ylim)

## determine sig threshold
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- new_pwas_dat %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # Using `match` to ensure we get the first occurrence if both fdr_p and max TWAS.P are repeated
  index1 <- match(TRUE, new_pwas_dat$fdr_p == result$fdr_p & new_pwas_dat$TWAS.P == result$TWAS.P)
  sig <- new_pwas_dat$TWAS.P[index1]
  print(paste("index row = ", index1))
  print(paste("p signficance threshold = ", sig))
  print(paste("-log10 sig threshold = ", -log10(sig)))
} else {
  index1 <- NA
  sig <- NA
  print("No valid significance threshold found")
}

# Make manhattan plot of PWAS results for females with male protein weights
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#FF7F00", size = 0.5) +
  geom_text_repel(
    data = subset(new_pwas_dat, TWAS.P <= sig), 
    aes(label = Gene, y = -log10(TWAS.P)),
    color = "#FF7F00",
    size = 3,
    nudge_y = 0.3,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_male-weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_male-weights_manhplot.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights_manhplot.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights_hg38_manhplot.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights_manhplot.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# Make manhattan plot of PWAS results for females with male protein weights - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#FF7F00", size = 0.5) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_male-weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_male-weights_manhplot2.pdf")
}  else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights_manhplot2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights_hg38_manhplot2.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights_manhplot2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

print("Females GWAS with MALE weights results completed for sex-strat W23 weights")

# remove variables
if (!is.na(opt$gene_file)) {
  rm(df, col, df_name)
}

rm(axis_set, col_to_move, data_cum, fdr_p, file_list, index1, manhplot, new_pwas_dat, 
new_pwas_file, pwas_dat, pwas_file, result, sig, ylim
); print("removed female GWAS with male weights variables")

#######################################################

# Males GWAS with FEMALE weights results
print("processing males GWAS with Wingo et al. 2023 female brain protein weights results")

#Make a list of files in directory to merge ".dat" files
file_list <- list.files(path = file.path(sex_dir, mf_dir), pattern = ".dat", full.names = TRUE)

#Exclude .dat.MHC files if --mhc is TRUE
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_list <- file_list[!grepl(".dat.MHC$", file_list)]
}

#make dataframe; remove PANEL and FILE columns; clean up df before saving
pwas_dat = ldply(file_list, read_delim, show_col_types = FALSE); print("loaded .dat files")
pwas_dat = select(pwas_dat, -PANEL, -FILE); print("removed PANEL and FILE column")
names(pwas_dat)[1] = "ENSG_ID"
pwas_dat$ID = gsub("ENSG00000\\d\\d\\d\\d\\d\\d.", "", as.character(pwas_dat$ENSG_ID))
if (!is.na(opt$gene_file)) {
  pwas_dat$Gene = pwas_dat$ID
  col_to_move = c("ID", "Gene")
  pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))

  # define the replacement rules as a named vector
  df <- list(pwas_dat = pwas_dat)
  col <- "Gene"

  # apply the function to the dataframe and column if exists
  for (df_name in names(df)) {
    for (col in col) {
      if (col %in% colnames(df[[df_name]])) {
        df[[df_name]] <- replace_genes(df[[df_name]], col, replacements)
      }
    }
  }

  # unpack modified data frame back to individual variable
  list2env(df, .GlobalEnv)
} else {
   col_to_move = "ID"
   pwas_dat = pwas_dat %>% relocate(all_of(col_to_move), .after = "ENSG_ID"); print("nrow of pwas_dat"); print(nrow(pwas_dat))
   names(pwas_dat)[2] = "Gene"
}

# if --lift option is not NA, read in the file and add hg38 coordinates for weights
if (!is.na(opt$lift)) {
  coord_dat <- fread(opt$lift, sep = '\t')

  # merge pwas_dat with the coordinate file by ENSG_ID
  pwas_dat <- pwas_dat %>%
    left_join(coord_dat %>% dplyr::select(ENSG_ID, Start_hg38, End_hg38), by = "ENSG_ID")
}

# save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  pwas_file = paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights-raw.txt")
} else {
  pwas_file = paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights-raw.txt")
}
fwrite(pwas_dat, pwas_file, row.names=F,quote=F, sep='\t'); print("saved pwas_dat file: "); print(pwas_file)

#remove NAs, convert p and P0 & P1 to numeric, calculate fdr_p & sort df on fdr_p
new_pwas_dat = na.omit(pwas_dat); print("nrow of new_pwas_dat after removing NAs"); print(nrow(new_pwas_dat))
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$P0 = as.numeric(new_pwas_dat$P0)

if (!is.na(opt$lift)) {
  new_pwas_dat$Start_hg38 = as.numeric(new_pwas_dat$Start_hg38)
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

# if used opt$ext_mhc remove the rows with a P0 and P1 (gene start and end) within the extended region
if (!is.na(opt$ext_mhc) && as.numeric(opt$ext_mhc) > 0) {
  print(paste("Extended MHC region:", ext_mhc_start, "to", ext_mhc_end))
  new_pwas_dat <- new_pwas_dat[!(new_pwas_dat$CHR == 6 &
                         (new_pwas_dat$P0 >= ext_mhc_start & new_pwas_dat$P0 <= ext_mhc_end)), ]
}
fdr_p = p.adjust(new_pwas_dat$TWAS.P, method = "BH")
new_pwas_dat$fdr_p = fdr_p
new_pwas_dat = new_pwas_dat[order(new_pwas_dat$fdr_p),]

#create a variable for the path for where to save the summarized PWAS results as a txt file and save
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  new_pwas_file = paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_female-weights.txt")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  new_pwas_file = paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights.txt")
} else {
   new_pwas_file = paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights.txt")
}
fwrite(new_pwas_dat, new_pwas_file, row.names=F,quote=F, sep='\t'); print("saved new_pwas_dat file: "); print(new_pwas_file)

# print sum of genes with fdr_p below 0.05
new_pwas_dat %>% filter(fdr_p <0.05)
print("Number of genes that pass fdr_p < 0.05"); sum(new_pwas_dat$fdr_p <0.05)

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Set y-axis limit above highest SNP -log10(p) by transforming the largest exponent to positive and add 1
ylim = abs(floor(log10(min(new_pwas_dat$TWAS.P)))) +1; print("ylim = "); print(ylim)

## determine sig threshold
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- new_pwas_dat %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # Using `match` to ensure we get the first occurrence if both fdr_p and max TWAS.P are repeated
  index1 <- match(TRUE, new_pwas_dat$fdr_p == result$fdr_p & new_pwas_dat$TWAS.P == result$TWAS.P)
  sig <- new_pwas_dat$TWAS.P[index1]
  print(paste("index row = ", index1))
  print(paste("p signficance threshold = ", sig))
  print(paste("-log10 sig threshold = ", -log10(sig)))
} else {
  index1 <- NA
  sig <- NA
  print("No valid significance threshold found")
}

# Make manhattan plot of PWAS results for males with female protein weights
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#FF7F00", size = 0.5) +
  geom_text_repel(
    data = subset(new_pwas_dat, TWAS.P <= sig), 
    aes(label = Gene, y = -log10(TWAS.P)),
    color = "#FF7F00",
    size = 3,
    nudge_y = 0.3,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_female-weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_female-weights_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights_hg38_manhplot.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights_manhplot.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights_hg38_manhplot.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights_manhplot.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# Make manhattan plot of PWAS results for males with female protein weights - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, color = as_factor(CHR))) +
  geom_point(aes(y = -log10(TWAS.P)), size = 0.3) +
  scale_color_manual(values = rep(c("grey30", "grey60"), unique(length(axis_set$CHR)))) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_point(data = subset(new_pwas_dat, TWAS.P <= sig), aes(y = -log10(TWAS.P)), color = "#FF7F00", size = 0.5) +
  geom_hline(
    yintercept = -log10(sig), color = "red",  linewidth = 0.2,
    linetype = "solid"
  ) +
  scale_size_continuous(range = c(1, 1)) +  
  labs(
    x = "Chromosome",
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_bw() +
  theme(  
    legend.position = "none",
    axis.line.x.bottom = element_line(),  
    axis.line.y.left = element_line(),  
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_female-weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_female-weights_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights_hg38_manhplot2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights_manhplot2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights_hg38_manhplot2.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights_manhplot2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

print("Males GWAS with FEMALE weights results completed for sex-strat W23 weights")

# remove variables
if (!is.na(opt$gene_file)) {
  rm(df, col, df_name)
}

rm(axis_set, col_to_move, data_cum, fdr_p, file_list, index1, manhplot, new_pwas_dat, 
new_pwas_file, pwas_dat, pwas_file, result, sig, ylim
); print("removed male GWAS with female weights variables")

###################### sex-specific filtering #########################

# load PWAS result files
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  f = fread(paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights.txt")); print("loaded F-F noMHC_ext W23 PWAS")
  m = fread(paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_weights.txt")); print("loaded M-M noMHC_ext W23 PWAS")
  f_non = fread(paste0(file.path(non_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights.txt")); print("loaded F-N noMHC_ext W23 PWAS")
  m_non = fread(paste0(file.path(non_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-sex-strat-W23_weights.txt")); print("loaded M-N noMHC_ext W23 PWAS")
  mf = fread(paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_female-weights.txt")); print("loaded M-F noMHC_ext W23 PWAS")
  fm = fread(paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_sex-strat-W23_male-weights.txt")); print("loaded F-M noMHC_ext W23 PWAS")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  f = fread(paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_sex-strat-W23_weights.txt")); print("loaded F-F noMHC W23 PWAS")
  m = fread(paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_sex-strat-W23_weights.txt")); print("loaded M-M noMHC W23 PWAS")
  f_non = fread(paste0(file.path(non_dir, f_dir, f_f), "_noMHC_non-sex-strat-W23_weights.txt")); print("loaded F-N noMHC W23 PWAS")
  m_non = fread(paste0(file.path(non_dir, m_dir, m_f), "_noMHC_non-sex-strat-W23_weights.txt")); print("loaded M-N noMHC W23 PWAS")
  mf = fread(paste0(file.path(sex_dir, mf_dir, m_f), "_noMHC_sex-strat-W23_female-weights.txt")); print("loaded M-F noMHC W23 PWAS")
  fm = fread(paste0(file.path(sex_dir, fm_dir, f_f), "_noMHC_sex-strat-W23_male-weights.txt")); print("loaded F-M noMHC W23 PWAS")
} else {
  f = fread(paste0(file.path(sex_dir, f_dir, f_f), "_sex-strat-W23_weights.txt")); print("loaded F-F W23 PWAS")
  m = fread(paste0(file.path(sex_dir, m_dir, m_f), "_sex-strat-W23_weights.txt")); print("loaded M-M W23 PWAS")
  f_non = fread(paste0(file.path(non_dir, f_dir, f_f), "_non-sex-strat-W23_weights.txt")); print("loaded F-N W23 PWAS")
  m_non = fread(paste0(file.path(non_dir, m_dir, m_f), "_non-sex-strat-W23_weights.txt")); print("loaded M-N W23 PWAS")
  mf = fread(paste0(file.path(sex_dir, mf_dir, m_f), "_sex-strat-W23_female-weights.txt")); print("loaded M-F W23 PWAS")
  fm = fread(paste0(file.path(sex_dir, fm_dir, f_f), "_sex-strat-W23_male-weights.txt")); print("loaded F-M W23 PWAS")
}
f = as.data.frame(f); m = as.data.frame(m); f_non = as.data.frame(f_non); m_non = as.data.frame(m_non)
mf = as.data.frame(mf); fm = as.data.frame(fm) 

## Female-specific results
print("processing female-specific results")

#sig genes F-F; remove rows where F-F fdr_p > 0.05
fs = f[which(f$fdr_p<0.05),]
print(paste("Number of sig genes in fs = ", nrow(fs)))

#filter 1- not sig in M-M W23 PWAS or if sig in M-M W23 PWAS z-score opposite sign
#not sig in males p > 0.05; remove rows where males p < 0.05 OR
#sig in males p< 0.05 & opposite z-score direction; remove rows where males p<0.05 and z-score has same direction
fs_1.1 <- fs[!fs$ENSG_ID %in% m[m$TWAS.P < 0.05, "ENSG_ID"] | 
             (fs$ENSG_ID %in% m[m$TWAS.P < 0.05, "ENSG_ID"] & sign(fs$TWAS.Z) != sign(m$TWAS.Z[match(fs$ENSG_ID, m$ENSG_ID)])), 
]; print(paste("number of sig genes in fs_1.1 = ", nrow(fs_1.1)))

#filter 2- not sig in M-N W23 PWAS or if sig in M-N W23 PWAS z-score opposite sign
#not sig in males non-sex-strat W23 p>0.05; remove rows where males non-sex-strat W23 p <0.05 OR
#sig in males non-sex-strat W23 & opposite z-score sign; remove rows where males non-sex-strat W23 p>0.05 and has same z-score sign
fs_2.1 <- fs_1.1[!fs_1.1$ENSG_ID %in% m_non[m_non$TWAS.P < 0.05, "ENSG_ID"] | 
             (fs_1.1$ENSG_ID %in% m_non[m_non$TWAS.P < 0.05, "ENSG_ID"] & sign(fs_1.1$TWAS.Z) != sign(m_non$TWAS.Z[match(fs_1.1$ENSG_ID, m_non$ENSG_ID)])), 
]; print(paste("number of sig genes in fs_2.1 = ", nrow(fs_2.1)))

#filter 3- not sig in GWAS M-F W23 PWAS or if sig in GWAS M-F W23 PWAS z-score opposite sign
#not sig in GWAS M-F W23 PWAS p>0.05; remove rows where GWAS M-F W23 p<0.05
#sig in GWAS M-F W23 PWAS p<0.05 & opposite z-score sign; remove rows where GWAS M-F W23 p<0.05 & z-score same sign
fs_3.1 <- fs_2.1[!fs_2.1$ENSG_ID %in% mf[mf$TWAS.P < 0.05, "ENSG_ID"] | 
             (fs_2.1$ENSG_ID %in% mf[mf$TWAS.P < 0.05, "ENSG_ID"] & sign(fs_2.1$TWAS.Z) != sign(mf$TWAS.Z[match(fs_2.1$ENSG_ID, mf$ENSG_ID)])), 
]; print(paste("number of sig genes in fs_3.1 = ", nrow(fs_3.1)))

#########################

#sig for females non-sex-strat W23; removes rows where fdr_p >0.05
fs_non = f_non[which(f_non$fdr_p<0.05),]
print(paste("number of sig genes in fs_non = ", nrow(fs_non)))

#filter 1- not sig in M-M W23 PWAS or if sig has opposite z-score
#not sig in M-M W23 PWAS p>0.05; remove rows where p < 0.05
#sig in  M-M W23 PWAS p<0.05 & opposite z-score direction; remove rows where p<0.05 & same z-score sign
fs_1.1_non <- fs_non[!fs_non$ENSG_ID %in% m[m$TWAS.P < 0.05, "ENSG_ID"] | 
             (fs_non$ENSG_ID %in% m[m$TWAS.P < 0.05, "ENSG_ID"] & sign(fs_non$TWAS.Z) != sign(m$TWAS.Z[match(fs_non$ENSG_ID, m$ENSG_ID)])), 
]; print(paste("number of sig genes in fs_1.1_non = ", nrow(fs_1.1_non)))

#filter 2- not sig in  M non-sex-strat W23 PWAS or if sig z score opposite sign
#not sig in  M non-sex-strat W23 p>0.05; remove rows where p<0.05
#sig in  M non-sex-strat W23 p<0.05 & opposite z score; remove rows where p<0.05 and z score same direction
fs_2.1_non <- fs_1.1_non[!fs_1.1_non$ENSG_ID %in% m_non[m_non$TWAS.P < 0.05, "ENSG_ID"] | 
             (fs_1.1_non$ENSG_ID %in% m_non[m_non$TWAS.P < 0.05, "ENSG_ID"] & sign(fs_1.1_non$TWAS.Z) != sign(m_non$TWAS.Z[match(fs_1.1_non$ENSG_ID, m_non$ENSG_ID)])), 
]; print(paste("number of sig genes in fs_2.1_non = ", nrow(fs_2.1_non)))

#filter 3- not sig in  M-F W23 PWAS or if sig opposite z-score
#not sig in  M-F W23 p>0.05; remove rows where p<0.05
#sig in  M-F W23 p<0.05 & opposite z score; remove rows where p<0.05 and z score same sign
fs_3.1_non <- fs_2.1_non[!fs_2.1_non$ENSG_ID %in% mf[mf$TWAS.P < 0.05, "ENSG_ID"] | 
             (fs_2.1_non$ENSG_ID %in% mf[mf$TWAS.P < 0.05, "ENSG_ID"] & sign(fs_2.1_non$TWAS.Z) != sign(mf$TWAS.Z[match(fs_2.1_non$ENSG_ID, mf$ENSG_ID)])), 
]; print(paste("number of sig genes in fs_3.1_non = ", nrow(fs_3.1_non)))

#########################
##prepping dataframe to create female-specific manhattan plot
##prep dataframe to create manhattanplot labeling only top female-specific SNPs between the two datasets, indicating by color which SNPs are unique to a dataset and which ones overlap

#make new subset of GWAS F-F W23 and GWAS F non-sex-strat W23 data with only necessary columns; add specific column to each with default value of 0 then update to matching ENSG_IDs in corresponding specific gene df
f_dat = as.data.frame(f)
f_non_dat = as.data.frame(f_non)
f_dat$specific = 0
f_non_dat$specific = 0
f_dat$specific[f_dat$ENSG_ID %in% fs_3.1$ENSG_ID & f_dat$TWAS.P %in% fs_3.1$TWAS.P] = 1; print(paste("number of f_dat$specific == 1 genes prior to determining common genes = ", sum(f_dat$specific==1)))
f_non_dat$specific[f_non_dat$ENSG_ID %in% fs_3.1_non$ENSG_ID & f_non_dat$TWAS.P %in% fs_3.1_non$TWAS.P] = 2; print(paste("number of f_non_dat$specific == 2 genes prior to determining common genes = ", sum(f_non_dat$specific==2)))

#update specific column with a value of 3 for ENSG_IDs that pass the filters in both dataset that match
common_ids = intersect(fs_3.1$ENSG_ID, fs_3.1_non$ENSG_ID); print(paste("number of genes in common between f_dat and f_non_dat = ", length(common_ids)))
f_dat$specific[f_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in f_dat that are shared with f_non_dat = ", sum(f_dat$specific==3)))
f_non_dat$specific[f_non_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in f_non_dat that are shared with f_dat = ", sum(f_non_dat$specific==3)))

#get top SNPs - identify sex-specific genes in both datasets that are unique to one dataset
f_non_only <- setdiff(fs_3.1_non$ENSG_ID, fs_3.1$ENSG_ID); print(paste("number of f_non_only genes = ", length(f_non_only)))
f_sex_only <- setdiff(fs_3.1$ENSG_ID, fs_3.1_non$ENSG_ID); print(paste("number of f_sex_genes only = ", length(f_sex_only)))

#add topSNP column to both datasets GWAS F-F W23 and GWAS F non-sex-strat W23, then update topSNP column based on genes unique to a dataset
f_dat$topSNP <- 0
f_non_dat$topSNP <- 0
f_dat$topSNP[f_dat$ENSG_ID %in% f_sex_only] <- 1
f_non_dat$topSNP[f_non_dat$ENSG_ID %in% f_non_only] <- 2

#get top SNPs - identify sex-specific genes that are in both datasets and select the gene with lowest p value for labeling
fs_merge <- merge(fs_3.1, fs_3.1_non, by = c("ENSG_ID", "ID", "Gene"), suffixes = c("_sex", "_non"))
fs_merge$topSNP <- ifelse(fs_merge$TWAS.P_sex < fs_merge$TWAS.P_non, 1,
                            ifelse(fs_merge$TWAS.P_sex > fs_merge$TWAS.P_non, 2, 3)
)
print(paste("number of genes fs_merge$topSNP == 1 prior to updating f_dat = ", sum(fs_merge$topSNP==1)))
print(paste("number of genes fs_merge$topSNP == 2 prior to updating f_non_dat = ", sum(fs_merge$topSNP==2)))
print(paste("number of genes fs_merge$topSNP == 3 prior to updating f_dat = ", sum(fs_merge$topSNP==3)))

#change value of topSNP in f_dat and f_non_dat for common genes with lowest p value for labeling
#genes that have the same p value in f_sex and f_non that are common genes, for labelling purposes add to topSNP in f_dat
f_dat$topSNP[f_dat$ENSG_ID %in% fs_merge$ENSG_ID & f_dat$ENSG_ID %in% fs_merge$ENSG_ID[fs_merge$topSNP ==1]] <- 3
print(paste("number of genes f_dat$topSNP == 1 after initial update = ", sum(f_dat$topSNP==1))); print(paste("number of genes f_dat$topSNP == 3 after initial update = ", sum(f_dat$topSNP==3)))
f_dat$topSNP[f_dat$ENSG_ID %in% fs_merge$ENSG_ID & f_dat$ENSG_ID %in% fs_merge$ENSG_ID[fs_merge$topSNP ==3]] <- 3
print(paste("final number of genes f_dat$topSNP == 1 after update = ", sum(f_dat$topSNP==1))); print(paste("final number of genes f_dat$topSNP == 3 after update = ", sum(f_dat$topSNP==3)))
f_non_dat$topSNP[f_non_dat$ENSG_ID %in% fs_merge$ENSG_ID & f_non_dat$ENSG_ID %in% fs_merge$ENSG_ID[fs_merge$topSNP ==2]] <- 3
print(paste("final number of genes f_non_dat$topSNP == 2 after update = ", sum(f_non_dat$topSNP==2))); print(paste("final number of genes f_non_dat$topSNP == 3 after update = ", sum(f_non_dat$topSNP==3)))
print(paste("total number of topSNPs = ", sum(sum(f_dat$topSNP==1)+sum(f_dat$topSNP==3)+sum(f_non_dat$topSNP==2)+sum(f_non_dat$topSNP==3))))

#create new data frame by concatenating f_sex and f_non dfs
f_dat$Discovery = "Primary"
f_non_dat$Discovery = "Secondary"
new_pwas_dat = rbind(f_dat, f_non_dat)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)

if (!is.na(opt$lift)) {
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

##determine which female-specific genes are novel vs known
#create new dataframe x from new_pwas_dat and order based on CHR and P1, then subset to fdr_p < 0.05 and where specific is greater than 0
#then add a TOP column and update it with the BEST.GWAS.ID for where topSNP is greater than 0 to identify unique topSNPs/genes
if (!is.na(opt$ad_ref)) {
  # create new variable for new_pwas_dat, filter to sex-specific genes
  x <- as.data.frame(new_pwas_dat)
  x <- x[order(x[,"CHR"], x[,"P1"]),]
  x <- x[which(x$fdr_p<0.05),]
  x <- x[which(x$specific!=0),]
  x$TOP = NA
  x$TOP <- ifelse(x$specific > 0, x$BEST.GWAS.ID, NA)

  # read in ad_ref
  ref_ad <- fread(opt$ad_ref); print("loaded AD risk loci reference file")
  ref_ad <- ref_ad[ref_ad$Sig_thresh == "gwas", ]
  x$NEW_hit <- NA
  x$LOC_con <- NA
  x$LOC_stu <- NA

  # find, if existent, the known locus next to top SNPs
  for (t in x[which(!is.na(x$TOP)),"TOP"]) {
    ra <- x[which(x$TOP==t),c("CHR","P1")]
    ref_ad_t <- ref_ad[which(ref_ad$CHR==ra$CHR),]
    if (length(ref_ad_t$CHR)>0) {
      st_dist <- 1000001; en_dist <- 1000001; 
      st_dist <- min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      st_dist_w <- which.min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      en_dist <- min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      en_dist_w <- which.min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      LOC_con <- ""; LOC_stu <- ""
      if (st_dist<=1000000 | en_dist<=1000000) {
        if (st_dist<=en_dist | st_dist==en_dist) {
          LOC_con <- ref_ad_t[st_dist_w[1],"Locus_consensus"]
          LOC_stu <- ref_ad_t[st_dist_w[1],"Study"]
        } else {
          LOC_con <- ref_ad_t[en_dist_w[1],"Locus_consensus"]
          LOC_stu <- ref_ad_t[en_dist_w[1],"Study"]
        }
        x[which(x$TOP==t),c("NEW_hit")] <- "N"
      } else {
        x[which(x$TOP==t),c("NEW_hit")] <- "Y"
      }
      x[which(x$TOP==t),c("LOC_con")] <- as.character(LOC_con)
      x[which(x$TOP==t),c("LOC_stu")] <- as.character(LOC_stu)
    } else {
      x[which(x$TOP==t),c("NEW_hit")] <- "Y"
    }
  }
  x[which(is.na(x$LOC_con)),"LOC_con"] <- ""
  x[which(is.na(x$LOC_stu)),"LOC_stu"] <- ""
  x[which(is.na(x$NEW_hit)),"NEW_hit"] <- ""

  # convert empty strings to NA then save
  x[x == ""] = NA

  # save
  if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(x, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female novelty file: ", file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"))
  } else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(x, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female novelty file: ", file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"))
  } else {
    fwrite(x, paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female novelty file: ", file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_top-female-specific-genes_novelty.txt"))
  }
}

# save combined dataframe
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(new_pwas_dat, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved combined stats dataframe: ", file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights.txt"))
    print("This dataframe is great for making additional manhattan plots")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(new_pwas_dat, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved combined stats dataframe: ", file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_weights.txt"))
    print("This dataframe is great for making additional manhattan plots")
} else {
    fwrite(new_pwas_dat, paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved combined stats dataframe: ", file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_weights.txt"))
    print("This dataframe is great for making additional manhattan plots")
}

#save gene list for LocusZoom and LocusCompare plots
if (!is.na(opt$gene_file)) {
  y = subset(x, select = c(ENSG_ID, ID, Gene, BEST.GWAS.ID, EQTL.ID, specific))
} else {
   y = subset(x, select = c(ENSG_ID, Gene, BEST.GWAS.ID, EQTL.ID, specific))
}
y$ENSG = gsub("\\..*$", "", as.character(y$ENSG_ID))
y = unique(y)
max_len <- max(sum(y$specific == 1), sum(y$specific == 2), sum(y$specific == 3))

non_ENSG.ID <- c(y$ENSG_ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_ENSG <- c(y$ENSG[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_GENE <- c(y$Gene[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_GWAS.ID <- c(y$BEST.GWAS.ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_QTL.ID <- c(y$EQTL.ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
sex_ENSG.ID <- c(y$ENSG_ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_ENSG <- c(y$ENSG[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_GENE <- c(y$Gene[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_GWAS.ID <- c(y$BEST.GWAS.ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_QTL.ID <- c(y$EQTL.ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
shared_ENSG.ID <- c(y$ENSG_ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_ENSG <- c(y$ENSG[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_GENE <- c(y$Gene[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_GWAS.ID <- c(y$BEST.GWAS.ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_QTL.ID <- c(y$EQTL.ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))

if (!is.na(opt$gene_file)) {
  non_ID <- c(y$ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
  sex_ID <- c(y$ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
  shared_ID <- c(y$ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))

  fs_genes = data.frame(
    non_ENSG.ID = non_ENSG.ID,
    non_ENSG = non_ENSG,
    non_ID = non_ID,
    non_GENE = non_GENE,
    non_GWAS.ID = non_GWAS.ID,
    non_QTL.ID = non_QTL.ID,
    sex_ENSG.ID = sex_ENSG.ID,
    sex_ENSG = sex_ENSG,
    sex_ID = sex_ID,
    sex_GENE = sex_GENE,
    sex_GWAS.ID = sex_GWAS.ID,
    sex_QTL.ID = sex_QTL.ID,
    shared_ENSG.ID = shared_ENSG.ID,
    shared_ENSG = shared_ENSG,
    shared_ID = shared_ID,
    shared_GENE = shared_GENE,
    shared_GWAS.ID = shared_GWAS.ID,
    shared_QTL.ID = shared_QTL.ID
  )
}  else {
  fs_genes = data.frame(
    non_ENSG.ID = non_ENSG.ID,
    non_ENSG = non_ENSG,
    non_GENE = non_GENE,
    non_GWAS.ID = non_GWAS.ID,
    non_QTL.ID = non_QTL.ID,
    sex_ENSG.ID = sex_ENSG.ID,
    sex_ENSG = sex_ENSG,
    sex_GENE = sex_GENE,
    sex_GWAS.ID = sex_GWAS.ID,
    sex_QTL.ID = sex_QTL.ID,
    shared_ENSG.ID = shared_ENSG.ID,
    shared_ENSG = shared_ENSG,
    shared_GENE = shared_GENE,
    shared_GWAS.ID = shared_GWAS.ID,
    shared_QTL.ID = shared_QTL.ID
  )
}
# save female-specific gene list
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(fs_genes, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_female-specific-genes-list.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female-specific gene list: ", file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_female-specific-genes-list.txt"))
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(fs_genes, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_female-specific-genes-list.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female-specific gene list: ", file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_female-specific-genes-list.txt"))
} else {
    fwrite(fs_genes, paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_female-specific-genes-list.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female-specific gene list: ", file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_female-specific-genes-list.txt"))
}

### create female-specific gene manhattan plot
## determine sig threshold for f
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- f %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  dplyr::slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # using match to ensure we get the first occurence if there are duplicates in fcr_p and TWAS.P
  index1.1 <- match(TRUE, f$fdr_p == result$fdr_p & f$TWAS.P == result$TWAS.P)
  sig1 <- f$TWAS.P[index1.1]
  print(paste("f index1.1 row = ", index1.1))
  print(paste("f p significance threshold = ", sig1))
  print(paste("f -log10 sig threshold = ", -log10(sig1)))
} else {
  index1.1 <- NA
  sig1 <- NA
  print("No valid significance threshold found")
}

## determine sig threshold for f_non
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- f_non %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  dplyr::slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # using match to ensure we get the first occurence if there are duplicates in fcr_p and TWAS.P
  index2.1 <- match(TRUE, f_non$fdr_p == result$fdr_p & f_non$TWAS.P == result$TWAS.P)
  sig2 <- f_non$TWAS.P[index2.1]
  print(paste("f_non index2.1 row = ", index2.1))
  print(paste("f_non p significance threshold = ", sig2))
  print(paste("f_non -log10 sig threshold = ", -log10(sig2)))
} else {
  index2.1 <- NA
  sig2 <- NA
  print("No valid significance threshold found")
}

#determine the lowest p value for each dataset to set the y axis limit
if (min(f$TWAS.P) < min(f_non$TWAS.P)) {
  ylim <- abs(floor(log10(min(f$TWAS.P)))) +1
  print(paste("ylim obtained from f = ", ylim))
} else {
  ylim <- abs(floor(log10(min(f_non$TWAS.P)))) +1
  print(paste("ylim obtained from f_non = ", ylim))
}

if (sig1 < sig2) {
  sig = sig2
  print("The lowest significance threshold was obtained from the secondary discovery for females")
} else if (sig2 < sig1) {
  sig = sig1
  print("The lowest significance threshold was obtained from the primary discovery for females")
} else {
   sig = sig1
   print("The lowest significance threshold is the same between both discoveries")
}
print(paste0("The lowest significance threshold for females is p = ", sig, " and -log10(p) = ", -log10(sig)))

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Create a new column for indicating color directly
new_pwas_dat$label_color <- ifelse(new_pwas_dat$topSNP == 1, "#953553",
                                   ifelse(new_pwas_dat$topSNP == 2, "#11a4a8",
                                          ifelse(new_pwas_dat$topSNP == 3, "#3154ba", NA))
)

#create new column for indicating nudge_y directly
new_pwas_dat$nudge_y <- ifelse(new_pwas_dat$topSNP == 1, 0.1,
                               ifelse(new_pwas_dat$topSNP == 2, 0.7, 
                                      ifelse(new_pwas_dat$topSNP == 3, 2, NA))
)

new_pwas_dat <- new_pwas_dat %>%
  mutate(specific = if_else(specific == 0, if_else(CHR %% 2 == 0, 4, 5), specific)
)

# with legend
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, y = -log10(TWAS.P))) +
  geom_point(data = subset(new_pwas_dat, specific == 1), aes(color = "#953553"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 2), aes(color = "#11a4a8"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 3), aes(color = "#3154ba"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 4), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey60", size = 0.3) +
  geom_point(data = subset(new_pwas_dat, specific == 5), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey30", size = 0.3) +
  geom_text_repel(
    data = new_pwas_dat[new_pwas_dat$topSNP %in% c(1, 2, 3), ],
    aes(label = Gene, y = -log10(TWAS.P)),
    color = new_pwas_dat$label_color[new_pwas_dat$topSNP %in% c(1, 2, 3)],
    size = 3,
    nudge_y = new_pwas_dat$nudge_y[new_pwas_dat$topSNP %in% c(1, 2, 3)],
    min.segment.length = 0.1,
    box.padding = unit(0.2, "lines")
  ) +
  scale_color_manual(
    breaks = c("#953553", "#11a4a8", "#3154ba"),
    limits = c("#953553", "#11a4a8", "#3154ba"),
    values = c("#953553", "#11a4a8", "#3154ba"),
    labels = c("Female-specific in primary discovery", "Female-specific in secondary discovery", "Female-specific in both discoveries"),
  ) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_hline(yintercept = -log10(sig), color = "black", linewidth = 0.2, linetype = "solid") +
  scale_size_continuous(range = c(1, 1)) +  
  labs(x = "Chromosome", y = "-log<sub>10</sub>(p)", color = "") +
  theme_bw() +
  theme(
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.background = element_rect(linewidth = 0.2),
    legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(size = 7.5),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1))
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_hg38_top-female-specific-genes.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-female-specific-genes.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_hg38_top-female-specific-genes.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_top-female-specific-genes.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_hg38_top-female-specific-genes.pdf")
} else {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_top-female-specific-genes.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# with legend - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, y = -log10(TWAS.P))) +
  geom_point(data = subset(new_pwas_dat, specific == 1), aes(color = "#953553"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 2), aes(color = "#11a4a8"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 3), aes(color = "#3154ba"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 4), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey60", size = 0.3) +
  geom_point(data = subset(new_pwas_dat, specific == 5), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey30", size = 0.3) +
  scale_color_manual(
    breaks = c("#953553", "#11a4a8", "#3154ba"),
    limits = c("#953553", "#11a4a8", "#3154ba"),
    values = c("#953553", "#11a4a8", "#3154ba"),
    labels = c("Female-specific in primary discovery", "Female-specific in secondary discovery", "Female-specific in both discoveries"),
  ) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_hline(yintercept = -log10(sig), color = "black", linewidth = 0.2, linetype = "solid") +
  scale_size_continuous(range = c(1, 1)) +  
  labs(x = "Chromosome", y = "-log<sub>10</sub>(p)", color = "") +
  theme_bw() +
  theme(
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.background = element_rect(linewidth = 0.2),
    legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(size = 7.5),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1))
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_hg38_top-female-specific-genes2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-female-specific-genes2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_hg38_top-female-specific-genes2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_top-female-specific-genes2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_hg38_top-female-specific-genes2.pdf")
} else {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_top-female-specific-genes2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

#########################
##prepping dataframe to create female-specific scatter plot

#make new subset of P2 F-F W23 and P2 F-N W23 data with only necessary columns
f_dat = dplyr::select(f, ENSG_ID, ID, Gene, CHR, P0, P1, TWAS.P, fdr_p); print("created new subset dataframe f_dat for scatterplot dataframe")
f_non_dat = dplyr::select(f_non, ENSG_ID, ID, Gene, CHR, P0, P1, TWAS.P, fdr_p); print("created new subset dataframe f_non_dat for scatterplot dataframe")

##make scatterplot with negative fake data, label all female-specific data points
##generate random numbers for missing data
#add specific column with value of 0 and create variables for identifying unique and common IDs per female-specific dataset
merged_dat = merge(f_dat, f_non_dat, by = c("ENSG_ID", "ID", "Gene", "CHR", "P0", "P1"), suffixes = c("_sex", "_non"), all = TRUE); print("created merged_dat - scatterplot dataframe")
merged_dat$specific = 0
sex_specific = setdiff(fs_3.1$ENSG_ID, fs_3.1_non$ENSG_ID); print("sex_specific genes identified for females")
non_specific = setdiff(fs_3.1_non$ENSG_ID, fs_3.1$ENSG_ID); print("non_specific genes identified for females")
common_ids = intersect(fs_3.1$ENSG_ID, fs_3.1_non$ENSG_ID); print("common_ids genes identified for females")

#add values to specific column depending on if the gene is found female-specific in sex-strat W23 only, non-sex-strat W23 only, or in both, as 1, 2, 3, respectively
merged_dat$specific[merged_dat$ENSG_ID %in% sex_specific] = 1; print("merged_dat$specific == 1 genes identified")
merged_dat$specific[merged_dat$ENSG_ID %in% non_specific] = 2; print("merged_dat$specific == 2 genes identified")
merged_dat$specific[merged_dat$ENSG_ID %in% common_ids] = 3; print("merged_dat$specific == 3 genes identified")

#indicate if NA in row for TWAS.P in sex-strat W23 or non-sex-strat W23 as 1
merged_dat$missing <- 0
merged_dat$missing <- ifelse(is.na(merged_dat$TWAS.P_sex) | is.na(merged_dat$TWAS.P_non), 1, 0); print("merged_dat$missing values for TWAS.P labeled")

#add column for calculated -log10 TWAS.P values
merged_dat$non_negLog10 <- -log10(merged_dat$TWAS.P_non); print("calculated merged_dat$non_negLog10 for TWAS.P_non")
merged_dat$sex_negLog10 <- -log10(merged_dat$TWAS.P_sex); print("calculated merged_dat$sex_negLog10 for TWAS.P_sex")

#add a random number between -1 and -0.5 for NA value in sex-strat W23 and non-sex-strat W23 TWAS.P 
merged_dat$sex_negLog10 <- ifelse(is.na(merged_dat$sex_negLog10), runif(1, min = -1, max = -0.5), merged_dat$sex_negLog10); print("random value between -1 and -0.5 given for sex_negLog10 for NA TWAS.P_sex")
merged_dat$non_negLog10 <- ifelse(is.na(merged_dat$non_negLog10), runif(1, min = -1, max = -0.5), merged_dat$non_negLog10); print("random value between -1 and -0.5 given for non_negLog10 for NA TWAS.P_non")

if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(merged_dat, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights-scatterplot.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved scatterplot stats dataframe (merged_dat): ", file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights-scatterplot.txt"))
    print("This dataframe is great for making additional scatterplots")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(merged_dat, paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_weights-scatterplot.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved scatterplot stats dataframe (merged_dat): ", file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_weights-scatterplot.txt"))
    print("This dataframe is great for making additional scatterplots")
} else {
    fwrite(merged_dat, paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_weights-scatterplot.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved scatterplot stats dataframe (merged_dat): ", file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_weights-scatterplot.txt"))
    print("This dataframe is great for making additional scatterplots")
}

#change sig values to -log10
sig1.1 <- -log10(sig1)
sig2.1 <- -log10(sig2)

#obtain axis limits; if the axis limit is odd add 1
xlim = abs(floor(log10(min(f_non_dat$TWAS.P)))) +1
ylim2 <- abs(floor(log10(min(f_dat$TWAS.P)))) +1

if (xlim %% 2 != 0) {
  xlim <- xlim + 1
}

if (ylim2 %% 2 != 0) {
  ylim2 <- ylim2 + 1
}

# convert non_negLog10 and sex_negLog10 to numeric
merged_dat$non_negLog10 = as.numeric(merged_dat$non_negLog10); merged_dat$sex_negLog10 = as.numeric(merged_dat$sex_negLog10)

#create a new column indicating color directly for labels
merged_dat$label_color <- ifelse(merged_dat$specific == 1, "#953553",
                                ifelse(merged_dat$specific == 2, "#11a4a8",
                                        ifelse(merged_dat$specific == 3, "#3154ba", NA))
)

#create new column indicating nudge_y directly for labels
merged_dat$nudge_y <- ifelse(merged_dat$specific == 1, 1,
                            ifelse(merged_dat$specific == 2, -0.4, 2)
)

#create new column indicating nudge_x directly for labels
merged_dat$nudge_x <- ifelse(merged_dat$specific == 1, -0.5,
                            ifelse(merged_dat$specific == 2, 0.8, 1.5)
)

# with legend 
scatterplot <- ggplot(merged_dat, aes(x = non_negLog10, y = sex_negLog10)) +
  geom_point(data = subset(merged_dat, missing == 1), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, missing == 0), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, specific == 1), aes(color = "#953553"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 2), aes(color = "#11a4a8"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 3), aes(color = "#3154ba"), size = 1) +
  geom_text_repel(
    data = merged_dat[merged_dat$specific %in% c(1, 2, 3), ],
    aes(label = Gene, y = sex_negLog10),
    color = merged_dat$label_color[merged_dat$specific %in% c(1, 2, 3)],
    size = 3,
    nudge_y = merged_dat$nudge_y[merged_dat$specific %in% c(1, 2, 3)],
    nudge_x = merged_dat$nudge_x[merged_dat$specific %in% c(1, 2, 3)],
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = sig1.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = sig2.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_abline(intercept = 0, slope = 1, color="#4D4D4D", linewidth=0.25) +
  labs(
    x = "Female secondary discovery AD PWAS -log<sub>10</sub>(p)",
    y = "Female primary discovery AD PWAS -log<sub>10</sub>(p)",
    color = ""
  )+
  scale_color_manual(
    values = c("#999999", "#953553", "#11a4a8", "#3154ba"),
    labels = c("Non female-specific", "Female-specific in primary discovery", "Female-specific in secondary discovery", "Female-specific in both discoveries"),
    breaks = c("#999999", "#953553", "#11a4a8", "#3154ba")
  ) +
  scale_y_continuous(limits = c(-2, ylim2), breaks = seq(0, ylim2, by = 2)) +
  scale_x_continuous(limits = c(-2, xlim), breaks = seq(0, xlim, by = 2)) +
 theme_minimal() +
  theme(  
    legend.box.background = element_rect(linewidth = 0.1),
    legend.key.width = unit(0.5, "cm"), # Adjust width of the legend box
    legend.box.margin = margin(0, 0.2, 0, 0, unit = "cm"),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    axis.ticks.y.left = element_line(),
    axis.title.y = element_markdown(size = 8),
    axis.title.x = element_markdown(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  guides(color = guide_legend()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-female-specific-genes_scatter.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_top-female-specific-genes_scatter.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_top-female-specific-genes_scatter.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(scatterplot)
dev.off()

# with legend - without gene labels if wanting to manually label afterwards
scatterplot <- ggplot(merged_dat, aes(x = non_negLog10, y = sex_negLog10)) +
  geom_point(data = subset(merged_dat, missing == 1), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, missing == 0), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, specific == 1), aes(color = "#953553"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 2), aes(color = "#11a4a8"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 3), aes(color = "#3154ba"), size = 1) +
  geom_hline(
    yintercept = sig1.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = sig2.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_abline(intercept = 0, slope = 1, color="#4D4D4D", linewidth=0.25) +
  labs(
    x = "Female secondary discovery AD PWAS -log<sub>10</sub>(p)",
    y = "Female primary discovery AD PWAS -log<sub>10</sub>(p)",
    color = ""
  )+
  scale_color_manual(
    values = c("#999999", "#953553", "#11a4a8", "#3154ba"),
    labels = c("Non female-specific", "Female-specific in primary discovery", "Female-specific in secondary discovery", "Female-specific in both discoveries"),
    breaks = c("#999999", "#953553", "#11a4a8", "#3154ba")
  ) +
  scale_y_continuous(limits = c(-2, ylim2), breaks = seq(0, ylim2, by = 2)) +
  scale_x_continuous(limits = c(-2, xlim), breaks = seq(0, xlim, by = 2)) +
 theme_minimal() +
  theme(  
    legend.box.background = element_rect(linewidth = 0.1),
    legend.key.width = unit(0.5, "cm"), # Adjust width of the legend box
    legend.box.margin = margin(0, 0.2, 0, 0, unit = "cm"),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    axis.ticks.y.left = element_line(),
    axis.title.y = element_markdown(size = 8),
    axis.title.x = element_markdown(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  guides(color = guide_legend()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-female-specific-genes_scatter2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_noMHC_non-strat_sex-strat_W23_top-female-specific-genes_scatter2.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, f_dir, f_f), "_non-strat_sex-strat_W23_top-female-specific-genes_scatter2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(scatterplot)
dev.off()

print("Female-specific gene results completed for non-sex-strat and sex-strat W23 weigthts")

rm(axis_set, common_ids, data_cum, en_dist, en_dist_w, f_dat, f_non_dat, f_non_only, f_sex_only, fs, fs_1.1, fs_1.1_non, fs_2.1, fs_2.1_non,
fs_3.1, fs_3.1_non, fs_genes, fs_merge, fs_non, index1.1, index2.1, LOC_con, LOC_stu, max_len,
merged_dat, new_pwas_dat, non_ENSG, non_ENSG.ID, non_ID, non_GENE, non_GWAS.ID, non_QTL.ID, non_specific, ra, ref_ad_t, result,
sex_ENSG, sex_ENSG.ID, sex_ID, sex_GENE, sex_GWAS.ID, sex_QTL.ID, shared_ENSG, shared_ENSG.ID, shared_ID, shared_GENE, shared_GWAS.ID, shared_QTL.ID, sex_specific,
sig, sig1, sig1.1, sig2, sig2.1, st_dist, st_dist_w, t, x, y, ylim, ylim2, xlim
); print("removed female-specific gene variables")

#######################################################################

## Male-specific results
print("processing male-specific results")

#sig genes M-M; remove rows where M-M fdr_p > 0.05
ms = m[which(m$fdr_p<0.05),]
print(paste("Number of sig genes in ms = ", nrow(ms)))

#filter 1- not sig in F-F W23 PWAS or if sig in F-F W23 PWAS z-score opposite sign
#not sig in females p > 0.05; remove rows where females p < 0.05
#sig in females p< 0.05 & opposite z-score direction; remove rows where females p<0.05 and z-score has same direction
ms_1.1 <- ms[!ms$ENSG_ID %in% f[f$TWAS.P < 0.05, "ENSG_ID"] |
             (ms$ENSG_ID %in% f[f$TWAS.P < 0.05, "ENSG_ID"]& sign(ms$TWAS.Z) != sign(f$TWAS.Z[match(ms$ENSG_ID, f$ENSG_ID)])), 
]; print(paste("number of sig genes in ms_1.1 = ", nrow(ms_1.1)))

#filter 2-
#not sig in females non-sex-strat W23 p>0.05; remove rows where females non-sex-strat W23 p <0.05
#sig in females non-sex-strat W23 & opposite z-score sign; remove rows where females non-sex-strat W23 p>0.05 and has same z-score sign
ms_2.1 <- ms_1.1[!ms_1.1$ENSG_ID %in% f_non[f_non$TWAS.P < 0.05, "ENSG_ID"] | 
             (ms_1.1$ENSG_ID %in% f_non[f_non$TWAS.P < 0.05, "ENSG_ID"] & sign(ms_1.1$TWAS.Z) != sign(f_non$TWAS.Z[match(ms_1.1$ENSG_ID, f_non$ENSG_ID)])), 
]; print(paste("number of sig genes in ms_2.1 = ", nrow(ms_2.1)))

#filter 3- not sig in F-M W23 PWAS or if sig in F-M W23 PWAS z-score opposite sign
#not sig in F-M W23 PWAS p>0.05; remove rows where F-M W23 p<0.05
#sig in F-M W23 PWAS p<0.05 & opposite z-score sign; remove rows where F-M W23 p<0.05 & z-score same sign
ms_3.1 <- ms_2.1[!ms_2.1$ENSG_ID %in% fm[fm$TWAS.P < 0.05, "ENSG_ID"] | 
             (ms_2.1$ENSG_ID %in% fm[fm$TWAS.P < 0.05, "ENSG_ID"] & sign(ms_2.1$TWAS.Z) != sign(fm$TWAS.Z[match(ms_2.1$ENSG_ID, fm$ENSG_ID)])), 
]; print(paste("number of sig genes in ms_3.1 = ", nrow(ms_3.1)))

#########################

#sig for males non-sex-strat W23; removes rows where fdr_p >0.05
ms_non = m_non[which(m_non$fdr_p<0.05),]
print(paste("number of sig genes in ms_non = ", nrow(ms_non)))

#filter 1- not sig in F-F W23 PWAS or if sig has opposite z-score
#not sig in F-F W23 PWAS p>0.05; remove rows where p < 0.05
#sig in F-F W23 PWAS p<0.05 & opposite z-score direction; remove rows where p<0.05 & same z-score sign
ms_1.1_non <- ms_non[!ms_non$ENSG_ID %in% f[f$TWAS.P < 0.05, "ENSG_ID"] | 
             (ms_non$ENSG_ID %in% f[f$TWAS.P < 0.05, "ENSG_ID"] & sign(ms_non$TWAS.Z) != sign(f$TWAS.Z[match(ms_non$ENSG_ID, f$ENSG_ID)])), 
]; print(paste("number of sig genes in ms_1.1_non = ", nrow(ms_1.1_non)))

#filter 2- not sig in F non-sex-strat W23 PWAS or if sig z score opposite sign
#not sig in F non-sex-strat W23 p>0.05; remove rows where p<0.05
#sig in F non-sex-strat W23 p<0.05 & opposite z score; remove rows where p<0.05 and z score same direction
ms_2.1_non <- ms_1.1_non[!ms_1.1_non$ENSG_ID %in% f_non[f_non$TWAS.P<0.05,"ENSG_ID"] |
              (ms_1.1_non$ENSG_ID %in% f_non[f_non$TWAS.P < 0.05,"ENSG_ID"] & sign(ms_1.1_non$TWAS.Z) != sign(f_non$TWAS.Z[match(ms_1.1_non$ENSG_ID, f_non$ENSG_ID)])),
]; print(paste("number of sig genes in ms_2.1_non = ", nrow(ms_2.1_non)))

#filter 3- not sig in F-M W23 PWAS or if sig opposite z-score
#move forward with ms_2.1_non since no genes removed in part 2 filter 2
#not sig in F-M W23 p>0.05; remove rows where p<0.05
#sig in F-M W23 p<0.05 & opposite z score; remove rows where p<0.05 and z score same sign
ms_3.1_non <- ms_2.1_non[!ms_2.1_non$ENSG_ID %in% fm[fm$TWAS.P<0.05,"ENSG_ID"] |
              (ms_2.1_non$ENSG_ID %in% fm[fm$TWAS.P < 0.05,"ENSG_ID"] & sign(ms_2.1_non$TWAS.Z) != sign(fm$TWAS.Z[match(ms_2.1_non$ENSG_ID, fm$ENSG_ID)])),
]; print(paste("number of sig genes in ms_3.1_non = ", nrow(ms_3.1_non)))

#########################
##prepping dataframe to create male-specific manhattan plot
##prep dataframe to create manhattanplot labeling only top male-specific SNPs between the two datasets, indicating by color which SNPs are unique to a dataset and which ones overlap

#make new subset of GWAS M-M W23 and GWAS M non-sex-strat W23 data with only necessary columns; add specific column to each with default value of 0 then update to matching ENSG_IDs in corresponding specific gene df
m_dat = as.data.frame(m)
m_non_dat = as.data.frame(m_non)
m_dat$specific = 0
m_non_dat$specific = 0
m_dat$specific[m_dat$ENSG_ID %in% ms_3.1$ENSG_ID & m_dat$TWAS.P %in% ms_3.1$TWAS.P] = 1;print(paste("number of m_dat$specific == 1 genes prior to determining common genes = ", sum(m_dat$specific==1)))
m_non_dat$specific[m_non_dat$ENSG_ID %in% ms_3.1_non$ENSG_ID & m_non_dat$TWAS.P %in% ms_3.1_non$TWAS.P] = 2; print(paste("number of m_non_dat$specific == 2 genes prior to determining common genes = ", sum(m_non_dat$specific==2)))

#update specific column with a value of 3 for ENSG_IDs that pass the filters in both dataset that match
common_ids = intersect(ms_3.1$ENSG_ID, ms_3.1_non$ENSG_ID); print(paste("number of genes in common between m_dat and m_non_dat = ", length(common_ids)))
m_dat$specific[m_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in m_dat that are shared with m_non_dat = ", sum(m_dat$specific==3)))
m_non_dat$specific[m_non_dat$ENSG_ID %in% common_ids] = 3; print(paste("number of genes in m_non_dat that are shaerd with m_dat = ", sum(m_non_dat$specific==3)))

#get top SNPs - identify sex-specific genes in both datasets that are unique to one dataset
m_non_only <- setdiff(ms_3.1_non$ENSG_ID, ms_3.1$ENSG_ID); print(paste("number of m_non_only genes = ", length(m_non_only)))
m_sex_only <- setdiff(ms_3.1$ENSG_ID, ms_3.1_non$ENSG_ID); print(paste("number of m_sex_only genes = ", length(m_sex_only)))

#add topSNP column to both datasets GWAS M-M W23 and GWAS M non-sex-strat W23, then update topSNP column based on genes unique to a dataset
m_dat$topSNP <- 0
m_non_dat$topSNP <- 0
m_dat$topSNP[m_dat$ENSG_ID %in% m_sex_only] <- 1
m_non_dat$topSNP[m_non_dat$ENSG_ID %in% m_non_only] <- 2

#get top SNPs - identify sex-specific genes that are in both datasets and select the gene with lowest p value for labeling
ms_non_dat <- dplyr::select(ms_3.1_non, ENSG_ID, ID, Gene, CHR, P0, P1, TWAS.P, fdr_p)
ms_sex_dat <- dplyr::select(ms_3.1, ENSG_ID, ID, Gene, CHR, P0, P1, TWAS.P, fdr_p)
ms_merge <- merge(ms_3.1, ms_3.1_non, by = c("ENSG_ID", "ID", "Gene"), suffixes = c("_sex", "_non"))
ms_merge$topSNP <- ifelse(ms_merge$TWAS.P_sex < ms_merge$TWAS.P_non, 1,
                            ifelse(ms_merge$TWAS.P_sex > ms_merge$TWAS.P_non, 2, 3)
)
print(paste("number of genes ms_merge$topSNP == 1 prior to updating m_dat = ", sum(ms_merge$topSNP==1)))
print(paste("number of genes ms_merge$topSNP == 2 prior to updating m_non_dat = ", sum(ms_merge$topSNP==2)))
print(paste("number of genes ms_merge$topSNP == 3 prior to updating m_dat = ", sum(ms_merge$topSNP==3)))

#change value of topSNP in m_dat and m_non_dat for common genes with lowest p value for labeling
#genes that have the same p value in m_sex and m_non that are common genes, for labelling purposes add to topSNP in m_dat 
m_dat$topSNP[m_dat$ENSG_ID %in% ms_merge$ENSG_ID & m_dat$ENSG_ID %in% ms_merge$ENSG_ID[ms_merge$topSNP ==1]] <- 3
print(paste("number of genes m_dat$topSNP == 1 after initial update = ", sum(m_dat$topSNP==1))); print(paste("number of genes m_dat$topSNP == 2 after initial update = ", sum(m_dat$topSNP==3)))
m_dat$topSNP[m_dat$ENSG_ID %in% ms_merge$ENSG_ID & m_dat$ENSG_ID %in% ms_merge$ENSG_ID[ms_merge$topSNP ==3]] <- 3
print(paste("final number of genes m_dat$topSNP == 1 after update = ", sum(m_dat$topSNP==1))); print(paste("final number of genes m_dat$topSNP == 3 after update = ", sum(m_dat$topSNP==3)))
m_non_dat$topSNP[m_non_dat$ENSG_ID %in% ms_merge$ENSG_ID & m_non_dat$ENSG_ID %in% ms_merge$ENSG_ID[ms_merge$topSNP ==2]] <- 3
print(paste("number of genes m_non_dat$topSNP == 2 after update = ", sum(m_non_dat$topSNP==2))); print(paste("final number of genes m_non_dat$topSNP == 3 after update = ", sum(m_non_dat$topSNP==3)))
print(paste("total number of topSNPs = ", sum(sum(m_dat$topSNP==1)+sum(m_dat$topSNP==3)+sum(m_non_dat$topSNP==2)+sum(m_non_dat$topSNP==3))))

#create new data frame by concatenating m_sex and m_non dfs
m_dat$Discovery = "Primary"
m_non_dat$Discovery = "Secondary"
new_pwas_dat = rbind(m_dat, m_non_dat)
new_pwas_dat$P1 = as.numeric(new_pwas_dat$P1)
new_pwas_dat$TWAS.P = as.numeric(new_pwas_dat$TWAS.P)

if (!is.na(opt$lift)) {
  new_pwas_dat$End_hg38 = as.numeric(new_pwas_dat$End_hg38)
}

##determine which male-specific genes are novel vs known
#create new dataframe x from new_pwas_dat and order based on CHR and P1, then subset to fdr_p < 0.05 and where specific is greater than 0
#then add a TOP column and update it with the BEST.GWAS.ID for where topSNP is greater than 0 to identify unique topSNPs/genes
if (!is.na(opt$ad_ref)) {
  # create new variable for new_pwas_dat, filter to sex-specific genes
  x <- as.data.frame(new_pwas_dat)
  x <- x[order(x[,"CHR"], x[,"P1"]),]
  x <- x[which(x$fdr_p<0.05),]
  x <- x[which(x$specific!=0),]
  x$TOP = NA
  x$TOP <- ifelse(x$specific > 0, x$BEST.GWAS.ID, NA)

  # read in ad_ref
  ref_ad <- fread(opt$ad_ref); print("loaded AD risk loci reference file")
  ref_ad <- ref_ad[ref_ad$Sig_thresh == "gwas", ]
  x$NEW_hit <- NA
  x$LOC_con <- NA
  x$LOC_stu <- NA

  # find, if existent, the known locus next to top SNPs
  for (t in x[which(!is.na(x$TOP)),"TOP"]) {
    ra <- x[which(x$TOP==t),c("CHR","P1")]
    ref_ad_t <- ref_ad[which(ref_ad$CHR==ra$CHR),]
    if (length(ref_ad_t$CHR)>0) {
      st_dist <- 1000001; en_dist <- 1000001; 
      st_dist <- min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      st_dist_w <- which.min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      en_dist <- min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      en_dist_w <- which.min(abs(ref_ad_t$GRCh37_POS - ra$P1))
      LOC_con <- ""; LOC_stu <- ""
      if (st_dist<=1000000 | en_dist<=1000000) {
        if (st_dist<=en_dist | st_dist==en_dist) {
          LOC_con <- ref_ad_t[st_dist_w[1],"Locus_consensus"]
          LOC_stu <- ref_ad_t[st_dist_w[1],"Study"]
        } else {
          LOC_con <- ref_ad_t[en_dist_w[1],"Locus_consensus"]
          LOC_stu <- ref_ad_t[en_dist_w[1],"Study"]
        }
        x[which(x$TOP==t),c("NEW_hit")] <- "N"
      } else {
        x[which(x$TOP==t),c("NEW_hit")] <- "Y"
      }
      x[which(x$TOP==t),c("LOC_con")] <- as.character(LOC_con)
      x[which(x$TOP==t),c("LOC_stu")] <- as.character(LOC_stu)
    } else {
      x[which(x$TOP==t),c("NEW_hit")] <- "Y"
    }
  }
  x[which(is.na(x$LOC_con)),"LOC_con"] <- ""
  x[which(is.na(x$LOC_stu)),"LOC_stu"] <- ""
  x[which(is.na(x$NEW_hit)),"NEW_hit"] <- ""

  # convert empty strings to NA then save
  x[x == ""] = NA

  # save
  if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(x, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female novelty file: ", file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"))
  } else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(x, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female novelty file: ", file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"))
  } else {
    fwrite(x, paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved female novelty file: ", file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_top-male-specific-genes_novelty.txt"))
  }
}

# save combined dataframe
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(new_pwas_dat, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved combined stats dataframe: ", file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights.txt"))
    print("This dataframe is great for making additional manhattan plots")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(new_pwas_dat, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved combined stats dataframe: ", file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_weights.txt"))
    print("This dataframe is great for making additional manhattan plots")
} else {
    fwrite(new_pwas_dat, paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_weights.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved combined stats dataframe: ", file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_weights.txt"))
    print("This dataframe is great for making additional manhattan plots")
}

#save gene list for LocusZoom and LocusCompare plots
if (!is.na(opt$gene_file)) {
  y = subset(x, select = c(ENSG_ID, ID, Gene, BEST.GWAS.ID, EQTL.ID, specific))
} else {
   y = subset(x, select = c(ENSG_ID, Gene, BEST.GWAS.ID, EQTL.ID, specific))
}
y$ENSG = gsub("\\..*$", "", as.character(y$ENSG_ID))
y = unique(y)
max_len <- max(sum(y$specific == 1), sum(y$specific == 2), sum(y$specific == 3))

non_ENSG.ID <- c(y$ENSG_ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_ENSG <- c(y$ENSG[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_GENE <- c(y$Gene[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_GWAS.ID <- c(y$BEST.GWAS.ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
non_QTL.ID <- c(y$EQTL.ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
sex_ENSG.ID <- c(y$ENSG_ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_ENSG <- c(y$ENSG[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_GENE <- c(y$Gene[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_GWAS.ID <- c(y$BEST.GWAS.ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
sex_QTL.ID <- c(y$EQTL.ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
shared_ENSG.ID <- c(y$ENSG_ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_ENSG <- c(y$ENSG[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_GENE <- c(y$Gene[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_GWAS.ID <- c(y$BEST.GWAS.ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))
shared_QTL.ID <- c(y$EQTL.ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))

if (!is.na(opt$gene_file)) {
  non_ID <- c(y$ID[y$specific == 2], rep(NA, max_len - sum(y$specific == 2)))
  sex_ID <- c(y$ID[y$specific == 1], rep(NA, max_len - sum(y$specific == 1)))
  shared_ID <- c(y$ID[y$specific == 3], rep(NA, max_len - sum(y$specific == 3)))

  ms_genes = data.frame(
    non_ENSG.ID = non_ENSG.ID,
    non_ENSG = non_ENSG,
    non_ID = non_ID,
    non_GENE = non_GENE,
    non_GWAS.ID = non_GWAS.ID,
    non_QTL.ID = non_QTL.ID,
    sex_ENSG.ID = sex_ENSG.ID,
    sex_ENSG = sex_ENSG,
    sex_ID = sex_ID,
    sex_GENE = sex_GENE,
    sex_GWAS.ID = sex_GWAS.ID,
    sex_QTL.ID = sex_QTL.ID,
    shared_ENSG.ID = shared_ENSG.ID,
    shared_ENSG = shared_ENSG,
    shared_ID = shared_ID,
    shared_GENE = shared_GENE,
    shared_GWAS.ID = shared_GWAS.ID,
    shared_QTL.ID = shared_QTL.ID
  )
}  else {
  ms_genes = data.frame(
    non_ENSG.ID = non_ENSG.ID,
    non_ENSG = non_ENSG,
    non_GENE = non_GENE,
    non_GWAS.ID = non_GWAS.ID,
    non_QTL.ID = non_QTL.ID,
    sex_ENSG.ID = sex_ENSG.ID,
    sex_ENSG = sex_ENSG,
    sex_GENE = sex_GENE,
    sex_GWAS.ID = sex_GWAS.ID,
    sex_QTL.ID = sex_QTL.ID,
    shared_ENSG.ID = shared_ENSG.ID,
    shared_ENSG = shared_ENSG,
    shared_GENE = shared_GENE,
    shared_GWAS.ID = shared_GWAS.ID,
    shared_QTL.ID = shared_QTL.ID
  )
}

#save male-specific gene list
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(ms_genes, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_male-specific-genes-list.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved male-specific gene list: ", file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_male-specific-genes-list.txt"))
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(ms_genes, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_male-specific-genes-list.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved male-specific gene list: ", file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_male-specific-genes-list.txt"))
} else {
    fwrite(ms_genes, paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_male-specific-genes-list.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved male-specific gene list: ", file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_male-specific-genes-list.txt"))
}

### create male-specific gene manhattan plot
## determine sig threshold for m
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- m %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  dplyr::slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # using match to ensure we get the first occurence if there are duplicates in fdr_p and TWAS.P
  index1.1 <- match(TRUE, m$fdr_p == result$fdr_p & m$TWAS.P == result$TWAS.P)
  sig1 <- m$TWAS.P[index1.1]
  print(paste("m index1.1 row = ", index1.1))
  print(paste("m p significance threshold = ", sig1))
  print(paste("m -log10 sig threshold = ", -log10(sig1)))
} else {
  index1.1 <- NA
  sig1 <- NA
  print("No valid significance threshold found")
}

## determine sig threshold for m_non
# Find the row with the highest fdr_p less than 0.05, and among ties, the highest TWAS.P
result <- m_non %>%
  filter(fdr_p < 0.05) %>%
  arrange(desc(fdr_p), desc(TWAS.P)) %>%
  dplyr::slice(1)

# Check if result is non-empty and obtain the index
if (nrow(result) > 0) {
  # using match to ensure we get the first occurence if there are duplicates in fcr_p and TWAS.P
  index2.1 <- match(TRUE, m_non$fdr_p == result$fdr_p & m_non$TWAS.P == result$TWAS.P)
  sig2 <- m_non$TWAS.P[index2.1]
  print(paste("m_non index2.1 row = ", index2.1))
  print(paste("m_non p significance threshold = ", sig2))
  print(paste("m_non -log10 sig threshold = ", -log10(sig2)))
} else {
  index2.1 <- NA
  sig2 <- NA
  print("No valid significance threshold found")
}

if (sig1 < sig2) {
  sig = sig2
  print("The lowest significance threshold was obtained from the secondary discovery for males")
} else if (sig2 < sig1) {
  sig = sig1
  print("The lowest significance threshold was obtained from the primary discovery for males")
} else {
   sig = sig1
   print("The lowest significance threshold is the same between both discoveries")
}
print(paste0("The lowest significance threshold for males is p = ", sig, " and -log10(p) = ", -log10(sig)))

#determine the lowest p value for each dataset to set the y axis limit
if (min(m$TWAS.P) < min(m_non$TWAS.P)) {
  ylim <- abs(floor(log10(min(m$TWAS.P)))) +1
  print(paste("ylim obtained from m = ", ylim))
} else {
  ylim <- abs(floor(log10(min(m_non$TWAS.P)))) +1
  print(paste("ylim obtained from m_non = ", ylim))
}

## prep for manhattan plot- use P0 and P1 if --lift option is not utilized otherwise use Start_hg38 and End_hg38
if (!is.na(opt$lift)) {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
  data_cum = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(max_P1 = max(End_hg38)) %>%
  dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
  dplyr::select(CHR, P1_add)

  # merge above data with original dataframe and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
  new_pwas_dat = new_pwas_dat %>%
  inner_join(data_cum, by = "CHR") %>%
  dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative hg38 bp position calculated for each row")

  # get center position of each chromosome for labeling x-axis in the middle of each chromosome
  axis_set = new_pwas_dat %>%
  group_by(CHR) %>%
  dplyr::summarize(center = mean(P1_cum)); print("calculated center hg38 position for each chromosome")
} else {
  # combine chromosome number and chromosome end position (bp) to create column with position to use for x-axis putting chromosome 1 first and 22 last
   data_cum = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(max_P1 = max(P1)) %>%
   dplyr::mutate(P1_add = lag(cumsum(max_P1), default = 0)) %>%
   dplyr::select(CHR, P1_add)

   # Merge above data with original data frame and calculate cumulative bp position for each row by adding relative position and adding factor together creating P1_cum column
   new_pwas_dat = new_pwas_dat %>%
   inner_join(data_cum, by = "CHR") %>%
   dplyr::mutate(P1_cum = P1 + P1_add); print("cumulative bp position calculated for each row")

   # Get center position of each chromosome for labeling x-axis in the middle of each chromosome
   axis_set = new_pwas_dat %>%
   group_by(CHR) %>%
   dplyr::summarize(center = mean(P1_cum)); print("calculated center position for each chromosome")
}

# Create a new column for indicating color directly
new_pwas_dat$label_color <- ifelse(new_pwas_dat$topSNP == 1, "#db7161",
                                   ifelse(new_pwas_dat$topSNP == 2, "#3680b6",
                                          ifelse(new_pwas_dat$topSNP == 3, "#152955", NA))
)

#create new column for indicating nudge_y directly
new_pwas_dat$nudge_y <- ifelse(new_pwas_dat$topSNP == 1, 0.5,
                               ifelse(new_pwas_dat$topSNP == 2, 1.5, 
                                      ifelse(new_pwas_dat$topSNP == 3, 0.1, NA))
)

new_pwas_dat <- new_pwas_dat %>%
  mutate(specific = if_else(specific == 0, if_else(CHR %% 2 == 0, 4, 5), specific)
)

# with legend
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, y = -log10(TWAS.P))) +
  geom_point(data = subset(new_pwas_dat, specific == 1), aes(color = "#db7161"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 2), aes(color = "#3680b6"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 3), aes(color = "#152955"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 4), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey60", size = 0.3) +
  geom_point(data = subset(new_pwas_dat, specific == 5), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey30", size = 0.3) +
  geom_text_repel(
    data = new_pwas_dat[new_pwas_dat$topSNP %in% c(1, 2, 3), ],
    aes(label = Gene, y = -log10(TWAS.P)),
    color = new_pwas_dat$label_color[new_pwas_dat$topSNP %in% c(1, 2, 3)],
    size = 3,
    nudge_y = new_pwas_dat$nudge_y[new_pwas_dat$topSNP %in% c(1, 2, 3)],
    min.segment.length = 0.1,
    box.padding = unit(0.2, "lines")
  ) +
  scale_color_manual(
    breaks = c("#db7161", "#3680b6", "#152955"),
    limits = c("#db7161", "#3680b6", "#152955"),
    values = c("#db7161", "#3680b6", "#152955"),
    labels = c("Male-specific in primary discovery", "Male-specific in secondary discovery", "Male-specific in both discoveries"),
  ) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_hline(yintercept = -log10(sig), color = "black", linewidth = 0.2, linetype = "solid") +  scale_size_continuous(range = c(1, 1)) +  
  labs(x = "Chromosome", y = "-log<sub>10</sub>(p)", color = "") +
  theme_bw() +
  theme(
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.background = element_rect(linewidth = 0.2),
    legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(size = 7.5),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1))
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_hg38_top-male-specific-genes.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-male-specific-genes.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_hg38_top-male-specific-genes.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_top-male-specific-genes.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_hg38_top-male-specific-genes.pdf")
} else {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_top-male-specific-genes.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

# with legend - without gene labels if wanting to manually label afterwards
manhplot <- ggplot(new_pwas_dat, aes(x = P1_cum, y = -log10(TWAS.P))) +
  geom_point(data = subset(new_pwas_dat, specific == 1), aes(color = "#db7161"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 2), aes(color = "#3680b6"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 3), aes(color = "#152955"), size = 0.5) +
  geom_point(data = subset(new_pwas_dat, specific == 4), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey60", size = 0.3) +
  geom_point(data = subset(new_pwas_dat, specific == 5), aes(x = P1_cum, y = -log10(TWAS.P)), color = "grey30", size = 0.3) +
  scale_color_manual(
    breaks = c("#db7161", "#3680b6", "#152955"),
    limits = c("#db7161", "#3680b6", "#152955"),
    values = c("#db7161", "#3680b6", "#152955"),
    labels = c("Male-specific in primary discovery", "Male-specific in secondary discovery", "Male-specific in both discoveries"),
  ) +
  scale_x_continuous(expand = c(0, 0.1), 
                     label = axis_set$CHR[!(axis_set$CHR %in% c("19", "21"))],
                     breaks = axis_set$center[!(axis_set$CHR %in% c("19", "21"))]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim), breaks = seq(0, ylim, by = 2)) +
  geom_hline(yintercept = -log10(sig), color = "black", linewidth = 0.2, linetype = "solid") +  scale_size_continuous(range = c(1, 1)) +  
  labs(x = "Chromosome", y = "-log<sub>10</sub>(p)", color = "") +
  theme_bw() +
  theme(
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.ticks.y.left = element_line(),
    axis.title = element_text(size = rel(0.8)),
    axis.title.y = element_markdown(),
    axis.title.x = element_markdown(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box.background = element_rect(linewidth = 0.2),
    legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(size = 7.5),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 1))
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0 && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_hg38_top-male-specific-genes2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-male-specific-genes2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_hg38_top-male-specific-genes2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_top-male-specific-genes2.pdf")
} else if (!is.na(opt$lift)) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_hg38_top-male-specific-genes2.pdf")
} else {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_top-male-specific-genes2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(manhplot)
dev.off()

#########################
##prepping dataframe to create male-specific scatter plot

# make new subset of M-M W23 and M-N W23 data with only necessary columns
m_dat = dplyr::select(m, ENSG_ID, ID, Gene, CHR, P0, P1, TWAS.P, fdr_p); print("created new subset dataframe m_dat for scatterplot dataframe")
m_non_dat = dplyr::select(m_non, ENSG_ID, ID, Gene, CHR, P0, P1, TWAS.P, fdr_p); print("created new subset dataframe m_non_dat for scatterplot dataframe")

##make scatterplot with negative fake data, label all male-specific data points
##generate random numbers for missing data
#add specific column with value of 0 and create variables for identifying unique and common IDs per male-specific dataset
merged_dat = merge(m_dat, m_non_dat, by = c("ENSG_ID", "ID", "Gene", "CHR", "P0", "P1"), suffixes = c("_sex", "_non"), all = TRUE); print("created merged_dat - scatterplot dataframe")
merged_dat$specific = 0
sex_specific = setdiff(ms_3.1$ENSG_ID, ms_3.1_non$ENSG_ID); print("sex_specific genes identified for males")
non_specific = setdiff(ms_3.1_non$ENSG_ID, ms_3.1$ENSG_ID); print("non_specific genes identified for males")
common_ids = intersect(ms_3.1$ENSG_ID, ms_3.1_non$ENSG_ID); print("common_ids genes identified for males")

#add values to specific column depending on if the gene is found male-specific in sex-strat W23 only, non-sex-strat W23 only, or in both, as 1, 2, 3, respectively
merged_dat$specific[merged_dat$ENSG_ID %in% sex_specific] = 1; print("merged_dat$specific == 1 genes identified")
merged_dat$specific[merged_dat$ENSG_ID %in% non_specific] = 2; print("merged_dat$specific == 2 genes identified")
merged_dat$specific[merged_dat$ENSG_ID %in% common_ids] = 3; print("merged_dat$specific == 3 genes identified")

#indicate if NA in row for TWAS.P in sex-strat W23 or non-sex-strat W23 as 1
merged_dat$missing <- 0
merged_dat$missing <- ifelse(is.na(merged_dat$TWAS.P_sex) | is.na(merged_dat$TWAS.P_non), 1, 0); print("merged_dat$missing values for TWAS.P labeled")

#add column for calculated -log10 TWAS.P values
merged_dat$non_negLog10 <- -log10(merged_dat$TWAS.P_non); print("calculated merged_dat$non_negLog10 for TWAS.P_non")
merged_dat$sex_negLog10 <- -log10(merged_dat$TWAS.P_sex); print("calculated merged_dat$sex_negLog10 for TWAS.P_sex")

#add a random number between -1 and -0.5 for NA value in sex-strat W23 and non-sex-strat W23 TWAS.P 
merged_dat$sex_negLog10 <- ifelse(is.na(merged_dat$sex_negLog10), runif(1, min = -1, max = -0.5), merged_dat$sex_negLog10); print("random value between -1 and -0.5 given for sex_negLog10 for NA TWAS.P_sex")
merged_dat$non_negLog10 <- ifelse(is.na(merged_dat$non_negLog10), runif(1, min = -1, max = -0.5), merged_dat$non_negLog10); print("random value between -1 and -0.5 given for non_negLog10 for NA TWAS.P_non")

if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
    fwrite(merged_dat, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights-scatterplot.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved scatterplot stats dataframe (merged_dat): ", file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_weights-scatterplot.txt"))
    print("This dataframe is great for making additional scatterplots")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
    fwrite(merged_dat, paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_weights-scatterplot.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved scatterplot stats dataframe (merged_dat): ", file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_weights-scatterplot.txt"))
    print("This dataframe is great for making additional scatterplots")
} else {
    fwrite(merged_dat, paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_weights-scatterplot.txt"), col.names = T, row.names = F, quote = F, na = "NA", sep = '\t')
    print(paste0("saved scatterplot stats dataframe (merged_dat): ", file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_weights-scatterplot.txt"))
    print("This dataframe is great for making additional scatterplots")
}

#change sig values to -log10
sig1.1 <- -log10(sig1)
sig2.1 <- -log10(sig2)

#obtain axis limits; if the axis limit is odd add 1
xlim = abs(floor(log10(min(m_non_dat$TWAS.P)))) +1
ylim2 <- abs(floor(log10(min(m_dat$TWAS.P)))) +1

if (xlim %% 2 != 0) {
  xlim <- xlim + 1
}

if (ylim2 %% 2 != 0) {
  ylim2 <- ylim2 + 1
}

# convert non_negLog10 and sex_negLog10 to numeric
merged_dat$non_negLog10 = as.numeric(merged_dat$non_negLog10); merged_dat$sex_negLog10 = as.numeric(merged_dat$sex_negLog10)

#create a new column indicating color directly for labels
merged_dat$label_color <- ifelse(merged_dat$specific == 1, "#db7161",
                                ifelse(merged_dat$specific == 2, "#3680b6",
                                        ifelse(merged_dat$specific == 3, "#152955", NA))
)

#create new column indicating nudge_y directly for labels
merged_dat$nudge_y <- ifelse(merged_dat$specific == 1, 0.8,
                            ifelse(merged_dat$specific == 2, -0.1, 
                            ifelse(merged_dat$specific ==3, 0.8, NA))
)

#create new column indicating nudge_x directly for labels
merged_dat$nudge_x <- ifelse(merged_dat$specific == 1, 0,
                            ifelse(merged_dat$specific == 2, 0.8, 
                            ifelse(merged_dat$specific == 3, 0.6, NA))
)

# with legend
scatterplot <- ggplot(merged_dat, aes(x = non_negLog10, y = sex_negLog10)) +
  geom_point(data = subset(merged_dat, missing == 1), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, missing == 0), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, specific == 1), aes(color = "#db7161"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 2), aes(color = "#3680b6"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 3), aes(color = "#152955"), size = 1) +
  geom_text_repel(
    data = merged_dat[merged_dat$specific %in% c(1, 2, 3), ],
    aes(label = Gene, y = sex_negLog10),
    color = merged_dat$label_color[merged_dat$specific %in% c(1, 2, 3)],
    size = 3,
    nudge_y = merged_dat$nudge_y[merged_dat$specific %in% c(1, 2, 3)],
    nudge_x = merged_dat$nudge_x[merged_dat$specific %in% c(1, 2, 3)],
    force = 10,
    min.segment.length = 0.1,
    box.padding = unit(0.1, "lines")
  ) +
  geom_hline(
    yintercept = sig1.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = sig2.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_abline(intercept = 0, slope = 1, color="#4D4D4D", linewidth=0.25) +
  labs(
    x = "Male secondary discovery AD PWAS -log<sub>10</sub>(p)",
    y = "Male primary discovery AD PWAS -log<sub>10</sub>(p)",
    color = ""
  )+
  scale_color_manual(
    values = c("#999999", "#db7161", "#3680b6", "#152955"),
    labels = c("Non male-specific", "Male-specific in primary discovery", "Male-specific in secondary discovery", "Male-specific in both discoveries"),
    breaks = c("#999999", "#db7161", "#3680b6", "#152955")
  ) +
  scale_y_continuous(limits = c(-2, ylim2), breaks = seq(0, ylim2, by = 2)) +
  scale_x_continuous(limits = c(-2, xlim), breaks = seq(0, xlim, by = 2)) +
 theme_minimal() +
  theme(  
    legend.box.background = element_rect(linewidth = 0.1),
    legend.key.width = unit(0.5, "cm"), # Adjust width of the legend box
    legend.box.margin = margin(0, 0.2, 0, 0, unit = "cm"),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    axis.ticks.y.left = element_line(),
    axis.title.y = element_markdown(size = 8),
    axis.title.x = element_markdown(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  guides(color = guide_legend()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-male-specific-genes_scatter.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_top-male-specific-genes_scatter.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_top-male-specific-genes_scatter.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(scatterplot)
dev.off()

# with legend - without gene labels if wanting to manually label afterwards
scatterplot <- ggplot(merged_dat, aes(x = non_negLog10, y = sex_negLog10)) +
  geom_point(data = subset(merged_dat, missing == 1), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, missing == 0), aes(color = "#999999"), size = 0.8) +
  geom_point(data = subset(merged_dat, specific == 1), aes(color = "#db7161"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 2), aes(color = "#3680b6"), size = 1) +
  geom_point(data = subset(merged_dat, specific == 3), aes(color = "#152955"), size = 1) +
  geom_hline(
    yintercept = sig1.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = sig2.1, color = "red",  linewidth=0.3,
    linetype = "dashed"
  ) +
  geom_abline(intercept = 0, slope = 1, color="#4D4D4D", linewidth=0.25) +
  labs(
    x = "Male secondary discovery AD PWAS -log<sub>10</sub>(p)",
    y = "Male primary discovery AD PWAS -log<sub>10</sub>(p)",
    color = ""
  )+
  scale_color_manual(
    values = c("#999999", "#db7161", "#3680b6", "#152955"),
    labels = c("Non male-specific", "Male-specific in primary discovery", "Male-specific in secondary discovery", "Male-specific in both discoveries"),
    breaks = c("#999999", "#db7161", "#3680b6", "#152955")
  ) +
  scale_y_continuous(limits = c(-2, ylim2), breaks = seq(0, ylim2, by = 2)) +
  scale_x_continuous(limits = c(-2, xlim), breaks = seq(0, xlim, by = 2)) +
 theme_minimal() +
  theme(  
    legend.box.background = element_rect(linewidth = 0.1),
    legend.key.width = unit(0.5, "cm"), # Adjust width of the legend box
    legend.box.margin = margin(0, 0.2, 0, 0, unit = "cm"),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    axis.ticks.y.left = element_line(),
    axis.title.y = element_markdown(size = 8),
    axis.title.x = element_markdown(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  guides(color = guide_legend()
)
if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE" && !is.na(opt$ext_mhc) && opt$ext_mhc > 0) {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_ext", opt$ext_mhc, "Mb_non-strat_sex-strat_W23_top-male-specific-genes_scatter2.pdf")
} else if (!is.na(opt$no_mhc) && opt$no_mhc == "TRUE") {
  file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_noMHC_non-strat_sex-strat_W23_top-male-specific-genes_scatter2.pdf")
} else {
   file_name <- paste0(file.path(sex_dir, m_dir, m_f), "_non-strat_sex-strat_W23_top-male-specific-genes_scatter2.pdf")
}
pdf(file = file_name, height = 4.2949, width = 7.08661)
plot(scatterplot)
dev.off()

