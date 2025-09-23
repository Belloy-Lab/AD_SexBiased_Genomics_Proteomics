library(data.table)
library(plyr)
library(dplyr)
library(openxlsx)
library(readr)
library(methods)
library(optparse)

# Define options
option_list = list(
  make_option("--file", action="store", type='character',
              help="Name of the GWAS summary statistics file [required]"),
  make_option("--dir", action="store", type='character',
              help="Directory path for summary statistics file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$file) || is.na(opt$dir)) {
  stop("Both --file and --dir options are required.")
}

input_path = file.path(opt$dir, opt$file)
df = fread(input_path, header = F, sep = '\t')
df_clean = na.omit(df)
fwrite(df_clean, input_path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
