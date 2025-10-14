#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(CMplot)
})

# ----------------------------- CLI --------------------------------------------
option_list <- list(
  make_option(c("--path"), type = "character",
              help = "Top folder containing Brain/CSF/<sex>/<discovery>/ subfolders (Primary, Secondary). REQUIRED."),
  make_option(c("--remove_chr"), type = "character", default = "",
              help = "Comma-separated chromosomes to remove, e.g. 6,19,21 (optional)."),
  make_option(c("--remove_genes"), type = "character", default = "",
              help = "Comma-separated gene symbols to remove (no spaces), e.g. G1,G2,G3 (optional).")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$path)) stop("--path is required", call. = FALSE)
top_path <- opt$path

# ----------------------------- Config -----------------------------------------
tissues <- c("Brain", "CSF")
sexes   <- c("Female", "Male")
discovs <- c("Primary", "Secondary")

female_col <- "#700000"
male_col   <- "#3780b5"

parse_remove_chr <- function(x) {
  if (!nzchar(x)) integer(0) else suppressWarnings(unique(as.integer(strsplit(x, ",", fixed=TRUE)[[1]])))
}
parse_remove_genes <- function(x) {
  if (!nzchar(x)) character(0) else toupper(strsplit(x, ",", fixed=TRUE)[[1]])
}
remove_chr_vec   <- parse_remove_chr(opt$remove_chr)
remove_genes_vec <- parse_remove_genes(opt$remove_genes)

# ------------------------- File discovery & readers ---------------------------
find_merged_file <- function(dir_td) {
  if (!dir.exists(dir_td)) stop("Directory not found: ", dir_td)
  cand <- list.files(dir_td, pattern = "_merged_with_FDR\\.txt$", full.names = TRUE)
  if (!length(cand)) stop("No '*_merged_with_FDR.txt' in: ", dir_td)
  cand[[1]]
}

read_merged <- function(fp) {
  df <- read_tsv(fp, show_col_types = FALSE)
  nm <- names(df); low <- tolower(nm)
  pick <- function(cands, prefer_exact = NULL) {
    idx <- which(low %in% tolower(cands))
    if (!length(idx)) return(NA_integer_)
    if (!is.null(prefer_exact)) {
      exact <- which(nm %in% prefer_exact)
      exact <- intersect(idx, exact)
      if (length(exact)) return(exact[1])
    }
    idx[1]
  }
  chr_i  <- pick(c("chr","chrom","chromosome"), prefer_exact = "CHR")
  gene_i <- pick(c("gene","GENE"), prefer_exact = "GENE")
  p_i    <- pick(c("pwas.p","pwas_p","p","pvalue","p_value","pwas.p.","pwas.pval","pwas.pvalue"), prefer_exact = "PWAS.P")
  fdr_i  <- pick(c("fdr","FDR"), prefer_exact = "FDR")
  p0_i   <- pick(c("p0","start","txstart","gene_start"), prefer_exact = "P0")
  p1_i   <- pick(c("p1","end","txend","gene_end"),     prefer_exact = "P1")
  miss <- c(CHR=chr_i, GENE=gene_i, P=p_i, FDR=fdr_i, P0=p0_i, P1=p1_i)
  miss <- names(miss)[is.na(unlist(miss))]
  if (length(miss)) stop("Missing columns in ", fp, ": ", paste(miss, collapse = ", "))

  tibble(
    CHR  = suppressWarnings(as.integer(df[[chr_i]])),
    GENE = as.character(df[[gene_i]]),
    P    = suppressWarnings(as.numeric(df[[p_i]])),
    FDR  = suppressWarnings(as.numeric(df[[fdr_i]])),
    P0   = suppressWarnings(as.numeric(df[[p0_i]])),
    P1   = suppressWarnings(as.numeric(df[[p1_i]]))
  ) |>
    mutate(POS = as.integer(round((P0 + P1)/2))) |>
    select(CHR, POS, GENE, P, FDR)
}

find_and_load_ss_genes <- function(dir_td) {
  if (is.null(dir_td) || !dir.exists(dir_td)) return(character(0))
  files <- list.files(dir_td, pattern = "(?i)sex.?specific.*\\.(tsv|txt|csv)$", full.names = TRUE)
  if (!length(files)) return(character(0))
  fp <- files[[1]]
  dat <- tryCatch({
    if (grepl("\\.csv$", fp, ignore.case = TRUE)) read_csv(fp, show_col_types = FALSE)
    else read_tsv(fp, show_col_types = FALSE)
  }, error = function(e) NULL)
  if (is.null(dat)) return(character(0))
  nm <- tolower(names(dat))
  if (!("gene" %in% nm)) return(character(0))
  gene_col <- names(dat)[match("gene", nm)]
  if ("sexspecific" %in% nm) {
    ss_col <- names(dat)[match("sexspecific", nm)]
    dat <- dat[!is.na(dat[[ss_col]]) & dat[[ss_col]] == TRUE, , drop = FALSE]
  }
  unique(toupper(na.omit(as.character(dat[[gene_col]]))))
}

get_ss_genes_for <- function(tissue, sex, file_map) {
  dir_prim <- dirname(file_map[[tissue]][[sex]][["Primary"]])
  dir_seco <- dirname(file_map[[tissue]][[sex]][["Secondary"]])
  unique(c(find_and_load_ss_genes(dir_prim), find_and_load_ss_genes(dir_seco)))
}

sanitize_p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p <= 0] <- 1e-300
  p[p >= 1] <- 1 - 1e-16
  p
}

collapse_minP <- function(df) {
  df %>% group_by(GENE, CHR, POS) %>%
    summarise(P = suppressWarnings(min(P, na.rm = TRUE)), .groups = "drop")
}

# ------------------------ Map files, read, prune early ------------------------
file_map <- list()
for (t in tissues) for (s in sexes) for (d in discovs) {
  dir_td <- file.path(top_path, t, s, d)
  file_map[[t]][[s]][[d]] <- find_merged_file(dir_td)
}
message("All required '*_merged_with_FDR.txt' files found.")

dat <- list()
for (t in tissues) for (s in sexes) for (d in discovs) {
  df <- read_merged(file_map[[t]][[s]][[d]])
  if (length(remove_chr_vec))   df <- filter(df, !(CHR %in% remove_chr_vec))
  if (length(remove_genes_vec)) df <- filter(df, !(toupper(GENE) %in% remove_genes_vec))
  dat[[t]][[s]][[d]] <- df
}

# ------------------------ Sex-specific unions per ring ------------------------
ss_FB <- get_ss_genes_for("Brain", "Female", file_map)
ss_MB <- get_ss_genes_for("Brain", "Male",   file_map)
ss_FC <- get_ss_genes_for("CSF",   "Female", file_map)
ss_MC <- get_ss_genes_for("CSF",   "Male",   file_map)

# Remove genes from sex-specific sets as well (if requested)
if (length(remove_genes_vec)) {
  ss_FB <- setdiff(ss_FB, remove_genes_vec)
  ss_MB <- setdiff(ss_MB, remove_genes_vec)
  ss_FC <- setdiff(ss_FC, remove_genes_vec)
  ss_MC <- setdiff(ss_MC, remove_genes_vec)
}

# ------------------------ Build FOUR ring blocks (ALL genes, grey base) -------
FB_df <- collapse_minP(bind_rows(dat$Brain$Female))
MB_df <- collapse_minP(bind_rows(dat$Brain$Male))
FC_df <- collapse_minP(bind_rows(dat$CSF$Female))
MC_df <- collapse_minP(bind_rows(dat$CSF$Male))

mk_block_all <- function(df, which_ring) {
  if (!nrow(df)) return(tibble(SNP=character(0), CHR=integer(0), BP=integer(0), p1=numeric(0), p2=numeric(0), p3=numeric(0), p4=numeric(0), GENE=character(0)))
  snp <- paste(df$GENE, df$CHR, df$POS, which_ring, sep=":")
  blk <- tibble(
    SNP = snp, CHR = df$CHR, BP = df$POS,
    p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, p4 = NA_real_,
    GENE = df$GENE
  )
  if (which_ring == "FB") blk$p4 <- sanitize_p(df$P) # FB to inner (p4) after inversion
  if (which_ring == "MB") blk$p3 <- sanitize_p(df$P) # MB to p3
  if (which_ring == "FC") blk$p2 <- sanitize_p(df$P) # FC to p2
  if (which_ring == "MC") blk$p1 <- sanitize_p(df$P) # MC to outer (p1)
  blk
}

block_FB <- mk_block_all(FB_df, "FB")
block_MB <- mk_block_all(MB_df, "MB")
block_FC <- mk_block_all(FC_df, "FC")
block_MC <- mk_block_all(MC_df, "MC")

dat_plot <- bind_rows(block_FB, block_MB, block_FC, block_MC) %>%
  filter(!is.na(CHR), CHR >= 1, CHR <= 22)

if (!nrow(dat_plot)) stop("No points to plot after filtering.", call. = FALSE)

# ------------------------ Thresholds from FDR<=0.05 (raw p) -------------------
fdr_p_cut <- function(df) {
  df <- df %>% filter(!is.na(FDR), !is.na(P))
  if (!nrow(df)) return(NA_real_)
  hit <- df %>% filter(FDR <= 0.05)
  if (!nrow(hit)) return(NA_real_)
  suppressWarnings(max(hit$P, na.rm = TRUE))
}
# Order thresholds to match inverted ring order: p1=MC, p2=FC, p3=MB, p4=FB
th_list <- list(
  fdr_p_cut(bind_rows(dat$CSF$Male)),    # MC -> outermost
  fdr_p_cut(bind_rows(dat$CSF$Female)),  # FC
  fdr_p_cut(bind_rows(dat$Brain$Male)),  # MB
  fdr_p_cut(bind_rows(dat$Brain$Female)) # FB -> innermost
)

# ------------------------ Colors & chr labels ---------------------------------
# Grey base; only highlights will be colored
base_cols <- matrix(rep("grey60", 4), nrow = 4, byrow = TRUE)

present_chr <- sort(unique(dat_plot$CHR))
chr_labels  <- as.character(present_chr)

# Ring-specific highlights (color only sex-specific genes)
pick_ring_highlights <- function(block_df, ss_upper, col) {
  idx <- which(toupper(block_df$GENE) %in% ss_upper)
  if (!length(idx)) return(list(genes = character(0), cols = character(0)))
  list(genes = block_df$SNP[idx], cols = rep(col, length(idx)))
}
h_FB <- pick_ring_highlights(block_FB, ss_FB, female_col)
h_MB <- pick_ring_highlights(block_MB, ss_MB, male_col)
h_FC <- pick_ring_highlights(block_FC, ss_FC, female_col)
h_MC <- pick_ring_highlights(block_MC, ss_MC, male_col)

hl_genes <- c(h_MC$genes, h_FC$genes, h_MB$genes, h_FB$genes) # order doesn't matter
hl_cols  <- c(h_MC$cols,  h_FC$cols,  h_MB$cols,  h_FB$cols)

# ------------------------ Plot helper (white background) ----------------------
plot_circos_to_device <- function(device = c("pdf","png"), filename, width, height, dpi = 300) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(paste0(filename, ".pdf"), width = width, height = height, useDingbats = FALSE, paper = "special", bg = "white")
  } else {
    png(paste0(filename, ".png"), width = width * dpi, height = height * dpi, res = dpi, bg = "white")
  }
  par(bg = "white")
  on.exit(dev.off(), add = TRUE)
  
  op <- par(no.readonly = TRUE)
  par(
    bg  = "white",
    mar = c(0, 0, 0, 0),   # no outer margins
    oma = c(0, 0, 0, 0),   # no outer margin area
    plt = c(0, 1, 0, 1),   # use full device area
    xaxs = "i", yaxs = "i" # no extra axis padding
  )
  on.exit(par(op), add = TRUE)

  CMplot(
    dat_plot %>% select(SNP, CHR, BP, p1, p2, p3, p4),
    type = "p",
    plot.type = "c",
    col = base_cols,
    pch = c(17, 17, 19, 19),
    band = 0.5,
    ylim = c(0, 15),
    r = 0.1,
    cir.axis = FALSE,
    outward = FALSE,                # still plotting inside; ring order inverted via p1..p4 mapping
    chr.labels = chr_labels,
    cir.axis.col = "black",
    cir.chr.h = 1,
    chr.den.col = "grey60",
    threshold = th_list,
    threshold.lwd = c(1.5),
    threshold.col = rep("black", 4),
    highlight = if (length(hl_genes)) hl_genes else NULL,
    highlight.col = if (length(hl_cols))  hl_cols  else NULL,
    highlight.cex = 1.5,            # larger colored dots
    amplify = FALSE,
    signal.line = 1,
    LOG10 = TRUE,
    file.output = FALSE,            # we control the device to force white background
    verbose = TRUE
  )
}

# ------------------------ COMBINED PLOT (PDF + PNG) ---------------------------
message("Rendering combined circos on white background (PDF + PNG)…")
plot_circos_to_device("pdf", "CircosePlot_combined", width = 8, height = 8, dpi = 300)
plot_circos_to_device("png", "CircosePlot_combined", width = 8, height = 8, dpi = 300)
message("Done: CircosePlot_combined.pdf/png")

# ------------------------ CSF-only (printing disabled) ------------------------
# Uncomment to generate CSF-only single-figure (white background)
# dat_csf <- block_FC %>% select(SNP, CHR, BP, p3)
# plot_circos_to_device("pdf", "CircosePlot_CSF", width = 6, height = 6, dpi = 300)
# plot_circos_to_device("png", "CircosePlot_CSF", width = 6, height = 6, dpi = 300)

# -------------

# Internally CercosPlot_V1.R