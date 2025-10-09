#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(grid)     # for unit()
  library(tools)    # for toTitleCase()
})

# ----------------------------- CLI: --path, --remove_chr ----------------------
option_list <- list(
  make_option(c("--path"),        type = "character", help = "Top folder containing Brain/CSF/<sex>/<discovery>/ subfolders."),
  make_option(c("--remove_chr"),  type = "character", default = "",
              help = "Comma-separated chromosomes to remove in plots, e.g. '6,19,21'. Optional.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$path)) stop("--path is required")
top_path <- opt$path

# Where plots go
plots_dir <- file.path(top_path, "Scatter_Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Globals ----------------------------------------
tissues  <- c("Brain", "CSF")
sexes    <- c("Female", "Male")
discovs  <- c("Primary", "Secondary")

parse_remove_chr <- function(x) {
  if (is.null(x) || x == "") return(integer(0))
  unique(as.integer(unlist(strsplit(x, "[,;\\s]+"))))
}
remove_chr_vec <- parse_remove_chr(opt$remove_chr)

# ----------------------------- File discovery ---------------------------------
find_merged_file <- function(dir_td) {
  if (!dir.exists(dir_td)) stop("Directory not found: ", dir_td)
  cand <- list.files(dir_td, pattern = "_merged_with_FDR\\.txt$", full.names = TRUE)
  if (!length(cand)) stop("No '*_merged_with_FDR.txt' in: ", dir_td)
  cand[[1]]
}

# Robust reader that returns only the columns needed for plotting,
# resolving duplicate 'Gene'/'GENE' issues cleanly.
read_merged <- function(fp) {
  df <- readr::read_tsv(fp, show_col_types = FALSE)

  nm <- names(df)
  low <- tolower(nm)

  pick_col <- function(cands, prefer_exact = NULL) {
    idx <- which(low %in% tolower(cands))
    if (!length(idx)) return(NA_integer_)
    if (!is.null(prefer_exact)) {
      exact <- which(nm %in% prefer_exact)
      exact <- intersect(idx, exact)
      if (length(exact)) return(exact[1])
    }
    idx[1]
  }

  chr_i  <- pick_col(c("chr", "chrom", "chromosome"), prefer_exact = "CHR")
  gene_i <- pick_col(c("gene", "GENE"), prefer_exact = "GENE")
  p_i    <- pick_col(c("PWAS.p", "PWAS_p", "p", "pvalue", "p_value", "PWAS.P"), prefer_exact = "PWAS.P")
  fdr_i  <- pick_col(c("fdr", "FDR"), prefer_exact = "FDR")

  miss <- c(CHR = chr_i, GENE = gene_i, `PWAS.P` = p_i, FDR = fdr_i)
  miss <- names(miss)[is.na(unlist(miss))]
  if (length(miss)) stop("Missing required columns in ", fp, ": ", paste(miss, collapse = ", "))

  tibble::tibble(
    CHR     = suppressWarnings(as.integer(df[[chr_i]])),
    GENE    = as.character(df[[gene_i]]),
    `PWAS.P` = suppressWarnings(as.numeric(df[[p_i]])),
    FDR     = suppressWarnings(as.numeric(df[[fdr_i]]))
  )
}

# Build map to Primary/Secondary merged files
file_map <- list()
for (t in tissues) for (s in sexes) for (d in discovs) {
  dir_td <- file.path(top_path, t, s, d)
  file_map[[t]][[s]][[d]] <- find_merged_file(dir_td)
}
message("All required '*_merged_with_FDR.txt' files found.")

# Read all data
dat <- list()
for (t in tissues) for (s in sexes) for (d in discovs) {
  dat[[t]][[s]][[d]] <- read_merged(file_map[[t]][[s]][[d]])
}

# ============================= Steps 1, 2, 3 ==================================
# (1) Find & load sex-specific genes from disk (Primary & Secondary folders) ----
find_and_load_ss_genes <- function(dir_td) {
  if (is.null(dir_td) || !dir.exists(dir_td)) return(character(0))

  # Flexible filename patterns; tighten if you use a fixed name
  files <- list.files(
    dir_td,
    pattern = "(?i)sex.?specific.*\\.(tsv|txt|csv)$",
    full.names = TRUE
  )
  if (!length(files)) return(character(0))
  fp <- files[[1]]

  # Read (TSV/TXT/CSV agnostic)
  dat <- tryCatch({
    if (grepl("\\.csv$", fp, ignore.case = TRUE)) {
      readr::read_csv(fp, show_col_types = FALSE)
    } else {
      readr::read_tsv(fp, show_col_types = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(dat)) return(character(0))

  nm <- tolower(names(dat))
  if (!("gene" %in% nm)) return(character(0))
  gene_col <- names(dat)[match("gene", nm)]

  # If SexSpecific exists, require TRUE; else, take all rows
  if ("sexspecific" %in% nm) {
    ss_col <- names(dat)[match("sexspecific", nm)]
    dat <- dat[!is.na(dat[[ss_col]]) & dat[[ss_col]] == TRUE, , drop = FALSE]
  }
  unique(stats::na.omit(as.character(dat[[gene_col]])))
}

# Union of sex-specific genes from Primary & Secondary folders
get_ss_genes_for <- function(tissue, sex, file_map) {
  dir_prim <- dirname(file_map[[tissue]][[sex]][["Primary"]])
  dir_seco <- dirname(file_map[[tissue]][[sex]][["Secondary"]])
  g1 <- find_and_load_ss_genes(dir_prim)
  g2 <- find_and_load_ss_genes(dir_seco)
  unique(c(g1, g2))
}

# (2) Make scatter accepts file_map and uses union list -------------------------
.fdr_p_cut <- function(df) {
  df <- df |> dplyr::filter(!is.na(FDR), !is.na(`PWAS.P`))
  if (!nrow(df)) return(NA_real_)
  hit <- df |> dplyr::filter(FDR <= 0.05)
  if (!nrow(hit)) return(NA_real_)
  max(hit$`PWAS.P`, na.rm = TRUE)
}

make_scatter_primary_vs_secondary <- function(tissue, sex,
                                              dat_list = dat,
                                              file_map  = file_map,
                                              remove_chr = remove_chr_vec) {
  stopifnot(tissue %in% names(dat_list))
  stopifnot(sex    %in% names(dat_list[[tissue]]))
  stopifnot(all(c("Primary","Secondary") %in% names(dat_list[[tissue]][[sex]])))

  df_prim <- dat_list[[tissue]][[sex]][["Primary"]]
  df_seco <- dat_list[[tissue]][[sex]][["Secondary"]]

  need <- c("CHR","GENE","PWAS.P","FDR")
  miss1 <- setdiff(need, names(df_prim))
  miss2 <- setdiff(need, names(df_seco))
  if (length(miss1) || length(miss2)) {
    stop("Scatter needs columns CHR, GENE, PWAS.P, FDR in both Primary and Secondary. ",
         "Missing: Primary{", paste(miss1, collapse=","), "} Secondary{", paste(miss2, collapse=","), "}")
  }

  # Apply --remove_chr
  if (length(remove_chr)) {
    df_prim <- df_prim |> dplyr::filter(!(CHR %in% remove_chr))
    df_seco <- df_seco |> dplyr::filter(!(CHR %in% remove_chr))
  }

  # Merge by GENE (keep CHR from Primary for reference)
  merged <- df_prim |>
  dplyr::select(CHR, GENE, P_primary = `PWAS.P`, FDR_primary = FDR) |>
  dplyr::full_join(
    df_seco |> dplyr::select(GENE, P_secondary = `PWAS.P`, FDR_secondary = FDR),
    by = "GENE"
  ) |>
  dplyr::mutate(
    was_missing_primary   = is.na(P_primary),
    was_missing_secondary = is.na(P_secondary),
    # Assign neutral dummy P for missing side so the point appears near the axis
    P_primary   = ifelse(was_missing_primary,   1, P_primary),
    P_secondary = ifelse(was_missing_secondary, 1, P_secondary),
    imputed_any = was_missing_primary | was_missing_secondary
  )

  if (!nrow(merged)) {
    warning(sprintf("No overlap after merge for %s %s; skipping plot.", tissue, sex))
    return(list(plot = NULL, out_base = NULL))
  }

  # Axis thresholds
  pcut_primary   <- .fdr_p_cut(df_prim)
  pcut_secondary <- .fdr_p_cut(df_seco)
  thr_x <- if (is.na(pcut_secondary)) NA_real_ else -log10(pcut_secondary)
  thr_y <- if (is.na(pcut_primary))   NA_real_ else -log10(pcut_primary)

  # (1+2) Sex-specific genes (union from Primary & Secondary)
  ss_genes <- get_ss_genes_for(tissue, sex, file_map)
  message(sprintf("[Scatter] %s %s — sex-specific genes loaded: %d",
                  tissue, toTitleCase(sex), length(ss_genes)))
  
  merged <- merged %>%
  dplyr::group_by(GENE) %>%
  dplyr::slice_min(order_by = pmin(P_primary, P_secondary, na.rm = TRUE),
                   n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()
  
  merged <- merged |>
    dplyr::mutate(
      x = -log10(P_secondary),
      y = -log10(P_primary),
      label_me = GENE %in% ss_genes
    )

  lab_col <- if (tolower(sex) == "female") "#700000" else "#3780b5"
  
  genes_to_exclude <- c("CEBPZOS", "PACSIN1")
  merged <- dplyr::filter(merged, !GENE %in% genes_to_exclude)

  g <- ggplot(merged, aes(x, y)) +
  geom_point(data = subset(merged, imputed_any & !label_me),
             inherit.aes = TRUE, shape = 17, size = 1.3, alpha = 0.7, colour = "grey50") +
  geom_point(data = subset(merged, !label_me & !imputed_any),
             inherit.aes = TRUE, alpha = 0.7, size = 1.3, colour = "grey50") +
  geom_point(data = subset(merged, label_me & !imputed_any),
             inherit.aes = TRUE, alpha = 0.95, size = 1.6, colour = lab_col) +
  geom_point(data = subset(merged, label_me & imputed_any),
             inherit.aes = TRUE, shape = 17, alpha = 0.95, size = 1.6, colour = lab_col)

  lab_df <- subset(merged, label_me & !is.na(GENE) & GENE != "")
  if (nrow(lab_df)) {
    g <- g + ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = GENE),
      size = 3.0,
      colour = lab_col,
      box.padding = grid::unit(0.15, "lines"),
      max.overlaps = Inf,
      min.segment.length = 0
    )
  }

  if (!is.na(thr_x)) g <- g + geom_vline(xintercept = thr_x, linetype = "dashed", alpha = 0.8, colour = "black", linewidth = 0.4)
  if (!is.na(thr_y)) g <- g + geom_hline(yintercept = thr_y, linetype = "dashed", alpha = 0.8, colour = "black", linewidth = 0.4)
  g <- g + geom_abline(slope = 1, intercept = 0, linewidth = 0.4, alpha = 0.8, colour = "grey40")

  xaxis <- if (tolower(sex) == "female") "Female secondary discovery AD PWAS −log10(P)" else "Male secondary discovery AD PWAS −log10(P)"
  yaxis <- if (tolower(sex) == "female") "Female primary discovery AD PWAS −log10(P)" else "Male primary discovery AD PWAS −log10(P)"
  g <- g + xlab(xaxis) + ylab(yaxis) + theme(axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6))
  
  # ---- Axes & theme polish ----
  # 1) pick a common max so x and y share the same range
  max_xy <- max(
    merged$x, merged$y,
    ifelse(is.na(thr_x), -Inf, thr_x),
    ifelse(is.na(thr_y), -Inf, thr_y),
    na.rm = TRUE
  )
  # round up to the next even number so breaks can be {0,2,4,...}
  lim_max       <- max(0, max_xy)
  lim_max_even  <- 2 * ceiling(lim_max / 2)
  brks_every_2  <- seq(0, lim_max_even, by = 2)
  
  g <- g +
    # 1) white background
    theme_bw(base_size = 12) +
    # 2) both axes start at 0, ticks every 2
    scale_x_continuous(limits = c(-1, lim_max_even), breaks = brks_every_2, expand = c(0,0)) +
    scale_y_continuous(limits = c(-1, lim_max_even), breaks = brks_every_2, expand = c(0,0)) +
    # 4) thin light-grey major gridlines (and no minor grid)
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.25),
      panel.grid.minor = element_blank()
    )
  
  out_base <- file.path(plots_dir, paste0("Scatter_", tissue, "_", sex, "_Primary_vs_Secondary"))
  list(plot = g, out_base = out_base)
}

# (3) Printing & saving loop passes file_map -----------------------------------
plots_to_write <- list()
for (t in tissues) for (s in sexes) {
  res <- make_scatter_primary_vs_secondary(tissue = t, sex = s, file_map = file_map)
  if (is.null(res$plot)) next
  print(res$plot)
  ggsave(filename = paste0(res$out_base, ".pdf"), plot = res$plot,
         width = 7, height = 6, units = "in", dpi = 300)
  ggsave(filename = paste0(res$out_base, ".png"), plot = res$plot,
         width = 7, height = 6, units = "in", dpi = 300)
  plots_to_write[[paste(t, s, sep = "_")]] <- res$out_base
}

message("Scatter plots written for: ", paste(names(plots_to_write), collapse = ", "))
message("All plots completed. Output in: ", plots_dir)

# Internal Ref PWAS_ScatterPlot_BrainCSF_V5.R
