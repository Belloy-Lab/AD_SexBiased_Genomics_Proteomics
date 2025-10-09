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
})

# ----------------------------- CLI: --path, --remove_chr ----------------------
option_list <- list(
  make_option(c("--path"), type = "character",
              help = "Top-level path containing Brain/ and CSF/ folders",
              metavar = "character"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Chromosomes to EXCLUDE from Manhattan (comma/space separated; e.g., '6,19,21'). Applies to Manhattan plots only.",
              metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$path)) stop("ERROR: --path is required.\n", call. = FALSE)
top_path <- normalizePath(opt$path, mustWork = TRUE)

# Parse --remove_chr into integer vector 1..22
parse_remove_chr <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(integer(0))
  toks <- unlist(strsplit(x, "[,;\\s]+"))
  toks <- gsub("^chr", "", toks, ignore.case = TRUE)
  toks[toks %in% c("X","x","Y","y","M","m")] <- NA
  suppressWarnings({
    vals <- as.integer(toks)
  })
  vals <- vals[!is.na(vals)]
  vals[vals >= 1 & vals <= 22]
}
remove_chr_vec <- parse_remove_chr(opt$remove_chr)

# ---------------------------- Helpers: filesystem -----------------------------
# Find exactly one *_merged_with_FDR.txt in a directory; otherwise stop with error
find_merged_file <- function(dir_path) {
  files <- list.files(dir_path, pattern = "_merged_with_FDR\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    stop(paste0("Missing merged file in: ", dir_path, "\n",
                "Please run the Sex_specific_gene_list.R script to make combined PWAS results before running this script."),
         call. = FALSE)
  }
  if (length(files) > 1) {
    stop(paste0("Found multiple '*_merged_with_FDR.txt' files in: ", dir_path, "\n",
                paste0(" - ", basename(files), collapse = "\n"),
                "\nPlease keep exactly ONE merged file per folder."),
         call. = FALSE)
  }
  files[[1]]
}

# ---- Sex-specific list helpers ----
find_ss_list <- function(folder) {
  # look for either _sex_specific_gene_list.txt or .tsv
  cand <- list.files(folder,
                     pattern = "_sex_specific_gene_list\\.(txt|tsv)$",
                     full.names = TRUE, ignore.case = TRUE)
  if (length(cand) == 0) return(NA_character_)
  if (length(cand) > 1) stop("Multiple *_sex_specific_gene_list.[txt|tsv] found in: ", folder)
  cand[[1]]
}

read_ss_genes <- function(filepath) {
  if (is.na(filepath)) return(character(0))
  ss <- suppressWarnings(
    readr::read_tsv(filepath, col_types = readr::cols(.default = readr::col_character()))
  )
  if (!("GENE" %in% names(ss))) {
    warning("Sex-specific list missing 'GENE' column: ", filepath)
    return(character(0))
  }
  # trim whitespace and drop empty/NA; return unique vector
  genes <- unique(stats::na.omit(trimws(ss$GENE)))
  genes[nzchar(genes)]
}

# Create Plots directory (mkdir -p)
plots_dir <- file.path(top_path, "Mahattan_Plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- Helpers: data parsing ---------------------------
# Normalize CHR to integer 1..22 (drop others later)
norm_chr <- function(chr) {
  x <- as.character(chr)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("X","x")] <- "23"
  x[x %in% c("Y","y")] <- "24"
  suppressWarnings(as.integer(x))
}

# Read a merged_with_FDR file and coerce key columns
read_merged <- function(path) {
  df <- suppressMessages(readr::read_tsv(path, guess_max = 100000))
  needed <- c("CHR", "P0", "PWAS.P", "FDR")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop(paste0("File lacks required columns: ", paste(missing, collapse = ", "),
                "\nFile: ", path), call. = FALSE)
  }
  df <- df %>%
    mutate(
      CHR = norm_chr(CHR),
      P0  = suppressWarnings(as.numeric(P0)),
      P1  = if ("P1" %in% names(.)) suppressWarnings(as.numeric(P1)) else NA_real_,
      PWAS.P = suppressWarnings(as.numeric(PWAS.P)),
      FDR    = suppressWarnings(as.numeric(FDR)),
      GENE = if ("GENE" %in% names(.)) as.character(GENE) else NA_character_,
      ID   = if ("ID"   %in% names(.)) as.character(ID)   else NA_character_
    ) %>%
    filter(!is.na(CHR), CHR >= 1, CHR <= 22, !is.na(P0), !is.na(PWAS.P), !is.na(FDR)) %>%
    mutate(pos = ifelse(!is.na(P1), floor((P0 + P1)/2), P0))
  return(df)
}

# Compute cumulative genomic position and chromosome centers for axis
add_cumulative_position <- function(df) {
  chr_lengths <- df %>%
    group_by(CHR) %>%
    summarize(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR)

  offsets <- c(0, cumsum(chr_lengths$chr_len)[-nrow(chr_lengths)])
  chr_lengths <- chr_lengths %>% mutate(offset = offsets)

  df2 <- df %>% left_join(chr_lengths, by = "CHR") %>% mutate(P_cum = pos + offset)

  axis_set <- df2 %>%
    group_by(CHR) %>%
    summarize(center = mean(range(P_cum, na.rm = TRUE)), .groups = "drop")

  list(df = df2, axis = axis_set)
}

# p* cutoff: largest PWAS.P among rows with FDR <= fdr_thr (NA if none)
largest_p_among_fdr <- function(df, fdr_thr = 0.05) {
  sig <- df %>% filter(!is.na(FDR), FDR <= fdr_thr)
  if (nrow(sig) == 0) return(NA_real_)
  suppressWarnings(max(sig$PWAS.P, na.rm = TRUE))
}

# Save both PDF and PNG
save_plot_dual <- function(plot, out_base, width = 7.2, height = 4.4, units = "in", dpi = 300) {
  ggsave(paste0(out_base, ".pdf"), plot = plot, width = width, height = height, units = units)
  ggsave(paste0(out_base, ".png"), plot = plot, width = width, height = height, units = units, dpi = dpi)
  message("Saved: ", out_base, ".pdf and .png")
}

# ------------------------------ Plot: Manhattan -------------------------------
plot_manhattan <- function(df, title, out_base, fdr_thr = 0.05, remove_chr = integer(0), y_max_override = NULL) {

  # 1) Apply chromosome removal first
  if (length(remove_chr)) {
    df <- dplyr::filter(df, !CHR %in% remove_chr)
  }
  if (nrow(df) == 0) { message("[manhattan] No data for: ", title); return(invisible(NULL)) }

  # 2) Ensure required columns and compute derived fields
  if (!("PWAS.P" %in% names(df))) stop("[manhattan] Missing PWAS.P to compute -log10(P).")
  if (!("P0" %in% names(df)))    stop("[manhattan] Missing P0 (and optional P1) to compute positions.")
  if (!("P1" %in% names(df)))    df$P1 <- df$P0
  if (!("GENE" %in% names(df)))  df$GENE <- NA_character_

  # --- EXCLUDE GENES FROM PLOTTING ONLY ---
  genes_to_exclude <- c("CEBPZOS", "PACSIN1")
  df <- dplyr::filter(df, !(GENE %in% genes_to_exclude))
  # ----------------------------------------

  if (!("ID"   %in% names(df)))  df$ID   <- NA_character_
  if (!("CHR"  %in% names(df)))  stop("[manhattan] Missing CHR column.")
  df <- dplyr::mutate(df,
                      CHR   = as.integer(CHR),
                      PWAS.P = as.numeric(PWAS.P),
                      P0    = as.numeric(P0),
                      P1    = as.numeric(P1),
                      pos   = ifelse(!is.na(P0) & !is.na(P1), (P0 + P1)/2, P0),
                      logp  = -log10(pmax(PWAS.P, .Machine$double.xmin)))

  # 3) Rebuild cumulative positions and x-axis centers fresh
  chr_sizes <- df %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(cumstart = dplyr::lag(cumsum(chr_len), default = 0))
  df <- df %>%
    dplyr::left_join(chr_sizes, by = "CHR") %>%
    dplyr::mutate(pos_cum = pos + cumstart) %>%
    dplyr::arrange(CHR, pos_cum)
  centers <- df %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(center = (min(pos_cum) + max(pos_cum))/2, .groups = "drop")
  attr(df, "axis_centers") <- centers

  # 4) Threshold and y-span
  p_star <- {
    ok <- !is.na(df$FDR) & df$FDR <= fdr_thr & !is.na(df$PWAS.P)
    if (!any(ok)) NA_real_ else max(df$PWAS.P[ok], na.rm = TRUE)
  }

  # Base on strongest signal: top = max(-log10(P)) + 2, but keep old safeguards
  max_logp <- max(df$logp, na.rm = TRUE)
  ymax_calc <- max(max_logp + 1,
                   if (is.finite(-log10(p_star))) -log10(p_star) else 0,
                   4,
                   na.rm = TRUE)

  ymax <- if (!is.null(y_max_override)) y_max_override else ymax_calc

  # 5) Which points to label
  lab_df <- df %>%
    dplyr::filter(!is.na(FDR) & FDR <= fdr_thr) %>%
    dplyr::mutate(label_col = dplyr::coalesce(GENE, ID))

  # 7) Plot
   p <- ggplot2::ggplot(df, ggplot2::aes(x = pos_cum, y = logp)) +
    ggplot2::geom_point(ggplot2::aes(color = as.factor(CHR %% 2)),
                        size = 0.8, alpha = 0.9, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("#aaaaaa", "#555555")) +
    ggplot2::geom_point(data = lab_df, color = "black", size = 1.1, alpha = 0.95) +
    ggplot2::geom_hline(yintercept = if (is.finite(-log10(p_star))) -log10(p_star) else Inf,
                        alpha = 1, linetype = "dashed", color = "#B22222", linewidth = 0.4) + # linetype = "dashed",
    ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = label_col),
      size = 3,
      color = "black",
      segment.color = "black",
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    ggplot2::scale_x_continuous(
      name   = "Chromosome",
      breaks = attr(df, "axis_centers")$center,
      labels = attr(df, "axis_centers")$CHR,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      name   = "-log10(P)",
      limits = c(0, ceiling(ymax)),
      breaks = seq(0, ceiling(ymax/2)*2, by = 2),
      expand = c(0, 0)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      # ↓ shrink title/axes/legend by 10%
      axis.title  = ggplot2::element_text(size = ggplot2::rel(0.7)),
      axis.text   = ggplot2::element_text(size = ggplot2::rel(0.5)),
      legend.title= ggplot2::element_text(size = ggplot2::rel(0.7)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.7))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ceiling(ymax)), clip = "on") #+

  save_plot_dual(p, paste0(out_base, "_manhattan"))
}

# --- Manhattan (second set) — label only genes in sex-specific list; labels+dots in sex color ---
plot_manhattan_sexcolor <- function(df, title, out_base, sex, label_genes,
                                    fdr_thr = 0.05, remove_chr = integer(0), y_max_override = NULL) {

  # require a non-empty gene list; otherwise skip
  if (missing(label_genes) || length(label_genes) == 0) {
    message("[manhattan-SS] No sex-specific label genes provided for: ", title)
    return(invisible(NULL))
  }

  # 1) Apply chromosome removal first
  if (length(remove_chr)) df <- dplyr::filter(df, !CHR %in% remove_chr)
  if (nrow(df) == 0) { message("[manhattan-SS] No data for: ", title); return(invisible(NULL)) }

  # 2) Ensure required columns and compute derived fields
  if (!("PWAS.P" %in% names(df))) stop("[manhattan-SS] Missing PWAS.P.")
  if (!("P0" %in% names(df)))    stop("[manhattan-SS] Missing P0.")
  if (!("P1" %in% names(df)))    df$P1 <- df$P0
  if (!("GENE" %in% names(df)))  df$GENE <- NA_character_
  
  # --- EXCLUDE GENES FROM PLOTTING ONLY ---
  genes_to_exclude <- c("CEBPZOS", "PACSIN1")
  df <- dplyr::filter(df, !(GENE %in% genes_to_exclude))
  # ----------------------------------------
  
  if (!("ID"   %in% names(df)))  df$ID   <- NA_character_
  if (!("CHR"  %in% names(df)))  stop("[manhattan-SS] Missing CHR.")
  df <- dplyr::mutate(df,
                      CHR    = as.integer(CHR),
                      PWAS.P = as.numeric(PWAS.P),
                      P0     = as.numeric(P0),
                      P1     = as.numeric(P1),
                      pos    = ifelse(!is.na(P0) & !is.na(P1), (P0 + P1)/2, P0),
                      logp   = -log10(pmax(PWAS.P, .Machine$double.xmin)))

  # 3) Rebuild cumulative positions and x-axis centers fresh
  chr_sizes <- df %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(cumstart = dplyr::lag(cumsum(chr_len), default = 0))
  df <- df %>%
    dplyr::left_join(chr_sizes, by = "CHR") %>%
    dplyr::mutate(pos_cum = pos + cumstart) %>%
    dplyr::arrange(CHR, pos_cum)
  centers <- df %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(center = (min(pos_cum) + max(pos_cum))/2, .groups = "drop")
  attr(df, "axis_centers") <- centers

  # 4) Threshold and y-span
  ok_sig <- !is.na(df$FDR) & df$FDR <= fdr_thr & !is.na(df$PWAS.P)
  p_star <- if (!any(ok_sig)) NA_real_ else max(df$PWAS.P[ok_sig], na.rm = TRUE)

  max_logp <- max(df$logp, na.rm = TRUE)
  ymax_calc <- max(max_logp + 1,
                   if (is.finite(-log10(p_star))) -log10(p_star) else 0,
                   4,
                   na.rm = TRUE)

  ymax <- if (!is.null(y_max_override)) y_max_override else ymax_calc

  # 5) LABEL ONLY: significant AND in the provided sex-specific GENE list
  lab_df <- df %>%
    dplyr::filter(ok_sig, !is.na(GENE), GENE %in% label_genes) %>%
    dplyr::mutate(label_col = GENE)

  # sex color for labels and labeled dots
  sex_col <- if (identical(sex, "Female")) "#700000" else "#3780b5"

  # 7) Plot (identical to normal, except sex_col used for label + overlay points)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pos_cum, y = logp)) +
    ggplot2::geom_point(ggplot2::aes(color = as.factor(CHR %% 2)),
                        size = 0.8, alpha = 0.9, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("#aaaaaa", "#555555")) +
    ggplot2::geom_point(data = lab_df, color = sex_col, size = 1.1, alpha = 0.95) +
    ggplot2::geom_hline(yintercept = if (is.finite(-log10(p_star))) -log10(p_star) else Inf,
                        alpha = 0.9, linetype = "dashed", color = "black", linewidth = 0.4) + # linetype = "dashed",
    ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = label_col),
      size = 3,
      color = sex_col,
      segment.color = sex_col,
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    ggplot2::scale_x_continuous(
      name   = "Chromosome",
      breaks = attr(df, "axis_centers")$center,
      labels = attr(df, "axis_centers")$CHR,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      name   = "-log10(P)",
      limits = c(0, ceiling(ymax)),
      breaks = seq(0, ceiling(ymax/2)*2, by = 2),
      expand = c(0, 0)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title  = ggplot2::element_text(size = ggplot2::rel(0.7)),
      axis.text   = ggplot2::element_text(size = ggplot2::rel(0.5)),
      legend.title= ggplot2::element_text(size = ggplot2::rel(0.7)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.7))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ceiling(ymax)), clip = "on") #+
  save_plot_dual(p, paste0(out_base, "_manhattan_SS"))
}

# ------------------------------ Driver logic ---------------------------------
tissues <- c("Brain", "CSF")
sexes   <- c("Female", "Male")
discovs <- c("Primary", "Secondary")

# Validate presence of merged files and build map
file_map <- list()
for (t in tissues) for (s in sexes) for (d in discovs) {
  dir_td <- file.path(top_path, t, s, d)
  if (!dir.exists(dir_td)) {
    stop(paste0("Missing folder: ", dir_td, "\n",
                "Please ensure the expected folder structure exists."), call. = FALSE)
  }
  file_map[[t]][[s]][[d]] <- find_merged_file(dir_td)
}
message("All required '*_merged_with_FDR.txt' files found.")

# Read all data
dat <- list()
for (t in tissues) for (s in sexes) for (d in discovs) {
  dat[[t]][[s]][[d]] <- read_merged(file_map[[t]][[s]][[d]])
}

# Manhattan plots (PDF & PNG)
for (t in tissues) for (s in sexes) for (d in discovs) {
  title_sub <- paste(t, s, d, "— Manhattan")
  out_base  <- file.path(plots_dir, paste0(t, "_", s, "_", d))
  df        <- dat[[t]][[s]][[d]]

  # Compute a shared y-limit: max(-log10(P)) + 2
  ymax_shared <- ceiling(max(-log10(pmax(df$PWAS.P, .Machine$double.xmin)), na.rm = TRUE) + 2)

  # normal plot (uses shared y-limit)
  plot_manhattan(df, title_sub, out_base,
               fdr_thr = 0.05, remove_chr = remove_chr_vec,
               y_max_override = ymax_shared)

  # sex-specific list for this folder
  folder   <- file.path(opt$path, t, s, d)
  ss_file  <- find_ss_list(folder)
  ss_genes <- read_ss_genes(ss_file)

  # second set: label only genes in the sex-specific list (same y-limit)
  if (length(ss_genes)) {
    plot_manhattan_sexcolor(df, title_sub, out_base, s, ss_genes,
                          fdr_thr = 0.05, remove_chr = remove_chr_vec,
                          y_max_override = ymax_shared)
  } else {
    message("[manhattan-SS] No sex-specific labels for: ", paste(t, s, d))
  }
}

message("All plots completed. Output in: ", plots_dir)

# Internal Ref PWAS_ManhattanPlot_BrainCSF_V3.R
