#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
option_list <- list(
  make_option(c("--work_dir"), type = "character", help = "Top-level working directory (e.g., /.../PWAS/EU_all)", metavar = "path"),
  make_option(c("--tissue"),   type = "character", help = "Tissue: Brain or CSF", metavar = "{Brain,CSF}"),
  make_option(c("--sex"),      type = "character", help = "Test sex: Male or Female", metavar = "{Male,Female}"),
  make_option(c("--prefix"),   type = "character", help = "Output filename prefix", metavar = "str"),
  make_option(c("--brain_map"),type = "character", help = "(Brain REQUIRED) TSV with columns: ENSG_ID, Start, End, Gene; extra columns allowed", metavar = "file"),
  make_option(c("--csf_map"),  type = "character", help = "(CSF REQUIRED) TSV with columns: Analytes, ENSG_ID, Gene, Start, End", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$work_dir), !is.null(opt$tissue), !is.null(opt$sex), !is.null(opt$prefix))
opt$tissue <- match.arg(opt$tissue, c("Brain", "CSF"))
opt$sex    <- match.arg(opt$sex, c("Male", "Female"))

if (opt$tissue == "Brain") stopifnot(!is.null(opt$brain_map))
if (opt$tissue == "CSF")   stopifnot(!is.null(opt$csf_map))

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
sex_test <- opt$sex
sex_ref  <- ifelse(sex_test == "Male", "Female", "Male")

path_test <- file.path(opt$work_dir, opt$tissue, sex_test)
path_ref  <- file.path(opt$work_dir, opt$tissue, sex_ref)

modes <- c("Primary", "Secondary")

# ------------------------------------------------------------
# IO helpers
# ------------------------------------------------------------
read_map_brain <- function(path) {
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    rename_with(~str_replace_all(., "\\s+", "")) %>%
    # normalize canonical names if variants appear
    rename(
      ENSG_ID = any_of(c("ENSG_ID", "ENSID", "ENSID", "ENSID")),
      Start   = any_of(c("Start")),
      End     = any_of(c("End")),
      Gene    = any_of(c("Gene", "GENE"))
    ) %>%
    select(any_of(c("ENSG_ID", "Start", "End", "Gene")))
}

read_map_csf <- function(path) {
  readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(
      Analytes = readr::col_character(),
      ENSG_ID  = readr::col_character(),
      Gene     = readr::col_character(),
      Start    = readr::col_double(),
      End      = readr::col_double(),
      CHR      = readr::col_double()
    )
  ) %>%
    rename_with(~stringr::str_replace_all(., "\\s+", "")) %>%
    dplyr::rename(
      Analytes = dplyr::any_of(c("Analytes")),
      ENSG_ID  = dplyr::any_of(c("ENSG_ID","ENSID")),
      Gene     = dplyr::any_of(c("Gene","GENE")),
      Start    = dplyr::any_of(c("Start")),
      End      = dplyr::any_of(c("End")),
      CHR      = dplyr::any_of(c("CHR"))
    ) %>%
    dplyr::select(dplyr::any_of(c("Analytes","ENSG_ID","Gene","Start","End","CHR")))
}

# Clean numeric helper
as_num_clean <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("NA", "Inf", ".", ""), NA, x)))

# Read *.dat folder
read_dat_folder <- function(dir_path, tissue) {
  files <- list.files(dir_path, pattern = "\\.dat$", full.names = TRUE)
  if (length(files) == 0) return(NULL)

  bind_rows(lapply(files, function(f) {
    df <- suppressMessages(read_tsv(f, show_col_types = FALSE, progress = FALSE))

    # Required columns
    req <- c("ID", "CHR", "P0", "P1", "TWAS.Z", "TWAS.P")
    miss <- setdiff(req, names(df))
    if (length(miss) > 0) stop(sprintf("%s is missing required columns: %s", basename(f), paste(miss, collapse=", ")))

    # Clean numerics
    df <- df %>% mutate(
      P0     = as_num_clean(P0),
      P1     = as_num_clean(P1),
      `TWAS.Z` = as_num_clean(`TWAS.Z`),
      `TWAS.P` = as_num_clean(`TWAS.P`)
    )

    # Brain-only: derive ENSG_ID & GENE from ID
    if (tissue == "Brain") {
      df <- df %>% mutate(
        ENSG_ID = if_else(str_detect(ID, "\\."), str_replace(ID, "\\..*$", ""), ID),
        GENE    = if_else(str_detect(ID, "\\."), str_replace(ID, ".*\\.", ""), ID)
      )
    }

    df
  }))
}

# Best per key
best_per_key <- function(df, key_col, extra_cols = character()) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>%
    dplyr::group_by(.data[[key_col]]) %>%
    dplyr::slice_min(order_by = `TWAS.P`, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::any_of(c(key_col, "ID", "ENSG_ID", "GENE", "CHR", "P0", "P1", "TWAS.Z", "TWAS.P", "FDR", extra_cols)))
}

# Build REF significance sign map across P/S/O
build_ref_sig_map <- function(dir_ref, key_col) {
  read_mode <- function(mode) read_dat_folder(file.path(dir_ref, mode), tissue = opt$tissue)
  ref_all <- bind_rows(read_mode("Primary"), read_mode("Secondary"), read_mode("Opposite"))
  if (is.null(ref_all) || nrow(ref_all) == 0) return(tibble(!!key_col := character(), signs = I(list())))

  ref_sig <- ref_all %>% filter(`TWAS.P` < 0.05)
  if (nrow(ref_sig) == 0) return(tibble(!!key_col := character(), signs = I(list())))

  ref_sig %>%
    group_by(.data[[key_col]]) %>%
    summarise(signs = list(sign(`TWAS.Z`)), .groups = "drop")
}

# Apply sex-specific rule
apply_sex_specific_rule <- function(test_df_best, ref_map, key_col) {
  if (is.null(test_df_best) || nrow(test_df_best) == 0) return(test_df_best)
  test_df_best %>%
    dplyr::left_join(ref_map, by = key_col) %>%
    dplyr::mutate(signs = ifelse(purrr::map_lgl(signs, is.null), list(integer()), signs)) %>%
    dplyr::mutate(keep = purrr::map2_lgl(signs, sign(`TWAS.Z`), ~{
      sigs <- .x
      if (length(sigs) == 0) return(TRUE)                  # no REF signal
      all(sigs != .y)                                      # all REF signs opposite to TEST sign
    })) %>%
    dplyr::filter(keep) %>%
    dplyr::select(-keep)
}

# --- CSF mapping: match by ID (Analyte) AND CHR to disambiguate duplicates ---
attach_csf <- function(df, map_csf) {
  if (is.null(df) || nrow(df) == 0) return(df)

  out <- df %>%
    dplyr::left_join(map_csf, by = c("ID" = "Analytes", "CHR" = "CHR"))

  # ---- Unify ENSG_ID regardless of suffixing after join ----
  if ("ENSG_ID.x" %in% names(out) && "ENSG_ID.y" %in% names(out)) {
    out <- out %>%
      dplyr::mutate(ENSG_ID = dplyr::coalesce(.data$ENSG_ID.x, .data$ENSG_ID.y)) %>%
      dplyr::select(-ENSG_ID.x, -ENSG_ID.y)
  } else if ("ENSG_ID.x" %in% names(out)) {
    out <- dplyr::rename(out, ENSG_ID = .data$ENSG_ID.x)
  } else if ("ENSG_ID.y" %in% names(out)) {
    out <- dplyr::rename(out, ENSG_ID = .data$ENSG_ID.y)
  }
  # If still missing (rare), create it so downstream select/order works
  if (!"ENSG_ID" %in% names(out)) out <- dplyr::mutate(out, ENSG_ID = NA_character_)

  # ---- Ensure GENE exists before coalescing, then push mapped fields ----
  if (!"GENE" %in% names(out)) out <- dplyr::mutate(out, GENE = NA_character_)

  out %>%
    dplyr::mutate(
      GENE = dplyr::coalesce(.data$GENE, .data$Gene),
      P0   = dplyr::coalesce(.data$Start, .data$P0),
      P1   = dplyr::coalesce(.data$End,   .data$P1)
    ) %>%
    # Drop mapping columns we’ve just consumed
    dplyr::select(-dplyr::any_of(c("Gene", "Start", "End")))
}

attach_brain <- function(df, map_brain) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>%
    left_join(map_brain, by = "ENSG_ID") %>%
    mutate(
      GENE = coalesce(Gene, GENE),
      P0   = coalesce(Start, P0),
      P1   = coalesce(End, P1)
    ) %>%
    select(-Gene, -Start, -End)
}

# Drop PANEL/FILE, standardize leading columns, and alias TWAS.* -> PWAS.* in the header only
write_tsv_header_alias <- function(df, path_out) {
  if (is.null(df)) {
    readr::write_lines("", path_out)
    return(invisible())
  }

  df2 <- df %>% dplyr::select(-dplyr::any_of(c("PANEL", "FILE", "EQTL.GWAS.Z")))

  # leading columns first
  lead <- c("ENSG_ID", "ID", "GENE", "CHR", "P0", "P1")
  cols <- c(intersect(lead, names(df2)), setdiff(names(df2), lead))
  df2  <- df2 %>% dplyr::select(dplyr::all_of(cols))

  # header alias WITHOUT regex
  nms <- names(df2)

  # 1) Direct TWAS.* → PWAS.*
  ix <- startsWith(nms, "TWAS.")
  nms[ix] <- paste0("PWAS.", substring(nms[ix], nchar("TWAS.") + 1))

  # 2) Reference ref_TWAS.* → ref_PWAS.*
  ix <- startsWith(nms, "ref_TWAS.")
  nms[ix] <- paste0("ref_PWAS.", substring(nms[ix], nchar("ref_TWAS.") + 1))

  names(df2) <- nms

  readr::write_tsv(df2, path_out)
}

# ------------------------------------------------------------
# Main mode processing
# ------------------------------------------------------------
process_mode <- function(mode, map_brain = NULL, map_csf = NULL) {
  message(sprintf("Processing %s — %s / %s", opt$tissue, sex_test, mode))

  # 1) TEST merged + FDR
  test_dir <- file.path(path_test, mode)
  test_df  <- read_dat_folder(test_dir, tissue = opt$tissue)
  if (is.null(test_df) || nrow(test_df) == 0) stop(sprintf("No .dat files found in %s", test_dir))
  test_df <- test_df %>% mutate(FDR = p.adjust(`TWAS.P`, method = "BH"))

  # 2) Save merged_with_FDR (after tissue mapping)
  merged_out <- file.path(path_test, mode, sprintf("%s_%s_merged_with_FDR.txt", opt$prefix, mode))
  test_df_mapped <- if (opt$tissue == "CSF") attach_csf(test_df, map_csf) else attach_brain(test_df, map_brain)
  write_tsv_header_alias(test_df_mapped, merged_out)

  # 3) REF same-mode in memory + FDR (not saved)
  ref_dir <- file.path(path_ref, mode)
  ref_df  <- read_dat_folder(ref_dir, tissue = opt$tissue)
  if (!is.null(ref_df) && nrow(ref_df) > 0) ref_df <- ref_df %>% mutate(FDR = p.adjust(`TWAS.P`, method = "BH"))

  # 4) Filter TEST to FDR<0.05 then best-per-key
  key_col <- ifelse(opt$tissue == "Brain", "GENE", "ID")
  test_filtered <- test_df %>% filter(FDR < 0.05)
  test_best     <- best_per_key(test_filtered, key_col = key_col)

  # 5) Join REF stats (same mode) onto TEST filtered distinct per key
  if (!is.null(ref_df) && nrow(ref_df) > 0) {
    # Join REF stats back to TEST by the analysis key (Brain=GENE, CSF=ID)
    ref_best <- best_per_key(ref_df, key_col = key_col) %>%
      dplyr::select(dplyr::all_of(c(key_col, "TWAS.Z", "TWAS.P", "FDR"))) %>%
      dplyr::rename(ref_TWAS.Z = `TWAS.Z`, ref_TWAS.P = `TWAS.P`, ref_TWAS.FDR = FDR)
    test_best <- dplyr::left_join(test_best, ref_best, by = key_col)
  }

  filtered_out <- file.path(path_test, mode, sprintf("%s_%s_filtered_list.txt", opt$prefix, mode))
  test_best_mapped <- if (opt$tissue == "CSF") attach_csf(test_best, map_csf) else attach_brain(test_best, map_brain)
  write_tsv_header_alias(test_best_mapped, filtered_out)

  # 7) Build REF sign map across P/S/O
  ref_sign_map <- build_ref_sig_map(path_ref, key_col = key_col)

  # 8) Apply sex-specific rule
  sex_spec <- apply_sex_specific_rule(test_best, ref_sign_map, key_col = key_col)

  # 9) Save sex-specific list
  sexspec_out <- file.path(path_test, mode, sprintf("%s_%s_sex_specific_gene_list.txt", opt$prefix, mode))
  sex_spec_mapped <- if (opt$tissue == "CSF") attach_csf(sex_spec, map_csf) else attach_brain(sex_spec, map_brain)
  write_tsv_header_alias(sex_spec_mapped, sexspec_out)
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------
message("Inputs:")
message(sprintf("  work_dir: %s", opt$work_dir))
message(sprintf("  tissue:   %s", opt$tissue))
message(sprintf("  sex:      %s", opt$sex))
message(sprintf("  prefix:   %s", opt$prefix))
if (opt$tissue == "Brain") message(sprintf("  brain_map: %s", opt$brain_map))
if (opt$tissue == "CSF")   message(sprintf("  csf_map:   %s", opt$csf_map))

# Validate folder structure
req_sub <- c("Primary", "Secondary", "Opposite")
check_subs <- function(base) all(dir.exists(file.path(base, req_sub)))
if (!check_subs(path_test)) stop(sprintf("Missing required subfolders under TEST: %s", paste(req_sub[!dir.exists(file.path(path_test, req_sub))], collapse=", ")))
if (!check_subs(path_ref))  stop(sprintf("Missing required subfolders under REF: %s",  paste(req_sub[!dir.exists(file.path(path_ref,  req_sub))], collapse=", ")))

# Load maps
map_brain <- if (opt$tissue == "Brain") read_map_brain(opt$brain_map) else NULL
map_csf   <- if (opt$tissue == "CSF")   read_map_csf(opt$csf_map)     else NULL

# Process modes
walk(modes, ~process_mode(.x, map_brain = map_brain, map_csf = map_csf))

message("Done.")
