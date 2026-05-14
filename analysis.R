# =============================================================================
# Head-to-head comparison of seven inflammation-nutrition composite indices
# for postmenopausal osteoporosis (NHANES 2005-2010)
# =============================================================================
#
# Indices evaluated (all computable from NHANES CBC + BIOPRO):
#   CALLY = 10 * Albumin * Lymphocyte / CRP            higher = better
#   SII   = Platelet * Neutrophil / Lymphocyte         higher = worse
#   NLR   = Neutrophil / Lymphocyte                    higher = worse
#   PLR   = Platelet / Lymphocyte                      higher = worse
#   MLR   = Monocyte / Lymphocyte                      higher = worse
#   PNI   = 10 * Albumin + 5 * Lymphocyte              higher = better
#   HALP  = Hgb * Albumin * Lymphocyte / Platelet      higher = better
#
# Reproducing this analysis:
#   1. Install R >= 4.2 and the packages listed in DEPS below
#   2. source("analysis.R") from the repository root
#   3. NHANES files are pulled live via the nhanesA package; outputs land in
#      ./output/figures and ./output/tables. The .RData snapshot at the end
#      lets you re-load all fitted models without re-running NHANES downloads.
#
# Outputs follow the manuscript hierarchy:
#   Main figures (5)  -> Figure_1 ... Figure_5
#   Main tables  (3)  -> Table_1  ... Table_3
#   Supp figures (11) -> Figure_S1 ... Figure_S11  (S1-S4 MLR; S5-S11 CALLY)
#   Supp tables  (7)  -> Table_S1  ... Table_S7    (S1-S5 CALLY; S6-S7 MLR)
#
# Each figure is saved as TIFF (600 dpi LZW) + PDF (vector) + PNG (300 dpi).
# Each table is saved as UTF-8 CSV and additionally bundled into a single
# multi-sheet xlsx if openxlsx is available.
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
# DEPS: install once with
#   install.packages(c("nhanesA","survey","tidyverse","rms","mediation",
#                      "pROC","tableone","splines","ggplot2","openxlsx"))

suppressPackageStartupMessages({
  library(nhanesA)
  library(survey)
  library(tidyverse)
  library(rms)
  library(mediation)
  library(pROC)
  library(tableone)
  library(splines)
  library(ggplot2)
})

# openxlsx is optional - degrade to CSV-only if missing
HAS_OPENXLSX <- tryCatch({
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    install.packages("openxlsx", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(openxlsx))
  TRUE
}, error = function(e) {
  message("openxlsx unavailable (", conditionMessage(e),
          "); CSV outputs only.")
  FALSE
})

# rms masks dplyr::select, force our preferred version
select <- dplyr::select

options(survey.lonely.psu = "adjust")
set.seed(20260511)


# ---- 0a. Output paths and helpers ------------------------------------------
# Outputs go to ./output relative to the working directory rather than the
# Desktop, so the repo stays self-contained.
OUT_ROOT <- file.path(getwd(), "output")
FIG_DIR  <- file.path(OUT_ROOT, "figures")
TBL_DIR  <- file.path(OUT_ROOT, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

message("Figures will be written to: ", FIG_DIR)
message("Tables  will be written to: ", TBL_DIR)

# Standard print column widths for journal submission (mm)
WIDTH_SINGLE  <-  89
WIDTH_ONEHALF <- 140
WIDTH_DOUBLE  <- 183
MM_TO_INCH    <- 1 / 25.4

# Plain ggplot theme matching most journal style guides
sci_theme <- function(base_size = 8) {
  theme_classic(base_family = "sans", base_size = base_size) +
    theme(
      axis.text        = element_text(color = "black", size = base_size),
      axis.title       = element_text(color = "black", size = base_size + 1),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                      hjust = 0.5),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(3, "mm"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin      = margin(3, 3, 3, 3, "mm")
    )
}

# ggplot saver: TIFF + PDF + PNG in one call
save_fig_ggplot <- function(plot, filename,
                            width_mm = WIDTH_SINGLE, height_mm = 75,
                            dpi = 600) {
  w <- width_mm  * MM_TO_INCH
  h <- height_mm * MM_TO_INCH
  base <- file.path(FIG_DIR, filename)
  ggsave(paste0(base, ".tiff"), plot = plot,
         width = w, height = h, units = "in", dpi = dpi,
         device = grDevices::tiff, compression = "lzw", bg = "white")
  ggsave(paste0(base, ".pdf"), plot = plot,
         width = w, height = h, units = "in",
         device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".png"), plot = plot,
         width = w, height = h, units = "in", dpi = 300, bg = "white")
  message("  saved: ", filename)
}

# Base R graphics saver (for pROC objects etc.)
save_fig_base <- function(plot_expr, filename,
                          width_mm = WIDTH_SINGLE, height_mm = 75,
                          dpi = 600) {
  w <- width_mm  * MM_TO_INCH
  h <- height_mm * MM_TO_INCH
  base <- file.path(FIG_DIR, filename)
  expr <- substitute(plot_expr)
  env  <- parent.frame()
  tiff(paste0(base, ".tiff"), width = w, height = h, units = "in",
       res = dpi, compression = "lzw", family = "sans", bg = "white")
  eval(expr, envir = env); dev.off()
  cairo_pdf(paste0(base, ".pdf"), width = w, height = h,
            family = "sans", bg = "white")
  eval(expr, envir = env); dev.off()
  png(paste0(base, ".png"), width = w, height = h, units = "in",
      res = 300, bg = "white", family = "sans")
  eval(expr, envir = env); dev.off()
  message("  saved: ", filename)
}

# Format helpers
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
fmt_or_ci_p <- function(or, l, u, p) {
  ifelse(is.na(or), "",
         sprintf("%.2f (%.2f-%.2f), P=%s", or, l, u, fmt_p(p)))
}

# Table saver: writes CSV and stashes the data.frame for the final xlsx bundle
.TABLE_BUFFER <- list()
save_table <- function(df, name, sheet_name = NULL) {
  csv_path <- file.path(TBL_DIR, paste0(name, ".csv"))
  readr::write_excel_csv(as.data.frame(df), csv_path, na = "")
  key <- if (is.null(sheet_name)) name else sheet_name
  .TABLE_BUFFER[[key]] <<- df
  message("  saved table: ", name, ".csv  (",
          nrow(df), " rows x ", ncol(df), " cols)")
  invisible(df)
}

bundle_tables_xlsx <- function(filename = "all_tables.xlsx") {
  if (!HAS_OPENXLSX) {
    message("  openxlsx not available; CSVs are in ", TBL_DIR)
    return(invisible(NULL))
  }
  if (length(.TABLE_BUFFER) == 0) {
    message("  no tables to bundle")
    return(invisible(NULL))
  }
  wb <- createWorkbook()
  hdr  <- createStyle(textDecoration = "bold", halign = "center",
                      border = "Bottom", borderStyle = "medium")
  body <- createStyle(halign = "left", valign = "center")
  for (nm in names(.TABLE_BUFFER)) {
    sheet <- substr(nm, 1, 31)
    addWorksheet(wb, sheet)
    df <- .TABLE_BUFFER[[nm]]
    writeData(wb, sheet, df, headerStyle = hdr)
    addStyle(wb, sheet, body,
             rows = 2:(nrow(df) + 1), cols = 1:ncol(df),
             gridExpand = TRUE, stack = TRUE)
    setColWidths(wb, sheet, cols = 1:ncol(df), widths = "auto")
    freezePane(wb, sheet, firstRow = TRUE)
  }
  out <- file.path(TBL_DIR, filename)
  saveWorkbook(wb, out, overwrite = TRUE)
  message("  bundled workbook: ", filename,
          " (", length(.TABLE_BUFFER), " sheets)")
  invisible(out)
}


# ---- 1. Pull NHANES 2005-2010 ----------------------------------------------
# Three two-year cycles (D, E, F). Everything is coerced to character on
# arrival to dodge type clashes when binding rows across cycles.

safe_nhanes <- function(table_name, max_attempts = 4,
                        base_wait = 3, verbose = TRUE) {
  # NCHS' wwwn.cdc.gov server occasionally returns transient 5xx errors or
  # empty responses, especially on the older cycles (D and E). Treat any
  # empty result as a failure and retry with exponential backoff before
  # giving up.
  for (attempt in seq_len(max_attempts)) {
    df <- tryCatch(
      nhanes(table_name),
      error = function(e) {
        if (verbose) message(sprintf("    attempt %d/%d failed: %s",
                                     attempt, max_attempts, e$message))
        NULL
      }
    )
    if (!is.null(df) && nrow(df) > 0) {
      return(df %>% mutate(across(everything(), as.character)))
    }
    if (attempt < max_attempts) {
      wait <- base_wait * 2^(attempt - 1)
      if (verbose) message(sprintf("    empty/failed; retrying in %ds", wait))
      Sys.sleep(wait)
    }
  }
  NULL
}

# Mandatory NHANES tables and their fallback names. The default name
# follows the convention TABLE_<cycle suffix>. A few tables have legacy or
# revised filenames; list those alternatives here so the retry logic can
# walk through them before declaring a table missing.
TABLE_FALLBACKS <- list(
  DEMO   = c("DEMO"),
  BMX    = c("BMX"),
  DXXFEM = c("DXXFEM"),
  CBC    = c("CBC"),
  BIOPRO = c("BIOPRO"),
  CRP    = c("CRP"),
  SMQ    = c("SMQ"),
  ALQ    = c("ALQ"),
  RHQ    = c("RHQ"),
  DIQ    = c("DIQ"),
  BPQ    = c("BPQ"),
  MCQ    = c("MCQ"),
  # NCHS released harmonized RIA data as VID2_C/D in 2010, then replaced
  # them with LC-MS/MS-standardized VID_C/D files in October 2015. Both
  # names exist on the server; try the standardized one first.
  VID    = c("VID", "VID2")
)

# Tables we will refuse to proceed without. If any of these fails for all
# three cycles, the script stops and tells the user to retry later.
REQUIRED_TABLES <- c("DEMO", "BMX", "DXXFEM", "CBC", "BIOPRO", "CRP")

fetch_table <- function(stem, suffix) {
  # Try every fallback name for this stem until one succeeds.
  candidates <- paste0(TABLE_FALLBACKS[[stem]], "_", suffix)
  for (tn in candidates) {
    message("    trying ", tn)
    d <- safe_nhanes(tn)
    if (!is.null(d)) return(d)
  }
  NULL
}

download_cycle <- function(suffix) {
  message("Downloading cycle ", suffix)
  out <- list()
  missing <- character()
  for (stem in names(TABLE_FALLBACKS)) {
    message("  ", stem, "_", suffix)
    d <- fetch_table(stem, suffix)
    if (!is.null(d)) {
      out[[stem]] <- d
    } else {
      missing <- c(missing, stem)
      message(sprintf("    NOT RETRIEVED after retries: %s_%s", stem, suffix))
    }
  }
  # Hard-stop on missing required tables - any one of these failing would
  # silently bias the analytic sample.
  req_missing <- intersect(missing, REQUIRED_TABLES)
  if (length(req_missing) > 0) {
    stop(sprintf(
      "Cycle %s is missing required tables (%s). NCHS server is likely having a transient issue; rerun in a few minutes.",
      suffix, paste(req_missing, collapse = ", ")
    ))
  }
  attr(out, "missing") <- missing
  out
}

merge_cycle <- function(lst) {
  Reduce(function(x, y) full_join(x, y, by = "SEQN"), lst)
}

message("\n=== NHANES download ===")
cyc_D <- download_cycle("D")
cyc_E <- download_cycle("E")
cyc_F <- download_cycle("F")

m_D <- merge_cycle(cyc_D)
m_E <- merge_cycle(cyc_E)
m_F <- merge_cycle(cyc_F)

# Report optional-table coverage so the analyst can see, before any models
# are fit, exactly which (non-required) tables fell through. The downstream
# code already tolerates these via is.na() handling, but the user now sees
# the picture rather than scrolling back through download logs.
message("\n=== Table coverage summary ===")
for (cyc in c("D", "E", "F")) {
  miss <- attr(get(paste0("cyc_", cyc)), "missing")
  if (length(miss) == 0) {
    message(sprintf("  Cycle %s: all tables retrieved", cyc))
  } else {
    message(sprintf("  Cycle %s: missing optional tables: %s",
                    cyc, paste(miss, collapse = ", ")))
  }
}

all_data <- bind_rows(m_D, m_E, m_F)
message("\nCombined sample: n = ", nrow(all_data))


# ---- 2. Variable extraction and typing -------------------------------------
yn_to_01 <- function(x) {
  case_when(
    x %in% c("1", "Yes") ~ 1,
    x %in% c("2", "No")  ~ 0,
    TRUE ~ NA_real_
  )
}

get_col <- function(df, var) {
  if (var %in% names(df)) df[[var]] else rep(NA_character_, nrow(df))
}

# coalesce numeric columns by name (for variables that changed code across
# cycles, e.g. LBDLYMNO vs LBXLYMNO for absolute lymphocyte count)
get_first <- function(df, ...) {
  vars <- c(...)
  cols <- lapply(vars, function(v) as.numeric(get_col(df, v)))
  Reduce(coalesce, cols)
}

dat <- tibble(
  seqn      = as.numeric(get_col(all_data, "SEQN")),
  age       = as.numeric(get_col(all_data, "RIDAGEYR")),
  gender_c  = get_col(all_data, "RIAGENDR"),
  race_c    = get_col(all_data, "RIDRETH1"),
  edu_c     = get_col(all_data, "DMDEDUC2"),
  pir       = as.numeric(get_col(all_data, "INDFMPIR")),
  psu       = as.numeric(get_col(all_data, "SDMVPSU")),
  strata    = as.numeric(get_col(all_data, "SDMVSTRA")),
  wt_mec    = as.numeric(get_col(all_data, "WTMEC2YR")),
  bmi       = as.numeric(get_col(all_data, "BMXBMI")),

  bmd_fn    = as.numeric(get_col(all_data, "DXXNKBMD")),
  bmd_tot   = as.numeric(get_col(all_data, "DXXOFBMD")),

  crp       = as.numeric(get_col(all_data, "LBXCRP")),
  alb       = as.numeric(get_col(all_data, "LBXSAL")),
  lymph_abs = get_first(all_data, "LBDLYMNO", "LBXLYMNO"),
  neut_abs  = get_first(all_data, "LBDNENO",  "LBXNENO"),
  plt       = as.numeric(get_col(all_data, "LBXPLTSI")),
  mono_abs  = get_first(all_data, "LBDMONO",  "LBXMONO"),
  hgb       = as.numeric(get_col(all_data, "LBXHGB")),

  smq020_c  = get_col(all_data, "SMQ020"),
  smq040_c  = get_col(all_data, "SMQ040"),
  alq101_c  = get_col(all_data, "ALQ101"),
  alq130    = as.numeric(get_col(all_data, "ALQ130")),

  rhq060    = as.numeric(get_col(all_data, "RHQ060")),
  rhq031_c  = get_col(all_data, "RHQ031"),
  rhq540_c  = get_col(all_data, "RHQ540"),

  diq010_c  = get_col(all_data, "DIQ010"),
  bpq020_c  = get_col(all_data, "BPQ020"),
  mcq220_c  = get_col(all_data, "MCQ220"),

  vit_d     = get_first(all_data, "LBXVIDMS", "LBDVIDMS",
                                  "LBXVIDLC", "LBDVIDLC")
)


# ---- 3. Categorical recoding -----------------------------------------------
dat <- dat %>%
  mutate(
    gender = as.numeric(case_when(
      gender_c %in% c("1","Male")   ~ "1",
      gender_c %in% c("2","Female") ~ "2",
      TRUE ~ NA_character_
    )),

    race = factor(case_when(
      race_c %in% c("1","Mexican American")              ~ "Mexican American",
      race_c %in% c("2","Other Hispanic")                ~ "Other Hispanic",
      race_c %in% c("3","Non-Hispanic White","NH White") ~ "NH White",
      race_c %in% c("4","Non-Hispanic Black","NH Black") ~ "NH Black",
      race_c %in% c("5","Other Race - Including Multi-Racial",
                    "Other","Other Race")                ~ "Other",
      TRUE ~ NA_character_
    ), levels = c("NH White","NH Black","Mexican American",
                  "Other Hispanic","Other")),

    # nhanesA returns label strings in some cycles, codes in others; match both
    edu_cat = factor(case_when(
      edu_c %in% c("1","2") |
        grepl("Less than 9th|9-11th", edu_c, ignore.case = TRUE) ~ "<High school",
      edu_c == "3" |
        grepl("High school|GED",      edu_c, ignore.case = TRUE) ~ "High school",
      edu_c %in% c("4","5") |
        grepl("Some college|AA|College grad", edu_c,
              ignore.case = TRUE)                                ~ ">High school",
      TRUE ~ NA_character_
    ), levels = c("<High school","High school",">High school")),

    diabetes     = yn_to_01(diq010_c),
    hypertension = yn_to_01(bpq020_c),
    cancer       = yn_to_01(mcq220_c),
    hrt          = yn_to_01(rhq540_c),
    rhq031_yn    = yn_to_01(rhq031_c),

    smoke_status = factor(case_when(
      smq020_c %in% c("2","No") ~ "Never",
      smq020_c %in% c("1","Yes") & smq040_c %in% c("3","Not at all") ~ "Former",
      smq020_c %in% c("1","Yes") & smq040_c %in% c("1","Every day",
                                                    "2","Some days") ~ "Current",
      TRUE ~ NA_character_
    ), levels = c("Never","Former","Current")),

    alcohol = factor(case_when(
      alq101_c %in% c("2","No") ~ "None",
      alq130 <  1                ~ "Light",
      alq130 >= 1 & alq130 < 3   ~ "Moderate",
      alq130 >= 3                ~ "Heavy",
      TRUE ~ NA_character_
    ), levels = c("None","Light","Moderate","Heavy"))
  )


# ---- 4. Compute the seven inflammation-nutrition indices -------------------
# Units follow NHANES defaults:
#   Alb (LBXSAL)     g/dL
#   Lymph (LBDLYMNO) 10^3 cells/uL
#   Neut (LBDNENO)   10^3 cells/uL
#   Mono (LBDMONO)   10^3 cells/uL
#   Plt (LBXPLTSI)   10^3/uL
#   Hgb (LBXHGB)     g/dL
#   CRP (LBXCRP)     mg/dL
#
# Pre-specified direction for ROC and effect-size signs:
#   higher = protective : CALLY, PNI, HALP
#   higher = risk       : SII, NLR, PLR, MLR
#
# All seven indices are strictly positive in valid cases, so we log-transform
# directly without a +1 offset. The +1 offset that earlier drafts of this
# script used artificially compressed the dynamic range of MLR (typical range
# 0.05-1.0), producing implausibly wide CIs on the per-unit log scale.
dat <- dat %>%
  mutate(
    cally = (lymph_abs * alb) / crp * 10,
    sii   = plt * neut_abs / lymph_abs,
    nlr   = neut_abs / lymph_abs,
    plr   = plt / lymph_abs,
    mlr   = mono_abs / lymph_abs,
    pni   = 10 * alb + 5 * lymph_abs,
    halp  = (hgb * alb * lymph_abs) / plt,

    log_cally = log(cally),
    log_sii   = log(sii),
    log_nlr   = log(nlr),
    log_plr   = log(plr),
    log_mlr   = log(mlr),
    log_pni   = log(pni),
    log_halp  = log(halp)
  )

INDEX_INFO <- data.frame(
  short      = c("cally", "sii", "nlr", "plr", "mlr", "pni", "halp"),
  label      = c("CALLY", "SII", "NLR", "PLR", "MLR", "PNI", "HALP"),
  raw_var    = c("cally", "sii", "nlr", "plr", "mlr", "pni", "halp"),
  log_var    = c("log_cally","log_sii","log_nlr","log_plr","log_mlr",
                 "log_pni","log_halp"),
  direction  = c(">", "<", "<", "<", "<", ">", ">"),  # for pROC
  protective = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)


# ---- 5. Outcome definition (WHO criterion at femoral neck) -----------------
# T-score reference from NHANES III healthy young non-Hispanic white women,
# femoral neck mean 0.849 g/cm^2, SD 0.111 g/cm^2 (Looker et al. 1998).
fn_mean <- 0.849; fn_sd <- 0.111

dat <- dat %>%
  mutate(
    t_score_fn   = (bmd_fn - fn_mean) / fn_sd,
    osteoporosis = case_when(
      t_score_fn <= -2.5 ~ 1,
      t_score_fn >  -2.5 ~ 0,
      TRUE ~ NA_real_
    ),
    osteopenia = case_when(
      t_score_fn <= -1.0 & t_score_fn > -2.5 ~ 1,
      t_score_fn >  -1.0 ~ 0,
      TRUE ~ NA_real_
    )
  )


# ---- 6. Postmenopausal flag -------------------------------------------------
# Self-reported menopause when available; otherwise fall back to age >= 55 with
# no current menstruation (RHQ060 = age at last period).
dat <- dat %>%
  mutate(
    postmenopausal = case_when(
      gender == 2 & rhq031_yn == 0 ~ 1,
      gender == 2 & !is.na(rhq060) & age >= rhq060 ~ 1,
      gender == 2 & age >= 55 ~ 1,
      TRUE ~ 0
    ),
    age_cat = ifelse(age >= 65, ">=65", "50-64"),
    bmi_cat = factor(case_when(
      bmi < 25 ~ "Normal",
      bmi < 30 ~ "Overweight",
      TRUE     ~ "Obese"
    ), levels = c("Normal","Overweight","Obese"))
  )


# ---- 7. Analytic sample ----------------------------------------------------
# Complete-case for all seven indices to keep the head-to-head comparison
# on the same patients.
study_dat <- dat %>%
  filter(
    gender == 2,
    age >= 50,
    postmenopausal == 1,
    !is.na(bmd_fn),
    !is.na(cally), is.finite(cally),
    !is.na(sii),   is.finite(sii),
    !is.na(nlr),   is.finite(nlr),
    !is.na(plr),   is.finite(plr),
    !is.na(mlr),   is.finite(mlr),
    !is.na(pni),   is.finite(pni),
    !is.na(halp),  is.finite(halp),
    !is.na(bmi),
    cancer == 0 | is.na(cancer),
    hrt    == 0 | is.na(hrt)
  )

message("\n=== Analytic sample ===")
message("Final n          : ", nrow(study_dat))
message("Osteoporosis (n) : ",
        sum(study_dat$osteoporosis == 1, na.rm = TRUE))

study_dat <- study_dat %>%
  mutate(
    cally_q_num = ntile(cally, 4),
    cally_q     = factor(cally_q_num, levels = 1:4,
                         labels = c("Q1","Q2","Q3","Q4")),
    # Combining three 2-year cycles: per NCHS guidance, divide MEC weight by 3
    wt_new      = wt_mec / 3
  )


# ---- 7b. Study flow chart (Figure 1) ---------------------------------------
# Numbers are recovered by re-applying the cumulative filter to dat. Done this
# way to avoid coupling the chart numbers to the filter pipeline above.

n0 <- nrow(dat)
n1 <- dat %>% filter(gender == 2)                                  %>% nrow()
n2 <- dat %>% filter(gender == 2, age >= 50)                       %>% nrow()
n3 <- dat %>% filter(gender == 2, age >= 50, postmenopausal == 1)  %>% nrow()
n4 <- dat %>% filter(gender == 2, age >= 50, postmenopausal == 1,
                     !is.na(bmd_fn))                               %>% nrow()
n5 <- dat %>% filter(gender == 2, age >= 50, postmenopausal == 1,
                     !is.na(bmd_fn),
                     !is.na(cally), is.finite(cally),
                     !is.na(sii),   is.finite(sii),
                     !is.na(nlr),   is.finite(nlr),
                     !is.na(plr),   is.finite(plr),
                     !is.na(mlr),   is.finite(mlr),
                     !is.na(pni),   is.finite(pni),
                     !is.na(halp),  is.finite(halp))               %>% nrow()
n6 <- dat %>% filter(gender == 2, age >= 50, postmenopausal == 1,
                     !is.na(bmd_fn),
                     !is.na(cally), is.finite(cally),
                     !is.na(sii),   is.finite(sii),
                     !is.na(nlr),   is.finite(nlr),
                     !is.na(plr),   is.finite(plr),
                     !is.na(mlr),   is.finite(mlr),
                     !is.na(pni),   is.finite(pni),
                     !is.na(halp),  is.finite(halp),
                     !is.na(bmi))                                  %>% nrow()
n_final <- nrow(study_dat)

boxes <- data.frame(
  y = 8:1,
  label = c(
    sprintf("NHANES 2005-2010 participants\n(n = %s)",
            formatC(n0, big.mark = ",")),
    sprintf("Female\n(n = %s)",                       formatC(n1, big.mark = ",")),
    sprintf("Aged >= 50 years\n(n = %s)",             formatC(n2, big.mark = ",")),
    sprintf("Postmenopausal\n(n = %s)",               formatC(n3, big.mark = ",")),
    sprintf("Available femoral neck BMD\n(n = %s)",   formatC(n4, big.mark = ",")),
    sprintf("Complete data for all 7 indices\n(n = %s)",
            formatC(n5, big.mark = ",")),
    sprintf("Available BMI\n(n = %s)",                formatC(n6, big.mark = ",")),
    sprintf("Final analytic sample\n(n = %s)\nOsteoporosis: %d",
            formatC(n_final, big.mark = ","),
            sum(study_dat$osteoporosis == 1, na.rm = TRUE))
  ),
  stringsAsFactors = FALSE
)

excl <- data.frame(
  y = c(7.5, 6.5, 5.5, 4.5, 3.5, 2.5, 1.5),
  label = c(
    sprintf("Excluded males\n(n = %s)",       formatC(n0 - n1,      big.mark = ",")),
    sprintf("Age < 50 y\n(n = %s)",           formatC(n1 - n2,      big.mark = ",")),
    sprintf("Premenopausal\n(n = %s)",        formatC(n2 - n3,      big.mark = ",")),
    sprintf("Missing BMD\n(n = %s)",          formatC(n3 - n4,      big.mark = ",")),
    sprintf("Missing any of\n7 indices\n(n = %s)",
            formatC(n4 - n5, big.mark = ",")),
    sprintf("Missing BMI\n(n = %s)",          formatC(n5 - n6,      big.mark = ",")),
    sprintf("Cancer/HRT history\n(n = %s)",   formatC(n6 - n_final, big.mark = ","))
  ),
  stringsAsFactors = FALSE
)

p_flow <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = 0.05, xmax = 0.55,
                ymin = y - 0.32, ymax = y + 0.32),
            fill = "white", color = "black", linewidth = 0.4) +
  geom_text(data = boxes,
            aes(x = 0.30, y = y, label = label),
            size = 2.2, color = "black", lineheight = 0.95) +
  geom_segment(data = boxes[1:7, ],
               aes(x = 0.30, xend = 0.30,
                   y = y - 0.32, yend = y - 0.68),
               arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
               linewidth = 0.4, color = "black") +
  geom_segment(data = excl,
               aes(x = 0.30, xend = 0.65, y = y, yend = y),
               linewidth = 0.4, color = "black") +
  geom_rect(data = excl,
            aes(xmin = 0.65, xmax = 0.98,
                ymin = y - 0.22, ymax = y + 0.22),
            fill = "grey95", color = "black", linewidth = 0.3) +
  geom_text(data = excl,
            aes(x = 0.815, y = y, label = label),
            size = 1.9, color = "black", lineheight = 0.95) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 8.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2, "mm"))

save_fig_ggplot(p_flow, "Figure_1_Study_Flowchart",
                width_mm = WIDTH_ONEHALF, height_mm = 140)


# ---- 8. Survey design ------------------------------------------------------
design    <- svydesign(id = ~psu, strata = ~strata, weights = ~wt_new,
                       nest = TRUE, data = study_dat)
design_op <- subset(design, !is.na(osteoporosis))


# ---- 9. Table 1: baseline by osteoporosis status ---------------------------
tab1_vars <- c("age","race","edu_cat","pir","bmi","smoke_status","alcohol",
               "diabetes","hypertension","vit_d","bmd_fn","t_score_fn",
               "cally","sii","nlr","plr","mlr","pni","halp")
fac_vars  <- c("race","edu_cat","smoke_status","alcohol",
               "diabetes","hypertension")

tab1 <- svyCreateTableOne(vars = tab1_vars, strata = "osteoporosis",
                          data = design_op, factorVars = fac_vars, test = TRUE)
print(tab1, showAllLevels = TRUE, smd = TRUE, contDigits = 2)

tab1_mat <- print(tab1, showAllLevels = TRUE, smd = TRUE, contDigits = 2,
                  printToggle = FALSE, noSpaces = TRUE)
tab1_df  <- as.data.frame(tab1_mat, stringsAsFactors = FALSE)
tab1_df  <- cbind(Variable = rownames(tab1_df), tab1_df)
rownames(tab1_df) <- NULL
save_table(tab1_df, "Table_1_Baseline_by_Osteoporosis",
           sheet_name = "Table 1 Baseline")


# ---- 10. Three-model regression on CALLY (the original primary index) ------
exp_ci <- function(model, term = NULL) {
  est <- coef(model); se <- sqrt(diag(vcov(model)))
  out <- data.frame(
    Variable = names(est),
    OR  = round(exp(est), 3),
    L95 = round(exp(est - 1.96 * se), 3),
    U95 = round(exp(est + 1.96 * se), 3),
    p   = signif(2 * (1 - pnorm(abs(est / se))), 3)
  )
  if (!is.null(term)) out <- out[out$Variable == term, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Continuous BMD model (informational, not in main results)
m1_lin <- svyglm(bmd_fn ~ log_cally, design = design_op)
m2_lin <- svyglm(bmd_fn ~ log_cally + age + race + bmi, design = design_op)
m3_lin <- svyglm(bmd_fn ~ log_cally + age + race + bmi + edu_cat + pir +
                          smoke_status + alcohol + diabetes + hypertension + vit_d,
                 design = design_op)

# Logistic models on osteoporosis status
m1_log <- svyglm(osteoporosis ~ log_cally,
                 design = design_op, family = quasibinomial())
m2_log <- svyglm(osteoporosis ~ log_cally + age + race + bmi,
                 design = design_op, family = quasibinomial())
m3_log <- svyglm(osteoporosis ~ log_cally + age + race + bmi + edu_cat + pir +
                                smoke_status + alcohol + diabetes + hypertension + vit_d,
                 design = design_op, family = quasibinomial())

tbl2a <- bind_rows(
  cbind(Model = "Model 1 (Crude)",            exp_ci(m1_log, "log_cally")),
  cbind(Model = "Model 2 (Age + Race + BMI)", exp_ci(m2_log, "log_cally")),
  cbind(Model = "Model 3 (Fully adjusted)",   exp_ci(m3_log, "log_cally"))
) %>%
  mutate(`OR (95% CI), P` = fmt_or_ci_p(OR, L95, U95, p)) %>%
  select(Model, Exposure = Variable, OR, `Lower 95% CI` = L95,
         `Upper 95% CI` = U95, `P value` = p, `OR (95% CI), P`)
save_table(tbl2a, "Table_S1_CALLY_ThreeModels",
           sheet_name = "Table S1 CALLY 3-models")

# Quartile model + trend test
m3_q     <- svyglm(osteoporosis ~ cally_q + age + race + bmi + edu_cat + pir +
                                  smoke_status + alcohol + diabetes + hypertension + vit_d,
                   design = design_op, family = quasibinomial())
m3_trend <- svyglm(osteoporosis ~ cally_q_num + age + race + bmi + edu_cat + pir +
                                  smoke_status + alcohol + diabetes + hypertension + vit_d,
                   design = design_op, family = quasibinomial())

q_all   <- exp_ci(m3_q)
q_rows  <- q_all[grepl("^cally_q", q_all$Variable), ]
p_trend <- summary(m3_trend)$coefficients["cally_q_num", "Pr(>|t|)"]

tbl2b <- bind_rows(
  data.frame(
    Quartile = "Q1 (reference)",
    `Range`  = sprintf("<= %.2f",
                       quantile(study_dat$cally, 0.25, na.rm = TRUE)),
    OR = 1.00, `Lower 95% CI` = NA_real_, `Upper 95% CI` = NA_real_,
    `P value` = NA_real_, check.names = FALSE
  ),
  data.frame(
    Quartile = c("Q2", "Q3", "Q4"),
    `Range`  = c(
      sprintf("%.2f - %.2f",
              quantile(study_dat$cally, 0.25, na.rm = TRUE),
              quantile(study_dat$cally, 0.50, na.rm = TRUE)),
      sprintf("%.2f - %.2f",
              quantile(study_dat$cally, 0.50, na.rm = TRUE),
              quantile(study_dat$cally, 0.75, na.rm = TRUE)),
      sprintf("> %.2f",
              quantile(study_dat$cally, 0.75, na.rm = TRUE))
    ),
    OR = q_rows$OR,
    `Lower 95% CI` = q_rows$L95,
    `Upper 95% CI` = q_rows$U95,
    `P value` = q_rows$p,
    check.names = FALSE
  )
) %>%
  mutate(`OR (95% CI), P` = ifelse(
    Quartile == "Q1 (reference)",
    "1.00 (reference)",
    fmt_or_ci_p(OR, `Lower 95% CI`, `Upper 95% CI`, `P value`)
  ))

tbl2b <- bind_rows(
  tbl2b,
  data.frame(
    Quartile = "P for trend", Range = "",
    OR = NA, `Lower 95% CI` = NA, `Upper 95% CI` = NA,
    `P value` = p_trend,
    `OR (95% CI), P` = sprintf("P-trend = %s", fmt_p(p_trend)),
    check.names = FALSE
  )
)
save_table(tbl2b, "Table_S2_CALLY_Quartile",
           sheet_name = "Table S2 CALLY Quartile")


# ---- 10b. CALLY three-model forest plot (Figure S5) ------------------------
forest_models <- bind_rows(
  cbind(Model = "Model 1 (Crude)",          exp_ci(m1_log, "log_cally")),
  cbind(Model = "Model 2 (Age/Race/BMI)",   exp_ci(m2_log, "log_cally")),
  cbind(Model = "Model 3 (Fully adjusted)", exp_ci(m3_log, "log_cally"))
) %>% mutate(Model = factor(Model, levels = rev(unique(Model))))

p_forest_models <- ggplot(forest_models, aes(x = OR, y = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.15, linewidth = 0.5,
                 color = "black") +
  geom_point(size = 2.2, color = "#2C7FB8") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f),  P=%.3f",
                                OR, L95, U95, p)),
            hjust = -0.12, size = 2.4) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.55))) +
  labs(x = "OR for osteoporosis per unit log(CALLY)", y = NULL) +
  sci_theme()

save_fig_ggplot(p_forest_models, "Figure_S5_CALLY_ThreeModels",
                width_mm = WIDTH_ONEHALF, height_mm = 55)


# ---- 10c. CALLY quartile forest (Figure S6) --------------------------------
q_forest_df <- data.frame(
  Quartile = factor(c("Q4","Q3","Q2","Q1 (reference)"),
                    levels = c("Q4","Q3","Q2","Q1 (reference)")),
  OR = c(q_rows$OR[3],  q_rows$OR[2],  q_rows$OR[1],  1.00),
  L  = c(q_rows$L95[3], q_rows$L95[2], q_rows$L95[1], NA),
  U  = c(q_rows$U95[3], q_rows$U95[2], q_rows$U95[1], NA),
  p  = c(q_rows$p[3],   q_rows$p[2],   q_rows$p[1],   NA),
  stringsAsFactors = FALSE
)

p_forest_q <- ggplot(q_forest_df, aes(x = OR, y = Quartile)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.15, linewidth = 0.5,
                 color = "black", na.rm = TRUE) +
  geom_point(size = 2.2, color = "#D7301F") +
  geom_text(aes(label = ifelse(is.na(L), "1.00 (reference)",
                               sprintf("%.2f (%.2f-%.2f), P=%s",
                                       OR, L, U, fmt_p(p)))),
            hjust = -0.10, size = 2.4) +
  annotate("text", x = max(q_forest_df$U, na.rm = TRUE) * 1.05, y = 0.4,
           label = sprintf("P for trend = %s", fmt_p(p_trend)),
           size = 2.4, fontface = "italic", hjust = 1) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.65))) +
  labs(x = "OR for osteoporosis (95% CI)", y = NULL) +
  sci_theme()

save_fig_ggplot(p_forest_q, "Figure_S6_CALLY_Quartile",
                width_mm = WIDTH_ONEHALF, height_mm = 65)


# ---- 11. CALLY restricted cubic spline on osteoporosis (Figure S7) ---------
rcs_fit <- svyglm(osteoporosis ~ ns(log_cally, df = 4) +
                                 age + race + bmi + edu_cat + pir +
                                 smoke_status + alcohol + diabetes +
                                 hypertension + vit_d,
                  design = design_op, family = quasibinomial())

new_x <- seq(quantile(study_dat$log_cally, 0.025, na.rm = TRUE),
             quantile(study_dat$log_cally, 0.975, na.rm = TRUE),
             length.out = 200)
ref_idx <- which.min(abs(new_x - median(study_dat$log_cally, na.rm = TRUE)))

new_dat <- data.frame(
  log_cally    = new_x,
  age          = median(study_dat$age, na.rm = TRUE),
  race         = factor("NH White",     levels = levels(study_dat$race)),
  bmi          = median(study_dat$bmi, na.rm = TRUE),
  edu_cat      = factor(">High school", levels = levels(study_dat$edu_cat)),
  pir          = median(study_dat$pir, na.rm = TRUE),
  smoke_status = factor("Never",        levels = levels(study_dat$smoke_status)),
  alcohol      = factor("None",         levels = levels(study_dat$alcohol)),
  diabetes     = 0,
  hypertension = 0,
  vit_d        = median(study_dat$vit_d, na.rm = TRUE)
)

pred <- predict(rcs_fit, newdata = new_dat, type = "link", se.fit = TRUE)
se   <- if (is.list(pred)) pred$se.fit else sqrt(attr(pred, "var"))
fit  <- if (is.list(pred)) pred$fit    else as.numeric(pred)

plot_df <- data.frame(
  cally = exp(new_x),
  or    = exp(fit - fit[ref_idx]),
  lo    = exp(fit - 1.96 * se - fit[ref_idx]),
  hi    = exp(fit + 1.96 * se - fit[ref_idx])
)

p_rcs <- ggplot(plot_df, aes(x = cally, y = or)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, fill = "#2C7FB8") +
  geom_line(color = "#2C7FB8", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_x_log10() +
  labs(x = "CALLY index", y = "OR for osteoporosis (95% CI)") +
  sci_theme()

save_fig_ggplot(p_rcs, "Figure_S7_CALLY_RCS_Osteoporosis",
                width_mm = WIDTH_SINGLE, height_mm = 75)

message("\nNon-linearity test for CALLY (RCS, df=4):")
print(regTermTest(rcs_fit, ~ns(log_cally, df = 4)))


# ---- 11b. CALLY RCS on continuous BMD (Figure S8) --------------------------
rcs_bmd_fit <- svyglm(bmd_fn ~ ns(log_cally, df = 4) +
                               age + race + bmi + edu_cat + pir +
                               smoke_status + alcohol + diabetes +
                               hypertension + vit_d,
                      design = design_op)

pred_bmd <- predict(rcs_bmd_fit, newdata = new_dat, type = "response",
                    se.fit = TRUE)
se_bmd  <- if (is.list(pred_bmd)) pred_bmd$se.fit else sqrt(attr(pred_bmd, "var"))
fit_bmd <- if (is.list(pred_bmd)) pred_bmd$fit    else as.numeric(pred_bmd)

bmd_plot_df <- data.frame(
  cally = exp(new_x),
  est   = fit_bmd,
  lo    = fit_bmd - 1.96 * se_bmd,
  hi    = fit_bmd + 1.96 * se_bmd
)

p_rcs_bmd <- ggplot(bmd_plot_df, aes(x = cally, y = est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, fill = "#2C7FB8") +
  geom_line(color = "#2C7FB8", linewidth = 0.8) +
  geom_hline(yintercept = median(bmd_plot_df$est),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  scale_x_log10() +
  labs(x = "CALLY index",
       y = expression("Adjusted femoral neck BMD (g/cm"^2*")")) +
  sci_theme()

save_fig_ggplot(p_rcs_bmd, "Figure_S8_CALLY_RCS_BMD",
                width_mm = WIDTH_SINGLE, height_mm = 75)


# ---- 12. CALLY subgroup analyses + interaction tests -----------------------
subgroup_one <- function(var_name, exposure = "log_cally") {
  groups <- na.omit(unique(study_dat[[var_name]]))
  do.call(rbind, lapply(groups, function(g) {
    sub_design <- subset(design_op, get(var_name) == g)
    f <- as.formula(paste0(
      "osteoporosis ~ ", exposure,
      " + age + race + bmi + edu_cat + pir + smoke_status + alcohol + ",
      "diabetes + hypertension + vit_d"))
    fit <- tryCatch(
      svyglm(f, design = sub_design, family = quasibinomial()),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    ci <- exp_ci(fit, exposure)
    data.frame(Variable = var_name, Subgroup = as.character(g),
               OR = ci$OR, L = ci$L95, U = ci$U95, p = ci$p)
  }))
}

sub_res <- bind_rows(
  subgroup_one("age_cat"),
  subgroup_one("bmi_cat"),
  subgroup_one("diabetes"),
  subgroup_one("smoke_status")
)

p_int_tbl <- data.frame(Variable = character(), P_int = numeric())
for (v in c("age_cat","bmi_cat","diabetes","smoke_status")) {
  f_int <- as.formula(paste0(
    "osteoporosis ~ log_cally * ", v,
    " + age + race + bmi + edu_cat + pir + smoke_status + alcohol + ",
    "diabetes + hypertension + vit_d"))
  fit_i <- tryCatch(
    svyglm(f_int, design = design_op, family = quasibinomial()),
    error = function(e) NULL
  )
  if (!is.null(fit_i)) {
    test <- regTermTest(fit_i, as.formula(paste0("~log_cally:", v)))
    p_int_tbl <- rbind(p_int_tbl,
                       data.frame(Variable = v,
                                  P_int    = as.numeric(test$p)))
  }
}


# ---- 12b. CALLY subgroup forest (Figure S9) --------------------------------
var_label <- c(age_cat = "Age, years",
               bmi_cat = "BMI category",
               diabetes = "Diabetes",
               smoke_status = "Smoking status")

sub_plot_df <- sub_res %>%
  mutate(VarLabel = var_label[Variable],
         SubLabel = case_when(
           Variable == "diabetes" & Subgroup == "1" ~ "Yes",
           Variable == "diabetes" & Subgroup == "0" ~ "No",
           TRUE ~ Subgroup
         )) %>%
  left_join(p_int_tbl, by = "Variable") %>%
  mutate(full_label = paste0("  ", SubLabel),
         row_id     = row_number())

header_rows <- sub_plot_df %>%
  distinct(Variable, VarLabel, P_int) %>%
  mutate(SubLabel = "", OR = NA, L = NA, U = NA, p = NA,
         full_label = sprintf("%s (P-int = %.3f)", VarLabel,
                              ifelse(is.na(P_int), NA, P_int))) %>%
  mutate(row_id = -row_number())

forest_df <- bind_rows(header_rows, sub_plot_df) %>%
  arrange(Variable, row_id) %>%
  mutate(yorder = rev(row_number()),
         full_label = factor(full_label, levels = full_label[order(yorder)]))

p_forest_sub <- ggplot(forest_df, aes(x = OR, y = full_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.2, linewidth = 0.4,
                 color = "black", na.rm = TRUE) +
  geom_point(size = 1.8, color = "#D7301F", na.rm = TRUE) +
  geom_text(aes(label = ifelse(is.na(OR), "",
                               sprintf("%.2f (%.2f-%.2f)", OR, L, U))),
            hjust = -0.12, size = 2.3, na.rm = TRUE) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.45))) +
  labs(x = "OR for osteoporosis per unit log(CALLY)", y = NULL) +
  sci_theme() +
  theme(axis.text.y = element_text(hjust = 0))

save_fig_ggplot(p_forest_sub, "Figure_S9_CALLY_Subgroup",
                width_mm = WIDTH_ONEHALF, height_mm = 110)

tbl3 <- sub_res %>%
  left_join(p_int_tbl, by = "Variable") %>%
  mutate(
    `Variable label` = var_label[Variable],
    `Subgroup label` = case_when(
      Variable == "diabetes" & Subgroup == "1" ~ "Yes",
      Variable == "diabetes" & Subgroup == "0" ~ "No",
      TRUE ~ Subgroup
    ),
    `OR (95% CI)`   = sprintf("%.2f (%.2f-%.2f)", OR, L, U),
    `P value`       = fmt_p(p),
    `P interaction` = fmt_p(P_int)
  ) %>%
  select(`Variable` = `Variable label`,
         `Subgroup` = `Subgroup label`,
         OR, `Lower 95% CI` = L, `Upper 95% CI` = U,
         `OR (95% CI)`, `P value`, `P interaction`) %>%
  group_by(Variable) %>%
  mutate(`P interaction` = ifelse(row_number() == 1, `P interaction`, "")) %>%
  ungroup()
save_table(tbl3, "Table_S3_CALLY_Subgroup",
           sheet_name = "Table S3 CALLY Subgroup")


# ---- 13. CALLY sensitivity analyses ----------------------------------------
m_s1 <- svyglm(osteoporosis ~ log_cally + age + race + bmi + edu_cat + pir +
                              smoke_status + alcohol + diabetes + hypertension + vit_d,
               design = subset(design_op, crp <= 10),
               family = quasibinomial())

m_s2 <- svyglm(osteoporosis ~ log_cally + age + race + bmi + edu_cat + pir +
                              smoke_status + alcohol + hypertension + vit_d,
               design = subset(design_op, diabetes == 0),
               family = quasibinomial())

study_dat2 <- study_dat %>%
  mutate(t_score_tot = (bmd_tot - 0.942) / 0.122,
         op_tot      = ifelse(t_score_tot <= -2.5, 1, 0))
design_tot <- svydesign(id = ~psu, strata = ~strata, weights = ~wt_new,
                        nest = TRUE,
                        data = filter(study_dat2, !is.na(op_tot)))
m_s3 <- svyglm(op_tot ~ log_cally + age + race + bmi + edu_cat + pir +
                        smoke_status + alcohol + diabetes + hypertension + vit_d,
               design = design_tot, family = quasibinomial())

tbl_s1 <- bind_rows(
  cbind(Analysis = "S1: Excluding participants with CRP > 10 mg/L",
        exp_ci(m_s1, "log_cally")),
  cbind(Analysis = "S2: Excluding participants with diabetes",
        exp_ci(m_s2, "log_cally")),
  cbind(Analysis = "S3: Total hip BMD-defined osteoporosis",
        exp_ci(m_s3, "log_cally"))
) %>%
  mutate(`OR (95% CI), P` = fmt_or_ci_p(OR, L95, U95, p)) %>%
  select(Analysis, Exposure = Variable, OR,
         `Lower 95% CI` = L95, `Upper 95% CI` = U95,
         `P value` = p, `OR (95% CI), P`)
save_table(tbl_s1, "Table_S4_CALLY_Sensitivity",
           sheet_name = "Table S4 CALLY Sensitivity")


# ---- 13b. CALLY sensitivity forest (Figure S10) ----------------------------
sens_forest_df <- bind_rows(
  cbind(Analysis = "Main analysis (Model 3)",
        exp_ci(m3_log, "log_cally")),
  cbind(Analysis = "S1: Excluding CRP > 10 mg/L",
        exp_ci(m_s1, "log_cally")),
  cbind(Analysis = "S2: Excluding diabetes",
        exp_ci(m_s2, "log_cally")),
  cbind(Analysis = "S3: Total hip BMD outcome",
        exp_ci(m_s3, "log_cally"))
) %>% mutate(Analysis = factor(Analysis, levels = rev(Analysis)))

p_forest_sens <- ggplot(sens_forest_df, aes(x = OR, y = Analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.15, linewidth = 0.5,
                 color = "black") +
  geom_point(size = 2.0, color = "#2C7FB8") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), P=%s",
                                OR, L95, U95, fmt_p(p))),
            hjust = -0.10, size = 2.3) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.65))) +
  labs(x = "OR per unit log(CALLY)", y = NULL) +
  sci_theme()

save_fig_ggplot(p_forest_sens, "Figure_S10_CALLY_Sensitivity",
                width_mm = WIDTH_ONEHALF, height_mm = 65)


# ---- 14. Mediation analysis on CALLY (BMI and 25(OH)D) ---------------------
med_dat <- as.data.frame(study_dat) %>%
  dplyr::select(osteoporosis, log_cally, bmi, vit_d, age, race, edu_cat, pir,
                smoke_status, alcohol, diabetes, hypertension) %>%
  na.omit()

message("\nMediation by BMI:")
m_bmi <- lm(bmi ~ log_cally + age + race + edu_cat + pir +
                  smoke_status + alcohol + diabetes + hypertension,
            data = med_dat)
y_bmi <- glm(osteoporosis ~ bmi + log_cally + age + race + edu_cat + pir +
                            smoke_status + alcohol + diabetes + hypertension,
             family = binomial(), data = med_dat)
res_bmi <- mediate(m_bmi, y_bmi, treat = "log_cally", mediator = "bmi",
                   boot = TRUE, sims = 500)
print(summary(res_bmi))

message("\nMediation by 25(OH)D:")
if (sum(!is.na(med_dat$vit_d)) > 100) {
  m_vd <- lm(vit_d ~ log_cally + age + race + edu_cat + pir +
                     smoke_status + alcohol + diabetes + hypertension,
             data = med_dat)
  y_vd <- glm(osteoporosis ~ vit_d + log_cally + age + race + edu_cat + pir +
                             smoke_status + alcohol + diabetes + hypertension,
              family = binomial(), data = med_dat)
  res_vd <- mediate(m_vd, y_vd, treat = "log_cally", mediator = "vit_d",
                    boot = TRUE, sims = 500)
  print(summary(res_vd))
} else {
  message("Vit D too sparse, skipping")
  res_vd <- NULL
}

extract_med <- function(res, mediator_name) {
  if (is.null(res)) return(NULL)
  s <- summary(res)
  data.frame(
    Mediator = mediator_name,
    Effect   = c("ACME (indirect effect)",
                 "ADE (direct effect)",
                 "Total effect",
                 "Proportion mediated"),
    Estimate = c(s$d.avg, s$z.avg, s$tau.coef, s$n.avg),
    `Lower 95% CI` = c(s$d.avg.ci[1], s$z.avg.ci[1],
                       s$tau.ci[1],   s$n.avg.ci[1]),
    `Upper 95% CI` = c(s$d.avg.ci[2], s$z.avg.ci[2],
                       s$tau.ci[2],   s$n.avg.ci[2]),
    `P value` = c(s$d.avg.p, s$z.avg.p, s$tau.p, s$n.avg.p),
    check.names = FALSE
  )
}
tbl_s2 <- bind_rows(
  extract_med(res_bmi, "BMI"),
  extract_med(res_vd,  "25(OH)D")
)
if (!is.null(tbl_s2) && nrow(tbl_s2) > 0) {
  tbl_s2 <- tbl_s2 %>%
    mutate(
      Estimate       = round(Estimate, 4),
      `Lower 95% CI` = round(`Lower 95% CI`, 4),
      `Upper 95% CI` = round(`Upper 95% CI`, 4),
      `Estimate (95% CI), P` = ifelse(
        Effect == "Proportion mediated",
        sprintf("%.1f%% (%.1f%% - %.1f%%), P=%s",
                Estimate * 100, `Lower 95% CI` * 100, `Upper 95% CI` * 100,
                fmt_p(`P value`)),
        sprintf("%.4f (%.4f - %.4f), P=%s",
                Estimate, `Lower 95% CI`, `Upper 95% CI`,
                fmt_p(`P value`))
      )
    )
  save_table(tbl_s2, "Table_S5_CALLY_Mediation",
             sheet_name = "Table S5 CALLY Mediation")
}


# ---- 14b. Mediation forest (Figure S11) ------------------------------------
if (exists("tbl_s2") && !is.null(tbl_s2) && nrow(tbl_s2) > 0) {

  med_plot_df <- tbl_s2 %>%
    filter(Effect %in% c("ACME (indirect effect)",
                          "ADE (direct effect)",
                          "Total effect")) %>%
    mutate(EffectLabel = factor(
      recode(Effect,
             "ACME (indirect effect)" = "ACME (indirect)",
             "ADE (direct effect)"    = "ADE (direct)",
             "Total effect"           = "Total"),
      levels = c("Total", "ADE (direct)", "ACME (indirect)")
    ))

  prop_anno <- tbl_s2 %>%
    filter(Effect == "Proportion mediated") %>%
    mutate(anno = sprintf(
      "Proportion mediated: %.1f%%\n(95%% CI: %.1f%% - %.1f%%, P=%s)",
      Estimate * 100, `Lower 95% CI` * 100, `Upper 95% CI` * 100,
      fmt_p(`P value`))) %>%
    select(Mediator, anno)

  p_med <- ggplot(med_plot_df, aes(x = Estimate, y = EffectLabel)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = `Lower 95% CI`, xmax = `Upper 95% CI`),
                   height = 0.18, linewidth = 0.5, color = "black") +
    geom_point(size = 2.0, color = "#D7301F") +
    geom_text(aes(label = sprintf("%.4f (%.4f, %.4f)",
                                   Estimate, `Lower 95% CI`, `Upper 95% CI`)),
              hjust = -0.08, size = 2.0) +
    geom_text(data = prop_anno,
              aes(x = -Inf, y = -Inf, label = anno),
              hjust = -0.05, vjust = -0.6,
              size = 1.9, color = "grey25",
              inherit.aes = FALSE) +
    facet_wrap(~ Mediator, ncol = 1, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.45))) +
    labs(x = "Effect size on log-odds scale (95% CI)", y = NULL) +
    sci_theme() +
    theme(strip.background = element_rect(fill = "grey92", color = NA),
          strip.text       = element_text(face = "bold", size = 8))

  save_fig_ggplot(p_med, "Figure_S11_CALLY_Mediation",
                  width_mm = WIDTH_ONEHALF, height_mm = 100)
}


# ---- 15. Head-to-head comparison of seven indices --------------------------
# 15a: per-SD adjusted ORs on the standardized log scale (Figure 2, Table 2)
# 15b: ROC + DeLong against the best-performing index   (Figures 3-4, Table 3)
# 15c: Spearman correlation heatmap                     (Figure 5)

message("\n=== Head-to-head comparison of 7 indices ===")

# z-scoring each log-index makes ORs comparable across indices that have
# very different natural scales (PNI ~40-100, CALLY ~100-2000+)
study_dat <- study_dat %>%
  mutate(across(all_of(INDEX_INFO$log_var),
                ~ as.numeric(scale(.x)),
                .names = "z_{.col}"))

design    <- svydesign(id = ~psu, strata = ~strata, weights = ~wt_new,
                       nest = TRUE, data = study_dat)
design_op <- subset(design, !is.na(osteoporosis))

# 15a. Multi-index Model 3 ORs --------------------------------------------
multi_or_results <- list()
for (i in seq_len(nrow(INDEX_INFO))) {
  idx_short <- INDEX_INFO$short[i]
  idx_label <- INDEX_INFO$label[i]
  z_var     <- paste0("z_", INDEX_INFO$log_var[i])

  f_str <- paste0("osteoporosis ~ ", z_var,
                  " + age + race + bmi + edu_cat + pir + smoke_status + ",
                  "alcohol + diabetes + hypertension + vit_d")
  fit <- tryCatch(
    svyglm(as.formula(f_str), design = design_op, family = quasibinomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  ci_row <- exp_ci(fit, z_var)
  multi_or_results[[idx_short]] <- data.frame(
    Index     = idx_label,
    Direction = ifelse(INDEX_INFO$protective[i], "Protective", "Risk"),
    OR = ci_row$OR, L = ci_row$L95, U = ci_row$U95, p = ci_row$p,
    stringsAsFactors = FALSE
  )
}
multi_or <- bind_rows(multi_or_results)
print(multi_or)

tbl2 <- multi_or %>%
  mutate(`OR (95% CI), P` = fmt_or_ci_p(OR, L, U, p)) %>%
  select(Index, Direction,
         `OR per SD` = OR, `Lower 95% CI` = L, `Upper 95% CI` = U,
         `P value` = p, `OR (95% CI), P`)
save_table(tbl2, "Table_2_MultiIndex_AdjustedOR",
           sheet_name = "Table 2 MultiIndex OR")

multi_or_plot <- multi_or %>%
  arrange(OR) %>%
  mutate(Index = factor(Index, levels = unique(Index)))

p_multi_or <- ggplot(multi_or_plot, aes(x = OR, y = Index, color = Direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.18, linewidth = 0.5) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), P=%s",
                                OR, L, U, fmt_p(p))),
            hjust = -0.10, size = 2.3, color = "black") +
  scale_color_manual(values = c("Protective" = "#2C7FB8",
                                "Risk"       = "#D7301F")) +
  scale_x_log10(expand = expansion(mult = c(0.08, 0.70))) +
  labs(x = "Adjusted OR for osteoporosis per SD of log-transformed index",
       y = NULL) +
  sci_theme() +
  theme(legend.position = "top")

save_fig_ggplot(p_multi_or, "Figure_2_MultiIndex_AdjustedOR",
                width_mm = WIDTH_ONEHALF, height_mm = 80)


# 15b. ROC + DeLong --------------------------------------------------------
roc_dat_multi <- study_dat %>%
  filter(!is.na(osteoporosis)) %>%
  filter(if_all(all_of(INDEX_INFO$raw_var), ~ !is.na(.x) & is.finite(.x)))

roc_list <- list()
auc_rows <- list()
for (i in seq_len(nrow(INDEX_INFO))) {
  idx_short <- INDEX_INFO$short[i]
  idx_label <- INDEX_INFO$label[i]
  direction <- INDEX_INFO$direction[i]
  r <- roc(roc_dat_multi$osteoporosis, roc_dat_multi[[idx_short]],
           direction = direction, quiet = TRUE)
  roc_list[[idx_short]] <- r
  ci_r <- as.numeric(ci.auc(r))
  auc_rows[[idx_short]] <- data.frame(
    Index = idx_label, AUC = as.numeric(auc(r)),
    L = ci_r[1], U = ci_r[3],
    stringsAsFactors = FALSE
  )
  message(sprintf("  %-6s AUC = %.3f (95%% CI: %.3f-%.3f)",
                  idx_label, auc(r), ci_r[1], ci_r[3]))
}
auc_tbl <- bind_rows(auc_rows) %>% arrange(desc(AUC))

# The manuscript reports MLR as the winner. If a different index ends up on
# top after a data update, this will be obvious in the printed AUC table
# above and downstream Section 17 will simply describe that index instead.
winner_label <- auc_tbl$Index[1]
winner_short <- INDEX_INFO$short[match(winner_label, INDEX_INFO$label)]
message(sprintf("Best AUC: %s (%.3f)", winner_label, auc_tbl$AUC[1]))

delong_rows <- list()
for (i in seq_len(nrow(auc_tbl))) {
  idx_label <- auc_tbl$Index[i]
  idx_short <- INDEX_INFO$short[match(idx_label, INDEX_INFO$label)]
  if (idx_short == winner_short) {
    delong_rows[[idx_short]] <- "- (reference, highest AUC)"
  } else {
    t <- roc.test(roc_list[[winner_short]], roc_list[[idx_short]])
    delong_rows[[idx_short]] <- sprintf("Z = %.3f, P = %s",
                                         as.numeric(t$statistic),
                                         fmt_p(t$p.value))
  }
}

tbl3 <- auc_tbl %>%
  mutate(`AUC (95% CI)` = sprintf("%.3f (%.3f - %.3f)", AUC, L, U)) %>%
  rowwise() %>%
  mutate(`DeLong test vs best` = {
    s <- INDEX_INFO$short[match(Index, INDEX_INFO$label)]
    delong_rows[[s]]
  }) %>%
  ungroup() %>%
  select(Index, AUC,
         `Lower 95% CI` = L, `Upper 95% CI` = U,
         `AUC (95% CI)`, `DeLong test vs best`)
save_table(tbl3, "Table_3_MultiIndex_AUC_DeLong",
           sheet_name = "Table 3 MultiIndex AUC")

roc_palette <- c("CALLY" = "#E41A1C", "SII" = "#377EB8", "NLR" = "#4DAF4A",
                 "PLR"   = "#984EA3", "MLR" = "#FF7F00", "PNI" = "#A65628",
                 "HALP"  = "#F781BF")

save_fig_base(
  {
    par(mar = c(4, 4, 1, 1), family = "sans", cex = 0.85,
        mgp = c(2.3, 0.6, 0), tcl = -0.3)
    plot(roc_list[[INDEX_INFO$short[1]]],
         col = roc_palette[INDEX_INFO$label[1]], lwd = 1.4,
         legacy.axes = TRUE,
         xlab = "1 - Specificity", ylab = "Sensitivity", main = "")
    for (i in 2:nrow(INDEX_INFO)) {
      lines(roc_list[[INDEX_INFO$short[i]]],
            col = roc_palette[INDEX_INFO$label[i]], lwd = 1.4)
    }
    abline(a = 0, b = 1, lty = 3, col = "grey60")
    leg_labs <- sprintf("%-5s (AUC = %.3f)", auc_tbl$Index, auc_tbl$AUC)
    legend("bottomright", legend = leg_labs,
           col = roc_palette[auc_tbl$Index],
           lwd = 1.4, bty = "n", cex = 0.72)
  },
  filename  = "Figure_3_MultiIndex_ROC_Curves",
  width_mm  = WIDTH_SINGLE, height_mm = 95
)

auc_bar_df <- auc_tbl %>%
  mutate(Index = factor(Index, levels = Index))

p_auc_bar <- ggplot(auc_bar_df, aes(x = Index, y = AUC,
                                     fill = Index == winner_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = L, ymax = U), width = 0.20, linewidth = 0.4) +
  geom_text(aes(y = U, label = sprintf("%.3f", AUC)),
            vjust = -0.8, size = 2.2) +
  scale_fill_manual(values = c(`FALSE` = "grey75", `TRUE` = "#2C7FB8"),
                    guide = "none") +
  coord_cartesian(ylim = c(0.35, max(auc_bar_df$U) + 0.06)) +
  labs(x = NULL, y = "AUC (95% CI)") +
  sci_theme()

save_fig_ggplot(p_auc_bar, "Figure_4_MultiIndex_AUC_Barplot",
                width_mm = WIDTH_SINGLE, height_mm = 70)


# 15c. Spearman correlation heatmap ---------------------------------------
cor_mat <- cor(study_dat[, INDEX_INFO$raw_var],
               use = "pairwise.complete.obs", method = "spearman")
colnames(cor_mat) <- INDEX_INFO$label
rownames(cor_mat) <- INDEX_INFO$label

cor_long <- as.data.frame(as.table(cor_mat))
names(cor_long) <- c("X", "Y", "rho")

p_corr <- ggplot(cor_long, aes(x = X, y = Y, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)),
            size = 2.4,
            color = ifelse(abs(cor_long$rho) > 0.6, "white", "black")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  sci_theme() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
        panel.grid  = element_blank(),
        axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        legend.position = "right",
        legend.key.size = unit(4, "mm"))

save_fig_ggplot(p_corr, "Figure_5_MultiIndex_Correlation",
                width_mm = WIDTH_SINGLE, height_mm = 85)


# ---- 17. Detailed analysis of the best-performing index (MLR) --------------
# Figures S1-S4 and Tables S6-S7. Parallel structure to the CALLY block
# (Sections 10-13) for visual consistency in the supplement.

message("\n=== Detailed analyses for best-performing index (", winner_label, ") ===")

study_dat <- study_dat %>%
  mutate(
    mlr_q_num = ntile(mlr, 4),
    mlr_q     = factor(mlr_q_num, levels = 1:4,
                       labels = c("Q1","Q2","Q3","Q4"))
  )
design    <- svydesign(id = ~psu, strata = ~strata, weights = ~wt_new,
                       nest = TRUE, data = study_dat)
design_op <- subset(design, !is.na(osteoporosis))


# 17a. MLR three-model and quartile (Table S6, Figure S1) --------------------
m1_mlr <- svyglm(osteoporosis ~ log_mlr,
                 design = design_op, family = quasibinomial())
m2_mlr <- svyglm(osteoporosis ~ log_mlr + age + race + bmi,
                 design = design_op, family = quasibinomial())
m3_mlr <- svyglm(osteoporosis ~ log_mlr + age + race + bmi + edu_cat + pir +
                                smoke_status + alcohol + diabetes +
                                hypertension + vit_d,
                 design = design_op, family = quasibinomial())

mlr_three <- bind_rows(
  cbind(Model = "Model 1 (Crude)",            exp_ci(m1_mlr, "log_mlr")),
  cbind(Model = "Model 2 (Age + Race + BMI)", exp_ci(m2_mlr, "log_mlr")),
  cbind(Model = "Model 3 (Fully adjusted)",   exp_ci(m3_mlr, "log_mlr"))
)
print(mlr_three)

# Also report the more interpretable per-SD log(MLR) and per-doubling MLR
# OR for the fully adjusted model. The per-SD figure matches Table 2 for
# cross-index comparison; the per-doubling figure is what readers want for a
# biological interpretation.
sd_log_mlr <- sd(study_dat$log_mlr, na.rm = TRUE)
b_mlr <- coef(m3_mlr)["log_mlr"]
v_mlr <- vcov(m3_mlr)["log_mlr", "log_mlr"]
or_per_sd  <- exp(b_mlr * sd_log_mlr)
ci_per_sd  <- exp(c(b_mlr - 1.96 * sqrt(v_mlr),
                    b_mlr + 1.96 * sqrt(v_mlr)) * sd_log_mlr)
or_per_dbl <- exp(b_mlr * log(2))
ci_per_dbl <- exp(c(b_mlr - 1.96 * sqrt(v_mlr),
                    b_mlr + 1.96 * sqrt(v_mlr)) * log(2))
message(sprintf("MLR per 1 SD of log(MLR) [SD = %.3f]: OR = %.2f (%.2f-%.2f)",
                sd_log_mlr, or_per_sd, ci_per_sd[1], ci_per_sd[2]))
message(sprintf("MLR per doubling of MLR             : OR = %.2f (%.2f-%.2f)",
                or_per_dbl, ci_per_dbl[1], ci_per_dbl[2]))

m3_mlr_q     <- svyglm(osteoporosis ~ mlr_q + age + race + bmi + edu_cat + pir +
                                      smoke_status + alcohol + diabetes +
                                      hypertension + vit_d,
                       design = design_op, family = quasibinomial())
m3_mlr_trend <- svyglm(osteoporosis ~ mlr_q_num + age + race + bmi + edu_cat + pir +
                                      smoke_status + alcohol + diabetes +
                                      hypertension + vit_d,
                       design = design_op, family = quasibinomial())
mlr_p_trend <- summary(m3_mlr_trend)$coefficients["mlr_q_num", "Pr(>|t|)"]
mlr_q_all   <- exp_ci(m3_mlr_q)
mlr_q_rows  <- mlr_q_all[grepl("^mlr_q", mlr_q_all$Variable), ]

tbl_s6 <- bind_rows(
  mlr_three %>%
    mutate(Section  = Model,
           Exposure = "log(MLR), continuous",
           `OR (95% CI), P` = fmt_or_ci_p(OR, L95, U95, p)) %>%
    select(Section, Exposure,
           OR, `Lower 95% CI` = L95, `Upper 95% CI` = U95,
           `P value` = p, `OR (95% CI), P`),
  data.frame(
    Section  = "Per SD of log(MLR) (Model 3)",
    Exposure = sprintf("log(MLR) per 1 SD (=%.3f)", sd_log_mlr),
    OR  = round(or_per_sd, 3),
    `Lower 95% CI` = round(ci_per_sd[1], 3),
    `Upper 95% CI` = round(ci_per_sd[2], 3),
    `P value` = NA_real_,
    `OR (95% CI), P` = sprintf("%.2f (%.2f-%.2f)",
                               or_per_sd, ci_per_sd[1], ci_per_sd[2]),
    check.names = FALSE
  ),
  data.frame(
    Section  = "Per doubling of MLR (Model 3)",
    Exposure = "MLR per 2-fold increase",
    OR  = round(or_per_dbl, 3),
    `Lower 95% CI` = round(ci_per_dbl[1], 3),
    `Upper 95% CI` = round(ci_per_dbl[2], 3),
    `P value` = NA_real_,
    `OR (95% CI), P` = sprintf("%.2f (%.2f-%.2f)",
                               or_per_dbl, ci_per_dbl[1], ci_per_dbl[2]),
    check.names = FALSE
  ),
  data.frame(
    Section  = "Quartile (Model 3)",
    Exposure = c("Q1 (reference)", "Q2", "Q3", "Q4", "P for trend"),
    OR  = c(1.00, mlr_q_rows$OR, NA),
    `Lower 95% CI` = c(NA, mlr_q_rows$L95, NA),
    `Upper 95% CI` = c(NA, mlr_q_rows$U95, NA),
    `P value` = c(NA, mlr_q_rows$p, mlr_p_trend),
    `OR (95% CI), P` = c(
      "1.00 (reference)",
      fmt_or_ci_p(mlr_q_rows$OR, mlr_q_rows$L95, mlr_q_rows$U95,
                  mlr_q_rows$p),
      sprintf("P-trend = %s", fmt_p(mlr_p_trend))
    ),
    check.names = FALSE
  )
)
save_table(tbl_s6, "Table_S6_MLR_OR_Detailed",
           sheet_name = "Table S6 MLR Detailed OR")

mlr_fdf_cont <- mlr_three %>%
  transmute(Label = factor(Model, levels = rev(Model)),
            OR, L95, U95, p,
            Panel = "Continuous (per unit log(MLR))")
mlr_fdf_quart <- data.frame(
  Label = factor(c("Q4","Q3","Q2","Q1 (reference)"),
                 levels = c("Q4","Q3","Q2","Q1 (reference)")),
  OR  = c(mlr_q_rows$OR[3],  mlr_q_rows$OR[2],  mlr_q_rows$OR[1],  1.00),
  L95 = c(mlr_q_rows$L95[3], mlr_q_rows$L95[2], mlr_q_rows$L95[1], NA),
  U95 = c(mlr_q_rows$U95[3], mlr_q_rows$U95[2], mlr_q_rows$U95[1], NA),
  p   = c(mlr_q_rows$p[3],   mlr_q_rows$p[2],   mlr_q_rows$p[1],   NA),
  Panel = "Quartile (vs Q1, Model 3)"
)
mlr_forest_df <- bind_rows(mlr_fdf_cont, mlr_fdf_quart) %>%
  mutate(Panel = factor(Panel, levels = c(
    "Continuous (per unit log(MLR))", "Quartile (vs Q1, Model 3)")))

p_mlr_forest <- ggplot(mlr_forest_df, aes(x = OR, y = Label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.18,
                 linewidth = 0.5, color = "black", na.rm = TRUE) +
  geom_point(size = 2.0, color = "#2C7FB8", na.rm = TRUE) +
  geom_text(aes(label = ifelse(is.na(L95),
                               "1.00 (reference)",
                               sprintf("%.2f (%.2f-%.2f), P=%s",
                                       OR, L95, U95, fmt_p(p)))),
            hjust = -0.10, size = 2.3) +
  facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.60))) +
  labs(x = "OR for osteoporosis (95% CI)", y = NULL) +
  sci_theme() +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold", size = 8))

save_fig_ggplot(p_mlr_forest, "Figure_S1_MLR_OR_Detailed",
                width_mm = WIDTH_ONEHALF, height_mm = 95)


# 17b. MLR restricted cubic spline (Figure S2) -------------------------------
rcs_mlr_fit <- svyglm(osteoporosis ~ ns(log_mlr, df = 4) +
                                     age + race + bmi + edu_cat + pir +
                                     smoke_status + alcohol + diabetes +
                                     hypertension + vit_d,
                      design = design_op, family = quasibinomial())

new_x_mlr <- seq(quantile(study_dat$log_mlr, 0.025, na.rm = TRUE),
                 quantile(study_dat$log_mlr, 0.975, na.rm = TRUE),
                 length.out = 200)
ref_idx_mlr <- which.min(abs(new_x_mlr - median(study_dat$log_mlr, na.rm = TRUE)))

new_dat_mlr <- data.frame(
  log_mlr      = new_x_mlr,
  age          = median(study_dat$age, na.rm = TRUE),
  race         = factor("NH White",     levels = levels(study_dat$race)),
  bmi          = median(study_dat$bmi, na.rm = TRUE),
  edu_cat      = factor(">High school", levels = levels(study_dat$edu_cat)),
  pir          = median(study_dat$pir, na.rm = TRUE),
  smoke_status = factor("Never",        levels = levels(study_dat$smoke_status)),
  alcohol      = factor("None",         levels = levels(study_dat$alcohol)),
  diabetes     = 0,
  hypertension = 0,
  vit_d        = median(study_dat$vit_d, na.rm = TRUE)
)

pred_mlr <- predict(rcs_mlr_fit, newdata = new_dat_mlr,
                    type = "link", se.fit = TRUE)
se_mlr  <- if (is.list(pred_mlr)) pred_mlr$se.fit else sqrt(attr(pred_mlr, "var"))
fit_mlr <- if (is.list(pred_mlr)) pred_mlr$fit    else as.numeric(pred_mlr)

rcs_mlr_df <- data.frame(
  mlr = exp(new_x_mlr),
  or  = exp(fit_mlr - fit_mlr[ref_idx_mlr]),
  lo  = exp(fit_mlr - 1.96 * se_mlr - fit_mlr[ref_idx_mlr]),
  hi  = exp(fit_mlr + 1.96 * se_mlr - fit_mlr[ref_idx_mlr])
)

p_mlr_rcs <- ggplot(rcs_mlr_df, aes(x = mlr, y = or)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, fill = "#2C7FB8") +
  geom_line(color = "#2C7FB8", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_x_log10() +
  labs(x = "MLR (Monocyte-to-Lymphocyte Ratio)",
       y = "OR for osteoporosis (95% CI)") +
  sci_theme()

save_fig_ggplot(p_mlr_rcs, "Figure_S2_MLR_RCS",
                width_mm = WIDTH_SINGLE, height_mm = 75)

message("Non-linearity test for MLR (RCS, df=4):")
print(regTermTest(rcs_mlr_fit, ~ns(log_mlr, df = 4)))


# 17c. MLR subgroup analyses (Figure S3) -------------------------------------
mlr_sub_res <- bind_rows(
  subgroup_one("age_cat",      exposure = "log_mlr"),
  subgroup_one("bmi_cat",      exposure = "log_mlr"),
  subgroup_one("diabetes",     exposure = "log_mlr"),
  subgroup_one("smoke_status", exposure = "log_mlr")
)

mlr_p_int_tbl <- data.frame(Variable = character(), P_int = numeric())
for (v in c("age_cat","bmi_cat","diabetes","smoke_status")) {
  f_int <- as.formula(paste0(
    "osteoporosis ~ log_mlr * ", v,
    " + age + race + bmi + edu_cat + pir + smoke_status + alcohol + ",
    "diabetes + hypertension + vit_d"))
  fit_i <- tryCatch(
    svyglm(f_int, design = design_op, family = quasibinomial()),
    error = function(e) NULL
  )
  if (!is.null(fit_i)) {
    test <- regTermTest(fit_i, as.formula(paste0("~log_mlr:", v)))
    mlr_p_int_tbl <- rbind(mlr_p_int_tbl,
                           data.frame(Variable = v,
                                      P_int    = as.numeric(test$p)))
  }
}

mlr_sub_plot <- mlr_sub_res %>%
  mutate(VarLabel = var_label[Variable],
         SubLabel = case_when(
           Variable == "diabetes" & Subgroup == "1" ~ "Yes",
           Variable == "diabetes" & Subgroup == "0" ~ "No",
           TRUE ~ Subgroup
         )) %>%
  left_join(mlr_p_int_tbl, by = "Variable") %>%
  mutate(full_label = paste0("  ", SubLabel),
         row_id     = row_number())

mlr_headers <- mlr_sub_plot %>%
  distinct(Variable, VarLabel, P_int) %>%
  mutate(SubLabel = "", OR = NA, L = NA, U = NA, p = NA,
         full_label = sprintf("%s (P-int = %s)", VarLabel,
                              ifelse(is.na(P_int), "NA", fmt_p(P_int)))) %>%
  mutate(row_id = -row_number())

mlr_forest <- bind_rows(mlr_headers, mlr_sub_plot) %>%
  arrange(Variable, row_id) %>%
  mutate(yorder = rev(row_number()),
         full_label = factor(full_label, levels = full_label[order(yorder)]))

p_mlr_subgroup <- ggplot(mlr_forest, aes(x = OR, y = full_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.2, linewidth = 0.4,
                 color = "black", na.rm = TRUE) +
  geom_point(size = 1.8, color = "#D7301F", na.rm = TRUE) +
  geom_text(aes(label = ifelse(is.na(OR), "",
                               sprintf("%.2f (%.2f-%.2f)", OR, L, U))),
            hjust = -0.12, size = 2.3, na.rm = TRUE) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.45))) +
  labs(x = "OR for osteoporosis per unit log(MLR)", y = NULL) +
  sci_theme() +
  theme(axis.text.y = element_text(hjust = 0))

save_fig_ggplot(p_mlr_subgroup, "Figure_S3_MLR_Subgroup",
                width_mm = WIDTH_ONEHALF, height_mm = 110)


# 17d. MLR sensitivity analyses (Figure S4, Table S7) ------------------------
ms1_mlr <- svyglm(osteoporosis ~ log_mlr + age + race + bmi + edu_cat + pir +
                                 smoke_status + alcohol + diabetes +
                                 hypertension + vit_d,
                  design = subset(design_op, crp <= 10),
                  family = quasibinomial())
ms2_mlr <- svyglm(osteoporosis ~ log_mlr + age + race + bmi + edu_cat + pir +
                                 smoke_status + alcohol + hypertension + vit_d,
                  design = subset(design_op, diabetes == 0),
                  family = quasibinomial())
ms3_mlr <- svyglm(op_tot ~ log_mlr + age + race + bmi + edu_cat + pir +
                          smoke_status + alcohol + diabetes +
                          hypertension + vit_d,
                  design = design_tot, family = quasibinomial())

mlr_sens <- bind_rows(
  cbind(Section = "Sensitivity",
        Analysis = "Main analysis (Model 3)",     exp_ci(m3_mlr,  "log_mlr")),
  cbind(Section = "Sensitivity",
        Analysis = "S1: Excluding CRP > 10 mg/L", exp_ci(ms1_mlr, "log_mlr")),
  cbind(Section = "Sensitivity",
        Analysis = "S2: Excluding diabetes",      exp_ci(ms2_mlr, "log_mlr")),
  cbind(Section = "Sensitivity",
        Analysis = "S3: Total hip BMD outcome",   exp_ci(ms3_mlr, "log_mlr"))
) %>%
  mutate(`OR (95% CI), P` = fmt_or_ci_p(OR, L95, U95, p)) %>%
  select(Section, Analysis, OR,
         `Lower 95% CI` = L95, `Upper 95% CI` = U95,
         `P value` = p, `OR (95% CI), P`)

mlr_sub_tbl <- mlr_sub_plot %>%
  left_join(mlr_p_int_tbl, by = "Variable") %>%
  mutate(Section  = "Subgroup",
         Analysis = sprintf("%s = %s", var_label[Variable], SubLabel),
         `OR (95% CI), P`  = fmt_or_ci_p(OR, L, U, p)) %>%
  select(Section, Analysis,
         OR, `Lower 95% CI` = L, `Upper 95% CI` = U,
         `P value` = p, `OR (95% CI), P`)

tbl_s7 <- bind_rows(mlr_sens, mlr_sub_tbl)
save_table(tbl_s7, "Table_S7_MLR_Subgroup_Sensitivity",
           sheet_name = "Table S7 MLR SubSens")

sens_mlr_df <- bind_rows(
  cbind(Analysis = "Main analysis (Model 3)",
        exp_ci(m3_mlr, "log_mlr")),
  cbind(Analysis = "S1: Excluding CRP > 10 mg/L",
        exp_ci(ms1_mlr, "log_mlr")),
  cbind(Analysis = "S2: Excluding diabetes",
        exp_ci(ms2_mlr, "log_mlr")),
  cbind(Analysis = "S3: Total hip BMD outcome",
        exp_ci(ms3_mlr, "log_mlr"))
) %>% mutate(Analysis = factor(Analysis, levels = rev(Analysis)))

p_mlr_sens <- ggplot(sens_mlr_df, aes(x = OR, y = Analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = L95, xmax = U95), height = 0.15,
                 linewidth = 0.5, color = "black") +
  geom_point(size = 2.0, color = "#2C7FB8") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), P=%s",
                                OR, L95, U95, fmt_p(p))),
            hjust = -0.10, size = 2.3) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.65))) +
  labs(x = "OR per unit log(MLR)", y = NULL) +
  sci_theme()

save_fig_ggplot(p_mlr_sens, "Figure_S4_MLR_Sensitivity",
                width_mm = WIDTH_ONEHALF, height_mm = 65)


# ---- 16. Save snapshot and bundle xlsx -------------------------------------
save(study_dat, design, design_op,
     m1_log, m2_log, m3_log, m3_q, m3_trend,
     rcs_fit, sub_res,
     res_bmi, res_vd,
     roc_list, auc_tbl, multi_or, cor_mat,
     INDEX_INFO, winner_label, winner_short,
     tab1,
     file = file.path(OUT_ROOT, "analysis_results.RData"))

message("\nBundling tables into a single xlsx workbook")
bundle_tables_xlsx("all_tables.xlsx")


# ---- 18. Session info ------------------------------------------------------
message("\n=== Session info ===")
print(sessionInfo())
message("\nDone. Outputs in: ", OUT_ROOT)
