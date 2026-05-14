# Inflammation–nutrition indices and postmenopausal osteoporosis

Analysis code for a cross-sectional comparison of seven inflammation–nutrition
composite indices (CALLY, SII, NLR, PLR, MLR, PNI, HALP) and postmenopausal
osteoporosis, using three cycles of NHANES (2005–2010).

## What this repo produces

Running `analysis.R` end-to-end yields:

- **Main figures** (`Figure_1` to `Figure_5`) and **main tables** (`Table_1` to `Table_3`)
- **Supplementary figures** (`Figure_S1` to `Figure_S11`) and **supplementary tables** (`Table_S1` to `Table_S7`)
- A single bundled `all_tables.xlsx` with every table as a separate sheet
- `analysis_results.RData` containing the fitted models for re-inspection

Each figure is written as TIFF (600 dpi, LZW), PDF (vector), and PNG (300 dpi),
sized to the standard single / 1.5 / double column widths.

## Reproducing the analysis

1. Install R ≥ 4.2.
2. Install the required packages once:
   ```r
   install.packages(c(
     "nhanesA", "survey", "tidyverse", "rms", "mediation",
     "pROC", "tableone", "splines", "ggplot2", "openxlsx", "readr"
   ))
   ```
3. From the repository root:
   ```r
   source("analysis.R")
   ```

The NHANES files are pulled live by the `nhanesA` package — no manual downloads
needed. Outputs land in `./output/figures` and `./output/tables`. End-to-end
runtime is roughly 5–10 minutes on a modern laptop, dominated by the NHANES
downloads and the 500-iteration bootstrap in the mediation analysis.

## Repository layout

```
analysis.R                 single-file analysis (run this)
output/                    created on first run
  figures/                 .tiff / .pdf / .png × 16 figures
  tables/                  .csv × 10 tables + all_tables.xlsx
  analysis_results.RData   saved models and summary objects
.gitignore
LICENSE
README.md
```

## Data

NHANES is a public, de-identified, IRB-approved survey conducted by the U.S.
National Center for Health Statistics. No additional approval is required to
analyse the released files. The `nhanesA` package is the official R interface;
versions of individual NHANES tables published after this script was written
should remain backward-compatible.

## Notes on a few methodological choices

- **Log-transformation without an offset.** All seven indices are strictly
  positive in valid cases, so we use `log(x)` rather than `log(x + 1)`. The
  `+1` offset that earlier drafts used artificially compressed the dynamic
  range of MLR (typical 0.05–1.0) and inflated its confidence interval on the
  per-unit log scale. The per-SD and per-doubling reporting in Section 17 also
  helps make the continuous effect size interpretable across indices that span
  several orders of magnitude.

- **Complex survey weights.** Three two-year NHANES cycles are combined, with
  MEC weights divided by 3 per NCHS guidance.
  `options(survey.lonely.psu = "adjust")` is set globally to handle the rare
  PSU-of-one stratum.

- **Same-sample comparison.** The analytic sample requires complete data for
  all seven indices, BMI, and femoral-neck BMD. This costs sample size but is
  essential for an honest head-to-head comparison — otherwise each index would
  be estimated on a different population.

- **Best-performing index, programmatically.** Section 17 runs detailed
  analyses for whichever index has the highest AUC, taken from the printed
  ROC summary. In the released dataset this is MLR; if a future NHANES
  refresh shifts the ranking, the script automatically follows the new winner.

- **Resilient NHANES downloads.** NCHS' wwwn.cdc.gov server occasionally
  returns transient errors or empty payloads on older cycles. `safe_nhanes()`
  retries each request up to four times with exponential backoff, and the
  download routine tries known filename variants (e.g. `VID` then `VID2`)
  before declaring a table missing. The script hard-stops if any
  *required* table (DEMO, BMX, DXXFEM, CBC, BIOPRO, CRP) cannot be
  retrieved for a cycle after retries — those tables drive sample
  definition and the analysis would be biased without them. Optional
  tables (MCQ for self-reported cancer history, VID for 25(OH)D, etc.)
  are tolerated and any per-cycle gaps are listed in the coverage summary
  printed before modelling begins.

- **Citation.** If you use this code, please cite the accompanying paper
  (under review).

## License

MIT — see `LICENSE`.
