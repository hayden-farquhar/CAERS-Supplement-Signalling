# =============================================================================
# 01_data_acquisition.R
# Download FDA CAERS/HFCS data and filter to dietary supplement reports
# =============================================================================
# NOTE: As of Oct 2024, CAERS was renamed to the Human Foods Complaint System
# (HFCS). The data content is identical; only the system name changed.
#
# Product-based CSV columns (after clean_names):
#   date_fda_first_received_report, report_id, date_event, product_type,
#   product, product_code, description, patient_age, age_units, sex,
#   case_meddra_preferred_terms, case_outcome
# =============================================================================

library(tidyverse)
library(data.table)
library(janitor)

# --- Configuration -----------------------------------------------------------

raw_dir    <- "data/raw"
proc_dir   <- "data/processed"

# FDA download URLs (from HFCS page, formerly CAERS)
url_product <- "https://www.fda.gov/media/161096/download?attachment"
url_case    <- "https://www.fda.gov/media/170793/download?attachment"

# Dietary supplement industry code
SUPPLEMENT_CODE <- 54

# --- 1. Download raw data ----------------------------------------------------

cat("=== Downloading CAERS/HFCS data ===\n")

product_file <- file.path(raw_dir, "CAERS_product_based.csv")
case_file    <- file.path(raw_dir, "CAERS_case_based.csv")

if (!file.exists(product_file)) {
  cat("Downloading product-based CSV...\n")
  download.file(url_product, product_file, mode = "wb", quiet = FALSE)
  cat("Done.\n")
} else {
  cat("Product-based CSV already exists, skipping download.\n")
}

if (!file.exists(case_file)) {
  cat("Downloading case-based CSV...\n")
  download.file(url_case, case_file, mode = "wb", quiet = FALSE)
  cat("Done.\n")
} else {
  cat("Case-based CSV already exists, skipping download.\n")
}

# --- 2. Load product-based data ----------------------------------------------

cat("\n=== Loading product-based data ===\n")
product_raw <- fread(product_file, encoding = "UTF-8") |> clean_names()

cat(sprintf("Total rows (product-level): %s\n", format(nrow(product_raw), big.mark = ",")))
cat(sprintf("Total unique reports: %s\n", format(n_distinct(product_raw$report_id), big.mark = ",")))

cat("\nColumn names:\n")
cat(paste(" -", names(product_raw)), sep = "\n")

# --- 3. Explore industry codes -----------------------------------------------

cat("\n\n=== Product industry codes ===\n")
industry_summary <- product_raw |>
  count(product_code, description, sort = TRUE)
print(as_tibble(industry_summary), n = 25)

# --- 4. Filter to dietary supplements (product_code 54) ----------------------

cat("\n=== Filtering to dietary supplements (code 54) ===\n")

supplements <- product_raw |>
  filter(product_code == SUPPLEMENT_CODE)

cat(sprintf("Supplement rows (all roles): %s\n", format(nrow(supplements), big.mark = ",")))
cat(sprintf("Supplement unique reports: %s\n", format(n_distinct(supplements$report_id), big.mark = ",")))

# --- 5. Filter to SUSPECT role only ------------------------------------------

cat("\n=== Product roles in supplement reports ===\n")
supplements |> count(product_type) |> print()

supplements_suspect <- supplements |>
  filter(str_to_upper(product_type) == "SUSPECT")

cat(sprintf("\nSuspect supplement rows: %s\n", format(nrow(supplements_suspect), big.mark = ",")))
cat(sprintf("Suspect supplement unique reports: %s\n",
            format(n_distinct(supplements_suspect$report_id), big.mark = ",")))

# --- 6. Initial data quality checks -----------------------------------------

cat("\n=== Data quality summary ===\n")

# Missing values
missing_pct <- supplements_suspect |>
  summarise(across(everything(), ~ mean(is.na(.) | . == "") * 100)) |>
  pivot_longer(everything(), names_to = "column", values_to = "pct_missing") |>
  arrange(desc(pct_missing))

cat("Missing data percentages:\n")
print(as_tibble(missing_pct), n = Inf)

# Date range
cat("\n=== Date range ===\n")
dates_parsed <- supplements_suspect |>
  mutate(report_date = lubridate::mdy(date_fda_first_received_report)) |>
  filter(!is.na(report_date))

cat(sprintf("Date range: %s to %s\n",
            min(dates_parsed$report_date),
            max(dates_parsed$report_date)))
cat(sprintf("Reports with parseable dates: %s / %s\n",
            format(nrow(dates_parsed), big.mark = ","),
            format(nrow(supplements_suspect), big.mark = ",")))

# Gender distribution
cat("\n=== Gender distribution ===\n")
supplements_suspect |> count(sex, sort = TRUE) |> print()

# Age distribution
cat("\n=== Age summary ===\n")
age_data <- supplements_suspect |>
  filter(!is.na(patient_age) & patient_age != "") |>
  mutate(age = as.numeric(patient_age))

cat(sprintf("Reports with age data: %s\n", format(nrow(age_data), big.mark = ",")))
if (nrow(age_data) > 0) {
  cat(sprintf("Age range: %s to %s\n", min(age_data$age, na.rm = TRUE),
              max(age_data$age, na.rm = TRUE)))
  cat(sprintf("Median age: %s\n", median(age_data$age, na.rm = TRUE)))
}

# Age units
cat("\n=== Age units ===\n")
supplements_suspect |> count(age_units, sort = TRUE) |> print()

# Outcomes
cat("\n=== Outcome categories ===\n")
outcomes_split <- supplements_suspect |>
  filter(!is.na(case_outcome) & case_outcome != "") |>
  separate_rows(case_outcome, sep = ",\\s*") |>
  count(case_outcome, sort = TRUE)
print(as_tibble(outcomes_split), n = 20)

# Top MedDRA terms
cat("\n=== Top 20 MedDRA Preferred Terms ===\n")
pt_split <- supplements_suspect |>
  filter(!is.na(case_meddra_preferred_terms) & case_meddra_preferred_terms != "") |>
  separate_rows(case_meddra_preferred_terms, sep = ",\\s*") |>
  count(case_meddra_preferred_terms, sort = TRUE)
print(as_tibble(pt_split |> head(20)), n = 20)
cat(sprintf("\nTotal unique MedDRA PTs: %s\n",
            format(n_distinct(pt_split$case_meddra_preferred_terms), big.mark = ",")))

# --- 7. Save filtered data ---------------------------------------------------

cat("\n=== Saving filtered supplement data ===\n")

fwrite(supplements_suspect,
       file.path(proc_dir, "supplements_suspect_raw.csv"))
cat(sprintf("Saved: %s (%s rows)\n",
            file.path(proc_dir, "supplements_suspect_raw.csv"),
            format(nrow(supplements_suspect), big.mark = ",")))

fwrite(supplements,
       file.path(proc_dir, "supplements_all_roles.csv"))
cat(sprintf("Saved: %s (%s rows)\n",
            file.path(proc_dir, "supplements_all_roles.csv"),
            format(nrow(supplements), big.mark = ",")))

# --- 8. Reporting volume over time -------------------------------------------

cat("\n=== Reporting volume by year ===\n")

annual_counts <- supplements_suspect |>
  mutate(report_date = lubridate::mdy(date_fda_first_received_report),
         year = lubridate::year(report_date)) |>
  filter(!is.na(year)) |>
  count(year, name = "n_reports") |>
  arrange(year)

print(as_tibble(annual_counts), n = Inf)
fwrite(annual_counts, file.path(proc_dir, "annual_report_counts.csv"))

# --- 9. Load case-based view for cross-reference ----------------------------

cat("\n=== Loading case-based view ===\n")
case_raw <- fread(case_file, encoding = "UTF-8") |> clean_names()

cat(sprintf("Total rows in case-based view: %s\n", format(nrow(case_raw), big.mark = ",")))
cat("Column names:\n")
cat(paste(" -", names(case_raw)), sep = "\n")

# Cross-check: how many case-based rows mention supplement products?
supp_cases <- case_raw |>
  filter(str_detect(suspect_products, fixed("54 Vit/Min/Prot")))
cat(sprintf("\nCase-based rows with supplement suspect products: %s\n",
            format(nrow(supp_cases), big.mark = ",")))

# --- Summary -----------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("PHASE 1a COMPLETE: Data Acquisition\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Total CAERS rows:             %s\n", format(nrow(product_raw), big.mark = ",")))
cat(sprintf("  Supplement rows (all roles):   %s\n", format(nrow(supplements), big.mark = ",")))
cat(sprintf("  Supplement rows (suspect):     %s\n", format(nrow(supplements_suspect), big.mark = ",")))
cat(sprintf("  Unique suspect report IDs:     %s\n",
            format(n_distinct(supplements_suspect$report_id), big.mark = ",")))
cat(sprintf("  Unique MedDRA PTs:             %s\n",
            format(n_distinct(pt_split$case_meddra_preferred_terms), big.mark = ",")))
cat(strrep("=", 60), "\n")
cat("\nNext: Run 02_data_cleaning.R\n")
