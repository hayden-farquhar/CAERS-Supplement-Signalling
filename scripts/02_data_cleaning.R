# =============================================================================
# 02_data_cleaning.R
# Standardise product names, classify categories, clean MedDRA terms,
# construct analytic variables, and deduplicate
# =============================================================================
#
# Input:  data/processed/supplements_suspect_raw.csv (from 01_data_acquisition.R)
# Output: data/processed/supplements_cleaned.csv     (product-level)
#         data/processed/symptoms_long.csv            (product-symptom pairs)
#         data/processed/case_level_summary.csv       (one row per case)
# =============================================================================

library(tidyverse)
library(data.table)
library(janitor)

# --- Configuration -----------------------------------------------------------

proc_dir <- "data/processed"
map_dir  <- "data/mappings"

# --- 1. Load filtered supplement data ----------------------------------------

cat("=== Loading filtered supplement data ===\n")
supps <- fread(file.path(proc_dir, "supplements_suspect_raw.csv"),
               encoding = "UTF-8") |>
  as_tibble()
cat(sprintf("Loaded %s rows\n", format(nrow(supps), big.mark = ",")))

# --- 2. Parse dates ----------------------------------------------------------

cat("\n=== Parsing dates ===\n")

supps <- supps |>
  mutate(
    report_date    = lubridate::mdy(date_fda_first_received_report),
    event_date_raw = lubridate::mdy(date_event),
    report_year    = lubridate::year(report_date),
    report_quarter = lubridate::quarter(report_date),
    report_yq      = paste0(report_year, "-Q", report_quarter)
  )

cat(sprintf("Report date range: %s to %s\n",
            min(supps$report_date, na.rm = TRUE),
            max(supps$report_date, na.rm = TRUE)))

# --- 3. Deduplicate ----------------------------------------------------------

cat("\n=== Deduplication ===\n")
n_before <- nrow(supps)

# A report_id can appear multiple times if multiple suspect products per case.
# Remove truly duplicated rows (same report + same product + same symptoms).
supps <- supps |>
  distinct(report_id, product, case_meddra_preferred_terms, .keep_all = TRUE)

n_after <- nrow(supps)
cat(sprintf("Rows before dedup: %s\n", format(n_before, big.mark = ",")))
cat(sprintf("Rows after dedup:  %s\n", format(n_after, big.mark = ",")))
cat(sprintf("Removed %s duplicate rows\n", format(n_before - n_after, big.mark = ",")))

# --- 4. Standardise product names --------------------------------------------

cat("\n=== Standardising product names ===\n")

supps <- supps |>
  mutate(
    product_clean = product |>
      str_to_upper() |>
      str_squish() |>
      str_remove_all("[®™©]") |>
      str_remove_all("\\(DIETARY SUPPLEMENT\\)") |>
      str_remove_all("\\(DS\\)") |>
      str_remove_all("^BRAND:\\s*") |>
      str_replace_all("\\s+", " ") |>
      str_trim()
  )

n_products_raw   <- n_distinct(supps$product)
n_products_clean <- n_distinct(supps$product_clean)
cat(sprintf("Unique raw product names:     %s\n", format(n_products_raw, big.mark = ",")))
cat(sprintf("Unique cleaned product names: %s\n", format(n_products_clean, big.mark = ",")))

cat("\nTop 30 products by report count:\n")
supps |> count(product_clean, sort = TRUE) |> head(30) |> print(n = 30)

# --- 5. Product category classification --------------------------------------

cat("\n=== Classifying product categories ===\n")

supps <- supps |>
  mutate(
    product_category = case_when(
      str_detect(product_clean,
        "WEIGHT|DIET|SLIM|LEAN|FAT BURN|GARCINIA|HYDROXYCUT|KETO|METABOL") ~
        "Weight Loss/Diet",
      str_detect(product_clean,
        "ENERGY|CAFFEIN|STIMULANT|5[- ]HOUR|RED BULL|MONSTER|BANG|G FUEL") ~
        "Energy/Stimulant",
      str_detect(product_clean,
        "SPORT|PROTEIN|CREATINE|BCAA|WHEY|PRE[- ]?WORK|MUSCLE|GYM|MASS GAIN") ~
        "Sports Nutrition",
      str_detect(product_clean,
        "VITAMIN|VIT |MULTIVIT|PRENATAL|FOLIC|BIOTIN|B[- ]?12|B[- ]?COMPLEX") ~
        "Vitamin/Mineral",
      str_detect(product_clean,
        "MINERAL|CALCIUM|IRON|ZINC|MAGNESIUM|SELENIUM|POTASSIUM|CHROMIUM") ~
        "Vitamin/Mineral",
      str_detect(product_clean,
        "FISH OIL|OMEGA|COD LIVER|FLAX|KRILL") ~
        "Omega/Fish Oil",
      str_detect(product_clean,
        "PROBIOTIC|PREBIOTIC|LACTOBACILL|BIFIDO") ~
        "Probiotic",
      str_detect(product_clean,
        "HERB|BOTANICAL|GINKGO|GINSENG|ECHINACEA|TURMERIC|CURCUM|ST\\.? JOHN|VALERIAN|ASHWAGAND|KRATOM|KAVA|BLACK COHOSH|SAW PALMETTO|MILK THISTLE|GINGER|ELDERBERRY") ~
        "Herbal/Botanical",
      str_detect(product_clean,
        "MELATONIN|SLEEP") ~
        "Sleep Aid",
      str_detect(product_clean,
        "CBD|CANNABID|HEMP") ~
        "CBD/Hemp",
      str_detect(product_clean,
        "SEXUAL|MALE ENHANCE|RHINO|VIGOR|LIBIDO") ~
        "Sexual Enhancement",
      TRUE ~ "Other"
    )
  )

cat("Product category distribution:\n")
supps |>
  count(product_category, sort = TRUE) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  print(n = Inf)

# --- 6. Construct demographic variables --------------------------------------

cat("\n=== Constructing demographic variables ===\n")

# Standardise age to years
supps <- supps |>
  mutate(
    age_numeric = suppressWarnings(as.numeric(patient_age)),
    age_years = case_when(
      is.na(age_numeric) ~ NA_real_,
      str_detect(str_to_lower(age_units), "year")   ~ age_numeric,
      str_detect(str_to_lower(age_units), "month")  ~ age_numeric / 12,
      str_detect(str_to_lower(age_units), "day")    ~ age_numeric / 365.25,
      str_detect(str_to_lower(age_units), "week")   ~ age_numeric / 52.18,
      str_detect(str_to_lower(age_units), "decade") ~ age_numeric * 10,
      is.na(age_units) | age_units == "" ~ age_numeric,  # assume years
      TRUE ~ age_numeric
    ),
    # Filter implausible ages
    age_years = if_else(age_years < 0 | age_years > 120, NA_real_, age_years),
    # Age bins
    age_group = case_when(
      is.na(age_years) ~ "Unknown",
      age_years < 18   ~ "<18",
      age_years < 40   ~ "18-39",
      age_years < 60   ~ "40-59",
      age_years >= 60  ~ "60+",
      TRUE             ~ "Unknown"
    ),
    # Standardise sex
    sex_clean = case_when(
      str_detect(str_to_upper(sex), "^F") ~ "Female",
      str_detect(str_to_upper(sex), "^M") ~ "Male",
      TRUE ~ "Unknown"
    )
  )

cat("Age group distribution:\n")
supps |> count(age_group, sort = TRUE) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |> print(n = Inf)

cat("\nSex distribution:\n")
supps |> count(sex_clean, sort = TRUE) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |> print(n = Inf)

cat("\nAge summary (years, valid only):\n")
valid_ages <- supps |> filter(!is.na(age_years))
cat(sprintf("  N with age: %s (%.1f%%)\n",
            format(nrow(valid_ages), big.mark = ","),
            nrow(valid_ages) / nrow(supps) * 100))
cat(sprintf("  Mean: %.1f, Median: %.0f, IQR: %.0f-%.0f\n",
            mean(valid_ages$age_years),
            median(valid_ages$age_years),
            quantile(valid_ages$age_years, 0.25),
            quantile(valid_ages$age_years, 0.75)))

# --- 7. Parse and classify outcomes ------------------------------------------

cat("\n=== Parsing outcomes ===\n")

supps <- supps |>
  mutate(
    outcome_death           = str_detect(str_to_lower(case_outcome), "death"),
    outcome_life_threat     = str_detect(str_to_lower(case_outcome), "life threaten"),
    outcome_hospitalisation = str_detect(str_to_lower(case_outcome), "hospitali"),
    outcome_disability      = str_detect(str_to_lower(case_outcome), "disabil"),
    outcome_er_visit        = str_detect(str_to_lower(case_outcome), "emergency"),
    outcome_serious = outcome_death | outcome_life_threat |
                      outcome_hospitalisation | outcome_disability
  )

cat("Outcome prevalence:\n")
tibble(
  outcome = c("Death", "Life-threatening", "Hospitalisation",
              "Disability", "ER visit", "Any serious"),
  n = c(sum(supps$outcome_death, na.rm = TRUE),
        sum(supps$outcome_life_threat, na.rm = TRUE),
        sum(supps$outcome_hospitalisation, na.rm = TRUE),
        sum(supps$outcome_disability, na.rm = TRUE),
        sum(supps$outcome_er_visit, na.rm = TRUE),
        sum(supps$outcome_serious, na.rm = TRUE)),
  pct = round(n / nrow(supps) * 100, 1)
) |> print(n = Inf)

# --- 8. Normalise and parse MedDRA symptoms ----------------------------------

cat("\n=== Parsing MedDRA symptoms ===\n")

# MedDRA PTs have inconsistent casing (e.g., "NAUSEA" vs "Nausea")
# Normalise to Title Case
symptoms_long <- supps |>
  select(report_id, product_clean, product_category, case_meddra_preferred_terms,
         age_group, sex_clean, report_year, report_yq, report_date,
         starts_with("outcome_")) |>
  filter(!is.na(case_meddra_preferred_terms) & case_meddra_preferred_terms != "") |>
  separate_rows(case_meddra_preferred_terms, sep = ",\\s*") |>
  rename(meddra_pt = case_meddra_preferred_terms) |>
  mutate(meddra_pt = str_to_title(str_trim(meddra_pt)))

# Remove empty terms
symptoms_long <- symptoms_long |>
  filter(meddra_pt != "" & !is.na(meddra_pt))

cat(sprintf("Total product-symptom pairs: %s\n",
            format(nrow(symptoms_long), big.mark = ",")))
cat(sprintf("Unique MedDRA PTs (after normalisation): %s\n",
            format(n_distinct(symptoms_long$meddra_pt), big.mark = ",")))

cat("\nTop 30 MedDRA Preferred Terms:\n")
symptoms_long |>
  count(meddra_pt, sort = TRUE) |>
  head(30) |>
  print(n = 30)

# --- 9. Approximate MedDRA SOC mapping ---------------------------------------

cat("\n=== Approximate MedDRA SOC mapping ===\n")

# Pattern-based SOC assignment for the most common PTs
# A full MedDRA licence provides the complete hierarchy
soc_patterns <- tribble(
  ~pattern, ~soc,
  "nausea|vomit|diarrh|abdominal|constipat|gastro|dyspepsia|flatulence|bloat|rectal|oesophag|dysphagia|chok",
    "Gastrointestinal Disorders",
  "headache|dizz|syncope|tremor|paraesth|seizure|convuls|somnolen|migraine|hypoaesthe|ataxia|dysgeusia",
    "Nervous System Disorders",
  "rash|prurit|urticaria|erythema|acne|alopecia|dermatit|skin|blister|hyperhidr",
    "Skin and Subcutaneous Tissue Disorders",
  "fatigue|asthenia|malaise|pyrexia|pain$|^pain|chills|oedema|swelling|feeling|chest discomfort|death|sudden",
    "General Disorders",
  "tachycard|palpitat|hypotens|hypertens|chest pain|cardiac|arrhythm|myocard|bradycard",
    "Cardiac Disorders",
  "dyspno|cough|asthma|wheez|respiratory|pulmonary|throat tight|pharyngeal",
    "Respiratory Disorders",
  "hepat|jaundice|liver|biliru|transaminase|alt incr|ast incr|cholest",
    "Hepatobiliary Disorders",
  "renal|kidney|creatinine|oliguria|urin",
    "Renal and Urinary Disorders",
  "anxi|depress|hallucin|agitat|confus|mood|panic|suicid|psycho|insomnia|sleep|nightmare",
    "Psychiatric Disorders",
  "anaphyla|allerg|hypersensit|angioedema",
    "Immune System Disorders",
  "myalgia|arthralg|muscle|musculoskel|back pain|joint|rhabdomyo",
    "Musculoskeletal Disorders",
  "anaemi|thrombocyt|leukocyt|neutrop|blood|haemorrh|coagul",
    "Blood Disorders",
  "vision|eye|blind|diplop|blur",
    "Eye Disorders",
  "weight|appeti|anorexia|body mass",
    "Metabolism and Nutrition Disorders"
)

assign_soc <- function(pts) {
  pt_lower <- str_to_lower(pts)
  soc_result <- rep("Other/Unclassified", length(pt_lower))
  for (i in seq_len(nrow(soc_patterns))) {
    matches <- str_detect(pt_lower, soc_patterns$pattern[i])
    soc_result[matches & soc_result == "Other/Unclassified"] <- soc_patterns$soc[i]
  }
  soc_result
}

symptoms_long <- symptoms_long |>
  mutate(soc = assign_soc(meddra_pt))

cat("SOC distribution:\n")
symptoms_long |>
  count(soc, sort = TRUE) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  print(n = Inf)

# --- 10. Build 2x2 contingency table inputs ----------------------------------

cat("\n=== Building contingency table inputs ===\n")

# For disproportionality analysis, we need:
#   a = cases with product X AND event Y
#   b = cases with product X AND NOT event Y
#   c = cases without product X AND event Y
#   d = cases without product X AND NOT event Y
#
# We'll precompute the product-event pair counts here.
# The full 2x2 table construction happens in 03_disproportionality_analysis.R

# Count product-PT pairs
product_pt_counts <- symptoms_long |>
  count(product_clean, meddra_pt, name = "n_pairs") |>
  arrange(desc(n_pairs))

cat(sprintf("Unique product-PT pairs: %s\n",
            format(nrow(product_pt_counts), big.mark = ",")))
cat(sprintf("Pairs with N >= 3: %s\n",
            format(sum(product_pt_counts$n_pairs >= 3), big.mark = ",")))

cat("\nTop 20 product-PT pairs:\n")
product_pt_counts |> head(20) |> print(n = 20)

# Save pair counts for Phase 2
fwrite(product_pt_counts, file.path(proc_dir, "product_pt_counts.csv"))

# --- 11. Save cleaned datasets -----------------------------------------------

cat("\n=== Saving cleaned datasets ===\n")

# Main cleaned dataset (product-level)
fwrite(supps, file.path(proc_dir, "supplements_cleaned.csv"))
cat(sprintf("Saved: supplements_cleaned.csv (%s rows)\n",
            format(nrow(supps), big.mark = ",")))

# Long-format symptom table
fwrite(symptoms_long, file.path(proc_dir, "symptoms_long.csv"))
cat(sprintf("Saved: symptoms_long.csv (%s rows)\n",
            format(nrow(symptoms_long), big.mark = ",")))

# Case-level summary (one row per report)
case_level <- supps |>
  group_by(report_id) |>
  summarise(
    report_date    = first(report_date),
    report_year    = first(report_year),
    report_yq      = first(report_yq),
    n_suspect_products = n(),
    products       = paste(product_clean, collapse = " | "),
    categories     = paste(unique(product_category), collapse = " | "),
    age_years      = first(age_years),
    age_group      = first(age_group),
    sex            = first(sex_clean),
    outcome_death  = any(outcome_death),
    outcome_life_threat = any(outcome_life_threat),
    outcome_hospitalisation = any(outcome_hospitalisation),
    outcome_serious = any(outcome_serious),
    symptoms       = first(case_meddra_preferred_terms),
    .groups = "drop"
  )

fwrite(case_level, file.path(proc_dir, "case_level_summary.csv"))
cat(sprintf("Saved: case_level_summary.csv (%s cases)\n",
            format(nrow(case_level), big.mark = ",")))

# --- 12. Mandatory vs voluntary reporting sensitivity flag -------------------

cat("\n=== Mandatory reporting flag ===\n")

# Post-2006 reports are subject to mandatory serious AE reporting
supps <- supps |>
  mutate(mandatory_era = report_year >= 2007)

cat(sprintf("Pre-mandatory era (before 2007): %s rows\n",
            format(sum(!supps$mandatory_era, na.rm = TRUE), big.mark = ",")))
cat(sprintf("Mandatory era (2007+):           %s rows\n",
            format(sum(supps$mandatory_era, na.rm = TRUE), big.mark = ",")))

# Resave with mandatory flag
fwrite(supps, file.path(proc_dir, "supplements_cleaned.csv"))

# --- Summary -----------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("PHASE 1 COMPLETE: Data Cleaning\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Cleaned product-level rows:  %s\n", format(nrow(supps), big.mark = ",")))
cat(sprintf("  Unique report IDs:           %s\n", format(nrow(case_level), big.mark = ",")))
cat(sprintf("  Product-symptom pairs:       %s\n", format(nrow(symptoms_long), big.mark = ",")))
cat(sprintf("  Unique MedDRA PTs:           %s\n",
            format(n_distinct(symptoms_long$meddra_pt), big.mark = ",")))
cat(sprintf("  Unique product-PT pairs:     %s\n",
            format(nrow(product_pt_counts), big.mark = ",")))
cat(sprintf("  Date range:                  %s to %s\n",
            min(supps$report_date, na.rm = TRUE),
            max(supps$report_date, na.rm = TRUE)))
cat(sprintf("  Serious outcome rate:        %.1f%%\n",
            mean(case_level$outcome_serious, na.rm = TRUE) * 100))
cat(sprintf("  Product categories:          %s\n",
            paste(sort(unique(supps$product_category)), collapse = ", ")))
cat(strrep("=", 60), "\n")
cat("\nNext: Run 03_disproportionality_analysis.R\n")
