# =============================================================================
# 08_validation.R
# Validation and sensitivity analysis
# =============================================================================
#
# Input:  data/processed/symptoms_long.csv
#         data/processed/disproportionality_results.csv
#         outputs/tables/robust_signals.csv
#
# Output: outputs/tables/validated_signal_catalogue.csv
#         outputs/tables/sensitivity_mandatory_era.csv
#         outputs/tables/sensitivity_threshold_variation.csv
#         outputs/tables/confounding_by_indication_flags.csv
#         outputs/tables/product_name_consolidation.csv
#
# Sections:
#   1. Exemption 4 exclusion and signal recount
#   2. Sensitivity analysis — mandatory reporting era only (post-2006)
#   3. Sensitivity analysis — threshold variation
#   4. Quantitative bias analysis for underreporting
#   5. Confounding-by-indication flagging
#   6. Product name consolidation audit
#   7. Final validated signal catalogue
# =============================================================================

library(tidyverse)
library(data.table)

# --- Configuration -----------------------------------------------------------

proc_dir  <- "data/processed"
table_dir <- "outputs/tables"

# Thresholds (same as script 03)
PRR_THRESHOLD  <- 2
CHI2_THRESHOLD <- 4
MIN_N          <- 3
GPS_EB05       <- 2
BCPNN_IC025    <- 0

# Underreporting rate estimate (Timbo et al. 2018)
REPORTING_RATE <- 0.02

# --- Helper: run full disproportionality on a data subset --------------------

run_disproportionality <- function(symptoms_data, label = "full") {
  N_total <- n_distinct(symptoms_data$report_id)
  if (N_total < 100) {
    cat(sprintf("  [%s] Too few reports (%d), skipping\n", label, N_total))
    return(tibble())
  }

  product_totals <- symptoms_data |>
    distinct(report_id, product_clean) |>
    count(product_clean, name = "n_product")

  pt_totals <- symptoms_data |>
    distinct(report_id, meddra_pt) |>
    count(meddra_pt, name = "n_pt")

  pair_counts <- symptoms_data |>
    distinct(report_id, product_clean, meddra_pt) |>
    count(product_clean, meddra_pt, name = "a") |>
    filter(a >= MIN_N)

  if (nrow(pair_counts) == 0) return(tibble())

  ct <- pair_counts |>
    left_join(product_totals, by = "product_clean") |>
    left_join(pt_totals, by = "meddra_pt") |>
    mutate(
      a = as.double(a),
      n_product = as.double(n_product),
      n_pt = as.double(n_pt),
      b = n_product - a,
      c = n_pt - a,
      d = N_total - a - b - c,
      N = as.double(N_total)
    )

  # PRR
  ct <- ct |>
    mutate(
      prr = (a / (a + b)) / (c / (c + d)),
      chi2 = ((abs(a * d - b * c) - N / 2)^2 * N) /
             ((a + b) * (c + d) * (a + c) * (b + d)),
      prr_signal = prr >= PRR_THRESHOLD & chi2 >= CHI2_THRESHOLD & a >= MIN_N
    )

  # ROR (0.5 continuity correction)
  ct <- ct |>
    mutate(
      ror = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)),
      ror_se = sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5)),
      ror_lower = exp(log(ror) - 1.96 * ror_se),
      ror_signal = ror_lower > 1
    )

  # GPS
  ct <- ct |> mutate(expected = (n_product * n_pt) / N)

  mean_rr <- mean(ct$a / ct$expected, na.rm = TRUE)
  var_rr  <- var(ct$a / ct$expected, na.rm = TRUE)
  alpha_prior <- mean_rr^2 / var_rr
  beta_prior  <- mean_rr / var_rr

  ct <- ct |>
    mutate(
      alpha_post = alpha_prior + a,
      beta_post  = beta_prior + expected,
      ebgm = log2(alpha_post / beta_post),
      eb05 = log2(qgamma(0.05, shape = alpha_post, rate = beta_post)),
      gps_signal = 2^eb05 >= GPS_EB05
    )

  # BCPNN
  ct <- ct |>
    mutate(
      ic = log2((a + 0.5) / (expected + 0.5)),
      ic_var = 1 / (a + 0.5) - 1 / N,
      ic025 = ic - 1.96 * sqrt(pmax(ic_var, 0)),
      bcpnn_signal = ic025 > BCPNN_IC025
    )

  # Concordance
  ct <- ct |>
    mutate(
      n_methods = as.integer(prr_signal) + as.integer(ror_signal) +
                  as.integer(gps_signal) + as.integer(bcpnn_signal),
      robust_signal = n_methods >= 3
    )

  ct
}

# --- 1. Load data ------------------------------------------------------------

cat("=== Loading data ===\n")
symptoms <- fread(file.path(proc_dir, "symptoms_long.csv")) |> as_tibble()
robust   <- fread(file.path(table_dir, "robust_signals.csv")) |> as_tibble()
disp     <- fread(file.path(proc_dir, "disproportionality_results.csv")) |> as_tibble()

cat(sprintf("Symptoms rows:     %s\n", format(nrow(symptoms), big.mark = ",")))
cat(sprintf("Robust signals:    %s\n", format(nrow(robust), big.mark = ",")))
cat(sprintf("Full disprop rows: %s\n", format(nrow(disp), big.mark = ",")))


# =============================================================================
# SECTION 1: Exemption 4 exclusion
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 1: Exemption 4 Exclusion\n")
cat(strrep("=", 60), "\n")

# Exemption 4 is an FDA FOIA redaction category, not a real product.
# Reports coded as "EXEMPTION 4" represent products whose identity was
# withheld under FOIA exemption 4 (trade secrets / confidential commercial
# information). These must be excluded from signal detection.

n_ex4_reports <- symptoms |>
  filter(product_clean == "EXEMPTION 4") |>
  pull(report_id) |>
  n_distinct()

n_ex4_signals <- robust |>
  filter(product_clean == "EXEMPTION 4") |>
  nrow()

cat(sprintf("Exemption 4 reports:        %d\n", n_ex4_reports))
cat(sprintf("Exemption 4 robust signals: %d\n", n_ex4_signals))

# Remove Exemption 4 from symptoms and recompute
symptoms_clean <- symptoms |>
  filter(product_clean != "EXEMPTION 4")

robust_clean <- robust |>
  filter(product_clean != "EXEMPTION 4")

disp_clean <- disp |>
  filter(product_clean != "EXEMPTION 4")

cat(sprintf("\nAfter exclusion:\n"))
cat(sprintf("  Symptom rows:     %s -> %s\n",
            format(nrow(symptoms), big.mark = ","),
            format(nrow(symptoms_clean), big.mark = ",")))
cat(sprintf("  Robust signals:   %d -> %d (-%d)\n",
            nrow(robust), nrow(robust_clean), n_ex4_signals))
cat(sprintf("  Unique reports:   %s -> %s\n",
            format(n_distinct(symptoms$report_id), big.mark = ","),
            format(n_distinct(symptoms_clean$report_id), big.mark = ",")))


# =============================================================================
# SECTION 2: Sensitivity analysis — mandatory reporting era (post-2006)
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 2: Sensitivity Analysis — Mandatory Era Only\n")
cat(strrep("=", 60), "\n")

# The Dietary Supplements and Nonprescription Drug Consumer Protection Act
# (2006) made serious AE reporting mandatory for manufacturers from 2007.
# Pre-2007 reports are voluntary only and may have different reporting patterns.

symptoms_mandatory <- symptoms_clean |>
  filter(report_year >= 2007)

cat(sprintf("Reports in mandatory era (>=2007): %s / %s (%.1f%%)\n",
            format(n_distinct(symptoms_mandatory$report_id), big.mark = ","),
            format(n_distinct(symptoms_clean$report_id), big.mark = ","),
            n_distinct(symptoms_mandatory$report_id) /
              n_distinct(symptoms_clean$report_id) * 100))

cat("\nRunning disproportionality on mandatory-era subset...\n")
mandatory_results <- run_disproportionality(symptoms_mandatory, "mandatory")

n_mandatory_robust <- sum(mandatory_results$robust_signal, na.rm = TRUE)
cat(sprintf("Mandatory-era robust signals: %d (vs %d full dataset)\n",
            n_mandatory_robust, nrow(robust_clean)))

# Compare: how many of the full-dataset robust signals are also detected
# in the mandatory-era analysis?
robust_pairs <- robust_clean |>
  select(product_clean, meddra_pt) |>
  mutate(in_full = TRUE)

mandatory_robust_pairs <- mandatory_results |>
  filter(robust_signal) |>
  select(product_clean, meddra_pt) |>
  mutate(in_mandatory = TRUE)

concordance_era <- robust_pairs |>
  full_join(mandatory_robust_pairs, by = c("product_clean", "meddra_pt")) |>
  replace_na(list(in_full = FALSE, in_mandatory = FALSE))

n_both     <- sum(concordance_era$in_full & concordance_era$in_mandatory)
n_full_only <- sum(concordance_era$in_full & !concordance_era$in_mandatory)
n_mand_only <- sum(!concordance_era$in_full & concordance_era$in_mandatory)

cat(sprintf("\nMandatory era concordance:\n"))
cat(sprintf("  In both:              %d\n", n_both))
cat(sprintf("  Full dataset only:    %d\n", n_full_only))
cat(sprintf("  Mandatory era only:   %d\n", n_mand_only))
cat(sprintf("  Concordance rate:     %.1f%% of full signals retained\n",
            n_both / nrow(robust_clean) * 100))

# Save mandatory era results
sensitivity_mandatory <- tibble(
  analysis = "Mandatory era (2007+)",
  n_reports = n_distinct(symptoms_mandatory$report_id),
  n_pairs_tested = nrow(mandatory_results),
  n_robust_signals = n_mandatory_robust,
  n_concordant_with_full = n_both,
  pct_concordance = round(n_both / nrow(robust_clean) * 100, 1)
)

# Signals lost when restricting to mandatory era
lost_in_mandatory <- concordance_era |>
  filter(in_full & !in_mandatory) |>
  left_join(robust_clean, by = c("product_clean", "meddra_pt")) |>
  arrange(desc(a))

cat(sprintf("\nTop 20 signals lost when restricting to mandatory era:\n"))
lost_in_mandatory |>
  select(product_clean, meddra_pt, a, prr, n_methods) |>
  head(20) |>
  print(n = 20, width = Inf)


# =============================================================================
# SECTION 3: Sensitivity analysis — threshold variation
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 3: Sensitivity Analysis — Threshold Variation\n")
cat(strrep("=", 60), "\n")

# Test how robust signal counts change under stricter thresholds
# Applied to the Exemption-4-excluded dataset

threshold_scenarios <- tribble(
  ~scenario,             ~prr_thresh, ~chi2_thresh, ~gps_eb05, ~bcpnn_ic025, ~min_n,
  "Baseline",            2,           4,            2,         0,            3,
  "Strict PRR (>=3)",    3,           4,            2,         0,            3,
  "Strict GPS (>=2.5)",  2,           4,            2.5,       0,            3,
  "Strict GPS (>=3)",    2,           4,            3,         0,            3,
  "Min N>=5",            2,           4,            2,         0,            5,
  "Min N>=10",           2,           4,            2,         0,            10,
  "All strict",          3,           4,            2.5,       0,            5
)

threshold_results <- threshold_scenarios |>
  rowwise() |>
  mutate(
    data = list({
      d <- disp_clean |>
        filter(a >= min_n) |>
        mutate(
          prr_sig = prr >= prr_thresh & chi2 >= chi2_thresh & a >= min_n,
          ror_sig = ror_lower > 1,
          gps_sig = 2^eb05 >= gps_eb05,
          bcpnn_sig = ic025 > bcpnn_ic025,
          n_meth = as.integer(prr_sig) + as.integer(ror_sig) +
                   as.integer(gps_sig) + as.integer(bcpnn_sig),
          robust = n_meth >= 3
        )
      tibble(
        n_eligible = nrow(d),
        n_prr = sum(d$prr_sig, na.rm = TRUE),
        n_ror = sum(d$ror_sig, na.rm = TRUE),
        n_gps = sum(d$gps_sig, na.rm = TRUE),
        n_bcpnn = sum(d$bcpnn_sig, na.rm = TRUE),
        n_robust = sum(d$robust, na.rm = TRUE),
        n_all4 = sum(d$n_meth == 4, na.rm = TRUE)
      )
    })
  ) |>
  unnest(data) |>
  ungroup()

cat("Threshold sensitivity results:\n")
threshold_results |>
  select(scenario, n_eligible, n_robust, n_all4) |>
  print(n = Inf, width = Inf)

fwrite(threshold_results,
       file.path(table_dir, "sensitivity_threshold_variation.csv"))
cat("Saved: sensitivity_threshold_variation.csv\n")


# =============================================================================
# SECTION 4: Quantitative bias analysis — underreporting
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 4: Quantitative Bias Analysis — Underreporting\n")
cat(strrep("=", 60), "\n")

# Timbo et al. (2018) estimated ~2% reporting rate for dietary supplement AEs.
# Disproportionality ratios (PRR, ROR, etc.) are RELATIVE measures and remain
# valid under uniform underreporting. However, if reporting rates differ by
# product or outcome severity, this introduces differential reporting bias.
#
# We quantify the potential true burden and assess whether differential
# reporting by severity could affect our findings.

n_unique_reports <- n_distinct(symptoms_clean$report_id)
n_serious_reports <- symptoms_clean |>
  distinct(report_id, .keep_all = TRUE) |>
  filter(outcome_serious == TRUE) |>
  nrow()

estimated_true_total <- n_unique_reports / REPORTING_RATE
estimated_true_serious <- n_serious_reports / REPORTING_RATE

cat(sprintf("Observed reports:           %s\n",
            format(n_unique_reports, big.mark = ",")))
cat(sprintf("Estimated reporting rate:   %.0f%%\n", REPORTING_RATE * 100))
cat(sprintf("Estimated true AE cases:    %s\n",
            format(round(estimated_true_total), big.mark = ",")))
cat(sprintf("\nObserved serious reports:    %s\n",
            format(n_serious_reports, big.mark = ",")))
cat(sprintf("Estimated true serious:     %s\n",
            format(round(estimated_true_serious), big.mark = ",")))

# Differential reporting sensitivity: what if serious events are reported
# at 5x the rate of non-serious events?
# If overall rate = 2%, and p_serious = 0.344 (observed proportion),
# then: 0.02 = p_serious * r_serious + (1 - p_serious) * r_nonserious
# If r_serious = 5 * r_nonserious:
# 0.02 = 0.344 * 5r + 0.656 * r = (1.72 + 0.656) * r = 2.376r
# r_nonserious = 0.02 / 2.376 = 0.0084, r_serious = 0.042

p_serious_obs <- n_serious_reports / n_unique_reports
differential_factor <- 5

r_nonserious <- REPORTING_RATE / (p_serious_obs * differential_factor +
                                    (1 - p_serious_obs))
r_serious <- differential_factor * r_nonserious

cat(sprintf("\nDifferential reporting sensitivity (serious = %dx non-serious):\n",
            differential_factor))
cat(sprintf("  Non-serious reporting rate: %.2f%%\n", r_nonserious * 100))
cat(sprintf("  Serious reporting rate:     %.2f%%\n", r_serious * 100))

# Under differential reporting, the true serious proportion would differ
true_n_serious <- n_serious_reports / r_serious
true_n_nonserious <- (n_unique_reports - n_serious_reports) / r_nonserious
true_p_serious <- true_n_serious / (true_n_serious + true_n_nonserious)

cat(sprintf("  True serious proportion:   %.1f%% (vs %.1f%% observed)\n",
            true_p_serious * 100, p_serious_obs * 100))
cat(sprintf("  Bias magnitude:            %.1f percentage points\n",
            abs(p_serious_obs - true_p_serious) * 100))

# Impact on disproportionality: PRR/ROR are ratios of proportions and remain
# valid if underreporting is non-differential (same rate for all product-PT
# combinations). We document this assumption.
cat("\nNote: Disproportionality ratios (PRR, ROR, GPS, BCPNN) are valid under\n")
cat("non-differential underreporting. The 2% rate does not affect relative\n")
cat("signal strength. Only if reporting rates differ systematically across\n")
cat("products or events (differential misclassification) would ratios be biased.\n")
cat("Well-known products (e.g., Kratom) may have stimulated reporting, which\n")
cat("could inflate absolute counts but is partially mitigated by GPS shrinkage.\n")

bias_summary <- tibble(
  metric = c("Observed reports", "Reporting rate", "Estimated true cases",
             "Observed serious", "Estimated true serious",
             "Observed serious %", "Adjusted serious % (5x differential)",
             "Bias magnitude (pp)"),
  value = c(n_unique_reports, paste0(REPORTING_RATE * 100, "%"),
            round(estimated_true_total),
            n_serious_reports, round(estimated_true_serious),
            paste0(round(p_serious_obs * 100, 1), "%"),
            paste0(round(true_p_serious * 100, 1), "%"),
            round(abs(p_serious_obs - true_p_serious) * 100, 1))
)

cat("\nBias analysis summary:\n")
print(bias_summary, n = Inf, width = Inf)


# =============================================================================
# SECTION 5: Confounding-by-indication flagging
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 5: Confounding-by-Indication Flagging\n")
cat(strrep("=", 60), "\n")

# Protopathic bias / confounding by indication: the adverse event reported
# may be the condition for which the supplement was taken, not a side effect.
# We flag signals where the product's expected indication matches the PT.

indication_patterns <- tribble(
  ~product_pattern,          ~pt_pattern,                          ~reason,
  "PRESERVISION|AREDS",      "Macular|Blindness|Visual|Vision",   "Eye supplement taken for macular degeneration",
  "BENEFIBER|METAMUCIL|FIBER","Diarrhoea|Constipation|Bowel|Faecal|Flatulence|Abdominal|Bloating",
                                                                   "Fibre supplement taken for bowel conditions",
  "MELATONIN|SLEEP",         "Insomnia|Sleep|Somnolence",          "Sleep supplement taken for insomnia",
  "GLUCOSAMINE|CHONDROITIN|JOINT", "Arthralgia|Arthritis|Joint",  "Joint supplement taken for arthritis",
  "PROBIOTIC|FLORASTOR|CULTURELLE", "Diarrhoea|Clostridium|Colitis",
                                                                   "Probiotic taken for GI conditions",
  "WEIGHT|DIET|SLIM|HYDROXYCUT", "Weight|Obesity|Overweight",     "Weight loss product taken for obesity",
  "IRON|FEOSOL|FERROUS",    "Anaemia|Iron Deficiency",             "Iron supplement taken for anaemia",
  "CALCIUM|CALTRATE|CITRACAL", "Osteo|Fracture|Bone",             "Calcium supplement taken for osteoporosis",
  "CRANBERRY",               "Urinary|Cystitis|Urine",             "Cranberry taken for urinary conditions",
  "PROSTATE|BETA PROSTATE",  "Prostat|Urinary|Micturition|Nocturia|Urine",
                                                                   "Prostate supplement taken for BPH",
  "GINKGO|PREVAGEN",        "Memory|Dementia|Cognit|Alzheimer",    "Cognitive supplement taken for memory loss",
  "FISH OIL|OMEGA|LOVAZA",  "Hyperlipid|Cholesterol|Triglyceride", "Fish oil taken for dyslipidaemia",
  "ST.? JOHN|HYPERICUM",    "Depress|Anxiety|Mood",                "St John's wort taken for depression",
  "VALERIAN",               "Insomnia|Anxiety|Sleep",              "Valerian taken for anxiety/insomnia",
  "SAW PALMETTO",           "Prostat|Urinary|Micturition",         "Saw palmetto taken for BPH"
)

# Apply flagging to robust signals
cbi_flags <- robust_clean |>
  mutate(cbi_flag = FALSE, cbi_reason = NA_character_)

for (i in seq_len(nrow(indication_patterns))) {
  prod_pat <- indication_patterns$product_pattern[i]
  pt_pat   <- indication_patterns$pt_pattern[i]
  reason   <- indication_patterns$reason[i]

  matches <- str_detect(cbi_flags$product_clean, regex(prod_pat, ignore_case = TRUE)) &
             str_detect(cbi_flags$meddra_pt, regex(pt_pat, ignore_case = TRUE))

  cbi_flags$cbi_flag[matches] <- TRUE
  cbi_flags$cbi_reason[matches] <- reason
}

cbi_flagged <- cbi_flags |>
  filter(cbi_flag) |>
  select(product_clean, meddra_pt, a, prr, n_methods, cbi_reason) |>
  arrange(desc(a))

cat(sprintf("Signals flagged for confounding by indication: %d / %d (%.1f%%)\n",
            nrow(cbi_flagged), nrow(robust_clean),
            nrow(cbi_flagged) / nrow(robust_clean) * 100))

cat("\nFlagged signals:\n")
print(cbi_flagged, n = Inf, width = Inf)

fwrite(cbi_flagged, file.path(table_dir, "confounding_by_indication_flags.csv"))
cat("Saved: confounding_by_indication_flags.csv\n")


# =============================================================================
# SECTION 6: Product name consolidation audit
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 6: Product Name Consolidation Audit\n")
cat(strrep("=", 60), "\n")

# Many products appear under multiple variant names in CAERS.
# Identify clusters of likely-duplicate product names that could be merged.

robust_products <- robust_clean |>
  distinct(product_clean) |>
  pull(product_clean) |>
  sort()

cat(sprintf("Unique products in robust signals: %d\n", length(robust_products)))

# Approach: extract a short "base name" and group products with the same base
# This catches variants like "HYDROXYCUT REGULAR RAPID RELEASE CAPLETS" vs
# "HYDROXYCUT HARDCORE"
consolidation <- tibble(product_clean = robust_products) |>
  mutate(
    # Extract first 2 words as base name
    base_name = str_extract(product_clean, "^\\S+\\s+\\S+") |>
                  str_replace_all("[^A-Z0-9 ]", "") |>
                  str_squish(),
    # Fallback for single-word products
    base_name = if_else(is.na(base_name), product_clean, base_name)
  )

# Count products per base name
base_groups <- consolidation |>
  count(base_name, name = "n_variants") |>
  filter(n_variants > 1) |>
  arrange(desc(n_variants))

cat(sprintf("Product name clusters (>1 variant): %d\n", nrow(base_groups)))

# Show the groups and their signal counts
consolidation_detail <- consolidation |>
  inner_join(base_groups |> select(base_name), by = "base_name") |>
  left_join(
    robust_clean |>
      count(product_clean, name = "n_signals"),
    by = "product_clean"
  ) |>
  left_join(
    robust_clean |>
      group_by(product_clean) |>
      summarise(total_cases = sum(a), .groups = "drop"),
    by = "product_clean"
  ) |>
  arrange(base_name, desc(n_signals))

cat("\nTop 30 consolidation clusters (by total variants):\n")
consolidation_detail |>
  group_by(base_name) |>
  summarise(
    n_variants = n(),
    variants = paste(product_clean, collapse = " | "),
    total_signals = sum(n_signals, na.rm = TRUE),
    total_cases = sum(total_cases, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_variants)) |>
  head(30) |>
  print(n = 30, width = Inf)

fwrite(consolidation_detail,
       file.path(table_dir, "product_name_consolidation.csv"))
cat("Saved: product_name_consolidation.csv\n")

# Quantify impact: how many signals would be affected by consolidation?
n_products_in_clusters <- nrow(consolidation_detail)
n_signals_in_clusters <- sum(consolidation_detail$n_signals, na.rm = TRUE)

cat(sprintf("\nProducts in multi-variant clusters: %d / %d (%.1f%%)\n",
            n_products_in_clusters, length(robust_products),
            n_products_in_clusters / length(robust_products) * 100))
cat(sprintf("Signals from clustered products:   %d / %d (%.1f%%)\n",
            n_signals_in_clusters, nrow(robust_clean),
            n_signals_in_clusters / nrow(robust_clean) * 100))


# =============================================================================
# SECTION 7: Final validated signal catalogue
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 7: Final Validated Signal Catalogue\n")
cat(strrep("=", 60), "\n")

# Annotate the robust signals with validation metadata

# Validation status from manual review of top products
# (Based on cross-referencing with FDA actions, published literature,
#  and pharmacological mechanisms)
validated_products <- tribble(
  ~product_pattern,             ~validation_status,  ~validation_note,
  "^KRATOM",                    "Validated",         "FDA import alert 54-15, DEA scheduling, extensive literature on mitragynine toxicity",
  "^HYDROXYCUT",                "Validated",         "2009 FDA recall, hepatotoxicity well-documented (Fong et al. 2010)",
  "^SUPER BETA PROSTATE",      "Validated",         "FDA warning letters, urologic symptoms consistent with beta-sitosterol effects",
  "^CENTRUM SILVER",           "Validated",         "Formulation-related choking/dysphagia, large tablet size, elderly population",
  "^CITRACAL",                 "Validated",         "Formulation-related choking, large calcium tablet size",
  "^HERBALIFE",                "Validated",         "2008-2017 FDA/regulatory scrutiny, hepatotoxicity case series",
  "^OXY ELITE PRO|^OXYELITE",  "Validated",         "2013 FDA recall, acute hepatitis outbreak (Hawaii cluster)",
  "^LIPODRENE",               "Validated",         "Sympathomimetic amine content, cardiovascular adverse events",
  "^PRESERVISION|^AREDS",     "Confounding",       "Confounding by indication: taken by AMD patients, reports AMD progression",
  "^BENEFIBER",               "Confounding",       "Confounding by indication: taken for bowel conditions, reports bowel symptoms",
  "^AG1|^ATHLETIC GREENS",    "Emerging",          "Hepatic enzyme signals emerging 2024; pre-regulatory detection",
  "^NUTRAFOL",                "Emerging",          "Hepatic enzyme signals emerging 2025; pre-regulatory detection",
  "^FLORASTOR",               "Partially validated","Saccharomyces boulardii fungaemia in immunocompromised; known rare risk",
  "^GRIPE WATER",             "Validated",         "Choking hazard in infants from liquid formulation; FDA consumer advisory",
  "^VITAMIN B6|^NATUREMADE VITAMIN B6", "Stimulated reporting",
                                                    "Well-known toxicity; 2023-2024 reporting spike driven by TGA/EFSA actions (Weber effect)"
)

# Apply validation annotations
validated_catalogue <- robust_clean |>
  mutate(
    validation_status = NA_character_,
    validation_note = NA_character_
  )

for (i in seq_len(nrow(validated_products))) {
  matches <- str_detect(validated_catalogue$product_clean,
                         regex(validated_products$product_pattern[i], ignore_case = TRUE))
  validated_catalogue$validation_status[matches & is.na(validated_catalogue$validation_status)] <-
    validated_products$validation_status[i]
  validated_catalogue$validation_note[matches & is.na(validated_catalogue$validation_note)] <-
    validated_products$validation_note[i]
}

# Add confounding-by-indication flags
cbi_pairs <- cbi_flagged |>
  select(product_clean, meddra_pt, cbi_reason)

validated_catalogue <- validated_catalogue |>
  left_join(cbi_pairs, by = c("product_clean", "meddra_pt")) |>
  mutate(
    # Override validation status for CBI-flagged pairs not already annotated
    validation_status = case_when(
      !is.na(validation_status) ~ validation_status,
      !is.na(cbi_reason) ~ "Confounding",
      TRUE ~ NA_character_
    ),
    validation_note = case_when(
      !is.na(validation_note) ~ validation_note,
      !is.na(cbi_reason) ~ cbi_reason,
      TRUE ~ NA_character_
    )
  ) |>
  select(-cbi_reason)

# Summary of validation coverage
validation_summary <- validated_catalogue |>
  count(validation_status, name = "n_signals") |>
  mutate(pct = round(n_signals / nrow(validated_catalogue) * 100, 1)) |>
  arrange(desc(n_signals))

cat("Validation status summary:\n")
print(validation_summary, n = Inf, width = Inf)

# Save final catalogue
fwrite(validated_catalogue,
       file.path(table_dir, "validated_signal_catalogue.csv"))
cat(sprintf("\nSaved: validated_signal_catalogue.csv (%s signals)\n",
            format(nrow(validated_catalogue), big.mark = ",")))

# Save sensitivity summary
sensitivity_summary <- bind_rows(
  sensitivity_mandatory,
  threshold_results |>
    transmute(
      analysis = scenario,
      n_reports = NA_integer_,
      n_pairs_tested = n_eligible,
      n_robust_signals = n_robust,
      n_concordant_with_full = NA_integer_,
      pct_concordance = NA_real_
    )
)

fwrite(sensitivity_summary,
       file.path(table_dir, "sensitivity_mandatory_era.csv"))
cat("Saved: sensitivity_mandatory_era.csv\n")


# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("PHASE 6 COMPLETE: Validation and Sensitivity Analysis\n")
cat(strrep("=", 60), "\n")

cat(sprintf("  Exemption 4 excluded:       %d signals removed\n", n_ex4_signals))
cat(sprintf("  Clean robust signals:       %d\n", nrow(robust_clean)))
cat(sprintf("  Mandatory-era concordance:  %.1f%%\n",
            n_both / nrow(robust_clean) * 100))
cat(sprintf("  CBI-flagged signals:        %d\n", nrow(cbi_flagged)))
cat(sprintf("  Product name clusters:      %d\n", nrow(base_groups)))

cat("\n  Validation coverage:\n")
validation_summary |>
  mutate(line = sprintf("    %-25s %d (%.1f%%)", validation_status, n_signals, pct)) |>
  pull(line) |>
  cat(sep = "\n")

cat(sprintf("\n    Unannotated:              %d (%.1f%%)\n",
            sum(is.na(validated_catalogue$validation_status)),
            sum(is.na(validated_catalogue$validation_status)) /
              nrow(validated_catalogue) * 100))

cat("\n  Threshold sensitivity:\n")
threshold_results |>
  mutate(line = sprintf("    %-25s %d robust signals", scenario, n_robust)) |>
  pull(line) |>
  cat(sep = "\n")

cat("\n")
cat(strrep("=", 60), "\n")
cat("Next: Run 09_visualisation.R to update figures, then begin manuscript\n")
cat(strrep("=", 60), "\n")
