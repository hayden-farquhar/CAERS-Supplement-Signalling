# =============================================================================
# 08b_extended_validation.R
# Extended validation and robustness analyses
# =============================================================================
#
# Sections:
#   1. Multiple testing / FDR (Benjamini-Hochberg)
#   2. Positive and negative control validation
#   3. Masking / competition bias diagnostic
#   4. Time-stratified disproportionality (2007-2015 vs 2016-2025)
#   5. OpenFDA enforcement cross-reference
#   6. Concomitant product characterisation
#   7. SOC-level signal aggregation
#   8. Bootstrap false discovery estimation
#   9. Regression-based signal detection (adjusted ORs)
# =============================================================================

library(tidyverse)
library(data.table)
library(jsonlite)

set.seed(42)

# --- Configuration -----------------------------------------------------------

proc_dir  <- "data/processed"
table_dir <- "outputs/tables"

PRR_THRESHOLD  <- 2
CHI2_THRESHOLD <- 4
MIN_N          <- 3
GPS_EB05       <- 2
BCPNN_IC025    <- 0

N_BOOTSTRAP    <- 50   # permutations for bootstrap FDR
N_REGRESSION   <- 300  # top pairs for regression analysis

# --- Load data ---------------------------------------------------------------

cat("=== Loading data ===\n")
symptoms <- fread(file.path(proc_dir, "symptoms_long.csv")) |> as_tibble()
supps    <- fread(file.path(proc_dir, "supplements_cleaned.csv")) |> as_tibble()
disp     <- fread(file.path(proc_dir, "disproportionality_results.csv")) |> as_tibble()
robust   <- fread(file.path(table_dir, "robust_signals.csv")) |> as_tibble()

# Exclude Exemption 4
symptoms <- symptoms |> filter(product_clean != "EXEMPTION 4")
supps    <- supps |> filter(product_clean != "EXEMPTION 4")
disp     <- disp |> filter(product_clean != "EXEMPTION 4")
robust   <- robust |> filter(product_clean != "EXEMPTION 4")

N_total <- n_distinct(symptoms$report_id)
cat(sprintf("Reports: %s, Eligible pairs: %s, Robust signals: %s\n",
            format(N_total, big.mark = ","),
            format(nrow(disp), big.mark = ","),
            format(nrow(robust), big.mark = ",")))


# =============================================================================
# SECTION 1: Multiple Testing / FDR
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 1: Multiple Testing / FDR (Benjamini-Hochberg)\n")
cat(strrep("=", 60), "\n")

# Compute p-values from the chi-squared statistics (PRR test)
disp <- disp |>
  mutate(
    chi2_p = pchisq(chi2, df = 1, lower.tail = FALSE),
    # BH-adjusted p-values
    chi2_p_bh = p.adjust(chi2_p, method = "BH"),
    # Bonferroni for reference
    chi2_p_bonf = p.adjust(chi2_p, method = "bonferroni"),
    # FDR-significant
    fdr_sig_05 = chi2_p_bh < 0.05,
    fdr_sig_01 = chi2_p_bh < 0.01,
    bonf_sig   = chi2_p_bonf < 0.05
  )

fdr_summary <- tibble(
  criterion = c(
    "Unadjusted p < 0.05",
    "BH-adjusted FDR < 0.05",
    "BH-adjusted FDR < 0.01",
    "Bonferroni p < 0.05",
    "Robust (>=3/4 methods)"
  ),
  n_signals = c(
    sum(disp$chi2_p < 0.05, na.rm = TRUE),
    sum(disp$fdr_sig_05, na.rm = TRUE),
    sum(disp$fdr_sig_01, na.rm = TRUE),
    sum(disp$bonf_sig, na.rm = TRUE),
    sum(disp$robust_signal, na.rm = TRUE)
  )
)

cat("Multiple testing comparison:\n")
print(fdr_summary, n = Inf, width = Inf)

# Compare FDR-significant signals with robust signals
fdr_robust_overlap <- disp |>
  summarise(
    both_fdr05_robust = sum(fdr_sig_05 & robust_signal, na.rm = TRUE),
    fdr05_only        = sum(fdr_sig_05 & !robust_signal, na.rm = TRUE),
    robust_only       = sum(!fdr_sig_05 & robust_signal, na.rm = TRUE),
    neither           = sum(!fdr_sig_05 & !robust_signal, na.rm = TRUE)
  )

cat("\nFDR 0.05 vs Robust signal concordance:\n")
print(fdr_robust_overlap)

# What proportion of robust signals survive FDR correction?
pct_robust_fdr <- fdr_robust_overlap$both_fdr05_robust /
  (fdr_robust_overlap$both_fdr05_robust + fdr_robust_overlap$robust_only) * 100
cat(sprintf("\nRobust signals also FDR<0.05: %.1f%%\n", pct_robust_fdr))

fwrite(fdr_summary, file.path(table_dir, "fdr_multiple_testing.csv"))


# =============================================================================
# SECTION 2: Positive and Negative Control Validation
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 2: Positive and Negative Control Validation\n")
cat(strrep("=", 60), "\n")

# Define positive controls: known causal product-AE relationships
positive_controls <- tribble(
  ~product_pattern,               ~pt_pattern,                        ~evidence,
  "^KRATOM",                      "Death",                            "FDA import alert, DEA, extensive literature",
  "^KRATOM",                      "Drug Dependence|Withdrawal",       "Known addictive properties of mitragynine",
  "^KRATOM",                      "Seizure",                          "Case series: Trakulsrichai et al. 2015",
  "^HYDROXYCUT",                  "Liver|Hepat|Jaundice",             "2009 FDA recall, Fong et al. 2010",
  "^HYDROXYCUT",                  "Rhabdomyolysis",                   "Case reports, sympathomimetic mechanism",
  "^OXY.?ELITE",                  "Hepatitis|Liver|Hepat",            "2013 FDA recall, Hawaii outbreak",
  "^EPHEDRA",                     "Myocardial|Cardiac Arrest|Stroke", "2004 FDA ban, Haller & Benowitz 2000",
  "^LIPODRENE",                   "Tachycardia|Palpitations|Cardiac", "Sympathomimetic amine content",
  "^HYDROXYCUT",                  "Nausea|Vomiting",                  "Common GI adverse effects, well-documented",
  "^KRATOM",                      "Nausea|Vomiting",                  "Common effects at standard doses",
  "^KRATOM",                      "Tachycardia|Palpitations",         "Sympathomimetic effects of mitragynine",
  "^RED YEAST",                   "Rhabdomyolysis|Myalgia",           "Contains monacolin K (lovastatin equivalent)",
  "^YOHIMBE",                     "Tachycardia|Hypertension|Anxiety", "Alpha-2 antagonist mechanism"
)

# Define negative controls: product-AE pairs with no plausible mechanism
negative_controls <- tribble(
  ~product_pattern,               ~pt_pattern,                        ~evidence,
  "^CENTRUM",                     "Rhabdomyolysis",                   "No myotoxic mechanism in standard multivitamin",
  "^CENTRUM",                     "Seizure",                          "No proconvulsant mechanism",
  "^FISH OIL|^OMEGA",             "Alopecia|Hair Loss",               "No mechanism for hair loss",
  "^FISH OIL|^OMEGA",             "Fracture",                         "No mechanism for bone fragility",
  "^VITAMIN D",                   "Tinnitus",                         "No ototoxic mechanism",
  "^VITAMIN D",                   "Rhabdomyolysis",                   "No myotoxic mechanism at standard doses",
  "^CALCIUM|^CALTRATE|^CITRACAL", "Seizure",                          "No proconvulsant mechanism",
  "^CALCIUM|^CALTRATE|^CITRACAL", "Alopecia|Hair Loss",               "No mechanism for hair loss",
  "^MELATONIN",                   "Renal Failure|Kidney",             "No nephrotoxic mechanism",
  "^MELATONIN",                   "Rhabdomyolysis",                   "No myotoxic mechanism",
  "^GLUCOSAMINE",                 "Blindness|Vision Loss",            "No ocular toxicity mechanism",
  "^GLUCOSAMINE",                 "Seizure",                          "No proconvulsant mechanism",
  "^PROBIOTIC|^CULTURELLE",       "Fracture",                         "No mechanism for bone fragility",
  "^BIOTIN",                      "Cardiac Arrest",                   "No cardiotoxic mechanism",
  "^FOLIC ACID",                  "Rhabdomyolysis",                   "No myotoxic mechanism"
)

# Function to match controls against disproportionality results
match_controls <- function(controls, disp_data, control_type) {
  results <- controls |>
    rowwise() |>
    mutate(
      matches = list({
        prod_match <- str_detect(disp_data$product_clean,
                                  regex(product_pattern, ignore_case = TRUE))
        pt_match   <- str_detect(disp_data$meddra_pt,
                                  regex(pt_pattern, ignore_case = TRUE))
        matched <- disp_data[prod_match & pt_match, ]
        if (nrow(matched) == 0) {
          tibble(found = FALSE, n_pairs = 0L, any_robust = FALSE,
                 max_prr = NA_real_, max_a = NA_integer_)
        } else {
          tibble(found = TRUE, n_pairs = nrow(matched),
                 any_robust = any(matched$robust_signal, na.rm = TRUE),
                 max_prr = max(matched$prr, na.rm = TRUE),
                 max_a = max(as.integer(matched$a), na.rm = TRUE))
        }
      })
    ) |>
    unnest(matches) |>
    ungroup() |>
    mutate(control_type = control_type)

  results
}

pos_results <- match_controls(positive_controls, disp, "positive")
neg_results <- match_controls(negative_controls, disp, "negative")

cat("Positive controls:\n")
pos_results |>
  select(product_pattern, pt_pattern, found, n_pairs, any_robust, max_prr, max_a) |>
  print(n = Inf, width = Inf)

cat("\nNegative controls:\n")
neg_results |>
  select(product_pattern, pt_pattern, found, n_pairs, any_robust, max_prr, max_a) |>
  print(n = Inf, width = Inf)

# Performance metrics
pos_found    <- pos_results |> filter(found)
neg_found    <- neg_results |> filter(found)
tp <- sum(pos_found$any_robust)
fn <- sum(!pos_found$any_robust) + sum(!pos_results$found)
fp <- sum(neg_found$any_robust)
tn <- sum(!neg_found$any_robust) + sum(!neg_results$found)

sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
ppv <- tp / max(tp + fp, 1)
npv <- tn / max(tn + fn, 1)

cat(sprintf("\nControl validation performance:\n"))
cat(sprintf("  True positives:  %d / %d positive controls detected as robust\n", tp, nrow(positive_controls)))
cat(sprintf("  False negatives: %d positive controls missed\n", fn))
cat(sprintf("  False positives: %d / %d negative controls detected as robust\n", fp, nrow(negative_controls)))
cat(sprintf("  True negatives:  %d negative controls correctly not detected\n", tn))
cat(sprintf("  Sensitivity:     %.1f%%\n", sensitivity * 100))
cat(sprintf("  Specificity:     %.1f%%\n", specificity * 100))
cat(sprintf("  PPV:             %.1f%%\n", ppv * 100))
cat(sprintf("  NPV:             %.1f%%\n", npv * 100))

control_results <- bind_rows(pos_results, neg_results)
fwrite(control_results, file.path(table_dir, "control_validation.csv"))

control_performance <- tibble(
  metric = c("Sensitivity", "Specificity", "PPV", "NPV", "TP", "FN", "FP", "TN"),
  value = c(round(sensitivity, 3), round(specificity, 3),
            round(ppv, 3), round(npv, 3), tp, fn, fp, tn)
)
fwrite(control_performance, file.path(table_dir, "control_performance.csv"))


# =============================================================================
# SECTION 3: Masking / Competition Bias Diagnostic
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 3: Masking / Competition Bias Diagnostic\n")
cat(strrep("=", 60), "\n")

# For products with a dominant signal, removing that pair can unmask
# suppressed signals for other events.

# Identify products with a dominant pair (one pair accounts for >30% of reports)
product_report_counts <- symptoms |>
  distinct(report_id, product_clean) |>
  count(product_clean, name = "total_reports")

dominant_pairs <- disp |>
  filter(robust_signal) |>
  inner_join(product_report_counts, by = "product_clean") |>
  mutate(pct_of_product = a / total_reports * 100) |>
  filter(pct_of_product > 30, total_reports >= 50) |>
  arrange(desc(pct_of_product)) |>
  select(product_clean, meddra_pt, a, total_reports, pct_of_product, prr)

cat(sprintf("Dominant product-PT pairs (>30%% of product reports, n>=50): %d\n",
            nrow(dominant_pairs)))
print(dominant_pairs |> head(20), n = 20, width = Inf)

# For each dominant pair, remove it and check what new signals emerge
unmasking_results <- list()

for (i in seq_len(min(nrow(dominant_pairs), 30))) {
  prod <- dominant_pairs$product_clean[i]
  pt   <- dominant_pairs$meddra_pt[i]

  # Remove the dominant pair's reports from the symptom data
  reports_with_dominant <- symptoms |>
    filter(product_clean == prod, meddra_pt == pt) |>
    pull(report_id) |>
    unique()

  symptoms_reduced <- symptoms |>
    filter(!(report_id %in% reports_with_dominant & product_clean == prod))

  N_reduced <- n_distinct(symptoms_reduced$report_id)

  # Recompute for remaining pairs of this product
  product_pairs_reduced <- symptoms_reduced |>
    filter(product_clean == prod) |>
    distinct(report_id, meddra_pt) |>
    count(meddra_pt, name = "a_new") |>
    filter(a_new >= MIN_N)

  if (nrow(product_pairs_reduced) == 0) next

  n_product_reduced <- n_distinct(
    symptoms_reduced$report_id[symptoms_reduced$product_clean == prod]
  )

  pt_totals_reduced <- symptoms_reduced |>
    distinct(report_id, meddra_pt) |>
    count(meddra_pt, name = "n_pt_reduced")

  reduced_prr <- product_pairs_reduced |>
    left_join(pt_totals_reduced, by = "meddra_pt") |>
    mutate(
      a = as.double(a_new),
      b = n_product_reduced - a,
      c = n_pt_reduced - a,
      d = N_reduced - a - b - c,
      prr_new = (a / (a + b)) / (c / (c + d)),
      chi2_new = ((abs(a * d - b * c) - N_reduced / 2)^2 * N_reduced) /
                 ((a + b) * (c + d) * (a + c) * (b + d)),
      prr_signal_new = prr_new >= PRR_THRESHOLD & chi2_new >= CHI2_THRESHOLD & a >= MIN_N
    )

  # Compare with original
  original_pairs <- disp |>
    filter(product_clean == prod) |>
    select(meddra_pt, prr_orig = prr, prr_signal_orig = prr_signal, robust_orig = robust_signal)

  comparison <- reduced_prr |>
    left_join(original_pairs, by = "meddra_pt") |>
    mutate(
      newly_significant = prr_signal_new & !coalesce(prr_signal_orig, FALSE),
      prr_change = prr_new / coalesce(prr_orig, prr_new)
    )

  n_unmasked <- sum(comparison$newly_significant, na.rm = TRUE)

  if (n_unmasked > 0) {
    unmasked <- comparison |>
      filter(newly_significant) |>
      mutate(dominant_product = prod, dominant_pt = pt) |>
      select(dominant_product, dominant_pt, unmasked_pt = meddra_pt,
             a_new, prr_new, prr_orig, prr_change)
    unmasking_results[[length(unmasking_results) + 1]] <- unmasked
  }
}

if (length(unmasking_results) > 0) {
  unmasked_all <- bind_rows(unmasking_results)
  cat(sprintf("\nNewly unmasked signals: %d\n", nrow(unmasked_all)))
  print(unmasked_all |> head(30), n = 30, width = Inf)
  fwrite(unmasked_all, file.path(table_dir, "masking_unmasked_signals.csv"))
} else {
  cat("\nNo newly unmasked signals detected.\n")
  unmasked_all <- tibble()
}


# =============================================================================
# SECTION 4: Time-Stratified Disproportionality
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 4: Time-Stratified Disproportionality\n")
cat(strrep("=", 60), "\n")

# Helper to run disproportionality on a subset
run_disprop_quick <- function(sym_data, label) {
  N <- n_distinct(sym_data$report_id)
  if (N < 500) return(tibble())

  product_totals <- sym_data |>
    distinct(report_id, product_clean) |>
    count(product_clean, name = "n_product")

  pt_totals <- sym_data |>
    distinct(report_id, meddra_pt) |>
    count(meddra_pt, name = "n_pt")

  pairs <- sym_data |>
    distinct(report_id, product_clean, meddra_pt) |>
    count(product_clean, meddra_pt, name = "a") |>
    filter(a >= MIN_N)

  if (nrow(pairs) == 0) return(tibble())

  ct <- pairs |>
    left_join(product_totals, by = "product_clean") |>
    left_join(pt_totals, by = "meddra_pt") |>
    mutate(
      a = as.double(a), n_product = as.double(n_product),
      n_pt = as.double(n_pt),
      b = n_product - a, c = n_pt - a,
      d = N - a - b - c, N_total = as.double(N),
      prr = (a / (a + b)) / (c / (c + d)),
      chi2 = ((abs(a * d - b * c) - N_total / 2)^2 * N_total) /
             ((a + b) * (c + d) * (a + c) * (b + d)),
      prr_signal = prr >= PRR_THRESHOLD & chi2 >= CHI2_THRESHOLD & a >= MIN_N,
      ror = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)),
      ror_se = sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5)),
      ror_lower = exp(log(ror) - 1.96 * ror_se),
      ror_signal = ror_lower > 1,
      expected = (n_product * n_pt) / N_total
    )

  mean_rr <- mean(ct$a / ct$expected, na.rm = TRUE)
  var_rr  <- var(ct$a / ct$expected, na.rm = TRUE)
  a_pr <- mean_rr^2 / var_rr; b_pr <- mean_rr / var_rr

  ct <- ct |>
    mutate(
      eb05 = log2(qgamma(0.05, shape = a_pr + a, rate = b_pr + expected)),
      gps_signal = 2^eb05 >= GPS_EB05,
      ic = log2((a + 0.5) / (expected + 0.5)),
      ic_var = 1 / (a + 0.5) - 1 / N_total,
      ic025 = ic - 1.96 * sqrt(pmax(ic_var, 0)),
      bcpnn_signal = ic025 > BCPNN_IC025,
      n_methods = as.integer(prr_signal) + as.integer(ror_signal) +
                  as.integer(gps_signal) + as.integer(bcpnn_signal),
      robust_signal = n_methods >= 3,
      period = label
    )

  ct
}

# Split into two eras
era1 <- symptoms |> filter(report_year >= 2007, report_year <= 2015)
era2 <- symptoms |> filter(report_year >= 2016)

cat(sprintf("Era 1 (2007-2015): %s reports\n",
            format(n_distinct(era1$report_id), big.mark = ",")))
cat(sprintf("Era 2 (2016-2025): %s reports\n",
            format(n_distinct(era2$report_id), big.mark = ",")))

cat("\nRunning disproportionality for Era 1...\n")
era1_results <- run_disprop_quick(era1, "2007-2015")
cat(sprintf("Era 1 robust signals: %d\n", sum(era1_results$robust_signal)))

cat("Running disproportionality for Era 2...\n")
era2_results <- run_disprop_quick(era2, "2016-2025")
cat(sprintf("Era 2 robust signals: %d\n", sum(era2_results$robust_signal)))

# Compare concordance between eras
era1_robust <- era1_results |>
  filter(robust_signal) |>
  select(product_clean, meddra_pt) |>
  mutate(in_era1 = TRUE)

era2_robust <- era2_results |>
  filter(robust_signal) |>
  select(product_clean, meddra_pt) |>
  mutate(in_era2 = TRUE)

era_concordance <- era1_robust |>
  full_join(era2_robust, by = c("product_clean", "meddra_pt")) |>
  replace_na(list(in_era1 = FALSE, in_era2 = FALSE))

n_both_eras  <- sum(era_concordance$in_era1 & era_concordance$in_era2)
n_era1_only  <- sum(era_concordance$in_era1 & !era_concordance$in_era2)
n_era2_only  <- sum(!era_concordance$in_era1 & era_concordance$in_era2)

cat(sprintf("\nTime-stratified concordance:\n"))
cat(sprintf("  In both eras:    %d\n", n_both_eras))
cat(sprintf("  Era 1 only:      %d\n", n_era1_only))
cat(sprintf("  Era 2 only:      %d\n", n_era2_only))
cat(sprintf("  Jaccard index:   %.3f\n",
            n_both_eras / (n_both_eras + n_era1_only + n_era2_only)))

# Era-2-only signals are potentially emerging
era2_only_signals <- era_concordance |>
  filter(!in_era1 & in_era2) |>
  left_join(era2_results |> select(product_clean, meddra_pt, a, prr, n_methods),
            by = c("product_clean", "meddra_pt")) |>
  arrange(desc(a))

cat(sprintf("\nTop 20 era-2-only signals (potentially emerging):\n"))
era2_only_signals |> head(20) |> print(n = 20, width = Inf)

time_strat_summary <- tibble(
  metric = c("Era 1 robust", "Era 2 robust", "Both eras", "Era 1 only",
             "Era 2 only", "Jaccard index"),
  value = c(sum(era1_results$robust_signal), sum(era2_results$robust_signal),
            n_both_eras, n_era1_only, n_era2_only,
            round(n_both_eras / (n_both_eras + n_era1_only + n_era2_only), 3))
)
fwrite(time_strat_summary, file.path(table_dir, "time_stratified_concordance.csv"))
fwrite(era2_only_signals, file.path(table_dir, "era2_only_emerging_signals.csv"))


# =============================================================================
# SECTION 5: OpenFDA Enforcement Cross-Reference
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 5: OpenFDA Enforcement Cross-Reference\n")
cat(strrep("=", 60), "\n")

# Query the OpenFDA food enforcement (recall) database for top products
top_products <- robust |>
  count(product_clean, sort = TRUE) |>
  head(50) |>
  pull(product_clean)

# Extract a short search term from each product name (first 1-2 meaningful words)
search_terms <- tibble(product_clean = top_products) |>
  mutate(
    search_term = str_extract(product_clean, "^[A-Z]+(?:\\s+[A-Z]+)?") |>
                    str_to_lower() |>
                    str_squish()
  ) |>
  filter(!is.na(search_term), nchar(search_term) >= 3) |>
  distinct(search_term, .keep_all = TRUE)

cat(sprintf("Querying OpenFDA for %d products...\n", nrow(search_terms)))

fda_results <- list()

for (i in seq_len(nrow(search_terms))) {
  term <- search_terms$search_term[i]
  prod <- search_terms$product_clean[i]

  url <- paste0(
    "https://api.fda.gov/food/enforcement.json?search=product_description:",
    URLencode(paste0('"', term, '"')),
    "+AND+product_description:",
    URLencode('"dietary supplement"'),
    "&limit=5"
  )

  result <- tryCatch({
    resp <- fromJSON(url, flatten = TRUE)
    if (!is.null(resp$results)) {
      tibble(
        product_clean = prod,
        search_term = term,
        n_recalls = nrow(resp$results),
        classifications = paste(unique(resp$results$classification), collapse = "; "),
        reasons = paste(str_trunc(unique(resp$results$reason_for_recall), 100),
                        collapse = " | "),
        statuses = paste(unique(resp$results$status), collapse = "; ")
      )
    } else {
      tibble(product_clean = prod, search_term = term,
             n_recalls = 0L, classifications = NA_character_,
             reasons = NA_character_, statuses = NA_character_)
    }
  }, error = function(e) {
    tibble(product_clean = prod, search_term = term,
           n_recalls = 0L, classifications = NA_character_,
           reasons = NA_character_, statuses = NA_character_)
  })

  fda_results[[i]] <- result
  Sys.sleep(0.3)  # rate limiting
}

fda_enforcement <- bind_rows(fda_results)

products_with_recalls <- fda_enforcement |> filter(n_recalls > 0)
cat(sprintf("\nProducts with FDA enforcement actions: %d / %d queried\n",
            nrow(products_with_recalls), nrow(search_terms)))

if (nrow(products_with_recalls) > 0) {
  cat("\nProducts with recalls:\n")
  products_with_recalls |>
    select(product_clean, n_recalls, classifications, reasons) |>
    print(n = Inf, width = Inf)
}

fwrite(fda_enforcement, file.path(table_dir, "openfda_enforcement.csv"))


# =============================================================================
# SECTION 6: Concomitant Product Characterisation
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 6: Concomitant Product Characterisation\n")
cat(strrep("=", 60), "\n")

# How many suspect products per report?
products_per_report <- supps |>
  group_by(report_id) |>
  summarise(n_products = n_distinct(product_clean), .groups = "drop")

cat("Suspect products per report:\n")
print(summary(products_per_report$n_products))

prods_table <- products_per_report |>
  mutate(n_cat = case_when(
    n_products == 1 ~ "1 product",
    n_products == 2 ~ "2 products",
    n_products == 3 ~ "3 products",
    TRUE ~ "4+ products"
  )) |>
  count(n_cat, name = "n_reports") |>
  mutate(pct = round(n_reports / sum(n_reports) * 100, 1))

cat("\nDistribution:\n")
print(prods_table, n = Inf, width = Inf)

# Does signal strength differ by number of concomitant products?
# For each product-PT pair, compute the average number of concomitant products
# in reports contributing to that pair
pair_concomitant <- symptoms |>
  distinct(report_id, product_clean, meddra_pt) |>
  left_join(products_per_report, by = "report_id") |>
  group_by(product_clean, meddra_pt) |>
  summarise(
    mean_products_per_report = mean(n_products, na.rm = TRUE),
    pct_multi_product = mean(n_products > 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Join with robust signal status
concomitant_by_signal <- disp |>
  select(product_clean, meddra_pt, robust_signal, prr) |>
  left_join(pair_concomitant, by = c("product_clean", "meddra_pt"))

cat("\nConcomitant product rates by signal status:\n")
concomitant_by_signal |>
  group_by(robust_signal) |>
  summarise(
    n_pairs = n(),
    mean_products_per_report = mean(mean_products_per_report, na.rm = TRUE),
    mean_pct_multi_product = mean(pct_multi_product, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print(width = Inf)

# Sensitivity: rerun on single-product reports only
cat("\nSensitivity: single-product reports only\n")
single_product_reports <- products_per_report |>
  filter(n_products == 1) |>
  pull(report_id)

symptoms_single <- symptoms |>
  filter(report_id %in% single_product_reports)

cat(sprintf("Single-product reports: %s / %s (%.1f%%)\n",
            format(length(single_product_reports), big.mark = ","),
            format(N_total, big.mark = ","),
            length(single_product_reports) / N_total * 100))

single_results <- run_disprop_quick(symptoms_single, "single_product")
n_single_robust <- sum(single_results$robust_signal, na.rm = TRUE)

cat(sprintf("Single-product robust signals: %d (vs %d all reports)\n",
            n_single_robust, nrow(robust)))

# Concordance with full analysis
single_robust_pairs <- single_results |>
  filter(robust_signal) |>
  select(product_clean, meddra_pt) |>
  mutate(in_single = TRUE)

full_robust_pairs <- robust |>
  select(product_clean, meddra_pt) |>
  mutate(in_full = TRUE)

concom_concordance <- full_robust_pairs |>
  full_join(single_robust_pairs, by = c("product_clean", "meddra_pt")) |>
  replace_na(list(in_full = FALSE, in_single = FALSE))

n_concom_both <- sum(concom_concordance$in_full & concom_concordance$in_single)
cat(sprintf("Concordance: %d / %d full signals (%.1f%%) also detected in single-product\n",
            n_concom_both, nrow(robust), n_concom_both / nrow(robust) * 100))

fwrite(prods_table, file.path(table_dir, "concomitant_distribution.csv"))


# =============================================================================
# SECTION 7: SOC-Level Signal Aggregation
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 7: SOC-Level Signal Aggregation\n")
cat(strrep("=", 60), "\n")

# Map robust signals to SOC via the symptoms_long SOC column
pt_to_soc <- symptoms |>
  distinct(meddra_pt, soc) |>
  filter(!is.na(soc), soc != "")

robust_with_soc <- robust |>
  left_join(pt_to_soc, by = "meddra_pt")

soc_summary <- robust_with_soc |>
  filter(!is.na(soc)) |>
  group_by(soc) |>
  summarise(
    n_signals = n(),
    n_products = n_distinct(product_clean),
    mean_prr = round(mean(prr, na.rm = TRUE), 2),
    median_cases = median(a),
    total_cases = sum(a),
    .groups = "drop"
  ) |>
  arrange(desc(n_signals))

cat("Robust signals by System Organ Class:\n")
print(soc_summary, n = Inf, width = Inf)

# SOC by product category
soc_by_category <- robust_with_soc |>
  filter(!is.na(soc)) |>
  left_join(
    supps |> distinct(product_clean, product_category),
    by = "product_clean"
  ) |>
  count(product_category, soc, name = "n_signals") |>
  arrange(product_category, desc(n_signals))

# Create a cross-tabulation for the top categories and SOCs
top_socs <- soc_summary |> head(10) |> pull(soc)
top_cats <- c("Vitamin/Mineral", "Weight Loss/Diet", "Herbal/Botanical",
              "Energy/Stimulant", "Sports Nutrition")

soc_cross <- soc_by_category |>
  filter(soc %in% top_socs, product_category %in% top_cats) |>
  pivot_wider(names_from = product_category, values_from = n_signals,
              values_fill = 0)

cat("\nSOC x Product Category cross-tabulation (top SOCs, top categories):\n")
print(soc_cross, n = Inf, width = Inf)

fwrite(soc_summary, file.path(table_dir, "soc_signal_summary.csv"))
fwrite(soc_by_category, file.path(table_dir, "soc_by_category.csv"))


# =============================================================================
# SECTION 8: Bootstrap False Discovery Estimation
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 8: Bootstrap False Discovery Estimation\n")
cat(strrep("=", 60), "\n")

# Permute product-event associations and count spurious signals under the null.
# Approach: shuffle which product sets are assigned to which reports,
# preserving the within-report product count and event count.

# Pre-compute report-level product and event lists
report_products <- symptoms |>
  distinct(report_id, product_clean)

report_events <- symptoms |>
  distinct(report_id, meddra_pt)

report_demo <- symptoms |>
  distinct(report_id)

# Group products by report for efficient shuffling
product_sets <- report_products |>
  group_by(report_id) |>
  summarise(products = list(product_clean), .groups = "drop")

cat(sprintf("Running %d bootstrap permutations...\n", N_BOOTSTRAP))

null_signal_counts <- integer(N_BOOTSTRAP)

for (b in seq_len(N_BOOTSTRAP)) {
  if (b %% 10 == 0) cat(sprintf("  Permutation %d / %d\n", b, N_BOOTSTRAP))

  # Shuffle: randomly reassign product sets to different reports
  shuffled <- product_sets |>
    mutate(products = sample(products))

  # Unnest to get permuted report × product pairs
  perm_rp <- shuffled |>
    unnest(products) |>
    rename(product_clean = products)

  # Join with events to create permuted product-event pairs
  perm_pairs <- perm_rp |>
    inner_join(report_events, by = "report_id", relationship = "many-to-many") |>
    distinct(report_id, product_clean, meddra_pt) |>
    count(product_clean, meddra_pt, name = "a") |>
    filter(a >= MIN_N)

  if (nrow(perm_pairs) == 0) {
    null_signal_counts[b] <- 0
    next
  }

  # Compute PRR (quick version — just PRR signal)
  perm_product_totals <- perm_rp |>
    distinct(report_id, product_clean) |>
    count(product_clean, name = "n_product")

  perm_pt_totals <- report_events |>
    count(meddra_pt, name = "n_pt")

  perm_ct <- perm_pairs |>
    left_join(perm_product_totals, by = "product_clean") |>
    left_join(perm_pt_totals, by = "meddra_pt") |>
    mutate(
      a = as.double(a), n_product = as.double(n_product),
      n_pt = as.double(n_pt),
      b = n_product - a, c = n_pt - a,
      d = N_total - a - b - c,
      prr = (a / (a + b)) / (c / (c + d)),
      chi2 = ((abs(a * d - b * c) - N_total / 2)^2 * N_total) /
             ((a + b) * (c + d) * (a + c) * (b + d)),
      prr_signal = prr >= PRR_THRESHOLD & chi2 >= CHI2_THRESHOLD
    )

  null_signal_counts[b] <- sum(perm_ct$prr_signal, na.rm = TRUE)
}

# Compare observed vs null
observed_prr_signals <- sum(disp$prr_signal, na.rm = TRUE)

cat(sprintf("\nBootstrap null distribution of PRR signals:\n"))
cat(sprintf("  Observed PRR signals:    %d\n", observed_prr_signals))
cat(sprintf("  Null mean:               %.1f\n", mean(null_signal_counts)))
cat(sprintf("  Null SD:                 %.1f\n", sd(null_signal_counts)))
cat(sprintf("  Null max:                %d\n", max(null_signal_counts)))
cat(sprintf("  Null 95th percentile:    %.0f\n", quantile(null_signal_counts, 0.95)))

# Empirical FDR estimate: E[false positives] / observed positives
empirical_fdr <- mean(null_signal_counts) / observed_prr_signals
cat(sprintf("  Empirical FDR estimate:  %.4f (%.2f%%)\n",
            empirical_fdr, empirical_fdr * 100))
cat(sprintf("  Signal-to-noise ratio:   %.1f\n",
            observed_prr_signals / max(mean(null_signal_counts), 1)))

bootstrap_summary <- tibble(
  metric = c("Observed PRR signals", "Null mean", "Null SD",
             "Null 95th pctile", "Empirical FDR", "Signal-to-noise"),
  value = c(observed_prr_signals, round(mean(null_signal_counts), 1),
            round(sd(null_signal_counts), 1),
            round(quantile(null_signal_counts, 0.95)),
            round(empirical_fdr, 4),
            round(observed_prr_signals / max(mean(null_signal_counts), 1), 1))
)
fwrite(bootstrap_summary, file.path(table_dir, "bootstrap_fdr.csv"))


# =============================================================================
# SECTION 9: Regression-Based Signal Detection (Adjusted ORs)
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 9: Regression-Based Signal Detection\n")
cat(strrep("=", 60), "\n")

# For top product-PT pairs, fit logistic regression adjusting for age and sex
# and compare adjusted OR to crude ROR.

# Pre-compute report-level data
report_has_product <- report_products |>
  mutate(has = TRUE) |>
  group_by(report_id) |>
  summarise(products = list(product_clean), .groups = "drop")

report_has_event <- report_events |>
  group_by(report_id) |>
  summarise(events = list(meddra_pt), .groups = "drop")

report_covariates <- symptoms |>
  distinct(report_id, age_group, sex_clean)

# Select top pairs by case count
top_pairs <- disp |>
  filter(robust_signal) |>
  arrange(desc(a)) |>
  head(N_REGRESSION) |>
  select(product_clean, meddra_pt, a, ror, ror_lower, robust_signal)

cat(sprintf("Fitting adjusted logistic regression for top %d pairs...\n",
            nrow(top_pairs)))

regression_results <- list()

for (i in seq_len(nrow(top_pairs))) {
  if (i %% 50 == 0) cat(sprintf("  Pair %d / %d\n", i, nrow(top_pairs)))

  prod <- top_pairs$product_clean[i]
  pt   <- top_pairs$meddra_pt[i]

  # Build report-level binary indicators
  reports_with_product <- report_products |>
    filter(product_clean == prod) |>
    pull(report_id)

  reports_with_event <- report_events |>
    filter(meddra_pt == pt) |>
    pull(report_id)

  model_data <- report_covariates |>
    mutate(
      has_product = report_id %in% reports_with_product,
      has_event = as.integer(report_id %in% reports_with_event)
    ) |>
    filter(!is.na(age_group), !is.na(sex_clean),
           age_group != "Unknown", sex_clean != "Unknown")

  result <- tryCatch({
    fit <- glm(has_event ~ has_product + age_group + sex_clean,
               data = model_data, family = binomial)

    coefs <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    product_coef <- coefs |> filter(term == "has_productTRUE")

    if (nrow(product_coef) == 1) {
      tibble(
        product_clean = prod, meddra_pt = pt,
        crude_ror = top_pairs$ror[i],
        adjusted_or = product_coef$estimate,
        adj_or_lower = product_coef$conf.low,
        adj_or_upper = product_coef$conf.high,
        adj_p_value = product_coef$p.value,
        adj_significant = product_coef$conf.low > 1
      )
    } else {
      tibble(product_clean = prod, meddra_pt = pt,
             crude_ror = top_pairs$ror[i],
             adjusted_or = NA_real_, adj_or_lower = NA_real_,
             adj_or_upper = NA_real_, adj_p_value = NA_real_,
             adj_significant = NA)
    }
  }, error = function(e) {
    tibble(product_clean = prod, meddra_pt = pt,
           crude_ror = top_pairs$ror[i],
           adjusted_or = NA_real_, adj_or_lower = NA_real_,
           adj_or_upper = NA_real_, adj_p_value = NA_real_,
           adj_significant = NA)
  })

  regression_results[[i]] <- result
}

adj_results <- bind_rows(regression_results) |>
  filter(!is.na(adjusted_or))

cat(sprintf("\nSuccessful regressions: %d / %d\n", nrow(adj_results), nrow(top_pairs)))

# Compare crude vs adjusted
adj_results <- adj_results |>
  mutate(
    or_ratio = adjusted_or / crude_ror,
    direction_change = case_when(
      crude_ror > 1 & adjusted_or < 1 ~ "Signal lost after adjustment",
      crude_ror < 1 & adjusted_or > 1 ~ "Signal gained after adjustment",
      TRUE ~ "Consistent direction"
    ),
    magnitude_change = case_when(
      abs(log(or_ratio)) < log(1.2) ~ "Minimal (<20% change)",
      abs(log(or_ratio)) < log(1.5) ~ "Moderate (20-50% change)",
      TRUE ~ "Substantial (>50% change)"
    )
  )

cat("\nCrude vs adjusted OR comparison:\n")
adj_results |>
  count(direction_change, name = "n_pairs") |>
  print(width = Inf)

cat("\nMagnitude of confounding:\n")
adj_results |>
  count(magnitude_change, name = "n_pairs") |>
  print(width = Inf)

cat(sprintf("\nCorrelation (log crude ROR vs log adjusted OR): %.3f\n",
            cor(log(adj_results$crude_ror), log(adj_results$adjusted_or),
                use = "complete.obs")))

# Signals that change classification after adjustment
signals_lost <- adj_results |>
  filter(crude_ror > 1 & !adj_significant) |>
  arrange(adj_p_value)

cat(sprintf("\nSignals where adjusted OR lower CI crosses 1: %d / %d (%.1f%%)\n",
            nrow(signals_lost), nrow(adj_results),
            nrow(signals_lost) / nrow(adj_results) * 100))

if (nrow(signals_lost) > 0) {
  cat("Top 20 signals weakened by demographic adjustment:\n")
  signals_lost |>
    select(product_clean, meddra_pt, crude_ror, adjusted_or,
           adj_or_lower, adj_p_value) |>
    head(20) |>
    print(n = 20, width = Inf)
}

fwrite(adj_results, file.path(table_dir, "regression_adjusted_ors.csv"))


# =============================================================================
# Summary
# =============================================================================

cat("\n\n")
cat(strrep("=", 60), "\n")
cat("EXTENDED VALIDATION COMPLETE\n")
cat(strrep("=", 60), "\n")

cat("\n1. MULTIPLE TESTING\n")
cat(sprintf("   BH FDR<0.05 signals: %d; Robust signals also FDR<0.05: %.1f%%\n",
            sum(disp$fdr_sig_05, na.rm = TRUE), pct_robust_fdr))

cat("\n2. CONTROL VALIDATION\n")
cat(sprintf("   Sensitivity: %.1f%%; Specificity: %.1f%%; PPV: %.1f%%\n",
            sensitivity * 100, specificity * 100, ppv * 100))

cat("\n3. MASKING ANALYSIS\n")
cat(sprintf("   Dominant pairs tested: %d; Unmasked signals: %d\n",
            min(nrow(dominant_pairs), 30), nrow(unmasked_all)))

cat("\n4. TIME STRATIFICATION\n")
cat(sprintf("   Era 1 robust: %d; Era 2 robust: %d; Jaccard: %.3f\n",
            sum(era1_results$robust_signal), sum(era2_results$robust_signal),
            n_both_eras / (n_both_eras + n_era1_only + n_era2_only)))

cat("\n5. OPENFDA ENFORCEMENT\n")
cat(sprintf("   Products with recalls: %d / %d queried\n",
            nrow(products_with_recalls), nrow(search_terms)))

cat("\n6. CONCOMITANT PRODUCTS\n")
cat(sprintf("   Single-product concordance: %.1f%%\n",
            n_concom_both / nrow(robust) * 100))

cat("\n7. SOC AGGREGATION\n")
cat(sprintf("   SOCs with signals: %d\n", nrow(soc_summary)))

cat("\n8. BOOTSTRAP FDR\n")
cat(sprintf("   Empirical FDR: %.2f%%; Signal-to-noise: %.1f\n",
            empirical_fdr * 100,
            observed_prr_signals / max(mean(null_signal_counts), 1)))

cat("\n9. REGRESSION ADJUSTMENT\n")
cat(sprintf("   Cor(crude, adjusted): %.3f; Signals lost after adjustment: %d / %d\n",
            cor(log(adj_results$crude_ror), log(adj_results$adjusted_or),
                use = "complete.obs"),
            nrow(signals_lost), nrow(adj_results)))

cat("\n")
cat(strrep("=", 60), "\n")
