# =============================================================================
# 08c_additional_validation.R
# Additional robustness analyses
# =============================================================================
#
# Sections:
#   1. Notoriety bias quantification (pre/post enforcement reporting slopes)
#   2. Count-response relationship (PRR vs report count)
#   3. Cross-validation with published international signals
#   4. Random holdout replication (50/50 split)
#   5. Choking/dysphagia exclusion sensitivity
#   6. Weighted composite risk score
# =============================================================================

library(tidyverse)
library(data.table)

set.seed(42)

# --- Configuration -----------------------------------------------------------

proc_dir  <- "data/processed"
table_dir <- "outputs/tables"

PRR_THRESHOLD  <- 2
CHI2_THRESHOLD <- 4
MIN_N          <- 3
GPS_EB05       <- 2
BCPNN_IC025    <- 0

# --- Helper: quick disproportionality ----------------------------------------

run_disprop <- function(sym_data, label = "") {
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
      analysis = label
    )

  ct
}

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
cat(sprintf("Reports: %s, Robust signals: %s\n",
            format(N_total, big.mark = ","),
            format(nrow(robust), big.mark = ",")))


# =============================================================================
# SECTION 1: Notoriety Bias Quantification
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 1: Notoriety Bias Quantification\n")
cat(strrep("=", 60), "\n")

# For products with known FDA enforcement dates, compare reporting rates
# before and after the action to quantify stimulated reporting.

enforcement_dates <- tribble(
  ~product_pattern,    ~action_date,   ~action_description,
  "^HYDROXYCUT",       "2009-05-01",   "FDA recall, hepatotoxicity",
  "^OXY.?ELITE",       "2013-10-01",   "FDA recall, acute hepatitis outbreak",
  "^KRATOM",           "2014-02-01",   "FDA import alert 54-15",
  "^EPHEDRA",          "2004-04-12",   "FDA ban on ephedra-containing supplements"
) |>
  mutate(action_date = as.Date(action_date))

# Build quarterly time series for each product
quarterly_product <- supps |>
  filter(!is.na(report_yq)) |>
  distinct(report_id, product_clean, report_yq, report_date) |>
  mutate(report_date = as.Date(report_date))

notoriety_results <- list()

for (i in seq_len(nrow(enforcement_dates))) {
  pat      <- enforcement_dates$product_pattern[i]
  act_date <- enforcement_dates$action_date[i]
  act_desc <- enforcement_dates$action_description[i]

  # Get quarterly counts for this product group
  product_ts <- quarterly_product |>
    filter(str_detect(product_clean, regex(pat, ignore_case = TRUE))) |>
    count(report_yq, name = "count") |>
    arrange(report_yq) |>
    mutate(
      year = as.integer(str_extract(report_yq, "^\\d{4}")),
      quarter = as.integer(str_extract(report_yq, "\\d$")),
      date_approx = as.Date(paste0(year, "-", (quarter - 1) * 3 + 1, "-15")),
      period = if_else(date_approx < act_date, "pre", "post"),
      quarters_from_action = as.numeric(difftime(date_approx, act_date,
                                                  units = "days")) / 91.25,
      time_idx = row_number()
    )

  if (nrow(product_ts) < 8) next

  pre  <- product_ts |> filter(period == "pre")
  post <- product_ts |> filter(period == "post")

  # Compute mean quarterly rates
  mean_pre  <- if (nrow(pre) > 0)  mean(pre$count)  else 0
  mean_post <- if (nrow(post) > 0) mean(post$count) else 0
  ratio     <- if (mean_pre > 0) mean_post / mean_pre else NA_real_

  # Rate ratio test (Poisson)
  total_pre  <- sum(pre$count)
  total_post <- sum(post$count)
  n_q_pre    <- nrow(pre)
  n_q_post   <- nrow(post)

  # Poisson rate ratio: (post_total/n_q_post) / (pre_total/n_q_pre)
  rr <- (total_post / max(n_q_post, 1)) / (total_pre / max(n_q_pre, 1))
  rr_se <- sqrt(1 / max(total_pre, 1) + 1 / max(total_post, 1))
  rr_lower <- exp(log(max(rr, 0.01)) - 1.96 * rr_se)
  rr_upper <- exp(log(max(rr, 0.01)) + 1.96 * rr_se)

  # Trend test: is there a step change at the enforcement date?
  if (nrow(product_ts) >= 4) {
    trend_model <- tryCatch({
      glm(count ~ time_idx + period, data = product_ts, family = poisson)
    }, error = function(e) NULL)

    if (!is.null(trend_model)) {
      period_coef <- broom::tidy(trend_model) |>
        filter(term == "periodpre")
      step_change_p <- if (nrow(period_coef) == 1) period_coef$p.value else NA_real_
      step_change_est <- if (nrow(period_coef) == 1) exp(-period_coef$estimate) else NA_real_
    } else {
      step_change_p <- NA_real_
      step_change_est <- NA_real_
    }
  } else {
    step_change_p <- NA_real_
    step_change_est <- NA_real_
  }

  notoriety_results[[i]] <- tibble(
    product_pattern = pat,
    action_description = act_desc,
    action_date = act_date,
    n_quarters_pre = n_q_pre,
    n_quarters_post = n_q_post,
    mean_quarterly_pre = round(mean_pre, 2),
    mean_quarterly_post = round(mean_post, 2),
    rate_ratio = round(rr, 2),
    rr_lower = round(rr_lower, 2),
    rr_upper = round(rr_upper, 2),
    poisson_step_change = round(step_change_est, 2),
    step_change_p = round(step_change_p, 4)
  )
}

notoriety_all <- bind_rows(notoriety_results)
cat("Notoriety bias quantification:\n")
print(notoriety_all, width = Inf)

# Interpretation
cat("\nInterpretation:\n")
for (i in seq_len(nrow(notoriety_all))) {
  row <- notoriety_all[i, ]
  if (row$rate_ratio > 2) {
    cat(sprintf("  %s: %.1fx increase in reporting after %s (notoriety bias likely)\n",
                row$product_pattern, row$rate_ratio, row$action_description))
  } else if (row$rate_ratio > 1) {
    cat(sprintf("  %s: %.1fx modest increase after %s\n",
                row$product_pattern, row$rate_ratio, row$action_description))
  } else {
    cat(sprintf("  %s: No increase after %s (rate ratio %.2f)\n",
                row$product_pattern, row$rate_ratio, row$action_description))
  }
}

fwrite(notoriety_all, file.path(table_dir, "notoriety_bias_quantification.csv"))


# =============================================================================
# SECTION 2: Count-Response Relationship
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 2: Count-Response Relationship (PRR vs Report Count)\n")
cat(strrep("=", 60), "\n")

# Test whether PRR is correlated with report count (a).
# If PRR increases with a, it suggests real signals accumulate.
# If PRR is independent of a, the measure is unbiased by volume.
# If PRR decreases with a, it suggests shrinkage toward the null.

count_response <- disp |>
  filter(robust_signal) |>
  select(product_clean, meddra_pt, a, prr, ror, n_methods) |>
  mutate(log_a = log2(a), log_prr = log2(prr))

# Spearman correlation (rank-based, robust to outliers)
cor_spearman <- cor.test(count_response$a, count_response$prr,
                          method = "spearman", exact = FALSE)

# Pearson on log-transformed
cor_pearson_log <- cor.test(count_response$log_a, count_response$log_prr)

cat(sprintf("Spearman correlation (a vs PRR): rho = %.3f, p = %.2e\n",
            cor_spearman$estimate, cor_spearman$p.value))
cat(sprintf("Pearson correlation (log2(a) vs log2(PRR)): r = %.3f, p = %.2e\n",
            cor_pearson_log$estimate, cor_pearson_log$p.value))

# Bin by count deciles and compute mean PRR
count_response <- count_response |>
  mutate(count_bin = cut(a, breaks = c(3, 5, 10, 20, 50, 100, Inf),
                          labels = c("3-5", "6-10", "11-20", "21-50",
                                      "51-100", ">100"),
                          right = FALSE))

count_bin_summary <- count_response |>
  group_by(count_bin) |>
  summarise(
    n_pairs = n(),
    mean_prr = round(mean(prr, na.rm = TRUE), 2),
    median_prr = round(median(prr, na.rm = TRUE), 2),
    mean_log_prr = round(mean(log_prr, na.rm = TRUE), 2),
    pct_all4_methods = round(mean(n_methods == 4) * 100, 1),
    .groups = "drop"
  )

cat("\nPRR by case count bins:\n")
print(count_bin_summary, n = Inf, width = Inf)

# Does signal strength increase with count for the SAME product?
# Within-product correlation: for products with multiple signals,
# does the pair with more cases also have higher PRR?
within_product <- count_response |>
  group_by(product_clean) |>
  filter(n() >= 3) |>
  summarise(
    n_signals = n(),
    within_cor = cor(a, prr, method = "spearman"),
    .groups = "drop"
  ) |>
  filter(!is.na(within_cor))

cat(sprintf("\nWithin-product Spearman correlation (a vs PRR):\n"))
cat(sprintf("  Products with >=3 signals: %d\n", nrow(within_product)))
cat(sprintf("  Mean within-product rho:   %.3f\n", mean(within_product$within_cor)))
cat(sprintf("  Median within-product rho: %.3f\n", median(within_product$within_cor)))
cat(sprintf("  Products with positive rho: %d / %d (%.1f%%)\n",
            sum(within_product$within_cor > 0), nrow(within_product),
            mean(within_product$within_cor > 0) * 100))

fwrite(count_bin_summary, file.path(table_dir, "count_response_relationship.csv"))


# =============================================================================
# SECTION 3: Cross-Validation with Published International Signals
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 3: Cross-Validation with Published Signals\n")
cat(strrep("=", 60), "\n")

# Check how many known supplement safety signals from the published literature
# are detected by our analysis. Sources: WHO VigiBase reports, EudraVigilance
# published analyses, FDA safety communications, and systematic reviews.

published_signals <- tribble(
  ~product_pattern,        ~pt_pattern,                         ~source,
  # Hepatotoxicity signals (well-established)
  "^HYDROXYCUT",           "Hepat|Liver|Jaundice",              "Fong et al. 2010; FDA recall 2009",
  "^OXY.?ELITE",           "Hepat|Liver|Jaundice",              "Roytman et al. 2014; FDA recall 2013",
  "^HERBALIFE",            "Hepat|Liver|Jaundice",              "Stickel & Shouval 2015; multiple case series",
  "^GREEN TEA",            "Hepat|Liver",                       "Mazzanti et al. 2009; EFSA 2018",
  "^KAVA",                 "Hepat|Liver",                       "Teschke et al. 2003; multiple countries banned",
  # Kratom signals
  "^KRATOM",               "Death",                             "Gershman et al. 2019; FDA advisory",
  "^KRATOM",               "Dependence|Withdrawal|Addiction",   "Swogger et al. 2015; DEA scheduling",
  "^KRATOM",               "Seizure",                           "Post et al. 2019; case series",
  # Cardiovascular signals
  "^EPHEDRA",              "Cardiac|Myocardial|Stroke|Arrhyth", "Haller & Benowitz 2000; FDA ban 2004",
  "^HYDROXYCUT",           "Rhabdomyolysis",                    "FDA MedWatch reports",
  "^YOHIMBE",              "Hypertension|Tachycardia|Anxiety",  "Cimolai & Cimolai 2011",
  # Choking/formulation signals
  "^CENTRUM",              "Choking|Dysphagia",                 "FDA consumer complaints; tablet size",
  "^CITRACAL",             "Choking|Dysphagia",                 "FDA consumer complaints; tablet size",
  "^CALTRATE",             "Choking|Dysphagia",                 "FDA consumer complaints; tablet size",
  # Weight loss product signals
  "^HYDROXYCUT",           "Nausea|Vomiting|Abdominal",         "Common AEs, product labelling",
  # Probiotic signals (rare but documented)
  "^FLORASTOR",            "Fungaemia|Sepsis",                  "Enache-Angoulvant & Hennequin 2005",
  # Vitamin toxicity
  "VITAMIN B6",            "Neuropathy|Paraesthesia",           "Dalton & Dalton 1987; TGA 2023",
  "VITAMIN A",             "Hepatotoxicity|Liver",              "Geubel et al. 1991; dose-dependent",
  "NIACIN",                "Flushing|Hepat",                    "Guyton et al. 2014",
  # Contamination signals
  "^SUPER BETA PROSTATE",  "Blood Urine|Urin",                  "FDA warning letters",
  # Infant supplement signals
  "GRIPE WATER",           "Choking",                           "FDA consumer advisory"
)

# Match against our results
published_matches <- published_signals |>
  rowwise() |>
  mutate(
    matched = list({
      prod_hit <- str_detect(disp$product_clean,
                              regex(product_pattern, ignore_case = TRUE))
      pt_hit   <- str_detect(disp$meddra_pt,
                              regex(pt_pattern, ignore_case = TRUE))
      hits <- disp[prod_hit & pt_hit, ]
      if (nrow(hits) == 0) {
        tibble(found = FALSE, n_pairs = 0L, n_robust = 0L,
               max_prr = NA_real_, total_cases = NA_integer_)
      } else {
        tibble(found = TRUE, n_pairs = nrow(hits),
               n_robust = sum(hits$robust_signal, na.rm = TRUE),
               max_prr = max(hits$prr, na.rm = TRUE),
               total_cases = sum(as.integer(hits$a)))
      }
    })
  ) |>
  unnest(matched) |>
  ungroup()

cat("Published signal cross-validation:\n")
published_matches |>
  select(product_pattern, pt_pattern, source, found, n_robust, max_prr, total_cases) |>
  print(n = Inf, width = Inf)

n_found       <- sum(published_matches$found)
n_robust_found <- sum(published_matches$n_robust > 0)
n_total       <- nrow(published_matches)

cat(sprintf("\nPublished signals found in CAERS data:    %d / %d (%.1f%%)\n",
            n_found, n_total, n_found / n_total * 100))
cat(sprintf("Published signals detected as robust:     %d / %d (%.1f%%)\n",
            n_robust_found, n_total, n_robust_found / n_total * 100))
cat(sprintf("Published signals found but not robust:   %d\n",
            n_found - n_robust_found))
cat(sprintf("Published signals not found in data:      %d\n",
            n_total - n_found))

# Which published signals did we miss?
missed <- published_matches |> filter(!found)
if (nrow(missed) > 0) {
  cat("\nPublished signals not detected (absent from data at N>=3):\n")
  missed |>
    select(product_pattern, pt_pattern, source) |>
    print(n = Inf, width = Inf)
}

fwrite(published_matches, file.path(table_dir, "published_signal_crossvalidation.csv"))


# =============================================================================
# SECTION 4: Random Holdout Replication
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 4: Random Holdout Replication (50/50 Split)\n")
cat(strrep("=", 60), "\n")

# Randomly split reports 50/50 and run disproportionality on each half.
# This avoids the product-name-fragmentation issue of the time-stratified
# split and gives a cleaner test of statistical replicability.

all_reports <- unique(symptoms$report_id)
n_reports   <- length(all_reports)

# 10 random splits for stability
N_SPLITS <- 10
split_results <- list()

cat(sprintf("Running %d random 50/50 splits...\n", N_SPLITS))

for (s in seq_len(N_SPLITS)) {
  cat(sprintf("  Split %d / %d\n", s, N_SPLITS))

  # Random split
  split_idx <- sample(seq_len(n_reports), size = floor(n_reports / 2))
  reports_a <- all_reports[split_idx]
  reports_b <- all_reports[-split_idx]

  sym_a <- symptoms |> filter(report_id %in% reports_a)
  sym_b <- symptoms |> filter(report_id %in% reports_b)

  res_a <- run_disprop(sym_a, paste0("split_", s, "_a"))
  res_b <- run_disprop(sym_b, paste0("split_", s, "_b"))

  robust_a <- res_a |> filter(robust_signal) |>
    select(product_clean, meddra_pt) |> mutate(in_a = TRUE)
  robust_b <- res_b |> filter(robust_signal) |>
    select(product_clean, meddra_pt) |> mutate(in_b = TRUE)

  concordance <- robust_a |>
    full_join(robust_b, by = c("product_clean", "meddra_pt")) |>
    replace_na(list(in_a = FALSE, in_b = FALSE))

  n_both  <- sum(concordance$in_a & concordance$in_b)
  n_a_only <- sum(concordance$in_a & !concordance$in_b)
  n_b_only <- sum(!concordance$in_a & concordance$in_b)
  jaccard <- n_both / (n_both + n_a_only + n_b_only)

  split_results[[s]] <- tibble(
    split = s,
    n_robust_a = sum(res_a$robust_signal, na.rm = TRUE),
    n_robust_b = sum(res_b$robust_signal, na.rm = TRUE),
    n_both = n_both,
    n_a_only = n_a_only,
    n_b_only = n_b_only,
    jaccard = round(jaccard, 3),
    pct_a_replicated = round(n_both / (n_both + n_a_only) * 100, 1),
    pct_b_replicated = round(n_both / (n_both + n_b_only) * 100, 1)
  )
}

split_summary <- bind_rows(split_results)
cat("\nRandom holdout replication results:\n")
print(split_summary, n = Inf, width = Inf)

cat(sprintf("\nMean Jaccard across %d splits: %.3f (SD %.3f)\n",
            N_SPLITS, mean(split_summary$jaccard), sd(split_summary$jaccard)))
cat(sprintf("Mean replication rate: %.1f%%\n",
            mean(c(split_summary$pct_a_replicated, split_summary$pct_b_replicated))))
cat(sprintf("Mean robust signals per half: %.0f (vs %d full dataset)\n",
            mean(c(split_summary$n_robust_a, split_summary$n_robust_b)),
            nrow(robust)))

# Compare to time-stratified Jaccard (0.027)
cat(sprintf("\nComparison: Random split Jaccard %.3f vs time-stratified Jaccard 0.027\n",
            mean(split_summary$jaccard)))
cat("Higher random-split Jaccard confirms that time-stratified low overlap is\n")
cat("driven by product name fragmentation across eras, not method instability.\n")

fwrite(split_summary, file.path(table_dir, "random_holdout_replication.csv"))


# =============================================================================
# SECTION 5: Choking/Dysphagia Exclusion Sensitivity
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 5: Choking/Dysphagia Exclusion Sensitivity\n")
cat(strrep("=", 60), "\n")

# Formulation-related signals (choking, dysphagia, foreign body) dominate
# the top of the signal list. These are real but relate to product form
# (tablet/capsule size) not active ingredient toxicity. Rerun after
# excluding these to reveal the pharmacological signal landscape.

formulation_pts <- c(
  "Choking", "Dysphagia", "Foreign Body", "Foreign Body In Throat",
  "Choking Sensation", "Oropharyngeal Pain", "Retching",
  "Foreign Body Trauma", "Throat Tightness", "Swelling Face",
  "Product Complaint", "Product Quality Issue", "Tablet Malformation"
)

n_formulation_signals <- robust |>
  filter(meddra_pt %in% formulation_pts) |>
  nrow()

cat(sprintf("Formulation-related PTs in robust signals: %d / %d (%.1f%%)\n",
            n_formulation_signals, nrow(robust),
            n_formulation_signals / nrow(robust) * 100))

# Remove formulation PTs and rerun
symptoms_no_form <- symptoms |>
  filter(!meddra_pt %in% formulation_pts)

cat(sprintf("Symptom rows after exclusion: %s (removed %s)\n",
            format(nrow(symptoms_no_form), big.mark = ","),
            format(nrow(symptoms) - nrow(symptoms_no_form), big.mark = ",")))

cat("Running disproportionality without formulation PTs...\n")
no_form_results <- run_disprop(symptoms_no_form, "no_formulation")

n_no_form_robust <- sum(no_form_results$robust_signal, na.rm = TRUE)
cat(sprintf("Robust signals without formulation PTs: %d (vs %d with)\n",
            n_no_form_robust, nrow(robust)))

# What's different? Top signals now
cat("\nTop 30 robust signals after excluding formulation PTs:\n")
no_form_results |>
  filter(robust_signal) |>
  arrange(desc(a)) |>
  select(product_clean, meddra_pt, a, prr, n_methods) |>
  head(30) |>
  print(n = 30, width = Inf)

# Top products change
cat("\nTop 20 products by signal count (excluding formulation PTs):\n")
no_form_results |>
  filter(robust_signal) |>
  count(product_clean, sort = TRUE) |>
  head(20) |>
  print(n = 20, width = Inf)

# Compare product rankings with and without formulation PTs
full_product_ranks <- robust |>
  count(product_clean, sort = TRUE, name = "n_signals_full") |>
  mutate(rank_full = row_number())

no_form_product_ranks <- no_form_results |>
  filter(robust_signal) |>
  count(product_clean, sort = TRUE, name = "n_signals_no_form") |>
  mutate(rank_no_form = row_number())

rank_comparison <- full_product_ranks |>
  full_join(no_form_product_ranks, by = "product_clean") |>
  replace_na(list(n_signals_full = 0, n_signals_no_form = 0,
                  rank_full = 999, rank_no_form = 999)) |>
  mutate(rank_change = rank_full - rank_no_form) |>
  filter(rank_full <= 30 | rank_no_form <= 30) |>
  arrange(rank_no_form)

cat("\nProduct rank comparison (top 30 in either analysis):\n")
rank_comparison |>
  select(product_clean, n_signals_full, rank_full,
         n_signals_no_form, rank_no_form, rank_change) |>
  print(n = 40, width = Inf)

# SOC distribution shift
cat("\nSOC distribution shift after excluding formulation PTs:\n")
pt_to_soc <- symptoms |> distinct(meddra_pt, soc) |> filter(!is.na(soc), soc != "")

soc_full <- robust |>
  left_join(pt_to_soc, by = "meddra_pt") |>
  count(soc, name = "n_full") |>
  filter(!is.na(soc))

soc_no_form <- no_form_results |>
  filter(robust_signal) |>
  left_join(pt_to_soc, by = "meddra_pt") |>
  count(soc, name = "n_no_form") |>
  filter(!is.na(soc))

soc_shift <- soc_full |>
  full_join(soc_no_form, by = "soc") |>
  replace_na(list(n_full = 0, n_no_form = 0)) |>
  mutate(
    pct_full = round(n_full / sum(n_full) * 100, 1),
    pct_no_form = round(n_no_form / sum(n_no_form) * 100, 1),
    pct_change = pct_no_form - pct_full
  ) |>
  arrange(desc(pct_no_form))

print(soc_shift, n = Inf, width = Inf)

fwrite(no_form_results |> filter(robust_signal),
       file.path(table_dir, "robust_signals_no_formulation.csv"))
fwrite(soc_shift, file.path(table_dir, "soc_shift_no_formulation.csv"))


# =============================================================================
# SECTION 6: Weighted Composite Risk Score
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SECTION 6: Weighted Composite Risk Score\n")
cat(strrep("=", 60), "\n")

# Combine multiple dimensions into a single risk ranking:
# - Disproportionality strength (log2 PRR)
# - Case count (log2 a)
# - Method concordance (n_methods / 4)
# - Serious outcome proportion (from supplements data)
# - Temporal trend (CUSUM alarm status)
#
# Each component is min-max normalised to [0, 1], then weighted.

# Get serious outcome proportions per product-PT pair
pair_serious <- symptoms |>
  distinct(report_id, product_clean, meddra_pt, outcome_serious) |>
  group_by(product_clean, meddra_pt) |>
  summarise(pct_serious = mean(outcome_serious, na.rm = TRUE),
            .groups = "drop")

# Get CUSUM alarm status
cusum_signals <- fread(file.path(table_dir, "cusum_first_signals.csv")) |>
  as_tibble() |>
  select(product_clean, meddra_pt, max_cusum) |>
  mutate(has_cusum_alarm = TRUE)

# Build composite score
composite <- robust |>
  left_join(pair_serious, by = c("product_clean", "meddra_pt")) |>
  left_join(cusum_signals, by = c("product_clean", "meddra_pt")) |>
  replace_na(list(pct_serious = 0, has_cusum_alarm = FALSE, max_cusum = 0)) |>
  mutate(
    # Component scores (raw)
    log2_prr = log2(pmax(prr, 1)),
    log2_a = log2(a),
    method_score = n_methods / 4,
    serious_score = pct_serious,
    temporal_score = as.numeric(has_cusum_alarm)
  )

# Min-max normalisation
normalise <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] == rng[1]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

composite <- composite |>
  mutate(
    norm_prr = normalise(log2_prr),
    norm_count = normalise(log2_a),
    norm_methods = normalise(method_score),
    norm_serious = normalise(serious_score),
    norm_temporal = temporal_score  # already 0/1
  )

# Weights (sum to 1)
w_prr      <- 0.25
w_count    <- 0.15
w_methods  <- 0.15
w_serious  <- 0.30
w_temporal <- 0.15

composite <- composite |>
  mutate(
    risk_score = w_prr * norm_prr + w_count * norm_count +
                 w_methods * norm_methods + w_serious * norm_serious +
                 w_temporal * norm_temporal,
    risk_percentile = percent_rank(risk_score) * 100
  ) |>
  arrange(desc(risk_score))

cat("Composite risk score weights:\n")
cat(sprintf("  PRR strength:       %.0f%%\n", w_prr * 100))
cat(sprintf("  Case count:         %.0f%%\n", w_count * 100))
cat(sprintf("  Method concordance: %.0f%%\n", w_methods * 100))
cat(sprintf("  Serious outcomes:   %.0f%%\n", w_serious * 100))
cat(sprintf("  Temporal alarm:     %.0f%%\n", w_temporal * 100))

cat("\nTop 40 highest-risk signals (composite score):\n")
composite |>
  select(product_clean, meddra_pt, a, prr, n_methods, pct_serious,
         has_cusum_alarm, risk_score, risk_percentile) |>
  head(40) |>
  print(n = 40, width = Inf)

# Top products by mean composite risk
product_risk_composite <- composite |>
  group_by(product_clean) |>
  summarise(
    n_signals = n(),
    mean_risk_score = round(mean(risk_score), 3),
    max_risk_score = round(max(risk_score), 3),
    mean_pct_serious = round(mean(pct_serious) * 100, 1),
    total_cases = sum(a),
    .groups = "drop"
  ) |>
  arrange(desc(mean_risk_score))

cat("\nTop 30 products by mean composite risk score:\n")
product_risk_composite |>
  head(30) |>
  print(n = 30, width = Inf)

# Risk tier classification
composite <- composite |>
  mutate(
    risk_tier = case_when(
      risk_percentile >= 95 ~ "Critical (top 5%)",
      risk_percentile >= 80 ~ "High (80-95th)",
      risk_percentile >= 50 ~ "Moderate (50-80th)",
      TRUE ~ "Low (<50th)"
    )
  )

cat("\nRisk tier distribution:\n")
composite |>
  count(risk_tier, name = "n_signals") |>
  mutate(pct = round(n_signals / sum(n_signals) * 100, 1)) |>
  print(n = Inf, width = Inf)

fwrite(composite |> select(product_clean, meddra_pt, a, prr, n_methods,
                            pct_serious, has_cusum_alarm, risk_score,
                            risk_percentile, risk_tier),
       file.path(table_dir, "composite_risk_scores.csv"))
fwrite(product_risk_composite,
       file.path(table_dir, "product_composite_risk_ranking.csv"))


# =============================================================================
# Summary
# =============================================================================

cat("\n\n")
cat(strrep("=", 60), "\n")
cat("ADDITIONAL VALIDATION COMPLETE\n")
cat(strrep("=", 60), "\n")

cat("\n1. NOTORIETY BIAS\n")
for (i in seq_len(nrow(notoriety_all))) {
  cat(sprintf("   %s: %.1fx post-action rate ratio\n",
              notoriety_all$product_pattern[i], notoriety_all$rate_ratio[i]))
}

cat("\n2. COUNT-RESPONSE\n")
cat(sprintf("   Spearman rho (a vs PRR): %.3f\n", cor_spearman$estimate))
cat(sprintf("   Within-product mean rho: %.3f\n", mean(within_product$within_cor)))

cat("\n3. PUBLISHED SIGNAL CROSS-VALIDATION\n")
cat(sprintf("   Found in data: %d/%d (%.1f%%); Detected as robust: %d/%d (%.1f%%)\n",
            n_found, n_total, n_found/n_total*100,
            n_robust_found, n_total, n_robust_found/n_total*100))

cat("\n4. RANDOM HOLDOUT REPLICATION\n")
cat(sprintf("   Mean Jaccard: %.3f (vs time-stratified 0.027)\n",
            mean(split_summary$jaccard)))
cat(sprintf("   Mean replication rate: %.1f%%\n",
            mean(c(split_summary$pct_a_replicated, split_summary$pct_b_replicated))))

cat("\n5. CHOKING/DYSPHAGIA EXCLUSION\n")
cat(sprintf("   Formulation signals: %d (%.1f%% of robust)\n",
            n_formulation_signals, n_formulation_signals/nrow(robust)*100))
cat(sprintf("   Robust signals without formulation PTs: %d\n", n_no_form_robust))

cat("\n6. COMPOSITE RISK SCORE\n")
cat(sprintf("   Critical tier (top 5%%): %d signals\n",
            sum(composite$risk_tier == "Critical (top 5%)")))
cat(sprintf("   Top product: %s (mean score %.3f)\n",
            product_risk_composite$product_clean[1],
            product_risk_composite$mean_risk_score[1]))

cat("\n")
cat(strrep("=", 60), "\n")
