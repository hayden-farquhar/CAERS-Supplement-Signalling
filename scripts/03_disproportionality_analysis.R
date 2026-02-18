# =============================================================================
# 03_disproportionality_analysis.R
# Classical pharmacovigilance signal detection: PRR, ROR, GPS, BCPNN
# =============================================================================
#
# Input:  data/processed/symptoms_long.csv
#         data/processed/supplements_cleaned.csv
# Output: outputs/tables/signal_catalogue.csv
#         outputs/tables/concordance_matrix.csv
#         outputs/tables/robust_signals.csv
#         data/processed/disproportionality_results.csv
#
# Methods:
#   PRR  — Proportional Reporting Ratio (Evans et al. 2001)
#   ROR  — Reporting Odds Ratio (van Puijenbroek et al. 2002)
#   GPS  — Gamma-Poisson Shrinker (DuMouchel 1999)
#   BCPNN — Bayesian Confidence Propagation Neural Network (Bate et al. 1998)
#
# All use 2x2 contingency tables:
#   a = reports with product X AND event Y
#   b = reports with product X AND NOT event Y
#   c = reports without product X AND event Y
#   d = reports without product X AND NOT event Y
# =============================================================================

library(tidyverse)
library(data.table)

# --- Configuration -----------------------------------------------------------

proc_dir   <- "data/processed"
table_dir  <- "outputs/tables"

# Signal thresholds (standard pharmacovigilance criteria)
PRR_THRESHOLD   <- 2     # PRR >= 2
CHI2_THRESHOLD  <- 4     # chi-squared >= 4
MIN_N           <- 3     # minimum case count
ROR_CI_LOWER    <- 1     # lower 95% CI of ROR > 1
GPS_EB05        <- 2     # EB05 (5th percentile of posterior) >= 2
BCPNN_IC025     <- 0     # IC025 (lower 95% CI of IC) > 0

# --- 1. Load data ------------------------------------------------------------

cat("=== Loading data ===\n")
symptoms <- fread(file.path(proc_dir, "symptoms_long.csv")) |> as_tibble()
supps    <- fread(file.path(proc_dir, "supplements_cleaned.csv")) |> as_tibble()

cat(sprintf("Product-symptom pairs: %s\n", format(nrow(symptoms), big.mark = ",")))
cat(sprintf("Unique reports: %s\n", format(n_distinct(symptoms$report_id), big.mark = ",")))

# --- 2. Build 2x2 contingency tables ----------------------------------------

cat("\n=== Constructing 2x2 contingency tables ===\n")

# Total number of unique reports in the supplement dataset
N_total <- n_distinct(symptoms$report_id)

# Count reports per product (regardless of PT)
product_totals <- symptoms |>
  distinct(report_id, product_clean) |>
  count(product_clean, name = "n_product")

# Count reports per PT (regardless of product)
pt_totals <- symptoms |>
  distinct(report_id, meddra_pt) |>
  count(meddra_pt, name = "n_pt")

# Count reports per product-PT pair
pair_counts <- symptoms |>
  distinct(report_id, product_clean, meddra_pt) |>
  count(product_clean, meddra_pt, name = "a")

# Filter to pairs with minimum case count
pair_counts <- pair_counts |>
  filter(a >= MIN_N)

cat(sprintf("Product-PT pairs with a >= %d: %s\n",
            MIN_N, format(nrow(pair_counts), big.mark = ",")))

# Join totals to construct 2x2 table
contingency <- pair_counts |>
  left_join(product_totals, by = "product_clean") |>
  left_join(pt_totals, by = "meddra_pt") |>
  mutate(
    # Convert to double to prevent integer overflow in chi-squared
    a = as.double(a),
    n_product = as.double(n_product),
    n_pt = as.double(n_pt),
    b = n_product - a,       # product X, not event Y
    c = n_pt - a,            # not product X, event Y
    d = N_total - a - b - c, # not product X, not event Y
    N = as.double(N_total)
  )

# Sanity check
stopifnot(all(contingency$b >= 0))
stopifnot(all(contingency$c >= 0))
stopifnot(all(contingency$d >= 0))

cat(sprintf("Contingency tables constructed: %s\n", format(nrow(contingency), big.mark = ",")))

# --- 3. PRR (Proportional Reporting Ratio) -----------------------------------

cat("\n=== Computing PRR ===\n")

contingency <- contingency |>
  mutate(
    # PRR = (a/(a+b)) / (c/(c+d))
    prr = (a / (a + b)) / (c / (c + d)),
    prr_log = log2(prr),
    # Standard error of log(PRR)
    prr_se = sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d)),
    prr_lower = exp(log(prr) - 1.96 * prr_se),
    prr_upper = exp(log(prr) + 1.96 * prr_se),
    # Chi-squared test (Yates corrected)
    chi2 = ((abs(a * d - b * c) - N / 2)^2 * N) /
           ((a + b) * (c + d) * (a + c) * (b + d)),
    # Signal flag
    prr_signal = prr >= PRR_THRESHOLD & chi2 >= CHI2_THRESHOLD & a >= MIN_N
  )

cat(sprintf("PRR signals: %s / %s pairs (%.1f%%)\n",
            sum(contingency$prr_signal, na.rm = TRUE),
            nrow(contingency),
            mean(contingency$prr_signal, na.rm = TRUE) * 100))

# --- 4. ROR (Reporting Odds Ratio) -------------------------------------------

cat("\n=== Computing ROR ===\n")

contingency <- contingency |>
  mutate(
    # ROR = (a * d) / (b * c), with 0.5 continuity correction for zero cells
    a_c = a + 0.5, b_c = b + 0.5, c_c = c + 0.5, d_c = d + 0.5,
    ror = (a_c * d_c) / (b_c * c_c),
    ror_log = log2(ror),
    ror_se = sqrt(1/a_c + 1/b_c + 1/c_c + 1/d_c),
    ror_lower = exp(log(ror) - 1.96 * ror_se),
    ror_upper = exp(log(ror) + 1.96 * ror_se),
    # Signal flag
    ror_signal = ror_lower > ROR_CI_LOWER
  ) |>
  select(-a_c, -b_c, -c_c, -d_c)

cat(sprintf("ROR signals: %s / %s pairs (%.1f%%)\n",
            sum(contingency$ror_signal, na.rm = TRUE),
            nrow(contingency),
            mean(contingency$ror_signal, na.rm = TRUE) * 100))

# --- 5. GPS (Gamma-Poisson Shrinker / MGPS) ---------------------------------

cat("\n=== Computing GPS (MGPS) ===\n")

# The GPS uses an empirical Bayes approach with a mixture of two gamma priors.
# Expected count: E = (n_product * n_pt) / N_total
# The posterior is a mixture of two gamma distributions.
# We use the method of moments to estimate the mixture parameters.

contingency <- contingency |>
  mutate(
    expected = (n_product * n_pt) / N_total,
    rr = a / expected  # relative reporting ratio
  )

# Fit mixture of two gammas using method of moments on log2(RR)
# Following DuMouchel (1999), fit to all pairs including those with a < 3
# For simplicity, we use a parametric approximation

# Estimate mixture parameters from the data
log2rr <- log2(contingency$rr + 0.5)  # add 0.5 to avoid log(0)

# Simple two-component mixture estimation using quantiles
# Component 1: background (centred near 0)
# Component 2: signal (shifted right)
mu1 <- median(log2rr)
sd1 <- mad(log2rr)

# Use empirical Bayes shrinkage formula
# EB estimate = (a + alpha1) / (E + beta1)
# where alpha1 and beta1 are estimated from the prior

# Simplified MGPS: use gamma(alpha, beta) prior estimated from data
# Method: match mean and variance of observed counts
mean_rr <- mean(contingency$rr, na.rm = TRUE)
var_rr  <- var(contingency$rr, na.rm = TRUE)

# Gamma prior parameters
alpha_prior <- mean_rr^2 / var_rr
beta_prior  <- mean_rr / var_rr

contingency <- contingency |>
  mutate(
    # Posterior parameters under gamma-Poisson model
    alpha_post = alpha_prior + a,
    beta_post  = beta_prior + expected,
    # EB estimate (posterior mean on log2 scale)
    ebgm = log2(alpha_post / beta_post),
    # EB05 = 5th percentile of the posterior (log2 scale)
    # For gamma distribution, use qgamma
    eb05 = log2(qgamma(0.05, shape = alpha_post, rate = beta_post)),
    eb95 = log2(qgamma(0.95, shape = alpha_post, rate = beta_post)),
    # Signal flag (on natural scale: EB05 >= 2 means log2(EB05) >= 1)
    gps_signal = 2^eb05 >= GPS_EB05
  )

cat(sprintf("GPS signals (EB05 >= 2): %s / %s pairs (%.1f%%)\n",
            sum(contingency$gps_signal, na.rm = TRUE),
            nrow(contingency),
            mean(contingency$gps_signal, na.rm = TRUE) * 100))

# --- 6. BCPNN (Information Component) ----------------------------------------

cat("\n=== Computing BCPNN ===\n")

# IC = log2(observed / expected) with Bayesian smoothing
# Using beta-binomial model (Bate et al. 1998; Noren et al. 2006)
#
# Prior: p_ij ~ Beta(alpha_ij, beta_ij)
# Posterior: p_ij | data ~ Beta(alpha_ij + a, beta_ij + N - a)
#
# IC = log2(p_ij / (p_i * p_j)) where p's are posterior estimates

contingency <- contingency |>
  mutate(
    # Marginal probabilities
    p_product = n_product / N,
    p_pt      = n_pt / N,

    # Expected probability under independence
    p_expected = p_product * p_pt,

    # Observed probability
    p_observed = a / N,

    # IC with Bayesian smoothing (add 0.5 for continuity correction)
    ic = log2((a + 0.5) / ((n_product * n_pt / N) + 0.5)),

    # Variance of IC (approximate, following Bate et al.)
    # V(IC) ≈ 1/(a + 0.5) - 1/N
    ic_var = 1 / (a + 0.5) - 1 / N,
    ic_se  = sqrt(pmax(ic_var, 0)),

    # IC025 = lower 2.5% credible interval
    ic025 = ic - 1.96 * ic_se,
    ic975 = ic + 1.96 * ic_se,

    # Signal flag
    bcpnn_signal = ic025 > BCPNN_IC025
  )

cat(sprintf("BCPNN signals (IC025 > 0): %s / %s pairs (%.1f%%)\n",
            sum(contingency$bcpnn_signal, na.rm = TRUE),
            nrow(contingency),
            mean(contingency$bcpnn_signal, na.rm = TRUE) * 100))

# --- 7. Concordance and robust signals --------------------------------------

cat("\n=== Signal concordance ===\n")

contingency <- contingency |>
  mutate(
    n_methods = as.integer(prr_signal) + as.integer(ror_signal) +
                as.integer(gps_signal) + as.integer(bcpnn_signal),
    robust_signal = n_methods >= 3
  )

# Concordance matrix
concordance <- contingency |>
  summarise(
    total_pairs    = n(),
    prr_signals    = sum(prr_signal),
    ror_signals    = sum(ror_signal),
    gps_signals    = sum(gps_signal),
    bcpnn_signals  = sum(bcpnn_signal),
    any_1_method   = sum(n_methods >= 1),
    any_2_methods  = sum(n_methods >= 2),
    any_3_methods  = sum(n_methods >= 3),
    all_4_methods  = sum(n_methods == 4)
  )

cat("Concordance summary:\n")
concordance |> pivot_longer(everything()) |> print(n = Inf)

# Method overlap matrix
cat("\nPairwise method overlap:\n")
methods <- c("prr_signal", "ror_signal", "gps_signal", "bcpnn_signal")
overlap <- matrix(0, 4, 4, dimnames = list(
  c("PRR", "ROR", "GPS", "BCPNN"),
  c("PRR", "ROR", "GPS", "BCPNN")
))
for (i in 1:4) {
  for (j in 1:4) {
    overlap[i, j] <- sum(contingency[[methods[i]]] & contingency[[methods[j]]])
  }
}
print(overlap)

# --- 8. Generate signal catalogue --------------------------------------------

cat("\n=== Generating signal catalogue ===\n")

# Robust signals (detected by >= 3 methods)
robust <- contingency |>
  filter(robust_signal) |>
  select(product_clean, meddra_pt, a, expected,
         prr, prr_lower, prr_upper, chi2, prr_signal,
         ror, ror_lower, ror_upper, ror_signal,
         ebgm, eb05, eb95, gps_signal,
         ic, ic025, ic975, bcpnn_signal,
         n_methods) |>
  arrange(desc(n_methods), desc(a))

cat(sprintf("Robust signals (>= 3 methods): %s\n", nrow(robust)))

cat("\nTop 30 robust signals by case count:\n")
robust |>
  select(product_clean, meddra_pt, a, prr, ror, ebgm, ic, n_methods) |>
  head(30) |>
  print(n = 30, width = Inf)

# All signals (any method)
all_signals <- contingency |>
  filter(n_methods >= 1) |>
  select(product_clean, meddra_pt, a, expected,
         prr, prr_lower, prr_upper, chi2, prr_signal,
         ror, ror_lower, ror_upper, ror_signal,
         ebgm, eb05, eb95, gps_signal,
         ic, ic025, ic975, bcpnn_signal,
         n_methods, robust_signal) |>
  arrange(desc(n_methods), desc(a))

cat(sprintf("Total signals (any method): %s\n", nrow(all_signals)))

# --- 9. Summary by product category -----------------------------------------

cat("\n=== Signals by product (top 20) ===\n")
robust |>
  count(product_clean, sort = TRUE, name = "n_signals") |>
  head(20) |>
  print(n = 20)

cat("\n=== Signals by MedDRA PT (top 20) ===\n")
robust |>
  count(meddra_pt, sort = TRUE, name = "n_signals") |>
  head(20) |>
  print(n = 20)

# --- 10. Save results --------------------------------------------------------

cat("\n=== Saving results ===\n")

# Full disproportionality results
fwrite(contingency, file.path(proc_dir, "disproportionality_results.csv"))
cat(sprintf("Saved: disproportionality_results.csv (%s pairs)\n",
            format(nrow(contingency), big.mark = ",")))

# Signal catalogue (all signals)
fwrite(all_signals, file.path(table_dir, "signal_catalogue.csv"))
cat(sprintf("Saved: signal_catalogue.csv (%s signals)\n",
            format(nrow(all_signals), big.mark = ",")))

# Robust signals only
fwrite(robust, file.path(table_dir, "robust_signals.csv"))
cat(sprintf("Saved: robust_signals.csv (%s signals)\n",
            format(nrow(robust), big.mark = ",")))

# Concordance
fwrite(as.data.frame(overlap), file.path(table_dir, "concordance_matrix.csv"),
       row.names = TRUE)

# --- Summary -----------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("PHASE 2 COMPLETE: Disproportionality Analysis\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Total product-PT pairs analysed: %s\n",
            format(nrow(contingency), big.mark = ",")))
cat(sprintf("  PRR signals:                     %s\n",
            format(sum(contingency$prr_signal, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  ROR signals:                     %s\n",
            format(sum(contingency$ror_signal, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  GPS signals:                     %s\n",
            format(sum(contingency$gps_signal, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  BCPNN signals:                   %s\n",
            format(sum(contingency$bcpnn_signal, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Robust signals (>= 3 methods):   %s\n",
            format(nrow(robust), big.mark = ",")))
cat(sprintf("  Signals by all 4 methods:        %s\n",
            format(sum(contingency$n_methods == 4, na.rm = TRUE), big.mark = ",")))
cat(strrep("=", 60), "\n")
cat("\nNext: Run 04_demographic_stratification.R\n")
