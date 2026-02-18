# =============================================================================
# 05_temporal_detection.R
# CUSUM control charts for sequential signal detection on quarterly time series
# =============================================================================

library(tidyverse)
library(data.table)

# --- Configuration -----------------------------------------------------------

proc_dir  <- "data/processed"
table_dir <- "outputs/tables"
fig_dir   <- "outputs/figures"

# CUSUM parameters
CUSUM_K <- 0.5  # allowance (reference value)
CUSUM_H <- 5    # decision interval (threshold)

# --- 1. Load data ------------------------------------------------------------

cat("=== Loading data ===\n")
symptoms <- fread(file.path(proc_dir, "symptoms_long.csv")) |> as_tibble()
robust   <- fread(file.path(table_dir, "robust_signals.csv")) |> as_tibble()

cat(sprintf("Loaded %s product-symptom pairs\n", format(nrow(symptoms), big.mark = ",")))
cat(sprintf("Robust signals to track: %s\n", format(nrow(robust), big.mark = ",")))

# --- 2. Build quarterly time series ------------------------------------------

cat("\n=== Building quarterly time series ===\n")

# Get all year-quarter combinations
all_yq <- symptoms |>
  distinct(report_yq) |>
  filter(!is.na(report_yq)) |>
  arrange(report_yq)

# Quarterly counts for each product-PT pair
quarterly <- symptoms |>
  distinct(report_id, product_clean, meddra_pt, report_yq) |>
  filter(!is.na(report_yq)) |>
  count(product_clean, meddra_pt, report_yq, name = "count")

# Focus on robust signals for CUSUM
robust_pairs <- robust |>
  select(product_clean, meddra_pt) |>
  distinct()

quarterly_robust <- quarterly |>
  inner_join(robust_pairs, by = c("product_clean", "meddra_pt"))

# Fill in missing quarters with 0
quarterly_complete <- quarterly_robust |>
  complete(nesting(product_clean, meddra_pt), report_yq = all_yq$report_yq,
           fill = list(count = 0)) |>
  arrange(product_clean, meddra_pt, report_yq)

cat(sprintf("Quarterly time series for %s product-PT pairs\n",
            format(n_distinct(paste(quarterly_complete$product_clean,
                                    quarterly_complete$meddra_pt)), big.mark = ",")))

# --- 3. Apply CUSUM ----------------------------------------------------------

cat("\n=== Applying CUSUM control charts ===\n")

# CUSUM: detect upward shifts in reporting rate
# C+ = max(0, C+_prev + (x_t - mu_0 - k))
# Signal when C+ > h

cusum_analysis <- quarterly_complete |>
  group_by(product_clean, meddra_pt) |>
  arrange(report_yq) |>
  mutate(
    # Baseline: overall mean count per quarter for this pair
    mu_0 = mean(count),
    # Standardised observations
    z = (count - mu_0) / pmax(sd(count), 1),
    # One-sided upward CUSUM
    cusum_pos = accumulate(z, ~ max(0, .x + .y - CUSUM_K), .init = 0)[-1],
    # Signal indicator
    cusum_signal = cusum_pos > CUSUM_H,
    # First signal quarter
    quarter_idx = row_number()
  ) |>
  ungroup()

# Identify first CUSUM signal for each pair
first_signals <- cusum_analysis |>
  filter(cusum_signal) |>
  group_by(product_clean, meddra_pt) |>
  summarise(
    first_signal_yq = first(report_yq),
    first_signal_idx = first(quarter_idx),
    max_cusum = max(cusum_pos),
    .groups = "drop"
  )

cat(sprintf("Pairs with CUSUM signal: %s / %s (%.1f%%)\n",
            nrow(first_signals),
            n_distinct(paste(quarterly_complete$product_clean,
                             quarterly_complete$meddra_pt)),
            nrow(first_signals) /
              n_distinct(paste(quarterly_complete$product_clean,
                               quarterly_complete$meddra_pt)) * 100))

cat(sprintf("Total CUSUM alarm quarters: %s\n",
            format(sum(cusum_analysis$cusum_signal), big.mark = ",")))

# --- 4. Compare CUSUM timing vs aggregate ------------------------------------

cat("\n=== CUSUM timing analysis ===\n")

# For signals detected by CUSUM, when was the first alarm relative to
# when the signal would have been detectable by aggregate disproportionality?

# The aggregate analysis uses all data (retrospective).
# CUSUM detects prospectively — earlier detection = more value.

# Distribution of first signal quarters
cat("First CUSUM signal distribution by year-quarter:\n")
first_signals |>
  mutate(year = str_extract(first_signal_yq, "^\\d{4}") |> as.integer()) |>
  count(year, sort = FALSE) |>
  print(n = Inf)

# Top 20 earliest signals (most potentially useful for early warning)
cat("\nEarliest CUSUM signals (most useful for prospective detection):\n")
first_signals |>
  arrange(first_signal_yq) |>
  head(20) |>
  print(n = 20, width = Inf)

# Latest signals (most recent emerging signals)
cat("\nMost recent CUSUM signals (potentially emerging):\n")
first_signals |>
  arrange(desc(first_signal_yq)) |>
  head(20) |>
  print(n = 20, width = Inf)

# --- 5. Identify temporal anomalies (surge detection) ------------------------

cat("\n=== Temporal anomaly detection ===\n")

# Identify pairs where recent quarters show unusual spikes
# Use a rolling baseline (trailing 8 quarters) vs current quarter

anomalies <- cusum_analysis |>
  group_by(product_clean, meddra_pt) |>
  arrange(report_yq) |>
  mutate(
    # Rolling 8-quarter baseline mean and SD
    roll_mean = slider::slide_dbl(count, mean, .before = 8, .after = 0,
                                   .complete = TRUE),
    roll_sd = slider::slide_dbl(count, sd, .before = 8, .after = 0,
                                 .complete = TRUE),
    # Z-score relative to rolling baseline
    z_rolling = (count - lag(roll_mean)) / pmax(lag(roll_sd), 1),
    is_spike = z_rolling > 3  # > 3 SD above rolling baseline
  ) |>
  ungroup()

# Recent spikes (last 2 years)
recent_spikes <- anomalies |>
  filter(is_spike & str_detect(report_yq, "^202[3-5]")) |>
  select(product_clean, meddra_pt, report_yq, count, roll_mean, z_rolling) |>
  arrange(desc(z_rolling))

cat(sprintf("Recent spikes (2023-2025, z > 3): %d\n", nrow(recent_spikes)))
if (nrow(recent_spikes) > 0) {
  cat("\nTop 20 recent spikes:\n")
  recent_spikes |> head(20) |> print(n = 20, width = Inf)
}

# --- 6. Save results ---------------------------------------------------------

cat("\n=== Saving results ===\n")

fwrite(cusum_analysis, file.path(proc_dir, "cusum_results.csv"))
cat(sprintf("Saved: cusum_results.csv (%s rows)\n",
            format(nrow(cusum_analysis), big.mark = ",")))

fwrite(first_signals, file.path(table_dir, "cusum_first_signals.csv"))
cat(sprintf("Saved: cusum_first_signals.csv (%s signals)\n",
            format(nrow(first_signals), big.mark = ",")))

if (nrow(recent_spikes) > 0) {
  fwrite(recent_spikes, file.path(table_dir, "recent_temporal_spikes.csv"))
}

# --- Summary -----------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("PHASE 4 COMPLETE: Temporal Signal Detection (CUSUM)\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Pairs tracked:            %s\n",
            format(n_distinct(paste(quarterly_complete$product_clean,
                                    quarterly_complete$meddra_pt)), big.mark = ",")))
cat(sprintf("  Pairs with CUSUM signal:  %s\n",
            format(nrow(first_signals), big.mark = ",")))
cat(sprintf("  Recent temporal spikes:   %d\n", nrow(recent_spikes)))
cat(strrep("=", 60), "\n")
