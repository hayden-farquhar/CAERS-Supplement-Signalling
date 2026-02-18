# =============================================================================
# 07_product_risk_profiling.R
# Serious outcome modelling and network visualisation
# =============================================================================

library(tidyverse)
library(data.table)
library(igraph)

# --- Configuration -----------------------------------------------------------

proc_dir    <- "data/processed"
table_dir   <- "outputs/tables"
fig_dir     <- "outputs/figures"
network_dir <- "outputs/networks"

# --- 1. Load data ------------------------------------------------------------

cat("=== Loading data ===\n")
cases <- fread(file.path(proc_dir, "case_level_summary.csv")) |> as_tibble()
supps <- fread(file.path(proc_dir, "supplements_cleaned.csv")) |> as_tibble()
robust <- fread(file.path(table_dir, "robust_signals.csv")) |> as_tibble()

cat(sprintf("Cases: %s, Products: %s, Robust signals: %s\n",
            format(nrow(cases), big.mark = ","),
            format(nrow(supps), big.mark = ","),
            format(nrow(robust), big.mark = ",")))

# --- 2. Serious outcome proportions by product category ----------------------

cat("\n=== Serious outcome proportions by product category ===\n")

category_outcomes <- supps |>
  group_by(product_category) |>
  summarise(
    n_reports = n(),
    n_cases = n_distinct(report_id),
    n_death = sum(outcome_death, na.rm = TRUE),
    n_life_threat = sum(outcome_life_threat, na.rm = TRUE),
    n_hosp = sum(outcome_hospitalisation, na.rm = TRUE),
    n_serious = sum(outcome_serious, na.rm = TRUE),
    pct_death = round(n_death / n_cases * 100, 2),
    pct_life_threat = round(n_life_threat / n_cases * 100, 2),
    pct_hosp = round(n_hosp / n_cases * 100, 2),
    pct_serious = round(n_serious / n_cases * 100, 2),
    .groups = "drop"
  ) |>
  arrange(desc(pct_serious))

cat("Serious outcome proportions by product category:\n")
print(category_outcomes, n = Inf, width = Inf)

fwrite(category_outcomes, file.path(table_dir, "category_outcome_proportions.csv"))

# --- 3. Logistic regression for serious outcomes -----------------------------

cat("\n=== Logistic regression: serious outcome modelling ===\n")

# Prepare case-level data with product category
model_data <- supps |>
  distinct(report_id, .keep_all = TRUE) |>
  filter(product_category != "Other") |>  # exclude non-specific category
  mutate(
    outcome_serious = as.integer(outcome_serious),
    product_category = factor(product_category),
    age_group = factor(age_group, levels = c("40-59", "<18", "18-39", "60+", "Unknown")),
    sex_clean = factor(sex_clean, levels = c("Female", "Male", "Unknown"))
  )

# Fit model
model <- glm(outcome_serious ~ product_category + age_group + sex_clean,
             data = model_data, family = binomial)

cat("\nLogistic regression results:\n")
summary(model) |> print()

# Odds ratios with CIs
or_table <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
  filter(term != "(Intercept)") |>
  mutate(across(where(is.numeric), ~ round(., 3))) |>
  arrange(desc(estimate))

cat("\nOdds ratios for serious outcomes:\n")
print(or_table, n = Inf, width = Inf)

fwrite(or_table, file.path(table_dir, "logistic_regression_ors.csv"))

# --- 4. Network visualisation: product-AE co-occurrence ----------------------

cat("\n=== Building product-AE network ===\n")

# Use robust signals for the network
# Filter to top signals to keep network readable
top_network <- robust |>
  filter(a >= 10) |>
  select(product_clean, meddra_pt, a, prr, n_methods)

cat(sprintf("Network edges (product-PT pairs, a >= 10): %s\n", nrow(top_network)))

# Create igraph bipartite network
products_in_net <- unique(top_network$product_clean)
pts_in_net <- unique(top_network$meddra_pt)

# Create edge list
edges <- top_network |>
  select(from = product_clean, to = meddra_pt, weight = a)

g <- graph_from_data_frame(edges, directed = FALSE)

# Set vertex attributes
V(g)$type <- V(g)$name %in% pts_in_net  # TRUE for PTs, FALSE for products
V(g)$node_type <- ifelse(V(g)$type, "MedDRA PT", "Product")

# Set vertex size by degree
V(g)$size <- degree(g)

cat(sprintf("Network: %d nodes (%d products, %d PTs), %d edges\n",
            vcount(g), sum(!V(g)$type), sum(V(g)$type), ecount(g)))

# Community detection
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)
cat(sprintf("Detected %d communities\n", length(unique(V(g)$community))))

# Save network
saveRDS(g, file.path(network_dir, "product_pt_network.rds"))

# Network statistics
cat("\nTop 20 products by network degree (number of connected PTs):\n")
product_degree <- tibble(
  product = V(g)$name[!V(g)$type],
  degree = degree(g)[!V(g)$type]
) |> arrange(desc(degree))
print(product_degree |> head(20), n = 20)

cat("\nTop 20 MedDRA PTs by network degree (number of connected products):\n")
pt_degree <- tibble(
  meddra_pt = V(g)$name[V(g)$type],
  degree = degree(g)[V(g)$type]
) |> arrange(desc(degree))
print(pt_degree |> head(20), n = 20)

# --- 5. Product risk ranking ------------------------------------------------

cat("\n=== Product risk ranking ===\n")

# Combine signal strength with serious outcome proportion
product_risk <- supps |>
  group_by(product_clean) |>
  summarise(
    n_reports = n(),
    pct_serious = mean(outcome_serious, na.rm = TRUE) * 100,
    pct_death = mean(outcome_death, na.rm = TRUE) * 100,
    pct_hosp = mean(outcome_hospitalisation, na.rm = TRUE) * 100,
    product_category = first(product_category),
    .groups = "drop"
  ) |>
  left_join(
    robust |>
      count(product_clean, name = "n_robust_signals"),
    by = "product_clean"
  ) |>
  replace_na(list(n_robust_signals = 0)) |>
  filter(n_reports >= 10) |>  # minimum report threshold
  arrange(desc(n_robust_signals), desc(pct_serious))

cat("Top 30 highest-risk products (by robust signal count):\n")
product_risk |>
  head(30) |>
  select(product_clean, product_category, n_reports, n_robust_signals,
         pct_serious, pct_death) |>
  print(n = 30, width = Inf)

fwrite(product_risk, file.path(table_dir, "product_risk_ranking.csv"))

# --- Summary -----------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("PHASE 5 COMPLETE: Product Risk Profiling\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Product categories ranked:  %d\n", nrow(category_outcomes)))
cat(sprintf("  Network nodes:              %d\n", vcount(g)))
cat(sprintf("  Network edges:              %d\n", ecount(g)))
cat(sprintf("  Network communities:        %d\n", length(unique(V(g)$community))))
cat(sprintf("  Products risk-ranked:       %s\n",
            format(nrow(product_risk), big.mark = ",")))
cat(strrep("=", 60), "\n")
