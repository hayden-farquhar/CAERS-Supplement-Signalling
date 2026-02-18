# =============================================================================
# 09_visualisation.R
# Publication-ready figures for the manuscript
# =============================================================================

library(tidyverse)
library(data.table)
library(igraph)
library(scales)
library(patchwork)

# --- Configuration -----------------------------------------------------------

proc_dir    <- "data/processed"
table_dir   <- "outputs/tables"
fig_dir     <- "outputs/figures"
network_dir <- "outputs/networks"

# Theme for publication figures
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# --- 1. Load data ------------------------------------------------------------

cat("=== Loading data ===\n")
supps    <- fread(file.path(proc_dir, "supplements_cleaned.csv")) |> as_tibble()
symptoms <- fread(file.path(proc_dir, "symptoms_long.csv")) |> as_tibble()
robust   <- fread(file.path(table_dir, "robust_signals.csv")) |> as_tibble()
disp     <- fread(file.path(proc_dir, "disproportionality_results.csv")) |> as_tibble()

# --- 2. Figure 1: Reporting volume over time ---------------------------------

cat("\n=== Figure 1: Reporting volume over time ===\n")

annual <- fread(file.path(proc_dir, "annual_report_counts.csv")) |> as_tibble()

p1 <- ggplot(annual |> filter(year >= 2004), aes(x = year, y = n_reports)) +
  geom_col(fill = "#2171B5", alpha = 0.8) +
  geom_vline(xintercept = 2006.5, linetype = "dashed", colour = "red", linewidth = 0.5) +
  annotate("text", x = 2007.5, y = max(annual$n_reports) * 0.95,
           label = "Mandatory\nreporting", hjust = 0, size = 3, colour = "red") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(2004, 2025, 2)) +
  labs(x = "Year", y = "Number of suspect supplement reports",
       title = "Annual reporting volume of dietary supplement adverse events (CAERS)") +
  theme_pub

ggsave(file.path(fig_dir, "fig1_reporting_volume.png"), p1,
       width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig1_reporting_volume.pdf"), p1,
       width = 8, height = 5, bg = "white")
cat("Saved: fig1_reporting_volume\n")

# --- 3. Figure 2: Method concordance Venn/upset ------------------------------

cat("\n=== Figure 2: Signal concordance ===\n")

concordance_data <- tibble(
  category = c("PRR only", "ROR only", "GPS only", "BCPNN only",
               "2 methods", "3 methods", "All 4 methods"),
  count = c(
    sum(disp$n_methods == 1 & disp$prr_signal, na.rm = TRUE),
    sum(disp$n_methods == 1 & disp$ror_signal, na.rm = TRUE),
    sum(disp$n_methods == 1 & disp$gps_signal, na.rm = TRUE),
    sum(disp$n_methods == 1 & disp$bcpnn_signal, na.rm = TRUE),
    sum(disp$n_methods == 2, na.rm = TRUE),
    sum(disp$n_methods == 3, na.rm = TRUE),
    sum(disp$n_methods == 4, na.rm = TRUE)
  )
)

# Bar chart of method agreement
method_counts <- tibble(
  method = c("BCPNN\n(IC025 > 0)", "ROR\n(LCI > 1)",
             "PRR\n(PRR≥2, χ²≥4)", "GPS\n(EB05 ≥ 2)"),
  n_signals = c(
    sum(disp$bcpnn_signal, na.rm = TRUE),
    sum(disp$ror_signal, na.rm = TRUE),
    sum(disp$prr_signal, na.rm = TRUE),
    sum(disp$gps_signal, na.rm = TRUE)
  )
) |>
  mutate(method = fct_reorder(method, n_signals))

p2a <- ggplot(method_counts, aes(x = method, y = n_signals)) +
  geom_col(fill = "#2171B5", alpha = 0.8) +
  geom_text(aes(label = comma(n_signals)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of signals detected",
       title = "A. Signals by detection method") +
  theme_pub

# Agreement levels
agreement_data <- tibble(
  n_methods = 1:4,
  count = c(
    sum(disp$n_methods == 1, na.rm = TRUE),
    sum(disp$n_methods == 2, na.rm = TRUE),
    sum(disp$n_methods == 3, na.rm = TRUE),
    sum(disp$n_methods == 4, na.rm = TRUE)
  ),
  label = c("1 method", "2 methods", "3 methods\n(robust)", "All 4 methods")
)

p2b <- ggplot(agreement_data, aes(x = factor(n_methods), y = count)) +
  geom_col(aes(fill = factor(n_methods)), alpha = 0.8, show.legend = FALSE) +
  geom_text(aes(label = comma(count)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Number of methods detecting signal", y = "Number of product-PT pairs",
       title = "B. Multi-method concordance") +
  theme_pub

p2 <- p2a + p2b
ggsave(file.path(fig_dir, "fig2_method_concordance.png"), p2,
       width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig2_method_concordance.pdf"), p2,
       width = 12, height = 5, bg = "white")
cat("Saved: fig2_method_concordance\n")

# --- 4. Figure 3: Top signals by composite risk score -----------------------

cat("\n=== Figure 3: Top signals (composite risk score) ===\n")

# Load composite risk scores
composite <- fread(file.path(table_dir, "composite_risk_scores.csv")) |> as_tibble()

# Top 25 robust signals by composite risk score
top_signals <- composite |>
  arrange(desc(risk_score)) |>
  head(25) |>
  mutate(
    label = paste0(str_trunc(product_clean, 35), " \u2014 ", meddra_pt),
    label = fct_reorder(label, risk_score),
    pct_serious_pct = pct_serious * 100
  )

p3 <- ggplot(top_signals, aes(x = label, y = prr)) +
  geom_point(aes(size = a, colour = pct_serious_pct), alpha = 0.8) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  scale_size_continuous(name = "Case count", range = c(2, 10)) +
  scale_colour_gradient(low = "#FDB863", high = "#C6233C",
                        name = "% Serious outcomes") +
  scale_y_log10() +
  labs(x = NULL, y = "Proportional Reporting Ratio (log scale)",
       title = "Top 25 signals by composite risk score",
       subtitle = "Size = case count, colour = % serious outcomes") +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig3_top_signals.png"), p3,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig3_top_signals.pdf"), p3,
       width = 10, height = 8, bg = "white")
cat("Saved: fig3_top_signals\n")

# --- 5. Figure 4: Product category outcomes ----------------------------------

cat("\n=== Figure 4: Category outcomes ===\n")

cat_outcomes <- fread(file.path(table_dir, "category_outcome_proportions.csv")) |>
  as_tibble()

# Only if this file exists
if (nrow(cat_outcomes) > 0) {
  cat_long <- cat_outcomes |>
    select(product_category, pct_death, pct_life_threat, pct_hosp) |>
    pivot_longer(-product_category, names_to = "outcome", values_to = "pct") |>
    mutate(
      outcome = recode(outcome,
                       pct_death = "Death",
                       pct_life_threat = "Life-threatening",
                       pct_hosp = "Hospitalisation"),
      product_category = fct_reorder(product_category, pct,
                                      .fun = sum, .desc = FALSE)
    )

  p4 <- ggplot(cat_long, aes(x = product_category, y = pct, fill = outcome)) +
    geom_col(position = "stack", alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c("Death" = "#C6233C",
                                 "Life-threatening" = "#E6873E",
                                 "Hospitalisation" = "#4A90D9"),
                      name = "Outcome") +
    labs(x = NULL, y = "Percentage of reports",
         title = "Serious outcome proportions by product category") +
    theme_pub

  ggsave(file.path(fig_dir, "fig4_category_outcomes.png"), p4,
         width = 9, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "fig4_category_outcomes.pdf"), p4,
         width = 9, height = 6, bg = "white")
  cat("Saved: fig4_category_outcomes\n")
}

# --- 6. Figure 5: Demographic stratification ---------------------------------

cat("\n=== Figure 5: Demographic signal counts ===\n")

strat_results <- fread(file.path(proc_dir, "stratified_results.csv")) |> as_tibble()

if (nrow(strat_results) > 0) {
  strat_summary <- strat_results |>
    filter(robust_signal == TRUE) |>
    count(stratum, name = "robust_signals") |>
    mutate(stratum = fct_reorder(stratum, robust_signals))

  p5 <- ggplot(strat_summary, aes(x = stratum, y = robust_signals)) +
    geom_col(fill = "#2171B5", alpha = 0.8) +
    geom_text(aes(label = comma(robust_signals)), hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Number of robust signals",
         title = "Robust signals by demographic stratum") +
    theme_pub

  ggsave(file.path(fig_dir, "fig5_stratification.png"), p5,
         width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "fig5_stratification.pdf"), p5,
         width = 8, height = 6, bg = "white")
  cat("Saved: fig5_stratification\n")
}

# --- 7. Figure 6: Network visualisation --------------------------------------

cat("\n=== Figure 6: Product-PT network ===\n")

network_file <- file.path(network_dir, "product_pt_network.rds")
if (file.exists(network_file)) {
  g <- readRDS(network_file)

  # Subset to top connected nodes for readability
  top_nodes <- names(sort(degree(g), decreasing = TRUE))[1:80]
  g_sub <- induced_subgraph(g, top_nodes)

  # Colour by node type
  colours <- ifelse(V(g_sub)$type, "#E6873E", "#2171B5")
  sizes <- pmin(degree(g_sub) * 0.8 + 3, 15)

  png(file.path(fig_dir, "fig6_network.png"), width = 2400, height = 2400,
      res = 300, bg = "white")
  plot(g_sub,
       vertex.color = colours,
       vertex.size = sizes,
       vertex.label = ifelse(degree(g_sub) >= 5,
                             str_trunc(V(g_sub)$name, 20), NA),
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.color = "grey80",
       edge.width = 0.5,
       layout = layout_with_fr(g_sub),
       main = "Product-Adverse Event Network (Top 80 Nodes)")
  legend("bottomright",
         legend = c("Product", "MedDRA PT"),
         fill = c("#2171B5", "#E6873E"),
         cex = 0.8, bty = "n")
  dev.off()

  # PDF version
  pdf(file.path(fig_dir, "fig6_network.pdf"), width = 8, height = 8, bg = "white")
  set.seed(42)
  plot(g_sub,
       vertex.color = colours,
       vertex.size = sizes,
       vertex.label = ifelse(degree(g_sub) >= 5,
                             str_trunc(V(g_sub)$name, 20), NA),
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.color = "grey80",
       edge.width = 0.5,
       layout = layout_with_fr(g_sub),
       main = "Product-Adverse Event Network (Top 80 Nodes)")
  legend("bottomright",
         legend = c("Product", "MedDRA PT"),
         fill = c("#2171B5", "#E6873E"),
         cex = 0.8, bty = "n")
  dev.off()

  cat("Saved: fig6_network.png, fig6_network.pdf\n")
}

# --- 8. Figure 7: CUSUM example plots ---------------------------------------

cat("\n=== Figure 7: CUSUM example plots ===\n")

cusum_file <- file.path(proc_dir, "cusum_results.csv")
if (file.exists(cusum_file)) {
  cusum <- fread(cusum_file) |> as_tibble()

  # Select notable signals for CUSUM illustration
  example_pairs <- tribble(
    ~product_clean, ~meddra_pt,
    "KRATOM", "Death",
    "HYDROXYCUT REGULAR RAPID RELEASE CAPLETS", "Liver Function Test Abnormal"
  )

  cusum_examples <- cusum |>
    inner_join(example_pairs, by = c("product_clean", "meddra_pt"))

  if (nrow(cusum_examples) > 0) {
    cusum_examples <- cusum_examples |>
      mutate(label = paste0(str_trunc(product_clean, 30), " — ", meddra_pt))

    p7 <- ggplot(cusum_examples, aes(x = report_yq, y = cusum_pos, group = label)) +
      geom_line(colour = "#2171B5", linewidth = 0.8) +
      geom_hline(yintercept = 5, linetype = "dashed", colour = "red") +
      facet_wrap(~label, scales = "free_y", ncol = 1) +
      scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 8)]) +
      labs(x = "Quarter", y = "CUSUM statistic",
           title = "CUSUM control charts for selected signals") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

    ggsave(file.path(fig_dir, "fig7_cusum_examples.png"), p7,
           width = 10, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(fig_dir, "fig7_cusum_examples.pdf"), p7,
           width = 10, height = 8, bg = "white")
    cat("Saved: fig7_cusum_examples\n")
  }
}

# --- Summary -----------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("VISUALISATION COMPLETE\n")
cat(strrep("=", 60), "\n")
fig_files <- list.files(fig_dir, pattern = "\\.(png|pdf)$")
cat(sprintf("Figures generated: %d\n", length(fig_files)))
cat(paste(" -", fig_files), sep = "\n")
cat(strrep("=", 60), "\n")
