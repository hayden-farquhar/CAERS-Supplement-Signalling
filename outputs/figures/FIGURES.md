# Figure Output → Manuscript Figure Mapping

The figure output filenames in this directory follow the analysis-script output
convention (i.e., the names produced by `scripts/09_visualisation.R`). These
file names predate the final manuscript figure numbering. The mapping below
allows readers to match each output file to its manuscript figure number.

| Output file | Manuscript Figure | Subject |
|---|---|---|
| `fig1_flow_diagram.png` / `.pdf` | **Figure 1** | Study flow diagram (231,897 → 3,017 robust signals) |
| `fig1_reporting_volume.png` / `.pdf` | **Figure 2** | Annual reporting volume 2004–2024 with mandatory-reporting threshold |
| `fig2_method_concordance.png` / `.pdf` | **Figure 3** | Per-method signal counts (Panel A) and multi-method concordance distribution (Panel B) |
| `fig3_top_signals.png` / `.pdf` | **Figure 4** | Top 25 product–PT signals by composite risk score (bubble plot) |
| `fig4_category_outcomes.png` / `.pdf` | **Figure 5** | Serious outcome proportions by dietary supplement product category |
| `fig6_network.png` / `.pdf` | **Figure 6** | Top products by adverse event hub connectivity (degree centrality in the bipartite network) |
| `fig5_stratification.png` / `.pdf` | **Figure 7** | Robust signal counts by demographic and product-category stratum |
| `fig7_cusum_examples.png` / `.pdf` | **Figure 8** | Illustrative CUSUM control-chart trajectories (Hydroxycut LFT; Kratom – Death) |

## Notes

- The two filenames with the `fig1_` prefix (`fig1_flow_diagram.png` and
  `fig1_reporting_volume.png`) are an artefact of the visualisation script —
  the flow diagram and the reporting-volume chart are produced together in the
  same script block. They correspond to manuscript Figures 1 and 2 respectively.
- The output filenames `fig5_stratification.png` and `fig6_network.png` map to
  manuscript Figures 7 and 6 respectively (i.e., their order differs between the
  analysis-script output and the manuscript). This reflects how the analysis was
  developed (network analysis ran first, before demographic stratification was
  added) versus how the results are presented in the manuscript (categorical
  product-level findings before demographic stratification).
- Both PDF and PNG versions are provided for each figure. The PNG versions
  matching the manuscript figure numbering 1–8 are also bundled with each
  manuscript submission package separately.
