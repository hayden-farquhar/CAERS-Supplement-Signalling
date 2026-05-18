# Changelog

All notable changes to this reproducibility archive are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-05-18

First tagged release. Snapshot of the reproducibility archive at the point of
submitting a peer-reviewed manuscript revision.

### Added

- `outputs/figures/FIGURES.md` — mapping document linking the analysis-script
  output filenames to the manuscript figure numbers (Figures 1–8). Resolves
  the visible mismatch between the script-output naming (where two files
  share the `fig1_` prefix and the network/stratification numbering reflects
  development order rather than manuscript order) and the final manuscript
  figure numbering.

### Changed

- `README.md` — corrected validation analyses count from 22 to 24 (reflects
  the additional analyses performed during pre-submission revision: missing-
  data IPW sensitivity, expanded confounding-by-indication treatment, and
  product-name consolidation outcome).
- `README.md` — added a "Key Findings" bullet describing the variant-name
  consolidation result: 117 clusters covering 62.6% of signal-producing
  products reduces the 3,017 product-name-level signals to an estimated
  1,800–2,200 distinct product-level signals; targeted manual consolidation
  of top products retained 12 of 13 testable Table 2 signals (+6.4% net).

### Unchanged

- Analysis scripts (`scripts/`) — no functional changes; outputs are
  reproducible from the same code that generated the originally archived
  results.
- Output tables (`outputs/tables/`) — preserved exactly as produced by the
  analysis pipeline.
- Output figures (`outputs/figures/`) — figure PDF/PNG files retain their
  original analysis-script-output filenames; manuscript-figure-numbering
  mapping documented in the new `FIGURES.md`.
