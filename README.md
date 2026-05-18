# Disproportionality Analysis of Dietary Supplement Adverse Events in the FDA CAERS Database, 2004-2025

Reproducibility repository for:

> Farquhar H. Disproportionality Analysis of Dietary Supplement Adverse Events in the FDA CAERS Database, 2004-2025. *Preprint (Authorea)*. https://doi.org/10.22541/au.177383914.49860394/v1. Manuscript under consideration at a peer-reviewed journal.

## Overview

This repository contains the analysis code and output files for a computational pharmacovigilance study of dietary supplement adverse events reported to the FDA CFSAN Adverse Event Reporting System (CAERS). The study applies four disproportionality methods (PRR, ROR, GPS, BCPNN), CUSUM temporal detection, and demographic stratification to 48,840 unique adverse event reports.

## Key Findings

- 3,017 consensus product-name-level safety signals detected by 3 or more of 4 methods
- 2,146 signals detected by all four methods
- 148 critical-tier signals by composite risk scoring
- 451 temporally emerging signals via CUSUM control charts
- 24 validation and sensitivity analyses supporting signal robustness
- Variant-name consolidation of the most heavily fragmented products (117 clusters covering 62.6% of signal-producing products) reduces the headline 3,017 product-name-level signals to an estimated 1,800–2,200 distinct product-level signals; targeted manual consolidation of top products retained 12 of 13 testable Table 2 signals as robust (+6.4% net signal count from pooled statistical power)

## Repository Structure

```
├── scripts/                  # Analysis scripts (R)
│   ├── 01_data_acquisition.R       # Download and filter CAERS data
│   ├── 02_data_cleaning.R          # Standardise products, clean MedDRA terms
│   ├── 03_disproportionality_analysis.R  # PRR, ROR, GPS, BCPNN
│   ├── 04_demographic_stratification.R   # Age, sex, category subgroups
│   ├── 05_temporal_detection.R     # CUSUM control charts
│   ├── 07_product_risk_profiling.R # Logistic regression, network analysis
│   ├── 08_validation.R             # Core validation analyses
│   ├── 08b_extended_validation.R   # Extended sensitivity analyses
│   ├── 08c_additional_validation.R # Additional robustness checks
│   └── 09_visualisation.R         # Publication figures
├── outputs/
│   ├── figures/              # Publication figures (PDF + PNG)
│   └── tables/               # Signal catalogues and result tables (CSV)
├── data/
│   └── README.md             # Instructions for downloading CAERS data
├── LICENSE                   # MIT License
└── README.md                 # This file
```

## Reproducing the Analysis

### Prerequisites

- R >= 4.4.0
- Required R packages: `tidyverse`, `data.table`, `janitor`, `lubridate`, `slider`, `igraph`, `scales`

Install all dependencies:

```r
install.packages(c("tidyverse", "data.table", "janitor", "lubridate", "slider", "igraph", "scales"))
```

### Running the Pipeline

1. Download the CAERS data (see `data/README.md`)
2. Run scripts sequentially from the project root:

```bash
Rscript scripts/01_data_acquisition.R
Rscript scripts/02_data_cleaning.R
Rscript scripts/03_disproportionality_analysis.R
Rscript scripts/04_demographic_stratification.R
Rscript scripts/05_temporal_detection.R
Rscript scripts/07_product_risk_profiling.R
Rscript scripts/08_validation.R
Rscript scripts/08b_extended_validation.R
Rscript scripts/08c_additional_validation.R
Rscript scripts/09_visualisation.R
```

Scripts are numbered sequentially. Script 06 (LSTM temporal detection) was deferred; CUSUM (script 05) is the primary temporal method.

### Output

- **`outputs/figures/`** — 8 publication figures in PDF and PNG formats
- **`outputs/tables/`** — Full signal catalogues, validation results, and summary statistics in CSV format

## Data Source

This study uses the FDA CFSAN Adverse Event Reporting System (CAERS), renamed the Human Foods Complaint System (HFCS) in October 2024. The database is publicly available from [fda.gov](https://www.fda.gov/food/compliance-enforcement-food/cfsan-adverse-event-reporting-system-caers). The analysis used data downloaded on 30 September 2025.

## Citation

If you use this code or data, please cite:

```
Farquhar H. Disproportionality Analysis of Dietary Supplement Adverse Events
in the FDA CAERS Database, 2004-2025. Preprint (Authorea).
https://doi.org/10.22541/au.177383914.49860394/v1
Manuscript under consideration at a peer-reviewed journal.
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Hayden Farquhar MBBS MPHTM — hayden.farquhar@icloud.com
ORCID: [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)
