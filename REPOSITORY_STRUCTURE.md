# File Organization for GitHub Repository

## Susceptibility Distribution Misspecification in Epidemic Forecasting

This document outlines the recommended organisation of the GitHub repository accompanying the manuscript **"An exploration into how susceptibility distribution misspecifications impact epidemic forecasting"**.

## Repository Structure

```text
susceptibility-distribution-misspecification/
├── README.md                              # Repository documentation
├── LICENSE.txt                            # MIT License
├── .gitignore                             # Git ignore rules for R projects
├── install_packages.R                     # Script to install required packages
├── 0_run_all.R                            # Driver script
├── R/                                     # Core functions
│   └── 1_helpers.R                        # Model, discretisation, likelihood, and fitting functions
├── scripts/                               # Analysis scripts
│   ├── 2_single_epidemic_misspec_hom.R
│   ├── 3_two_epidemics_misspec_hom.R
│   ├── 4_distribution_level_validation.R
│   ├── 5_single_epidemic_CV_sweep.R
│   ├── 6_two_epidemics_CV_sweep.R
│   └── 7_gamma_reduced_trajectory_validation.R
├── results/                               # Generated numerical results
│   ├── inference/
│   ├── cv_sensitivity/
│   ├── distribution_validation/
│   └── trajectory_validation/
├── figures/                               # Generated plots and visualisations
│   ├── inference/
│   ├── cv_sensitivity/
│   ├── distribution_validation/
│   └── trajectory_validation/
└── paper/                                 # Manuscript files
    ├── misspec_main.pdf
    └── misspec_p2_sm.pdf
```

## File Mapping from Working Directory to Repository

### Core Function Files

- `Generalfun_mispec_distribution.R` → `R/1_helpers.R`
- `MLE_functions_paper.R` → `R/1_helpers.R`

These two working files are combined into one helper file so that all analysis scripts call a single source file.

### Analysis Scripts

- `Single_Epidemic_misspec_hom.R` → `scripts/2_single_epidemic_misspec_hom.R`
- `Two_Epidemics_misspec_hom.R` → `scripts/3_two_epidemics_misspec_hom.R`
- `distribution_level_analysis_LA_SM.R` → `scripts/4_distribution_level_validation.R`
- `Single_Epidemic_CV_sweep.R` → `scripts/5_single_epidemic_CV_sweep.R`
- `Two_Epidemics_CV_sweep.R` → `scripts/6_two_epidemics_CV_sweep.R`
- `Match_hetsus_vs_reduced_ch5.R` → `scripts/7_gamma_reduced_trajectory_validation.R`

### Manuscript Files

- `misspec_main.pdf` → `paper/misspec_main.pdf`
- `misspec_p2_sm.pdf` → `paper/misspec_p2_sm.pdf`

## File Content Guidelines

### `R/1_helpers.R`

This file should contain only reusable functions. It should not run the full simulation analysis automatically when sourced.

### Analysis scripts

Each analysis script should contain:

1. A header explaining the purpose of the script.
2. Required library calls.
3. A source statement for `R/1_helpers.R` where needed.
4. Parameter definitions.
5. Simulation, fitting, analysis, and export sections.
6. Clear output filenames.

### Output directories

The `results/` and `figures/` directories are included for organisation, but many scripts currently write output files to the working directory. This is intentional to preserve the original analysis logic. Outputs can be moved into subfolders after a run if desired.

## Recommended Running Order

For the baseline paper results:

```r
source("scripts/2_single_epidemic_misspec_hom.R")
source("scripts/3_two_epidemics_misspec_hom.R")
```

For full supporting analyses:

```r
source("scripts/4_distribution_level_validation.R")
source("scripts/5_single_epidemic_CV_sweep.R")
source("scripts/6_two_epidemics_CV_sweep.R")
source("scripts/7_gamma_reduced_trajectory_validation.R")
```

Or use:

```bash
RUN_DISTRIBUTION_VALIDATION=true RUN_TRAJECTORY_VALIDATION=true RUN_CV_SWEEPS=true Rscript 0_run_all.R
```
