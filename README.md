# Susceptibility Distribution Misspecification in Epidemic Forecasting

This repository contains the R code and reproducibility materials for the manuscript **"An exploration into how susceptibility distribution misspecifications impact epidemic forecasting"** by Ibrahim Mohammed, Chris Robertson, and M. Gabriela M. Gomes.

## Overview

This study investigates how assumptions about the distribution of individual susceptibility affect epidemic inference and forecasting in heterogeneous SEIR models. The analysis compares gamma and lognormal susceptibility distributions when both are calibrated to have the same mean and coefficient of variation.

The repository focuses on four connected questions:

1. **Distributional shape**: How gamma and lognormal susceptibility distributions differ when their first two moments are matched.
2. **Discretisation**: How continuous susceptibility distributions are represented by finite susceptibility classes using log-affine and split-multiplicative calibration.
3. **Misspecification in inference**: How fitting the wrong susceptibility family affects maximum-likelihood estimates of \(R_0\), \(\nu\), \(t_0\), and \(c^*\).
4. **Forecasting impact**: How distributional misspecification affects future epidemic trajectories compared with the larger error introduced by ignoring heterogeneity altogether.

The main analyses compare correctly specified, misspecified, and homogeneous models under single-epidemic and two-epidemic fitting designs.

## Repository Structure

```text
susceptibility-distribution-misspecification/
├── README.md                              # Repository documentation
├── LICENSE.txt                            # MIT License
├── .gitignore                             # Git ignore rules
├── install_packages.R                     # Script to install required R packages
├── 0_run_all.R                            # Driver script for baseline and optional analyses
├── R/                                     # Core model and estimation functions
│   └── 1_helpers.R                        # Discretisation, ODE, simulation, likelihood, and fitting functions
├── scripts/                               # Analysis scripts
│   ├── 2_single_epidemic_misspec_hom.R     # Single-epidemic misspecification analysis
│   ├── 3_two_epidemics_misspec_hom.R       # Two-epidemic misspecification analysis
│   ├── 4_distribution_level_validation.R  # Distribution/discretisation validation
│   ├── 5_single_epidemic_CV_sweep.R        # Single-epidemic CV-sensitivity and forecast analysis
│   ├── 6_two_epidemics_CV_sweep.R          # Two-epidemic CV-sensitivity and forecast analysis
│   └── 7_gamma_reduced_trajectory_validation.R  # Gamma reduced-model trajectory validation
├── results/                               # Generated numerical outputs
│   ├── inference/
│   ├── cv_sensitivity/
│   ├── distribution_validation/
│   └── trajectory_validation/
├── figures/                               # Generated plots and visualisations
│   ├── inference/
│   ├── cv_sensitivity/
│   ├── distribution_validation/
│   └── trajectory_validation/
└── paper/                                 # Manuscript files included for reference
    ├── misspec_main.pdf                   # Main manuscript
    └── misspec_p2_sm.pdf                  # Supporting information
```

## Key Models and Methods

### Heterogeneous SEIR model

The code implements a heterogeneous SEIR of Gomes et al (2022) in which susceptibility varies across individuals. The continuous susceptibility distribution is discretised into finite groups with weights \(q_i\) and representative susceptibilities \(x_i\).

### Susceptibility distributions

The analyses compare:

- **Gamma susceptibility**, parameterised to have mean 1 and coefficient of variation \(\nu\).
- **Lognormal susceptibility**, also parameterised to have mean 1 and coefficient of variation \(\nu\).
- **Homogeneous susceptibility**, used as a reference model with no susceptibility variation.

### Discretisation methods

The supporting information analyses two calibration methods:

- **Log-affine calibration**, which transforms representatives using \(x_i' = \exp(A + B\log x_i)\).
- **Split-multiplicative calibration**, which rescales lower and upper halves of the discretised distribution.

### Statistical inference

The inference scripts use:

- Poisson observation models for daily incidence.
- Maximum likelihood estimation.
- Hessian-based Wald intervals.
- Parameter-bias, confidence-interval width, coverage, condition-number, and correlation diagnostics.

## Requirements

To run this code, the following are required:

- **R version**: 4.0.0 or higher recommended.
- **Required packages**:
  - `deSolve`
  - `tidyverse`
  - `dplyr`
  - `tidyr`
  - `ggplot2`
  - `MASS`
  - `gridExtra`
  - `grid`
  - `cowplot`
  - `scales`

Install the required packages by running:

```r
source("install_packages.R")
```

## Usage

### Quick start

1. Clone the repository:

```bash
git clone https://github.com/<your-github-username>/susceptibility-misspecification-forecasting.git
cd susceptibility-misspecification-forecasting
```

2. Install R dependencies:

```r
source("install_packages.R")
```

3. Run the two baseline analyses:

```r
source("scripts/2_single_epidemic_misspec_hom.R")
source("scripts/3_two_epidemics_misspec_hom.R")
```

Alternatively, run the driver script from a terminal:

```bash
Rscript 0_run_all.R
```

By default, the driver script runs the single-epidemic and two-epidemic baseline analyses. The distribution-level validation and CV-sensitivity sweeps are optional because they can be computationally expensive.

To run all analyses:

```bash
RUN_DISTRIBUTION_VALIDATION=true RUN_TRAJECTORY_VALIDATION=true RUN_CV_SWEEPS=true Rscript 0_run_all.R
```

## Script Descriptions

### `R/1_helpers.R`

Contains the shared functions used by the analysis scripts, including:

- logit and inverse-logit transforms;
- initial-condition construction;
- heterogeneous SEIR ODE system;
- gamma and lognormal discretisation functions;
- simulation functions;
- Poisson likelihood functions;
- maximum-likelihood fitting functions for gamma, lognormal, and homogeneous models.

### `scripts/2_single_epidemic_misspec_hom.R`

Runs the baseline single-epidemic misspecification analysis. It simulates incidence data under a chosen truth distribution, fits gamma, lognormal, and homogeneous models, and exports parameter summaries, confidence-interval diagnostics, and fitted-model outputs.

### `scripts/3_two_epidemics_misspec_hom.R`

Runs the corresponding two-epidemic analysis. It simulates two epidemics with different initial conditions and fits a shared parameter vector jointly across both epidemics.

### `scripts/4_distribution_level_validation.R`

Reproduces the distribution-level validation of the discretisation scheme. It evaluates gamma and lognormal discretisations under log-affine and split-multiplicative calibration.

### `scripts/5_single_epidemic_CV_sweep.R`

Performs the single-epidemic CV-sensitivity and forecast analysis across \(\nu \in \{0.5, 1, \sqrt{2}, 2\}\). It computes forecast-window differences in peak height, peak timing, and final size.

### `scripts/6_two_epidemics_CV_sweep.R`

Performs the two-epidemic version of the CV-sensitivity and forecast analysis. It jointly fits two epidemics and evaluates forecast metrics for each epidemic.

### `scripts/7_gamma_reduced_trajectory_validation.R`

Validates the discretised gamma heterogeneous SEIR model against the analytical reduced gamma SEIR model. It reproduces the supporting-information check that the discretised gamma system preserves the reduced-model trajectory, especially once the number of susceptibility groups is moderately large.

## Main Outputs

Typical output files include:

- `combined_valid_single_epi_sim_<truth>.csv`
- `combined_valid_two_epi_sim_<truth>.csv`
- `summary_single_epi_sim_<truth>.csv`
- `summary_2epi_sim_<truth>.csv`
- `peak_diff_single_sim_<truth>_allCV.csv`
- `peak_diff_two_epi_sim_<truth>_allCV.csv`
- `CV_sensitivity_single_epidemic_truth_<truth>.csv`
- `CV_sensitivity_two_epidemic_truth_<truth>.csv`
- `gamma_reduced_trajectory_validation_summary.csv`
- `gamma_reduced_trajectory_validation_trajectories.csv`
- covariance objects saved as `.rds` files.

Output filenames depend on the value of `SIMULATE_WITH` inside each script.

##  Some Findings Reproduced by the Code

1. Gamma and lognormal susceptibility distributions can produce distinct epidemic trajectories even when they share the same mean and coefficient of variation.
2. These differences are small when heterogeneity is low, but become more visible when \(\nu\) is moderate or high.
3. Single-epidemic inference can mask susceptibility-family misspecification because changes in heterogeneity and intervention parameters compensate for each other.
4. Jointly fitting two epidemics can reduce this compensation and expose bias under misspecification.
5. The forecasting cost of using the wrong heterogeneous family is generally smaller than the cost of ignoring heterogeneity entirely.

## File Mapping

Original working files are organised as follows:

- `Generalfun_mispec_distribution.R` and `MLE_functions_paper.R` → `R/1_helpers.R`
- `Single_Epidemic_misspec_hom.R` → `scripts/2_single_epidemic_misspec_hom.R`
- `Two_Epidemics_misspec_hom.R` → `scripts/3_two_epidemics_misspec_hom.R`
- `distribution_level_analysis_LA_SM.R` → `scripts/4_distribution_level_validation.R`
- `Single_Epidemic_CV_sweep.R` → `scripts/5_single_epidemic_CV_sweep.R`
- `Two_Epidemics_CV_sweep.R` → `scripts/6_two_epidemics_CV_sweep.R`
- `Match_hetsus_vs_reduced_ch5.R` → `scripts/7_gamma_reduced_trajectory_validation.R`
- `misspec_main.pdf` → `paper/misspec_main.pdf`
- `misspec_p2_sm.pdf` → `paper/misspec_p2_sm.pdf`

## Reproducibility Notes

The scripts set fixed random seeds where simulations are generated. Full runs with 200 replicates, especially the CV sweeps, may take several hours depending on hardware. For quick test runs, edit the relevant script and set `QUICK_TEST <- TRUE` where available.

A `sessionInfo.txt` file is written automatically when `0_run_all.R` is used.

## Citation

If you use this code or build on the analysis, please cite:

```bibtex
@article{mohammed2026misspecification,
  title={An exploration into how susceptibility distribution misspecifications impact epidemic forecasting},
  author={Mohammed, Ibrahim and Robertson, Chris and Gomes, M. Gabriela M.},
  journal={Submitted manuscript},
  year={2026}
}
```

## License

This project is licensed under the MIT License. See [`LICENSE.txt`](LICENSE.txt) for details.

## Authors and Affiliations

- **Ibrahim Mohammed**
  - Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK
  - Department of Mathematical Sciences, Abubakar Tafawa Balewa University, Bauchi, Nigeria

- **Chris Robertson**
  - Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK
  - Public Health Scotland, Glasgow, UK

- **M. Gabriela M. Gomes**
  - Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK
  - Centre for Mathematics and Applications (NOVA MATH), NOVA School of Science and Technology, Caparica, Portugal

## Support

For questions about the code or methods, please open an issue on the repository or contact the authors through their institutional affiliations.
