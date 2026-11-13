# Individual-level uncertainty in PRS
Simulations on the individual-level uncertainty in Polygenic Risk Scores (PRS), and how this uncertainty affects PRS-based disease risk prediction and classification of individuals into high-risk groups.

Code accompanying my thesis: **The Impact of Individual-Level Uncertainty in Genetic Risk on Disease
Prediction**, Aalto University, November 2026.

## Overview
 
The simulations compare two ways of predicting disease risk under the liability threshold model (LTM):
 
- **Method P (point estimate):** risk from the PRS alone, liability | PRS ~ N(PRS, 1 − h²).
- **Method U (uncertainty-aware):** risk also accounts for the posterior variance of the genetic value, liability | GWAS data ~ N(PRS, 1 − h² + posterior variance).
  
Both methods use the same simulated target population and the same GWAS draw.

 ### Simulation model
 
1. **Genetic architecture:** allele frequencies and true causal effects for `M` SNPs, given heritability `h2`.
2. **Target population:** genotypes of `n.target` individuals; invariant SNPs are removed and genotypes standardized.
3. **Disease status:** liability threshold model with prevalence `K`.
4. **GWAS:** one draw of effect estimates for GWAS sample size `N`, converted to posterior means and variances per SNP (normal prior, normal likelihood).
5. **PRS and uncertainty:** the PRS is the posterior mean of the genetic value and its posterior variance is the individual-level uncertainty.
6. **Classification:** individuals are classified as high risk if their predicted risk exceeds `C_FP / (C_FP + C_FN)`, for a range of false-positive and false-negative costs.

## Repository structure
 
```
R/
  thesis-functions.R    simulation, PRS, risk and classification functions
  plot-functions.R      helper functions for plotting
analysis/               scripts, run in numeric order
results/                saved simulation output
figures/                generated figures (not tracked by git)
PGSsimulations.Rproj    RStudio project
```

### Scripts
 
| Script | Description | Output |
|---|---|---|
| `01_prs_uncertainty_expl.R` | Exploratory: PRS vs posterior variance, PRS by disease status | `figures/` |
| `02_uncertainty_vs_params.R` | Effect of `N`, `M` and `h2` on posterior variance | `figures/` |
| `03_variance_decomposition.R` | Between- vs within-individual genetic variance | `figures/` |
| `04_disease_risk_uncertainty.R` | Disease risk with 95% credible interval as a function of PRS | `figures/` |
| `05_run_risk_comparison.R` | Simulation: disease risks from Methods P and U | `results/risk_comparison_*.rds` |
| `05a_plot_risk_comparison.R` | Scatter plots of disease risks from Method P vs Method U | `figures/` |
| `06_run_classification_simulation.R` | Simulation: classification by Methods P and U | `results/classification_sim_*.rds` |
| `06a_plot_groupsize_cost.R` | Proportion classified high risk and cost per individual | `figures/` |
 
Scripts `05` and `06` are the simulation loops. Their saved output is in `results/`, so the plotting scripts `05a` and `06a` can be run without rerunning them.

## Usage
 
Requirements: R and the packages `here`, `ggplot2`, `dplyr`, `tidyr` and `patchwork`.
 
```r
install.packages(c("here", "ggplot2", "dplyr", "tidyr", "patchwork"))
```
 
```bash
git clone https://github.com/janemms/thesis-prs-uncertainty
cd thesis-prs-uncertainty
mkdir figures # ignored by git, but the scripts save plots here
```
 
Open `PGSsimulations.Rproj` in RStudio (or run R from the repository root, since paths are resolved with `here`), then run the scripts in `analysis/`.
 
