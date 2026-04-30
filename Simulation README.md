# Out-of-Control Simulation Study (EWMA and Variants)

This repository contains R code for an out-of-control simulation study comparing several EWMA-based control chart methods under non-normal data conditions and mean shifts.

## Objective

This simulation study code evaluates how different control charts respond to process shifts under skewed and heavy-tailed distributions. Performance is measured using **Average Run Length (ARL1)**.

## Control Chart Methods

The following methods are implemented:

- EWMA (Exponentially Weighted Moving Average)
- WV-EWMA (Weighted Variance EWMA)
- WSD-EWMA (Weighted Standard Deviation EWMA)
- SC-EWMA (Skewness-Corrected EWMA)

## Distribution Settings

Two distributional frameworks are used:

- Skew-Normal distributions (varying skewness levels)
- Skew-Student-t distributions (varying skewness and excess kurtosis)

## Simulation Design

The simulation structure includes:

- Phase I sample size: 100 observations
- Phase II monitoring using EWMA-type charts
- 10,000 Monte Carlo replications per configuration
- Mean shifts: 0.25σ, 0.50σ, 0.75σ, 1.00σ, 1.50σ
- Maximum run length: 2,000 observations

## Output Metrics

For each configuration, the simulation study code computes:

- ARL1 (Average Run Length under shift)
- SDRL1 (Standard deviation of run length)
- MRL1 (Median run length)

## Output Files

The script generates the following outputs:

### CSV Files
- `simulation_results_skewness_oc.csv`
- `simulation_results_kurtosis_oc.csv`
- `simulation_results_all_oc.csv`

### Plots
- ARL1 vs skewness (full scale and zoomed versions)
- ARL1 vs excess kurtosis (full scale and zoomed versions)
- Separate plots for different shift sizes and λ values

## Libraries Required

Install the required R packages:

```r
install.packages(c(
  "sn",
  "moments",
  "ggplot2",
  "dplyr",
  "tidyr",
  "param2moment",
  "parallel"
))
