# Dissertation Code: Individuals-Based Modified EWMA Control Charts for Non-Normality conditions

This repository contains the R code used in my dissertation. It includes both an empirical study based on financial data and a simulation study evaluating the performance of EWMA-type control charts under non-normal distributions.

## Repository Structure
- **empirical-study/**  
  Analysis of USD/EUR exchange rate data. Includes data transformation, distribution fitting, stationarity testing, and implementation of EWMA-based control charts.

- **simulation-study/**  
  Monte Carlo simulation study assessing the performance of different EWMA-type control charts using ARL and false alarm rate metrics.

## Methods Implemented

- EWMA control charts  
- Weighted Variance EWMA (WV-EWMA)  
- Weighted Standard Deviation EWMA (WSD-EWMA)  
- Skewness-Corrected EWMA (SC-EWMA)  

## Requirements

Run the following in R to install required packages:

```r
install.packages(c(
  "moments",
  "tseries",
  "fitdistrplus",
  "sn",
  "forecast",
  "MASS",
  "ggplot2",
  "dplyr",
  "tidyr",
  "parallel"
))
