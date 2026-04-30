# ===============================================================================
# EMPIRICAL STUDY - USD/EUR EXCHANGE RATE MONITORING
# Dissertation: Comparing WV-EWMA, WSD-EWMA, SC-EWMA
# STAGE 1: Real Data Application with ARL and FAR Computation
# ===============================================================================

# Clear workspace
rm(list = ls())

# ==== Load Required Libraries ====
library(tseries)
library(forecast)
library(rugarch)
library(MASS)
library(nortest)
library(moments)
library(fitdistrplus)
library(sn)           # For skewed distributions

# ===============================================================================
# SECTION 1: LOAD AND PREPARE DATA (LOG-RETURNS)
# ===============================================================================

# Load data
merged_data <- read.csv("merged_data.csv")

# Extract USD/EUR spot rates
usd_eur_spot <- na.omit(merged_data$USD_EUR)

cat("==========================================================\n")
cat("USD/EUR Exchange Rate Monitoring - Empirical Study\n")
cat("==========================================================\n")
cat(paste("Data points (spot rates):", length(usd_eur_spot), "\n"))
cat(paste("Date range: 2015-01-02 to", tail(merged_data$DATE, 1), "\n\n"))

# ===============================================================================
# STEP 1: COMPUTE LOG-RETURNS
# ===============================================================================
# Compute LOG-RETURNS (not raw prices) 
# Financial returns are typically log-returns for nice statistical properties
usd_eur_returns <- diff(log(usd_eur_spot))
usd_eur_returns <- usd_eur_returns[is.finite(usd_eur_returns)]

cat("==========================================================\n")
cat("STEP 1: Descriptive Statistics (Log-Returns)\n")
cat("==========================================================\n")

# ===============================================================================
# BASIC STATISTICS
# ===============================================================================

mu_empirical <- mean(usd_eur_returns)
sigma_empirical <- sd(usd_eur_returns)

# Moving range estimator (SPC-consistent)
MR <- abs(diff(usd_eur_returns))
MRbar <- mean(MR, na.rm = TRUE)
d2 <- 1.128
sigma_mr <- MRbar / d2

# Skewness and kurtosis
gamma1_empirical <- skewness(usd_eur_returns)

kurt_raw <- kurtosis(usd_eur_returns)
kurt_excess <- kurt_raw - 3

# Sample size
n_obs <- length(usd_eur_returns)

# ===============================================================================
# PRINT RESULTS
# ===============================================================================

cat(paste("Mean (μ₀):", round(mu_empirical, 6), "\n"))
cat(paste("Std Dev (classical):", round(sigma_empirical, 6), "\n"))
cat(paste("Std Dev (MR-based):", round(sigma_mr, 6), "\n"))
cat(paste("Skewness (γ₁):", round(gamma1_empirical, 3), "\n"))
cat(paste("Excess Kurtosis:", round(kurt_excess, 3), "\n"))
cat(paste("Sample size (n):", n_obs, "\n\n"))

# ===============================================================================
# INTERPRETATION (IMPORTANT FOR THESIS)
# ===============================================================================

cat("Interpretation:\n")

# Skewness
if(abs(gamma1_empirical) < 0.1){
  cat("• Distribution is approximately symmetric\n")
} else if(gamma1_empirical > 0){
  cat("• Right-skewed distribution (positive skewness)\n")
} else {
  cat("• Left-skewed distribution (negative skewness)\n")
}

# Kurtosis
if(kurt_excess > 0){
  cat("• Heavy-tailed behaviour detected (leptokurtic)\n")
} else {
  cat("• No heavy tails (platykurtic or normal-like)\n")
}

# Variance comparison
cat("• MR-based sigma used for control chart estimation (SPC-consistent)\n\n")

# ===============================================================================
# VISUALISATION
# ===============================================================================

par(mfrow = c(2,2))

# Time series plot
plot(usd_eur_returns, type='l', col='blue',
     main='USD/EUR Log-Returns Over Time',
     ylab='Log-Return', xlab='Time')
abline(h=0, col='red', lty=2)

# Histogram with density
hist(usd_eur_returns, breaks=50, probability=TRUE,
     main='Histogram of Log-Returns',
     xlab='Log-Return', col='lightblue', border='white')
lines(density(usd_eur_returns), col='red', lwd=2)

# Q-Q plot
qqnorm(usd_eur_returns, main='Q-Q Plot: Log-Returns vs Normal')
qqline(usd_eur_returns, col='red', lwd=2)

# ACF
acf(usd_eur_returns, main='ACF: Log-Returns')

par(mfrow=c(1,1))

# ===============================================================================
# SECTION 2: STATIONARITY TESTS
# ===============================================================================

cat("==========================================================\n")
cat("STEP 2: Stationarity Testing\n")
cat("==========================================================\n")

# ADF test
adf_result <- adf.test(usd_eur_returns)
cat(paste("ADF Test p-value:", round(adf_result$p.value, 4), 
          ifelse(adf_result$p.value < 0.05, "✓ Stationary", "✗ Non-stationary"), "\n"))

# KPSS test
kpss_result <- kpss.test(usd_eur_returns, null="Level")
cat(paste("KPSS Test p-value:", round(kpss_result$p.value, 4),
          ifelse(kpss_result$p.value > 0.05, "✓ Stationary", "✗ Non-stationary"), "\n\n"))

# ===============================================================================
# SECTION 3: AUTOCORRELATION ANALYSIS
# ===============================================================================

cat("==========================================================\n")
cat("STEP 3: Autocorrelation Analysis\n")
cat("==========================================================\n")

# Ljung-Box test
lb_result <- Box.test(usd_eur_returns, lag=20, type="Ljung-Box")
cat(paste("Ljung-Box Test (lag=20) p-value:", round(lb_result$p.value, 4), "\n"))

if(lb_result$p.value < 0.05){
  cat("⚠️  Significant autocorrelation detected. Fitting ARIMA model...\n\n")
  
  # Fit ARIMA
  arima_model <- auto.arima(usd_eur_returns, seasonal=FALSE, trace=FALSE)
  cat("Best ARIMA model:", arimaorder(arima_model), "\n")
  cat(paste("AIC:", round(arima_model$aic, 2), 
            "BIC:", round(BIC(arima_model), 2), "\n"))
  
  # Extract residuals
  residuals_clean <- residuals(arima_model)
  
  # Test residuals for autocorrelation
  lb_resid <- Box.test(residuals_clean, lag=20, type="Ljung-Box")
  cat(paste("Ljung-Box on residuals p-value:", round(lb_resid$p.value, 4), "\n"))
  
  # Use residuals for control chart monitoring
  monitoring_series <- as.numeric(residuals_clean)
  cat("\n✓ Using ARIMA residuals for control chart monitoring\n\n")
  
} else {
  cat("✓ No significant autocorrelation. Using raw returns.\n\n")
  monitoring_series <- usd_eur_returns
}

# ===============================================================================
# SECTION 4: DISTRIBUTION FITTING
# ===============================================================================

cat("==========================================================\n")
cat("STEP 4: Distribution Fitting\n")
cat("==========================================================\n")

fit_results <- list()
aic_table <- data.frame(
  Distribution = character(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

x_fit <- as.numeric(monitoring_series)
x_fit <- x_fit[is.finite(x_fit)]

# ---------- Normal ----------
tryCatch({
  fit_norm <- fitdist(x_fit, "norm")
  fit_results[["norm"]] <- fit_norm
  aic_table <- rbind(aic_table, data.frame(
    Distribution = "Normal",
    AIC = fit_norm$aic,
    BIC = fit_norm$bic
  ))
  cat(paste("Normal      - AIC:", round(fit_norm$aic, 2),
            "BIC:", round(fit_norm$bic, 2), "\n"))
}, error = function(e) {
  cat("Normal      - Fitting failed\n")
})

# ---------- Student-t ----------
# fitdistr from MASS uses location m, scale s, df
tryCatch({
  fit_t <- fitdistr(
    x_fit,
    densfun = function(x, m, s, df) dt((x - m) / s, df = df) / s,
    start = list(m = mean(x_fit), s = sd(x_fit), df = 6)
  )
  
  loglik_t <- fit_t$loglik
  k_t <- length(fit_t$estimate)
  n_t <- length(x_fit)
  aic_t <- -2 * loglik_t + 2 * k_t
  bic_t <- -2 * loglik_t + log(n_t) * k_t
  
  fit_results[["t"]] <- list(
    estimate = fit_t$estimate,
    loglik = loglik_t,
    aic = aic_t,
    bic = bic_t
  )
  
  aic_table <- rbind(aic_table, data.frame(
    Distribution = "Student-t",
    AIC = aic_t,
    BIC = bic_t
  ))
  
  cat(paste("Student-t   - AIC:", round(aic_t, 2),
            "BIC:", round(bic_t, 2), "\n"))
}, error = function(e) {
  cat("Student-t   - Fitting failed\n")
})

# ---------- Skew-Normal ----------
tryCatch({
  fit_sn <- selm(x_fit ~ 1, family = "SN")
  aic_sn <- AIC(fit_sn)
  bic_sn <- BIC(fit_sn)
  
  fit_results[["sn"]] <- fit_sn
  aic_table <- rbind(aic_table, data.frame(
    Distribution = "Skew-Normal",
    AIC = aic_sn,
    BIC = bic_sn
  ))
  
  cat(paste("Skew-Normal - AIC:", round(aic_sn, 2),
            "BIC:", round(bic_sn, 2), "\n"))
}, error = function(e) {
  cat("Skew-Normal - Fitting failed\n")
})

# ---------- Skew-Student-t ----------
tryCatch({
  fit_st <- selm(x_fit ~ 1, family = "ST")
  aic_st <- AIC(fit_st)
  bic_st <- BIC(fit_st)
  
  fit_results[["st"]] <- fit_st
  aic_table <- rbind(aic_table, data.frame(
    Distribution = "Skew-Student-t",
    AIC = aic_st,
    BIC = bic_st
  ))
  
  cat(paste("Skew-Student-t - AIC:", round(aic_st, 2),
            "BIC:", round(bic_st, 2), "\n"))
}, error = function(e) {
  cat("Skew-Student-t - Fitting failed\n")
})

# ---------- Choose best ----------
aic_table <- aic_table[order(aic_table$AIC), ]
best_dist <- aic_table$Distribution[1]

cat("\nDistribution ranking by AIC:\n")
print(aic_table)

cat(paste("\n✓ Best distribution:", best_dist, "\n\n"))

# ===============================================================================
# SECTION 5: COMPUTE CONTROL LIMITS FOR ALL THREE METHODS
# INDIVIDUALS-BASED ADAPTATION
# ===============================================================================

cat("==========================================================\n")
cat("STEP 5: Computing Control Limits (λ = 0.10)\n")
cat("==========================================================\n")

# EWMA parameters
lambda <- 0.10
L <- 2.6952

# Process parameters from empirical data
mu0 <- mean(monitoring_series)
n <- 1   # individual observations

# ----- Individuals-chart estimator of sigma -----
MR <- abs(diff(monitoring_series))
MRbar <- mean(MR, na.rm = TRUE)
d2 <- 1.128
sigma_hat <- MRbar / d2

cat(paste("Process mean (μ₀):", round(mu0, 6), "\n"))
cat(paste("MR-bar:", round(MRbar, 6), "\n"))
cat(paste("d2:", d2, "\n"))
cat(paste("Estimated sigma (MR/d2):", round(sigma_hat, 6), "\n"))
cat(paste("Sample size (n):", n, "\n\n"))

# ----- Common EWMA scale -----
# Steady-state limits so your plotting section stays unchanged
ewma_sd <- sigma_hat * sqrt(lambda / (2 - lambda))

cat(paste("EWMA steady-state SD:", round(ewma_sd, 6), "\n\n"))

# ----- Common probability estimator -----
# P_hat_X = proportion of observations less than or equal to the mean
P_L <- mean(monitoring_series <= mu0, na.rm = TRUE)
P_U <- 1 - P_L

cat(paste("P_L (below or equal to mean):", round(P_L, 4), "\n"))
cat(paste("P_U (above mean):", round(P_U, 4), "\n\n"))

# ----- WV-EWMA LIMITS -----
cat("--- WV-EWMA (Weighted Variance adaptation) ---\n")

UCL_WV <- mu0 + L * ewma_sd * sqrt(2 * P_U)
LCL_WV <- mu0 - L * ewma_sd * sqrt(2 * P_L)
CL_WV <- mu0

cat(paste("UCL_WV:", round(UCL_WV, 6), "\n"))
cat(paste("CL_WV :", round(CL_WV, 6), "\n"))
cat(paste("LCL_WV:", round(LCL_WV, 6), "\n\n"))

# ----- WSD-EWMA LIMITS -----
cat("--- WSD-EWMA (Weighted Standard Deviation adaptation) ---\n")

UCL_WSD <- mu0 + L * ewma_sd * (2 * P_U)
LCL_WSD <- mu0 - L * ewma_sd * (2 * P_L)
CL_WSD <- mu0

cat(paste("UCL_WSD:", round(UCL_WSD, 6), "\n"))
cat(paste("CL_WSD :", round(CL_WSD, 6), "\n"))
cat(paste("LCL_WSD:", round(LCL_WSD, 6), "\n\n"))

# ----- SC-EWMA LIMITS -----
cat("--- SC-EWMA (Skewness Correction adaptation) ---\n")

gamma1 <- skewness(monitoring_series)

C4_star <- (4 * gamma1 / 3) / (1 + 0.2 * gamma1^2)

UCL_SC <- mu0 + (L + C4_star) * ewma_sd
LCL_SC <- mu0 + (-L + C4_star) * ewma_sd
CL_SC <- mu0

cat(paste("Skewness (γ₁):", round(gamma1, 6), "\n"))
cat(paste("C4*:", round(C4_star, 6), "\n"))
cat(paste("UCL_SC:", round(UCL_SC, 6), "\n"))
cat(paste("CL_SC :", round(CL_SC, 6), "\n"))
cat(paste("LCL_SC:", round(LCL_SC, 6), "\n\n"))

# ===============================================================================
# SECTION 6: APPLY EWMA AND COMPUTE STATISTICS
# ===============================================================================

cat("==========================================================\n")
cat("STEP 6: Applying EWMA Charts to Data\n")
cat("==========================================================\n")

Z <- numeric(length(monitoring_series))
Z[1] <- mu0

for(i in 2:length(monitoring_series)){
  Z[i] <- lambda * monitoring_series[i] + (1 - lambda) * Z[i-1]
}

OOC_WV  <- which(Z > UCL_WV  | Z < LCL_WV)
OOC_WSD <- which(Z > UCL_WSD | Z < LCL_WSD)
OOC_SC  <- which(Z > UCL_SC  | Z < LCL_SC)

cat(paste("WV-EWMA: ", length(OOC_WV), "out-of-control points\n"))
cat(paste("WSD-EWMA:", length(OOC_WSD), "out-of-control points\n"))
cat(paste("SC-EWMA: ", length(OOC_SC), "out-of-control points\n\n"))

# ===============================================================================
# SECTION 7: COMPUTE ARL AND FAR
# ===============================================================================

cat("==========================================================\n")
cat("STEP 7: Computing ARL and False Alarm Rate (FAR)\n")
cat("==========================================================\n")

# Function to compute ARL from run lengths
compute_ARL_from_OOC <- function(OOC_points, total_length){
  if(length(OOC_points) == 0){
    return(list(ARL = total_length, FAR = 0, n_signals = 0))
  }
  
  # Run lengths = differences between consecutive signals
  if(length(OOC_points) == 1){
    run_lengths <- OOC_points[1]
  } else {
    run_lengths <- diff(c(0, OOC_points))
  }
  
  ARL <- mean(run_lengths)
  FAR <- length(OOC_points) / total_length  # Proportion of signals
  
  list(ARL = ARL, FAR = FAR, n_signals = length(OOC_points), 
       run_lengths = run_lengths)
}

# Compute for each method
ARL_WV <- compute_ARL_from_OOC(OOC_WV, length(monitoring_series))
ARL_WSD <- compute_ARL_from_OOC(OOC_WSD, length(monitoring_series))
ARL_SC <- compute_ARL_from_OOC(OOC_SC, length(monitoring_series))

# Display results
cat("\n--- WV-EWMA ---\n")
cat(paste("ARL:", round(ARL_WV$ARL, 2), "\n"))
cat(paste("FAR:", round(ARL_WV$FAR * 100, 2), "%\n"))
cat(paste("Number of signals:", ARL_WV$n_signals, "\n\n"))

cat("--- WSD-EWMA ---\n")
cat(paste("ARL:", round(ARL_WSD$ARL, 2), "\n"))
cat(paste("FAR:", round(ARL_WSD$FAR * 100, 2), "%\n"))
cat(paste("Number of signals:", ARL_WSD$n_signals, "\n\n"))

cat("--- SC-EWMA ---\n")
cat(paste("ARL:", round(ARL_SC$ARL, 2), "\n"))
cat(paste("FAR:", round(ARL_SC$FAR * 100, 2), "%\n"))
cat(paste("Number of signals:", ARL_SC$n_signals, "\n\n"))

# ===============================================================================
# SECTION 8: VISUALIZATION
# ===============================================================================

cat("==========================================================\n")
cat("STEP 8: Generating Control Chart Visualizations\n")
cat("==========================================================\n")

# Create plots
par(mfrow=c(3,1), mar=c(4, 4, 3, 2))

# WV-EWMA Chart
WVplot <- plot(Z, type='l', col='blue', lwd=1.5,
               main='WV-EWMA Control Chart (USD/EUR Returns)',
               ylab='EWMA Statistic', xlab='Time',
               ylim=range(c(Z, UCL_WV, LCL_WV)))
abline(h=c(LCL_WV, CL_WV, UCL_WV), col=c('red', 'black', 'red'), 
       lty=c(2, 1, 2), lwd=2)
points(OOC_WV, Z[OOC_WV], col='red', pch=19, cex=1.2)
legend('topleft', 
       legend=c(paste('ARL:', round(ARL_WV$ARL, 1)), 
                paste('FAR:', round(ARL_WV$FAR*100, 2), '%')),
       bty='n', cex=0.9)

# WSD-EWMA Chart
WSDplot <- plot(Z, type='l', col='darkgreen', lwd=1.5,
                main='WSD-EWMA Control Chart (USD/EUR Returns)',
                ylab='EWMA Statistic', xlab='Time',
                ylim=range(c(Z, UCL_WSD, LCL_WSD)))
abline(h=c(LCL_WSD, CL_WSD, UCL_WSD), col=c('red', 'black', 'red'), 
       lty=c(2, 1, 2), lwd=2)
points(OOC_WSD, Z[OOC_WSD], col='red', pch=19, cex=1.2)
legend('topleft', 
       legend=c(paste('ARL:', round(ARL_WSD$ARL, 1)), 
                paste('FAR:', round(ARL_WSD$FAR*100, 2), '%')),
       bty='n', cex=0.9)

# SC-EWMA Chart
SCplot <- plot(Z, type='l', col='purple', lwd=1.5,
               main='SC-EWMA Control Chart (USD/EUR Returns)',
               ylab='EWMA Statistic', xlab='Time',
               ylim=range(c(Z, UCL_SC, LCL_SC)))
abline(h=c(LCL_SC, CL_SC, UCL_SC), col=c('red', 'black', 'red'), 
       lty=c(2, 1, 2), lwd=2)
points(OOC_SC, Z[OOC_SC], col='red', pch=19, cex=1.2)
legend('topleft', 
       legend=c(paste('ARL:', round(ARL_SC$ARL, 1)), 
                paste('FAR:', round(ARL_SC$FAR*100, 2), '%')),
       bty='n', cex=0.9)

par(mfrow=c(1,1))

# ===============================================================================
# SECTION 9: EXPORT RESULTS FOR SIMULATION STUDY
# ===============================================================================

cat("==========================================================\n")
cat("STEP 9: Exporting Parameters for Simulation Study\n")
cat("==========================================================\n")

# Package empirical parameters for simulation
empirical_params_export <- list(
  # Distribution parameters
  mu0 = mu0,
  sigma_hat = sigma_hat,
  MRbar = MRbar,
  d2 = d2,
  gamma1 = gamma1_empirical,
  gamma2 = kurt_excess,  # Excess kurtosis
  
  # Best-fit distribution
  dist_type = best_dist,
  dist_params = if(!is.null(fit_results[[best_dist]])) 
    as.list(fit_results[[best_dist]]$estimate) else NULL,
  
  # Sample size
  n = length(monitoring_series),
  
  # EWMA parameters used
  lambda = lambda,
  L = L,
  
  # WV-EWMA empirical limits
  WV_limits = list(
    UCL = UCL_WV,
    CL = CL_WV,
    LCL = LCL_WV,
    P_L = P_L,
    P_U = P_U,
    ewma_sd = ewma_sd
  ),
  
  # WSD-EWMA empirical limits
  WSD_limits = list(
    UCL = UCL_WSD,
    CL = CL_WSD,
    LCL = LCL_WSD,
    P_L = P_L,
    P_U = P_U,
    ewma_sd = ewma_sd
  ),
  
  # SC-EWMA empirical limits
  SC_limits = list(
    UCL = UCL_SC,
    CL = CL_SC,
    LCL = LCL_SC,
    gamma1 = gamma1,
    C4_star = C4_star,
    ewma_sd = ewma_sd
  ),
  
  # Empirical ARL and FAR
  empirical_performance = data.frame(
    Method = c("WV-EWMA", "WSD-EWMA", "SC-EWMA"),
    ARL = c(ARL_WV$ARL, ARL_WSD$ARL, ARL_SC$ARL),
    FAR = c(ARL_WV$FAR, ARL_WSD$FAR, ARL_SC$FAR),
    N_Signals = c(ARL_WV$n_signals, ARL_WSD$n_signals, ARL_SC$n_signals)
  )
)

# Save to RData file
save(empirical_params_export, file="empirical_parameters.RData")

# Also save as CSV for easy viewing
write.csv(empirical_params_export$empirical_performance, 
          "empirical_performance.csv", row.names=FALSE)

cat("\n✓ Empirical parameters saved to 'empirical_parameters.RData'\n")
cat("✓ Performance metrics saved to 'empirical_performance.csv'\n\n")

# Print summary for dissertation
cat("==========================================================\n")
cat("SUMMARY FOR DISSERTATION (Chapter 3)\n")
cat("==========================================================\n")
cat("\nEmpirical Parameters (USD/EUR Log-Returns):\n")
cat(paste("  μ₀ =", round(mu0, 6), "\n"))
cat(paste("  sigma_hat =", round(sigma_hat, 6), "\n"))
cat(paste("  γ₁ =", round(gamma1_empirical, 3), "\n"))
cat(paste("  Fitted distribution:", best_dist, "\n"))
cat(paste("  Sample size:", length(monitoring_series), "\n"))

cat("\nControl Limits (λ = 0.10, L = 2.7):\n")
cat("\n  WV-EWMA:\n")
cat(paste("    UCL =", round(UCL_WV, 6), "\n"))
cat(paste("    LCL =", round(LCL_WV, 6), "\n"))

cat("\n  WSD-EWMA:\n")
cat(paste("    UCL =", round(UCL_WSD, 6), "\n"))
cat(paste("    LCL =", round(LCL_WSD, 6), "\n"))

cat("\n  SC-EWMA:\n")
cat(paste("    UCL =", round(UCL_SC, 6), "\n"))
cat(paste("    LCL =", round(LCL_SC, 6), "\n"))

cat("\nEmpirical Performance:\n")
print(empirical_params_export$empirical_performance)

cat("\n✅ Empirical study complete!\n")
cat("==========================================================\n\n")

