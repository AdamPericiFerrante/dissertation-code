# ==============================================================================
# OUT-OF-CONTROL SIMULATION STUDY
# EWMA, WV-EWMA, WSD-EWMA, SC-EWMA
# ARL1 under mean shifts
# ==============================================================================


# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# 0. LIBRARIES
# ------------------------------------------------------------------------------

library(sn)
library(moments)
library(ggplot2)
library(dplyr)
library(tidyr)
library(param2moment)
library(parallel)

set.seed(12345)

# ------------------------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------------------------

lambda_L_pairs <- data.frame(
  Lambda = c(0.10, 0.20, 0.30),
  L      = c(2.6952, 2.8537, 2.9286)
)

skewness_levels <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

fixed_skewness_for_kurtosis <- 0.50
excess_kurtosis_levels <- c(2, 4, 6, 8)

# Out-of-control mean shifts in SD units
shift_sizes <- c(0.25, 0.50, 0.75, 1.00, 1.50)

# Simulation controls
N_sim <- 10000
m_phase1 <- 100
max_steps <- 2000
reference_n <- 20000

methods <- c("EWMA", "WV", "WSD", "SC")
distribution_names_skew <- c("Skew-Normal", "Skew-Student-t")

mu_r <- 0.0000
sigma_r <- 0.0060

n_cores <- max(1, detectCores() - 1)

# ------------------------------------------------------------------------------
# 2. PLOT SCALE SETTINGS
# ------------------------------------------------------------------------------

# Uniform x-axis intervals
x_breaks_skew <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
x_breaks_kurt <- c(2, 4, 6, 8)

# Uniform y-axis scales for full plots
y_limits_skew_full <- c(0, 300)
y_breaks_skew_full <- seq(0, 300, by = 50)

y_limits_kurt_full <- c(0, 300)
y_breaks_kurt_full <- seq(0, 300, by = 50)

# Uniform y-axis scales for zoomed plots excluding SC
y_limits_skew_zoom <- c(0, 50)
y_breaks_skew_zoom <- seq(0, 50, by = 10)

y_limits_kurt_zoom <- c(0, 50)
y_breaks_kurt_zoom <- seq(0, 50, by = 10)

# ------------------------------------------------------------------------------
# 3. HELPERS
# ------------------------------------------------------------------------------

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(NA_real_)
  moments::skewness(x)
}

safe_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || sd(x) == 0) return(NA_real_)
  moments::kurtosis(x) - 3
}

estimate_sigma_individuals <- function(x) {
  MR <- abs(diff(x))
  MRbar <- mean(MR, na.rm = TRUE)
  d2 <- 1.128
  sigma_hat <- MRbar / d2
  
  list(
    MRbar = MRbar,
    d2 = d2,
    sigma_hat = sigma_hat
  )
}

ewma_sd_fun <- function(sigma_hat, lambda) {
  sigma_hat * sqrt(lambda / (2 - lambda))
}

# ------------------------------------------------------------------------------
# 4. PARAMETER BUILDERS
# ------------------------------------------------------------------------------

get_sn_params_from_skewness <- function(target_skewness) {
  pars <- param2moment::moment2sn(
    mean = 0,
    sd = 1,
    skewness = target_skewness
  )
  
  list(
    xi = unname(pars[1]),
    omega = unname(pars[2]),
    alpha = unname(pars[3])
  )
}

get_st_params_from_moments <- function(target_skewness, target_excess_kurtosis) {
  pars <- param2moment::moment2st(
    mean = 0,
    sd = 1,
    skewness = target_skewness,
    kurtosis = 3 + target_excess_kurtosis
  )
  
  list(
    xi = unname(pars[1]),
    omega = unname(pars[2]),
    alpha = unname(pars[3]),
    nu = unname(pars[4])
  )
}

# ------------------------------------------------------------------------------
# 5. RETURN GENERATION
# ------------------------------------------------------------------------------

generate_returns <- function(n, dist_label, params, mu_r = 0, sigma_r = 0.006, shift = 0) {
  if (dist_label == "Skew-Normal") {
    z <- rsn(n, xi = params$xi, omega = params$omega, alpha = params$alpha)
  } else if (dist_label == "Skew-Student-t") {
    z <- rst(n, xi = params$xi, omega = params$omega, alpha = params$alpha, nu = params$nu)
  } else {
    stop("Unsupported distribution.")
  }
  
  mu_r + shift + sigma_r * as.numeric(z)
}

# ------------------------------------------------------------------------------
# 6. PRECOMPUTED REFERENCE STATS
# ------------------------------------------------------------------------------

compute_reference_stats <- function(dist_label,
                                    params,
                                    reference_n = 20000,
                                    mu_r = 0,
                                    sigma_r = 0.006) {
  ref <- generate_returns(
    n = reference_n,
    dist_label = dist_label,
    params = params,
    mu_r = mu_r,
    sigma_r = sigma_r,
    shift = 0
  )
  
  ref_mean <- mean(ref)
  P_L <- mean(ref <= ref_mean)
  P_U <- 1 - P_L
  gamma1 <- safe_skewness(ref)
  gamma2 <- safe_excess_kurtosis(ref)
  C4_star <- (4 * gamma1 / 3) / (1 + 0.2 * gamma1^2)
  
  list(
    ref_mean = ref_mean,
    P_L = P_L,
    P_U = P_U,
    gamma1 = gamma1,
    gamma2 = gamma2,
    C4_star = C4_star
  )
}

# ------------------------------------------------------------------------------
# 7. PHASE I ESTIMATION
# ------------------------------------------------------------------------------

estimate_phase1_hybrid <- function(x_phase1, ref_stats) {
  x_phase1 <- as.numeric(x_phase1)
  x_phase1 <- x_phase1[is.finite(x_phase1)]
  
  mu0 <- mean(x_phase1, na.rm = TRUE)
  sigma_info <- estimate_sigma_individuals(x_phase1)
  
  list(
    mu0 = mu0,
    sigma_hat = sigma_info$sigma_hat,
    P_L = ref_stats$P_L,
    P_U = ref_stats$P_U,
    gamma1 = ref_stats$gamma1,
    gamma2 = ref_stats$gamma2,
    C4_star = ref_stats$C4_star
  )
}

# ------------------------------------------------------------------------------
# 8. CONTROL LIMITS
# ------------------------------------------------------------------------------

get_limits <- function(est, lambda, L, method) {
  ewma_sd <- ewma_sd_fun(est$sigma_hat, lambda)
  
  if (method == "EWMA") {
    return(list(
      UCL = est$mu0 + L * ewma_sd,
      CL  = est$mu0,
      LCL = est$mu0 - L * ewma_sd
    ))
  }
  
  if (method == "WV") {
    return(list(
      UCL = est$mu0 + L * ewma_sd * sqrt(2 * est$P_U),
      CL  = est$mu0,
      LCL = est$mu0 - L * ewma_sd * sqrt(2 * est$P_L)
    ))
  }
  
  if (method == "WSD") {
    return(list(
      UCL = est$mu0 + L * ewma_sd * (2 * est$P_U),
      CL  = est$mu0,
      LCL = est$mu0 - L * ewma_sd * (2 * est$P_L)
    ))
  }
  
  if (method == "SC") {
    return(list(
      UCL = est$mu0 + (L + est$C4_star) * ewma_sd,
      CL  = est$mu0,
      LCL = est$mu0 + (-L + est$C4_star) * ewma_sd
    ))
  }
  
  stop("Unknown method.")
}

# ------------------------------------------------------------------------------
# 9. OUT-OF-CONTROL RUN LENGTH
# ------------------------------------------------------------------------------

simulate_run_length_oc <- function(dist_label,
                                   params,
                                   lambda,
                                   limits,
                                   shift,
                                   max_steps = 2000,
                                   mu_r = 0,
                                   sigma_r = 0.006) {
  z_ewma <- limits$CL
  
  for (t in 1:max_steps) {
    x_new <- generate_returns(
      n = 1,
      dist_label = dist_label,
      params = params,
      mu_r = mu_r,
      sigma_r = sigma_r,
      shift = shift
    )
    
    z_ewma <- lambda * x_new + (1 - lambda) * z_ewma
    
    if (z_ewma > limits$UCL || z_ewma < limits$LCL) {
      return(t)
    }
  }
  
  max_steps
}

# ------------------------------------------------------------------------------
# 10. ONE REPLICATION
# ------------------------------------------------------------------------------

one_replication_oc <- function(dist_label,
                               params,
                               lambda,
                               L,
                               method,
                               ref_stats,
                               shift,
                               m_phase1 = 100,
                               max_steps = 2000,
                               mu_r = 0,
                               sigma_r = 0.006) {
  x_phase1 <- generate_returns(
    n = m_phase1,
    dist_label = dist_label,
    params = params,
    mu_r = mu_r,
    sigma_r = sigma_r,
    shift = 0
  )
  
  est <- estimate_phase1_hybrid(x_phase1, ref_stats)
  limits <- get_limits(est, lambda, L, method)
  
  rl <- simulate_run_length_oc(
    dist_label = dist_label,
    params = params,
    lambda = lambda,
    limits = limits,
    shift = shift,
    max_steps = max_steps,
    mu_r = mu_r,
    sigma_r = sigma_r
  )
  
  c(rl = rl)
}

# ------------------------------------------------------------------------------
# 11. RUN ONE CONFIGURATION
# ------------------------------------------------------------------------------

run_configuration_oc <- function(dist_label,
                                 params,
                                 target_sk,
                                 target_kurt,
                                 scenario_name,
                                 lambda,
                                 L,
                                 method,
                                 shift_size,
                                 N_sim = 1000,
                                 m_phase1 = 100,
                                 max_steps = 2000,
                                 reference_n = 20000,
                                 mu_r = 0,
                                 sigma_r = 0.006,
                                 n_cores = 1) {
  
  ref_stats <- compute_reference_stats(
    dist_label = dist_label,
    params = params,
    reference_n = reference_n,
    mu_r = mu_r,
    sigma_r = sigma_r
  )
  
  shift_value <- shift_size * sigma_r
  
  sim_fun <- function(i) {
    one_replication_oc(
      dist_label = dist_label,
      params = params,
      lambda = lambda,
      L = L,
      method = method,
      ref_stats = ref_stats,
      shift = shift_value,
      m_phase1 = m_phase1,
      max_steps = max_steps,
      mu_r = mu_r,
      sigma_r = sigma_r
    )
  }
  
  if (.Platform$OS.type == "unix") {
    sim_out <- mclapply(1:N_sim, sim_fun, mc.cores = n_cores)
  } else {
    sim_out <- lapply(1:N_sim, sim_fun)
  }
  
  rl_vec <- as.numeric(unlist(sim_out))
  
  data.frame(
    Scenario = scenario_name,
    Distribution = dist_label,
    Target_Skewness = target_sk,
    Target_ExcessKurtosis = target_kurt,
    Achieved_Skewness = ref_stats$gamma1,
    Achieved_ExcessKurtosis = ref_stats$gamma2,
    Lambda = lambda,
    L = L,
    Method = method,
    Shift_Size = shift_size,
    ARL1 = mean(rl_vec, na.rm = TRUE),
    SDRL1 = sd(rl_vec, na.rm = TRUE),
    MRL1 = median(rl_vec, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# 12. SKEWNESS STUDY
# ------------------------------------------------------------------------------

run_skewness_study_oc <- function() {
  results <- list()
  counter <- 1
  
  for (dist_label in distribution_names_skew) {
    cat("=====================================================\n")
    cat("Out-of-control skewness study | Distribution:", dist_label, "\n")
    cat("=====================================================\n")
    
    for (target_sk in skewness_levels) {
      
      if (dist_label == "Skew-Normal") {
        params <- get_sn_params_from_skewness(target_sk)
        target_kurt <- NA_real_
      } else {
        target_kurt <- 4
        params <- get_st_params_from_moments(
          target_skewness = target_sk,
          target_excess_kurtosis = target_kurt
        )
      }
      
      ref_stats <- compute_reference_stats(
        dist_label = dist_label,
        params = params,
        reference_n = reference_n,
        mu_r = mu_r,
        sigma_r = sigma_r
      )
      
      cat("Target skewness =", target_sk,
          "| achieved skewness =", round(ref_stats$gamma1, 4), "\n")
      
      for (shift_size in shift_sizes) {
        for (i in 1:nrow(lambda_L_pairs)) {
          lambda <- lambda_L_pairs$Lambda[i]
          L <- lambda_L_pairs$L[i]
          
          for (method in methods) {
            cat("  shift =", shift_size,
                "| lambda =", lambda,
                "| method =", method, "\n")
            
            results[[counter]] <- run_configuration_oc(
              dist_label = dist_label,
              params = params,
              target_sk = target_sk,
              target_kurt = target_kurt,
              scenario_name = "Skewness",
              lambda = lambda,
              L = L,
              method = method,
              shift_size = shift_size,
              N_sim = N_sim,
              m_phase1 = m_phase1,
              max_steps = max_steps,
              reference_n = reference_n,
              mu_r = mu_r,
              sigma_r = sigma_r,
              n_cores = n_cores
            )
            
            counter <- counter + 1
          }
        }
      }
    }
  }
  
  bind_rows(results)
}

# ------------------------------------------------------------------------------
# 13. KURTOSIS STUDY
# ------------------------------------------------------------------------------

run_kurtosis_study_oc <- function() {
  results <- list()
  counter <- 1
  
  dist_label <- "Skew-Student-t"
  
  cat("=====================================================\n")
  cat("Out-of-control kurtosis study | Distribution: Skew-Student-t\n")
  cat("=====================================================\n")
  
  for (target_kurt in excess_kurtosis_levels) {
    params <- get_st_params_from_moments(
      target_skewness = fixed_skewness_for_kurtosis,
      target_excess_kurtosis = target_kurt
    )
    
    ref_stats <- compute_reference_stats(
      dist_label = dist_label,
      params = params,
      reference_n = reference_n,
      mu_r = mu_r,
      sigma_r = sigma_r
    )
    
    cat("Target excess kurtosis =", target_kurt,
        "| achieved excess kurtosis =", round(ref_stats$gamma2, 4), "\n")
    
    for (shift_size in shift_sizes) {
      for (i in 1:nrow(lambda_L_pairs)) {
        lambda <- lambda_L_pairs$Lambda[i]
        L <- lambda_L_pairs$L[i]
        
        for (method in methods) {
          cat("  shift =", shift_size,
              "| lambda =", lambda,
              "| method =", method, "\n")
          
          results[[counter]] <- run_configuration_oc(
            dist_label = dist_label,
            params = params,
            target_sk = fixed_skewness_for_kurtosis,
            target_kurt = target_kurt,
            scenario_name = "Kurtosis",
            lambda = lambda,
            L = L,
            method = method,
            shift_size = shift_size,
            N_sim = N_sim,
            m_phase1 = m_phase1,
            max_steps = max_steps,
            reference_n = reference_n,
            mu_r = mu_r,
            sigma_r = sigma_r,
            n_cores = n_cores
          )
          
          counter <- counter + 1
        }
      }
    }
  }
  
  bind_rows(results)
}

# ------------------------------------------------------------------------------
# 14. PLOT SETTINGS
# ------------------------------------------------------------------------------

theme_dissertation <- function() {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

# Fixed legend order everywhere
method_levels <- c("EWMA", "WV", "WSD", "SC")

method_colors <- c(
  "EWMA" = "#1f77b4",
  "WV"   = "#d62728",
  "WSD"  = "#2ca02c",
  "SC"   = "#9467bd"
)

# ------------------------------------------------------------------------------
# 15. FULL PLOTS WITH UNIFORM SCALES
# ------------------------------------------------------------------------------

plot_arl1_vs_skewness_full <- function(results_df, shift_size,
                                       save_dir = "plots_arl1_skewness_full") {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  df_plot <- results_df %>%
    filter(Scenario == "Skewness", Shift_Size == shift_size) %>%
    mutate(Method = factor(Method, levels = method_levels))
  
  for (dist_name in unique(df_plot$Distribution)) {
    for (lambda_val in unique(df_plot$Lambda)) {
      df_sub <- df_plot %>%
        filter(Distribution == dist_name, Lambda == lambda_val)
      
      p <- ggplot(df_sub,
                  aes(x = Target_Skewness, y = ARL1,
                      color = Method, group = Method)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        scale_color_manual(values = method_colors, drop = FALSE) +
        scale_x_continuous(
          breaks = x_breaks_skew,
          limits = range(x_breaks_skew)
        ) +
        scale_y_continuous(
          breaks = y_breaks_skew_full,
          limits = y_limits_skew_full
        ) +
        labs(
          title = "ARL1 vs Skewness",
          subtitle = paste("Distribution:", dist_name,
                           "| Lambda:", lambda_val,
                           "| Shift:", shift_size, "sigma"),
          x = "Skewness level",
          y = "ARL1",
          color = "Method"
        ) +
        theme_dissertation()
      
      file_name <- paste0(
        save_dir, "/arl1_",
        gsub("[^A-Za-z0-9]", "_", dist_name),
        "_lambda_", gsub("\\.", "_", as.character(lambda_val)),
        "_shift_", gsub("\\.", "_", as.character(shift_size)),
        ".png"
      )
      
      ggsave(file_name, plot = p, width = 8, height = 5, dpi = 300)
    }
  }
}

plot_arl1_vs_kurtosis_full <- function(results_df, shift_size,
                                       save_dir = "plots_arl1_kurtosis_full") {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  df_plot <- results_df %>%
    filter(Scenario == "Kurtosis", Shift_Size == shift_size) %>%
    mutate(Method = factor(Method, levels = method_levels))
  
  for (lambda_val in unique(df_plot$Lambda)) {
    df_sub <- df_plot %>% filter(Lambda == lambda_val)
    
    p <- ggplot(df_sub,
                aes(x = Target_ExcessKurtosis, y = ARL1,
                    color = Method, group = Method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = method_colors, drop = FALSE) +
      scale_x_continuous(
        breaks = x_breaks_kurt,
        limits = range(x_breaks_kurt)
      ) +
      scale_y_continuous(
        breaks = y_breaks_kurt_full,
        limits = y_limits_kurt_full
      ) +
      labs(
        title = "ARL1 vs Excess Kurtosis",
        subtitle = paste("Distribution: Skew-Student-t",
                         "| Lambda:", lambda_val,
                         "| Shift:", shift_size, "sigma"),
        x = "Excess kurtosis level",
        y = "ARL1",
        color = "Method"
      ) +
      theme_dissertation()
    
    file_name <- paste0(
      save_dir, "/arl1_kurtosis_lambda_",
      gsub("\\.", "_", as.character(lambda_val)),
      "_shift_", gsub("\\.", "_", as.character(shift_size)),
      ".png"
    )
    
    ggsave(file_name, plot = p, width = 8, height = 5, dpi = 300)
  }
}

# ------------------------------------------------------------------------------
# 16. ZOOMED PLOTS EXCLUDING SC
# ------------------------------------------------------------------------------

plot_arl1_vs_skewness_zoom <- function(results_df, shift_size,
                                       save_dir = "plots_arl1_skewness_zoom") {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  df_plot <- results_df %>%
    filter(Scenario == "Skewness",
           Shift_Size == shift_size,
           Method != "SC") %>%
    mutate(Method = factor(Method, levels = c("EWMA", "WV", "WSD")))
  
  zoom_colors <- method_colors[c("EWMA", "WV", "WSD")]
  
  for (dist_name in unique(df_plot$Distribution)) {
    for (lambda_val in unique(df_plot$Lambda)) {
      df_sub <- df_plot %>%
        filter(Distribution == dist_name, Lambda == lambda_val)
      
      p <- ggplot(df_sub,
                  aes(x = Target_Skewness, y = ARL1,
                      color = Method, group = Method)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        scale_color_manual(values = zoom_colors, drop = FALSE) +
        scale_x_continuous(
          breaks = x_breaks_skew,
          limits = range(x_breaks_skew)
        ) +
        scale_y_continuous(
          breaks = y_breaks_skew_zoom,
          limits = y_limits_skew_zoom
        ) +
        labs(
          title = "ARL1 vs Skewness (Zoomed)",
          subtitle = paste("Distribution:", dist_name,
                           "| Lambda:", lambda_val,
                           "| Shift:", shift_size, "sigma"),
          x = "Skewness level",
          y = "ARL1",
          color = "Method"
        ) +
        theme_dissertation()
      
      file_name <- paste0(
        save_dir, "/arl1_zoom_",
        gsub("[^A-Za-z0-9]", "_", dist_name),
        "_lambda_", gsub("\\.", "_", as.character(lambda_val)),
        "_shift_", gsub("\\.", "_", as.character(shift_size)),
        ".png"
      )
      
      ggsave(file_name, plot = p, width = 8, height = 5, dpi = 300)
    }
  }
}

plot_arl1_vs_kurtosis_zoom <- function(results_df, shift_size,
                                       save_dir = "plots_arl1_kurtosis_zoom") {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  df_plot <- results_df %>%
    filter(Scenario == "Kurtosis",
           Shift_Size == shift_size,
           Method != "SC") %>%
    mutate(Method = factor(Method, levels = c("EWMA", "WV", "WSD")))
  
  zoom_colors <- method_colors[c("EWMA", "WV", "WSD")]
  
  for (lambda_val in unique(df_plot$Lambda)) {
    df_sub <- df_plot %>% filter(Lambda == lambda_val)
    
    p <- ggplot(df_sub,
                aes(x = Target_ExcessKurtosis, y = ARL1,
                    color = Method, group = Method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = zoom_colors, drop = FALSE) +
      scale_x_continuous(
        breaks = x_breaks_kurt,
        limits = range(x_breaks_kurt)
      ) +
      scale_y_continuous(
        breaks = y_breaks_kurt_zoom,
        limits = y_limits_kurt_zoom
      ) +
      labs(
        title = "ARL1 vs Excess Kurtosis (Zoomed)",
        subtitle = paste("Distribution: Skew-Student-t",
                         "| Lambda:", lambda_val,
                         "| Shift:", shift_size, "sigma"),
        x = "Excess kurtosis level",
        y = "ARL1",
        color = "Method"
      ) +
      theme_dissertation()
    
    file_name <- paste0(
      save_dir, "/arl1_kurtosis_zoom_lambda_",
      gsub("\\.", "_", as.character(lambda_val)),
      "_shift_", gsub("\\.", "_", as.character(shift_size)),
      ".png"
    )
    
    ggsave(file_name, plot = p, width = 8, height = 5, dpi = 300)
  }
}

# ------------------------------------------------------------------------------
# 17. RUN STUDY
# ------------------------------------------------------------------------------

results_skewness_oc <- run_skewness_study_oc()
results_kurtosis_oc <- run_kurtosis_study_oc()
results_all_oc <- bind_rows(results_skewness_oc, results_kurtosis_oc)

write.csv(results_skewness_oc, "simulation_results_skewness_oc.csv", row.names = FALSE)
write.csv(results_kurtosis_oc, "simulation_results_kurtosis_oc.csv", row.names = FALSE)
write.csv(results_all_oc, "simulation_results_all_oc.csv", row.names = FALSE)

# Full plots with uniform scales
plot_arl1_vs_skewness_full(results_all_oc, shift_size = 0.50)
plot_arl1_vs_skewness_full(results_all_oc, shift_size = 1.00)

plot_arl1_vs_kurtosis_full(results_all_oc, shift_size = 0.50)
plot_arl1_vs_kurtosis_full(results_all_oc, shift_size = 1.00)

# Zoomed plots excluding SC
plot_arl1_vs_skewness_zoom(results_all_oc, shift_size = 0.50)
plot_arl1_vs_skewness_zoom(results_all_oc, shift_size = 1.00)

plot_arl1_vs_kurtosis_zoom(results_all_oc, shift_size = 0.50)
plot_arl1_vs_kurtosis_zoom(results_all_oc, shift_size = 1.00)

cat("=====================================================\n")
cat("OUT-OF-CONTROL SIMULATION COMPLETE\n")
cat("=====================================================\n")
cat("Saved files:\n")
cat(" - simulation_results_skewness_oc.csv\n")
cat(" - simulation_results_kurtosis_oc.csv\n")
cat(" - simulation_results_all_oc.csv\n")
cat(" - plots_arl1_skewness_full/*.png\n")
cat(" - plots_arl1_kurtosis_full/*.png\n")
cat(" - plots_arl1_skewness_zoom/*.png\n")
cat(" - plots_arl1_kurtosis_zoom/*.png\n")
cat("=====================================================\n")