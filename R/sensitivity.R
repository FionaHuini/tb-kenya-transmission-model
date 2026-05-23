# TB Model Sensitivity Analysis
# Partial Rank Correlation Coefficients (PRCC)
# Author: Fiona Huini

library(tidyverse)
library(deSolve)

source("R/model.R")
source("R/calibration.R")

# ── 1. Define Parameter Ranges ─────────────────────────────────────────────
# Each parameter varied ±30% around calibrated value

n_samples <- 500

param_ranges <- tibble(
  beta    = runif(n_samples, cal_params["beta"] * 0.7, cal_params["beta"] * 1.3),
  k       = runif(n_samples, params["k"] * 0.7,        params["k"] * 1.3),
  mu      = runif(n_samples, params["mu"] * 0.7,       params["mu"] * 1.3),
  tau     = runif(n_samples, params["tau"] * 0.7,      params["tau"] * 1.3),
  delta   = runif(n_samples, params["delta"] * 0.7,    params["delta"] * 1.3),
  gamma_r = runif(n_samples, params["gamma_r"] * 0.7,  params["gamma_r"] * 1.3),
  rho     = runif(n_samples, params["rho"] * 0.7,      params["rho"] * 1.3)
)

# ── 2. Run Model for Each Sample ───────────────────────────────────────────

run_sample <- function(i) {
  
  p <- cal_params
  p["beta"]    <- param_ranges$beta[i]
  p["k"]       <- param_ranges$k[i]
  p["mu"]      <- param_ranges$mu[i]
  p["tau"]     <- param_ranges$tau[i]
  p["delta"]   <- param_ranges$delta[i]
  p["gamma_r"] <- param_ranges$gamma_r[i]
  p["rho"]     <- param_ranges$rho[i]
  
  burn <- run_tb_model(p, c(S = N0*0.95, E = N0*0.04, I = N0*0.005,
                            T = N0*0.004, R = N0*0.001), seq(0, 50, by=1))
  
  eq <- c(S = tail(burn$S,1), E = tail(burn$E,1), I = tail(burn$I,1),
          T = tail(burn$T,1), R = tail(burn$R,1))
  
  result <- run_tb_model(p, eq, seq(0, 1, by=1))
  N      <- sum(result[1, c("S","E","I","T","R")])
  inc    <- (p["beta"] * result$S[1] * result$I[1] / N) / N * 100000
  
  return(inc)
}

cat("Running sensitivity analysis — this may take a minute...\n")
incidence_samples <- map_dbl(1:n_samples, run_sample)

# ── 3. PRCC Calculation ────────────────────────────────────────────────────

prcc_df <- param_ranges |>
  mutate(incidence = incidence_samples) |>
  pivot_longer(-incidence, names_to = "parameter", values_to = "value") |>
  group_by(parameter) |>
  summarise(
    prcc = cor(rank(value), rank(incidence), method = "pearson"),
    .groups = "drop"
  ) |>
  arrange(desc(abs(prcc)))

print(prcc_df)

# ── 4. Tornado Plot ────────────────────────────────────────────────────────

prcc_df |>
  mutate(
    parameter = factor(parameter, levels = parameter[order(abs(prcc))]),
    direction = ifelse(prcc > 0, "Positive", "Negative")
  ) |>
  ggplot(aes(x = prcc, y = parameter, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Positive" = "#E63946", 
                               "Negative" = "#457B9D")) +
  labs(
    x       = "Partial Rank Correlation Coefficient (PRCC)",
    y       = "Parameter",
    title   = "Sensitivity Analysis — Drivers of TB Incidence",
    caption = "PRCC > 0: parameter increases incidence. PRCC < 0: parameter reduces incidence.",
    fill    = "Effect"
  ) +
  theme_minimal()