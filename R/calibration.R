# TB Model Calibration
# Fitting beta to WHO Kenya 2000 baseline incidence
# Author: Fiona Huini
# Parameters: Kirimi et al. (2024)

library(deSolve)
library(tidyverse)

source("R/model.R")

# ── 1. Load WHO Kenya Data ─────────────────────────────────────────────────

kenya_tb <- read_csv("data/kenya_tb_who.csv")

observed <- kenya_tb |>
  select(year, e_inc_100k) |>
  filter(year >= 2000, year <= 2022) |>
  arrange(year)

observed_2000 <- observed |>
  filter(year == 2000) |>
  pull(e_inc_100k)

# ── 2. Calibration Function ────────────────────────────────────────────────

model_incidence_2000 <- function(beta_val) {
  
  cal_params       <- params
  cal_params["beta"] <- beta_val
  
  burn_state <- c(
    S = N0 * 0.95,
    E = N0 * 0.04,
    I = N0 * 0.005,
    T = N0 * 0.004,
    R = N0 * 0.001
  )
  
  burn_result <- run_tb_model(cal_params, burn_state, seq(0, 50, by = 1))
  
  eq_state <- c(
    S = tail(burn_result$S, 1),
    E = tail(burn_result$E, 1),
    I = tail(burn_result$I, 1),
    T = tail(burn_result$T, 1),
    R = tail(burn_result$R, 1)
  )
  
  result <- run_tb_model(cal_params, eq_state, seq(0, 1, by = 1))
  N      <- sum(result[1, c("S","E","I","T","R")])
  inc    <- (beta_val * result$S[1] * result$I[1] / N) / N * 100000
  return(inc)
}

# ── 3. Objective Function ──────────────────────────────────────────────────

objective <- function(beta_val) {
  (model_incidence_2000(beta_val) - observed_2000)^2
}

# ── 4. Optimise ────────────────────────────────────────────────────────────

fit        <- optimize(objective, interval = c(0.1, 50.0))
beta_fitted <- fit$minimum

cat("Fitted beta:", round(beta_fitted, 4), "\n")
cat("SSE:        ", round(fit$objective, 4), "\n")
cat("Calibrated R0:", round(R0_calibrated, 4), "\n")
cat("Modeled 2000 incidence:", round(model_incidence_2000(beta_fitted), 1), "\n")
cat("Observed 2000 incidence:", observed_2000, "\n")
# ── 5. R0 ──────────────────────────────────────────────────────────────────

cal_params         <- params
cal_params["beta"] <- beta_fitted
R0_calibrated      <- compute_R0(cal_params)

cat("Calibrated R0:", round(R0_calibrated, 4), "\n")

# ── 6. Verify Incidence Match ──────────────────────────────────────────────

cat("Modeled 2000 incidence:", round(model_incidence_2000(beta_fitted), 1), "\n")
cat("Observed 2000 incidence:", observed_2000, "\n")

# ── 7. Fit Plot ────────────────────────────────────────────────────────────

burn_state <- c(
  S = N0 * 0.95, E = N0 * 0.04, I = N0 * 0.005,
  T = N0 * 0.004, R = N0 * 0.001
)

burn_result <- run_tb_model(cal_params, burn_state, seq(0, 50, by = 1))

eq_state <- c(
  S = tail(burn_result$S, 1), E = tail(burn_result$E, 1),
  I = tail(burn_result$I, 1), T = tail(burn_result$T, 1),
  R = tail(burn_result$R, 1)
)

fit_result <- run_tb_model(cal_params, eq_state, seq(0, 22, by = 1))
N_fit      <- rowSums(fit_result[, c("S","E","I","T","R")])

modeled_incidence <- (cal_params["beta"] * fit_result$S *
                        fit_result$I / N_fit) / N_fit * 100000

fit_df <- tibble(
  year     = 2000:2022,
  modeled  = modeled_incidence,
  observed = observed$e_inc_100k
)

fit_df |>
  pivot_longer(-year, names_to = "Series", values_to = "Incidence") |>
  ggplot(aes(x = year, y = Incidence, color = Series, linetype = Series)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = function(x) filter(x, Series == "observed"), size = 2) +
  scale_color_manual(values = c("modeled" = "#E63946",
                                "observed" = "#1D3557")) +
  annotate("text", x = 2005, y = 300,
           label = "Model uses fixed parameters\n(declining trend reflects\ntreatment programme improvement)",
           size = 3, color = "grey40", hjust = 0) +
  labs(
    x       = "Year",
    y       = "TB Incidence per 100,000",
    title   = "Modelled vs Observed TB Incidence — Kenya 2000–2022",
    caption = "Observed: WHO Global TB Report. Modelled: fixed-parameter SETR model, beta calibrated to 2000 baseline."
  ) +
  theme_minimal()

# ── 8. Scenario Analysis 2022–2030 ────────────────────────────────────────

# End TB 2030 target: 80% reduction from 2015 baseline
baseline_2015 <- observed |>
  filter(year == 2015) |>
  pull(e_inc_100k)

end_tb_target <- baseline_2015 * 0.20  # 80% reduction

cat("2015 baseline incidence:", baseline_2015, "\n")
cat("End TB 2030 target:", round(end_tb_target, 1), "\n")

# Starting point — use 2022 equilibrium state
scenario_times <- seq(0, 8, by = 1)  # 2022–2030

burn_2022 <- run_tb_model(cal_params, eq_state, seq(0, 22, by = 1))

state_2022 <- c(
  S = tail(burn_2022$S, 1),
  E = tail(burn_2022$E, 1),
  I = tail(burn_2022$I, 1),
  T = tail(burn_2022$T, 1),
  R = tail(burn_2022$R, 1)
)

# Scenario 1 — Baseline (no change)
params_s1        <- cal_params
result_s1        <- run_tb_model(params_s1, state_2022, scenario_times)
N_s1             <- rowSums(result_s1[, c("S","E","I","T","R")])
inc_s1           <- (params_s1["beta"] * result_s1$S * result_s1$I / N_s1) / N_s1 * 100000

# Scenario 2 — Treatment scale-up (tau + 30%)
params_s2        <- cal_params
params_s2["tau"] <- cal_params["tau"] * 1.30
result_s2        <- run_tb_model(params_s2, state_2022, scenario_times)
N_s2             <- rowSums(result_s2[, c("S","E","I","T","R")])
inc_s2           <- (params_s2["beta"] * result_s2$S * result_s2$I / N_s2) / N_s2 * 100000

# Scenario 3 — Combined (tau +30%, beta -20%)
params_s3         <- cal_params
params_s3["tau"]  <- cal_params["tau"] * 1.30
params_s3["beta"] <- cal_params["beta"] * 0.80
result_s3         <- run_tb_model(params_s3, state_2022, scenario_times)
N_s3              <- rowSums(result_s3[, c("S","E","I","T","R")])
inc_s3            <- (params_s3["beta"] * result_s3$S * result_s3$I / N_s3) / N_s3 * 100000

# ── 9. Scenario Plot ───────────────────────────────────────────────────────

scenario_df <- tibble(
  year     = 2022:2030,
  Baseline                = inc_s1,
  `Treatment scale-up`    = inc_s2,
  `Combined intervention` = inc_s3
)

scenario_df |>
  pivot_longer(-year, names_to = "Scenario", values_to = "Incidence") |>
  ggplot(aes(x = year, y = Incidence, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = end_tb_target, linetype = "dashed",
             color = "darkgreen", linewidth = 0.8) +
  annotate("text", x = 2023, y = end_tb_target + 8,
           label = paste0("End TB 2030 target (", round(end_tb_target, 0), "/100k)"),
           color = "darkgreen", size = 3, hjust = 0) +
  scale_color_manual(values = c(
    "Baseline"                = "#E63946",
    "Treatment scale-up"      = "#457B9D",
    "Combined intervention"   = "#2D6A4F"
  )) +
  labs(
    x       = "Year",
    y       = "TB Incidence per 100,000",
    title   = "Projected TB Incidence Under Intervention Scenarios — Kenya 2022–2030",
    caption = "End TB target: 80% reduction from 2015 baseline. Parameters: Kirimi et al. (2024)."
  ) +
  theme_minimal()

cat("2015 baseline:", baseline_2015, "\n")
cat("End TB target:", round(end_tb_target, 1), "\n")
cat("2030 Baseline scenario:", round(tail(inc_s1, 1), 1), "\n")
cat("2030 Treatment scale-up:", round(tail(inc_s2, 1), 1), "\n")
cat("2030 Combined:", round(tail(inc_s3, 1), 1), "\n")