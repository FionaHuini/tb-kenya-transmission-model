# TB Transmission Model — SETR ODE System
# Kenya 2000–2022
# Author: Fiona Huini

library(deSolve)

#' TB Compartmental Model
#'
#' @param t time sequence
#' @param state named vector of initial compartment values
#' @param params named vector of model parameters
#'
#' @return list of derivatives

tb_model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    N <- S + E + I + T + R
    
    dS <- mu * N - beta * S * I / N - mu * S
    dE <- beta * S * I / N - (k + mu) * E + rho * R
    dI <- k * E - (tau + delta + mu) * I
    dT <- tau * I - (gamma_r + mu) * T
    dR <- gamma_r * T - (rho + mu) * R
    
    list(c(dS, dE, dI, dT, dR))
  })
}


#' Compute Basic Reproduction Number
#'
#' @param params named vector of model parameters
#' @return R0 scalar

compute_R0 <- function(params) {
  with(as.list(params), {
    R0 <- (beta * k) / ((k + mu) * (tau + delta + mu))
    return(R0)
  })
}


#' Run the TB model
#'
#' @param params named vector of parameters
#' @param initial_state named vector of initial conditions
#' @param times time sequence
#' @return data frame of compartment sizes over time

run_tb_model <- function(params, initial_state, times) {
  out <- ode(
    y = initial_state,
    times = times,
    func = tb_model,
    parms = params
  )
  as.data.frame(out)
}

# ── Baseline Parameters ────────────────────────────────────────────────────
# Source: Kirimi et al. (2024) doi:10.1155/2024/5883142
# delta corrected to annual rate (CFR 50% over 3yr duration)

params <- c(
  beta    = 0.15,   # transmission rate — calibrated separately
  k       = 0.002,   # progression latent to active
  mu      = 0.0147, # natural death rate
  tau     = 0.68,   # treatment initiation rate
  delta   = 0.15,   # TB mortality rate (annual)
  gamma_r = 0.75,   # treatment success rate
  rho     = 0.003   # relapse rate
)

# Kenya 2000 population
N0 <- 31000000

initial_state <- c(
  S = N0 * 0.643,
  E = N0 * 0.350,
  I = N0 * 0.005,
  T = N0 * 0.001,
  R = N0 * 0.001
)

times <- seq(0, 22, by = 1)

result <- run_tb_model(params, initial_state, times)

# R0 with literature parameters
R0_literature <- compute_R0(params)
cat("R0 with literature parameters:", round(R0_literature, 3), "\n")