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