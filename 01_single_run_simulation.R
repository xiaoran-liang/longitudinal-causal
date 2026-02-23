
# 01_single_run_simulation.R
# ------------------------------------------------------------
# Single simulation run: DGP + run 4 estimation methods
#   (1) IPTW via CBPS
#   (2) Basic NUC-based G-estimation
#   (3) Efficient NUC-based G-estimation
#   (4) IV decay model: estimate (beta, alpha) + implied total effects
# ------------------------------------------------------------

## ---- Packages ----
suppressPackageStartupMessages({
  library(MASS)   # used for ginv() inside IV code 
  library(CBPS)   # CBPS estimation for IPTW
})

## ---- Source method functions ----
source("functions_methods.R")

## ---- Parameter Specification ----
seed <- 284723
set.seed(seed)

n <- 5000
time_interval <- 1:3
T <- length(time_interval)

ar1_rho <- 0.98
v1_sd   <- 0.2

## ---- DGP PARAMETERS ----
params <- list(
  MAF   = 0.2,
  gamma = 0.5,
  
  eta_A = 0.2,
  eta_Y = 0.1,
  
  beta_y = 0.5,
  
  eta_f1 = 0.2,
  eta_f2 = 0.3,
  
  beta_f1 = 0.4,
  beta_f2 = -0.1,
  
  eta_v1 = 0.5,
  beta_v1 = -0.6,
  
  tau_v = 0,
  
  beta  = -1.1,
  alpha = 0.95
)

## ---- Ground Truth: total effects + implied direct effects ----
truth <- within(list(), {
  beta   <- params$beta
  alpha  <- params$alpha
  beta_y <- params$beta_y
  tau_v  <- params$tau_v
  beta_v1 <- params$beta_v1
  
  total_11 <- beta
  total_21 <- beta * alpha
  total_22 <- beta
  total_31 <- beta * alpha^2
  total_32 <- beta * alpha
  total_33 <- beta
  
  delta_11 <- total_11
  delta_21 <- total_21 - beta_y * delta_11 - tau_v * beta_v1
  delta_22 <- total_22
  delta_31 <- total_31 - delta_21 * beta_y - delta_11 * beta_y^2 - tau_v * beta_v1 * beta_y
  delta_32 <- total_32 - beta_y * delta_22 - tau_v * beta_v1
  delta_33 <- total_33
})

truth_total <- c(
  "beta_1(1)" = truth$total_11,
  "beta_2(1)" = truth$total_21,
  "beta_2(2)" = truth$total_22,
  "beta_3(1)" = truth$total_31,
  "beta_3(2)" = truth$total_32,
  "beta_3(3)" = truth$total_33
)

## ---- DGP: one dataset (T=3) ----
simulate_once <- function(n, T, params, truth, ar1_rho, v1_sd) {
  stopifnot(T == 3)
  
  G  <- rbinom(n = n, size = 2, prob = params$MAF)
  F1 <- rnorm(n, mean = 0.2, sd = 1)
  F2 <- rnorm(n, mean = -0.1, sd = 1)
  
  V1_mat <- matrix(NA_real_, nrow = n, ncol = T)
  for (i in seq_len(n)) {
    V1_mat[i, ] <- as.numeric(
      arima.sim(
        model = list(order = c(1, 0, 0), ar = ar1_rho),
        n = T,
        rand.gen = rnorm,
        sd = v1_sd
      )
    )
  }
  V11 <- V1_mat[, 1]
  V12 <- V1_mat[, 2]
  V13 <- V1_mat[, 3]
  
  A1 <- params$gamma * G +
    params$eta_f1 * F1 + params$eta_f2 * F2 +
    params$eta_v1 * V11 + rnorm(n)
  
  Y1 <- truth$delta_11 * A1 +
    params$beta_f1 * F1 + params$beta_f2 * F2 +
    params$beta_v1 * V11 + rnorm(n)
  
  V12 <- V12 + params$tau_v * A1
  
  A2 <- params$gamma * G +
    params$eta_A * A1 + params$eta_Y * Y1 +
    params$eta_f1 * F1 + params$eta_f2 * F2 +
    params$eta_v1 * V12 + rnorm(n)
  
  Y2 <- truth$delta_21 * A1 + truth$delta_22 * A2 +
    params$beta_y * Y1 +
    params$beta_f1 * F1 + params$beta_f2 * F2 +
    params$beta_v1 * V12 + rnorm(n)
  
  V13 <- V13 + params$tau_v * A2
  
  A3 <- params$gamma * G +
    params$eta_A * A2 + params$eta_Y * Y2 +
    params$eta_f1 * F1 + params$eta_f2 * F2 +
    params$eta_v1 * V13 + rnorm(n)
  
  Y3 <- truth$delta_31 * A1 + truth$delta_32 * A2 + truth$delta_33 * A3 +
    params$beta_y * Y2 +
    params$beta_f1 * F1 + params$beta_f2 * F2 +
    params$beta_v1 * V13 + rnorm(n)
  
  data.frame(
    id = seq_len(n),
    G = G,
    F1 = F1, F2 = F2,
    V11 = V11, V12 = V12, V13 = V13,
    A1 = A1, A2 = A2, A3 = A3,
    Y1 = Y1, Y2 = Y2, Y3 = Y3
  )
}

dat <- simulate_once(n, T, params, truth, ar1_rho, v1_sd)

## ---- Histories ----
H <- list(
  t1 = "F1 + F2 + V11",
  t2 = "A1 + Y1 + F1 + F2 + V11 + V12",
  t3 = "A2 + Y2 + A1 + Y1 + F1 + F2 + V11 + V12 + V13"
)

## ---- Print truth ----
cat(sprintf("n = %d, T = %d, seed = %d\n", n, T, seed))
cat("\nTruth (total effects):\n")
print(round(truth_total, 4))

## ---- (1) IPTW via CBPS ----
ps_forms <- list(
  paste0("A1 ~ ", H$t1),
  paste0("A2 ~ ", H$t2),
  paste0("A3 ~ ", H$t3)
)

iptw_out <- run_iptw_cbps_t3(
  dat = dat,
  ps_formulas = ps_forms,
  truth_total = truth_total,
  cbps_method = "exact",
  vcov_type = "HC1",
  ci_level = 0.95
)

cat("\nIPTW (CBPS exact): estimates + SEs\n")
print(round(cbind(
  Estimate = iptw_out$est,
  SE_model = iptw_out$se_model,
  SE_HC1   = iptw_out$se_sandwich
), 4))

cat("\nIPTW weight diagnostics (normalized cumulative weights):\n")
print(round(iptw_out$weight_diagnostics, 4))

## ---- (2) Basic G-estimation ----
g_basic_out <- run_gest_basic_t3(dat = dat, H = H, type = "HC1", truth_total = truth_total)

cat("\nBasic G-estimation (HC1):\n")
print(round(cbind(Estimate = g_basic_out$est, SE_HC1 = g_basic_out$se), 4))


## ---- (3) Efficient G-estimation ----
g_eff_out <- run_gest_efficient_t3(dat = dat, H = H, type = "HC1", truth_total = truth_total)

cat("\nEfficient G-estimation (HC1):\n")
print(round(cbind(Estimate = g_eff_out$est, SE_HC1 = g_eff_out$se), 4))


## ---- (4) IV decay model: (beta, alpha) + implied totals ----
iv_out <- run_iv_decay_t3(
  dat = dat,
  baseline_covars = c("F1","F2"),
  instrument = "G",
  init = c(-0.1, 1),
  type = "HC1",
  truth = c(beta = params$beta, alpha = params$alpha)
)

cat("\nIV decay parameters (HC1):\n")
print(round(cbind(Estimate = iv_out$est, SE_HC1 = iv_out$se), 4))

beta_hat  <- unname(iv_out$est["beta"])
alpha_hat <- unname(iv_out$est["alpha"])

iv_implied <- c(
  "beta_1(1)" = beta_hat,
  "beta_2(1)" = beta_hat * alpha_hat,
  "beta_2(2)" = beta_hat,
  "beta_3(1)" = beta_hat * alpha_hat^2,
  "beta_3(2)" = beta_hat * alpha_hat,
  "beta_3(3)" = beta_hat
)

cat("\nIV implied total effects (point estimates):\n")
print(round(iv_implied, 4))


## ---- Comparison table of the 6 total effects across four methods ----
rows <- list(truth = truth_total, IPTW_CBPS = iptw_out$est[names(truth_total)])

rows$G_basic <- g_basic_out$est[names(truth_total)]
rows$G_eff   <- g_eff_out$est[names(truth_total)]
rows$IV_implied <- iv_implied[names(truth_total)]

res_tbl <- do.call(rbind, rows)

cat("\nComparison (6 total effects):\n")
print(round(res_tbl, 4))
