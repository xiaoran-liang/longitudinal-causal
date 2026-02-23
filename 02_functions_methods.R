
# Function to produce the summary diagnostics for IPTW weights
diagnose_weights <- function(w,
                             normalize = TRUE,
                             probs = c(0.5, 0.9, 0.99, 0.999),
                             top_frac = c(0.01, 0.05)) {
  
  if (normalize) {
    w <- w / mean(w)
  }
  
  n <- length(w)
  
  w_q  <- quantile(w, probs = probs, names = FALSE)
  w_max <- max(w)
  ess <- (sum(w)^2) / sum(w^2)
  
  ord <- order(w, decreasing = TRUE)
  
  k1 <- ceiling(top_frac[1] * n)
  k5 <- ceiling(top_frac[2] * n)
  
  share_top1 <- sum(w[ord[1:k1]]) / sum(w)
  share_top5 <- sum(w[ord[1:k5]]) / sum(w)
  
  out <- c(
    q50 = w_q[1],
    q90 = w_q[2],
    q99 = w_q[3],
    q999 = w_q[4],
    max = w_max,
    ESS = ess,
    share_top1 = share_top1,
    share_top5 = share_top5
  )
  
  return(out)
}


# IPTW estimation for three time points, based on the CBPS method
run_iptw_cbps_t3 <- function(dat,
                             ps_formulas,
                             truth_total = NULL,
                             cbps_method = "exact",
                             vcov_type = "HC1",
                             ci_level = 0.95) {
  
  # ----------------------------
  # 1. CBPS propensity models
  # ----------------------------
  
  to_formula <- function(x) {
    if (inherits(x, "formula")) return(x)
    if (is.character(x) && length(x) == 1) return(stats::as.formula(x))
    stop("Each ps_formulas element must be a formula or single string.")
  }
  
  f1 <- to_formula(ps_formulas[[1]])
  f2 <- to_formula(ps_formulas[[2]])
  f3 <- to_formula(ps_formulas[[3]])
  
  cbps_1 <- CBPS::CBPS(f1, data = dat, method = cbps_method)
  cbps_2 <- CBPS::CBPS(f2, data = dat, method = cbps_method)
  cbps_3 <- CBPS::CBPS(f3, data = dat, method = cbps_method)
  
  w1 <- cbps_1$weights
  w2 <- cbps_2$weights
  w3 <- cbps_3$weights
  
  W1 <- w1
  W2 <- w1 * w2
  W3 <- w1 * w2 * w3
  
  # ----------------------------
  # 2. Weighted MSM regressions
  # ----------------------------
  
  reg1 <- stats::lm(Y1 ~ A1, data = dat, weights = W1)
  reg2 <- stats::lm(Y2 ~ A1 + A2, data = dat, weights = W2)
  reg3 <- stats::lm(Y3 ~ A1 + A2 + A3, data = dat, weights = W3)
  
  s1 <- summary(reg1)$coefficients
  s2 <- summary(reg2)$coefficients
  s3 <- summary(reg3)$coefficients
  
  est <- c(
    "beta_1(1)" = s1["A1", "Estimate"],
    "beta_2(1)" = s2["A1", "Estimate"],
    "beta_2(2)" = s2["A2", "Estimate"],
    "beta_3(1)" = s3["A1", "Estimate"],
    "beta_3(2)" = s3["A2", "Estimate"],
    "beta_3(3)" = s3["A3", "Estimate"]
  )
  
  se_model <- c(
    "beta_1(1)" = s1["A1", "Std. Error"],
    "beta_2(1)" = s2["A1", "Std. Error"],
    "beta_2(2)" = s2["A2", "Std. Error"],
    "beta_3(1)" = s3["A1", "Std. Error"],
    "beta_3(2)" = s3["A2", "Std. Error"],
    "beta_3(3)" = s3["A3", "Std. Error"]
  )
  
  # ----------------------------
  # 3. Sandwich SE (HC1 default)
  # ----------------------------
  
  get_sandwich_se <- function(reg) {
    vc <- sandwich::vcovHC(reg, type = vcov_type)
    lmtest::coeftest(reg, vcov. = vc)[, "Std. Error"]
  }
  
  se1_s <- get_sandwich_se(reg1)["A1"]
  se2_s <- get_sandwich_se(reg2)[c("A1","A2")]
  se3_s <- get_sandwich_se(reg3)[c("A1","A2","A3")]
  
  se_sandwich <- c(
    "beta_1(1)" = se1_s,
    "beta_2(1)" = se2_s[1],
    "beta_2(2)" = se2_s[2],
    "beta_3(1)" = se3_s[1],
    "beta_3(2)" = se3_s[2],
    "beta_3(3)" = se3_s[3]
  )
  
  # ----------------------------
  # 4. Confidence intervals
  # ----------------------------
  
  z <- stats::qnorm(1 - (1 - ci_level)/2)
  
  ci_model <- cbind(
    lower = est - z * se_model,
    upper = est + z * se_model
  )
  
  ci_sandwich <- cbind(
    lower = est - z * se_sandwich,
    upper = est + z * se_sandwich
  )
  
  cover_model <- NULL
  cover_sandwich <- NULL
  
  if (!is.null(truth_total)) {
    cover_model <- (truth_total > ci_model[,1]) &
      (truth_total < ci_model[,2])
    
    cover_sandwich <- (truth_total > ci_sandwich[,1]) &
      (truth_total < ci_sandwich[,2])
  }
  
  # ----------------------------
  # 5. Weight diagnostics
  # ----------------------------
  
  diag_W1 <- diagnose_weights(W1)
  diag_W2 <- diagnose_weights(W2)
  diag_W3 <- diagnose_weights(W3)
  
  diagnostics <- rbind(
    W1 = diag_W1,
    W2 = diag_W2,
    W3 = diag_W3
  )
  
  # ----------------------------
  # 6. Return
  # ----------------------------
  
  return(list(
    est = est,
    se_model = se_model,
    se_sandwich = se_sandwich,
    ci_model = ci_model,
    ci_sandwich = ci_sandwich,
    cover_model = cover_model,
    cover_sandwich = cover_sandwich,
    weight_diagnostics = diagnostics
  ))
}

# -----------------------------
# Basic G-estimation (NUC), T=3
# Residualize A only; just-identified GMM at each time point.
# -----------------------------

#' Residualize a continuous exposure A on history H via Gaussian GLM
#' @param dat data.frame
#' @param a_name exposure column name (e.g. "A1")
#' @param hist either a formula RHS string like "F1 + F2 + V11" or a formula "A1 ~ ..."
#' @return list(residual = A - E[A|H], fit = glm fit)
residualize_A_gaussian <- function(dat, a_name, hist) {
  f <- if (inherits(hist, "formula")) {
    hist
  } else if (is.character(hist) && length(hist) == 1) {
    stats::as.formula(paste(a_name, "~", hist))
  } else {
    stop("hist must be a formula or a single RHS string.")
  }
  
  fit <- stats::glm(f, family = stats::gaussian(), data = dat)
  p_hat <- stats::predict(fit, type = "response")
  R <- dat[[a_name]] - p_hat
  list(residual = R, fit = fit)
}

#' Compute HC0/HC1 just-identified sandwich SE for basic G-est at time t
#'
#' Moments at time 3: m = (E[R1*U1], E[R2*U2], E[R3*U3])
#' where U3=Y-b3*A3; U2=Y-b2*A2-b3*A3; U1=Y-b1*A1-b2*A2-b3*A3
#'
#' @param beta_hat numeric vector length p (p=1,2,3)
#' @param Y numeric
#' @param A_mat matrix n x p of (A1,...,Ap) with columns in order
#' @param R_mat matrix n x p of (R1,...,Rp) with columns in order
#' @param type "HC0" or "HC1"
#' @return list(se, V, J, Omega)
sandwich_basic_gest <- function(beta_hat, Y, A_mat, R_mat, type = c("HC0","HC1")) {
  type <- match.arg(type)
  n <- NROW(A_mat)
  p <- NCOL(A_mat)
  
  # Build U_k for k=1..p:
  # U_k = Y - sum_{j=k..p} beta_j * A_j
  # and moment k uses R_k * U_k
  U_list <- vector("list", p)
  for (k in 1:p) {
    Aj <- A_mat[, k:p, drop = FALSE]
    bj <- beta_hat[k:p]
    U_list[[k]] <- Y - as.numeric(Aj %*% bj)
  }
  
  g_mat <- matrix(NA_real_, nrow = n, ncol = p)
  for (k in 1:p) {
    g_mat[, k] <- R_mat[, k] * U_list[[k]]
  }
  
  # Omega_hat (centered)
  g_ctr <- scale(g_mat, center = TRUE, scale = FALSE)
  Omega_hat <- crossprod(g_ctr) / n
  if (type == "HC1") {
    Omega_hat <- (n/(n - p)) * Omega_hat
  }
  
  # Jacobian J = d E[g(beta)] / d beta^T
  # For this basic setup:
  # g_k(beta) = E[ R_k * (Y - sum_{j=k..p} beta_j A_j ) ]
  # so derivative wrt beta_j is: -E[ R_k * A_j ] for j>=k else 0
  D <- matrix(0, nrow = p, ncol = p)
  for (k in 1:p) {
    for (j in k:p) {
      D[k, j] <- mean(R_mat[, k] * A_mat[, j])
    }
  }
  J <- -D
  
  Jinv <- solve(J)
  V <- Jinv %*% Omega_hat %*% t(Jinv) / n
  se <- sqrt(diag(V))
  
  list(se = se, V = V, J = J, Omega = Omega_hat)
}

#' Basic G-estimation for T=3
#'
#' @param dat data.frame with A1,A2,A3,Y1,Y2,Y3 and covariates in histories
#' @param H list with elements t1,t2,t3 giving RHS strings for exposure models
#'          (e.g. H$t1="F1+F2+V11", H$t2="A1+Y1+...", H$t3="A2+Y2+...")
#' @param type "HC0" or "HC1" for sandwich SE
#' @param ci_level e.g. 0.95
#' @param truth_total optional named vector of length 6 for coverage
#' @return list(est, se, ci, cover, fits)
run_gest_basic_t3 <- function(dat, H, type = c("HC0","HC1"), ci_level = 0.95, truth_total = NULL) {
  type <- match.arg(type)
  
  # ---- residualize exposures ----
  r1 <- residualize_A_gaussian(dat, "A1", H$t1)
  r2 <- residualize_A_gaussian(dat, "A2", H$t2)
  r3 <- residualize_A_gaussian(dat, "A3", H$t3)
  
  R1 <- r1$residual
  R2 <- r2$residual
  R3 <- r3$residual
  
  # ---- helper to estimate beta at time t by minimizing Q ----
  fit_time <- function(t, Y, A_mat, R_mat) {
    p <- NCOL(A_mat)
    Q_obj <- function(beta, Y, A_mat, R_mat, W = diag(p)) {
      # build U_k and moments
      U_list <- vector("list", p)
      for (k in 1:p) {
        Aj <- A_mat[, k:p, drop = FALSE]
        bj <- beta[k:p]
        U_list[[k]] <- Y - as.numeric(Aj %*% bj)
      }
      m <- numeric(p)
      for (k in 1:p) {
        m[k] <- mean(R_mat[, k] * U_list[[k]])
      }
      drop(t(m) %*% W %*% m)
    }
    
    beta0 <- rep(0, p)
    opt <- stats::optim(
      par = beta0,
      fn  = Q_obj,
      method = "BFGS",
      Y = Y, A_mat = A_mat, R_mat = R_mat
    )
    
    beta_hat <- opt$par
    sw <- sandwich_basic_gest(beta_hat, Y = Y, A_mat = A_mat, R_mat = R_mat, type = type)
    
    list(beta = beta_hat, se = sw$se, opt = opt, sandwich = sw)
  }
  
  # ---- time 1 ----
  out1 <- fit_time(
    t = 1,
    Y = dat$Y1,
    A_mat = cbind(A1 = dat$A1),
    R_mat = cbind(R1 = R1)
  )
  
  # ---- time 2 ----
  out2 <- fit_time(
    t = 2,
    Y = dat$Y2,
    A_mat = cbind(A1 = dat$A1, A2 = dat$A2),
    R_mat = cbind(R1 = R1, R2 = R2)
  )
  
  # ---- time 3 ----
  out3 <- fit_time(
    t = 3,
    Y = dat$Y3,
    A_mat = cbind(A1 = dat$A1, A2 = dat$A2, A3 = dat$A3),
    R_mat = cbind(R1 = R1, R2 = R2, R3 = R3)
  )
  
  # ---- stack into 6 estimands ----
  est <- c(
    "beta_1(1)" = out1$beta[1],
    "beta_2(1)" = out2$beta[1],
    "beta_2(2)" = out2$beta[2],
    "beta_3(1)" = out3$beta[1],
    "beta_3(2)" = out3$beta[2],
    "beta_3(3)" = out3$beta[3]
  )
  
  se <- c(
    "beta_1(1)" = out1$se[1],
    "beta_2(1)" = out2$se[1],
    "beta_2(2)" = out2$se[2],
    "beta_3(1)" = out3$se[1],
    "beta_3(2)" = out3$se[2],
    "beta_3(3)" = out3$se[3]
  )
  
  # ---- CI + coverage ----
  z <- stats::qnorm(1 - (1 - ci_level)/2)
  ci <- cbind(lower = est - z * se, upper = est + z * se)
  
  cover <- NULL
  if (!is.null(truth_total)) {
    cover <- (truth_total > ci[, "lower"]) & (truth_total < ci[, "upper"])
  }
  
  list(
    est = est,
    se = se,
    ci = ci,
    cover = cover,
    residual_fits = list(A1 = r1$fit, A2 = r2$fit, A3 = r3$fit),
    optim = list(t1 = out1$opt, t2 = out2$opt, t3 = out3$opt),
    sandwich = list(t1 = out1$sandwich, t2 = out2$sandwich, t3 = out3$sandwich)
  )
}


# -----------------------------
# Efficient G-estimation (NUC), T=3
# "Double residualization": residualize both Y_t and A_j on H_t.
# Two-step GMM (W=I then Wopt = Omega^{-1}); HC0/HC1 sandwich.
# -----------------------------

# Helper: make formula from response + RHS string
.form_from_rhs <- function(resp, rhs) stats::as.formula(paste(resp, "~", rhs))

# Helper: residualize a variable on RHS using Gaussian model
residualize_gaussian <- function(dat, resp, rhs) {
  fit <- stats::glm(.form_from_rhs(resp, rhs), family = stats::gaussian(), data = dat)
  mu <- stats::predict(fit, type = "response")
  res <- dat[[resp]] - mu
  list(residual = res, fit = fit)
}

# Helper: safe inverse (use ginv if singular)
safe_invert <- function(M) {
  tryCatch(solve(M), error = function(e) MASS::ginv(M))
}

# Helper: build optimal weight from individual moments g_i (n x p)
optimal_weight <- function(g_mat) {
  n <- nrow(g_mat)
  g_ctr <- scale(g_mat, center = TRUE, scale = FALSE)
  Omega <- crossprod(g_ctr) / n
  list(Omega = Omega, Wopt = safe_invert(Omega))
}

# Helper: CI + optional coverage
ci_and_cover <- function(est, se, ci_level = 0.95, truth = NULL) {
  z <- stats::qnorm(1 - (1 - ci_level)/2)
  ci <- cbind(lower = est - z * se, upper = est + z * se)
  cover <- NULL
  if (!is.null(truth)) {
    cover <- (truth > ci[, "lower"]) & (truth < ci[, "upper"])
  }
  list(ci = ci, cover = cover)
}

# Sandwich for efficient (residualized) moments at time t with p params
# g_k = E[ rA_k,t * U_k ], with U_k defined on residual scale
sandwich_eff_gest <- function(beta_hat, g_builder, J_builder, type = c("HC0","HC1")) {
  type <- match.arg(type)
  
  g_mat <- g_builder(beta_hat)  # n x p
  n <- nrow(g_mat)
  p <- ncol(g_mat)
  
  g_ctr <- scale(g_mat, center = TRUE, scale = FALSE)
  Omega <- crossprod(g_ctr) / n
  if (type == "HC1") Omega <- (n/(n - p)) * Omega
  
  J <- J_builder()              # p x p
  Jinv <- solve(J)
  V <- Jinv %*% Omega %*% t(Jinv) / n
  se <- sqrt(diag(V))
  
  list(se = se, V = V, J = J, Omega = Omega, moments = colMeans(g_mat))
}

#' Efficient G-estimation wrapper for T=3
#'
#' @param dat data.frame with A1,A2,A3,Y1,Y2,Y3 and covariates
#' @param H list with elements t1,t2,t3 RHS strings for histories
#' @param type "HC0" or "HC1"
#' @param ci_level e.g. 0.95
#' @param truth_total optional named vector length 6 for coverage
#' @return list(est, se, ci, cover, details)
run_gest_efficient_t3 <- function(dat, H, type = c("HC0","HC1"), ci_level = 0.95, truth_total = NULL) {
  type <- match.arg(type)
  
  # -----------------------------
  # TIME 3 (Y3 with 3 parameters)
  # -----------------------------
  # Residuals for A's given H1/H2/H3
  rA1_H1 <- residualize_gaussian(dat, "A1", H$t1)$residual
  rA2_H1 <- residualize_gaussian(dat, "A2", H$t1)$residual
  rA3_H1 <- residualize_gaussian(dat, "A3", H$t1)$residual
  
  rA2_H2 <- residualize_gaussian(dat, "A2", H$t2)$residual
  rA3_H2 <- residualize_gaussian(dat, "A3", H$t2)$residual
  
  rA3_H3 <- residualize_gaussian(dat, "A3", H$t3)$residual
  
  # Residuals for Y3 given H1/H2/H3
  rY3_H1 <- residualize_gaussian(dat, "Y3", H$t1)$residual
  rY3_H2 <- residualize_gaussian(dat, "Y3", H$t2)$residual
  rY3_H3 <- residualize_gaussian(dat, "Y3", H$t3)$residual
  
  Q3 <- function(beta, W = diag(3)) {
    U3 <- rY3_H3 - beta[3] * rA3_H3
    U2 <- rY3_H2 - beta[2] * rA2_H2 - beta[3] * rA3_H2
    U1 <- rY3_H1 - beta[1] * rA1_H1 - beta[2] * rA2_H1 - beta[3] * rA3_H1
    m <- c(mean(rA1_H1 * U1), mean(rA2_H2 * U2), mean(rA3_H3 * U3))
    drop(t(m) %*% W %*% m)
  }
  
  beta0_3 <- c(0, 0, 0)
  step1_3 <- stats::optim(par = beta0_3, fn = function(b) Q3(b, W = diag(3)), method = "BFGS")
  beta_step1_3 <- step1_3$par
  
  # Build Omega and Wopt at step 1
  g_builder_3 <- function(beta) {
    U3 <- rY3_H3 - beta[3] * rA3_H3
    U2 <- rY3_H2 - beta[2] * rA2_H2 - beta[3] * rA3_H2
    U1 <- rY3_H1 - beta[1] * rA1_H1 - beta[2] * rA2_H1 - beta[3] * rA3_H1
    cbind(rA1_H1 * U1, rA2_H2 * U2, rA3_H3 * U3)
  }
  ow3 <- optimal_weight(g_builder_3(beta_step1_3))
  
  step2_3 <- stats::optim(
    par = beta_step1_3,
    fn  = function(b) Q3(b, W = ow3$Wopt),
    method = "BFGS"
  )
  beta3 <- step2_3$par
  
  # Jacobian for time 3
  J_builder_3 <- function() {
    D <- matrix(0, 3, 3)
    D[1, ] <- c(mean(rA1_H1 * rA1_H1),
                mean(rA1_H1 * rA2_H1),
                mean(rA1_H1 * rA3_H1))
    D[2, ] <- c(0,
                mean(rA2_H2 * rA2_H2),
                mean(rA2_H2 * rA3_H2))
    D[3, ] <- c(0, 0,
                mean(rA3_H3 * rA3_H3))
    -D
  }
  
  sw3 <- sandwich_eff_gest(beta3, g_builder = g_builder_3, J_builder = J_builder_3, type = type)
  se3 <- sw3$se
  
  # -----------------------------
  # TIME 2 (Y2 with 2 parameters)
  # -----------------------------
  # Residuals for A1,A2 given H1 (needed for U1), and A2 given H2 (needed for U2)
  rA1_H1_2 <- residualize_gaussian(dat, "A1", H$t1)$residual
  rA2_H1_2 <- residualize_gaussian(dat, "A2", H$t1)$residual
  rA2_H2_2 <- residualize_gaussian(dat, "A2", H$t2)$residual
  
  # Residuals for Y2 given H1 and H2
  rY2_H1 <- residualize_gaussian(dat, "Y2", H$t1)$residual
  rY2_H2 <- residualize_gaussian(dat, "Y2", H$t2)$residual
  
  Q2 <- function(beta, W = diag(2)) {
    U2 <- rY2_H2 - beta[2] * rA2_H2_2
    U1 <- rY2_H1 - beta[1] * rA1_H1_2 - beta[2] * rA2_H1_2
    m <- c(mean(rA1_H1_2 * U1), mean(rA2_H2_2 * U2))
    drop(t(m) %*% W %*% m)
  }
  
  beta0_2 <- c(0, 0)
  step1_2 <- stats::optim(par = beta0_2, fn = function(b) Q2(b, W = diag(2)), method = "BFGS")
  beta_step1_2 <- step1_2$par
  
  g_builder_2 <- function(beta) {
    U2 <- rY2_H2 - beta[2] * rA2_H2_2
    U1 <- rY2_H1 - beta[1] * rA1_H1_2 - beta[2] * rA2_H1_2
    cbind(rA1_H1_2 * U1, rA2_H2_2 * U2)
  }
  ow2 <- optimal_weight(g_builder_2(beta_step1_2))
  
  step2_2 <- stats::optim(
    par = beta_step1_2,
    fn  = function(b) Q2(b, W = ow2$Wopt),
    method = "BFGS"
  )
  beta2 <- step2_2$par
  
  J_builder_2 <- function() {
    D <- matrix(0, 2, 2)
    D[1, ] <- c(mean(rA1_H1_2 * rA1_H1_2),
                mean(rA1_H1_2 * rA2_H1_2))
    D[2, ] <- c(0,
                mean(rA2_H2_2 * rA2_H2_2))
    -D
  }
  
  sw2 <- sandwich_eff_gest(beta2, g_builder = g_builder_2, J_builder = J_builder_2, type = type)
  se2 <- sw2$se
  
  # -----------------------------
  # TIME 1 (Y1 with 1 parameter)
  # -----------------------------
  rA1_H1_1 <- residualize_gaussian(dat, "A1", H$t1)$residual
  rY1_H1   <- residualize_gaussian(dat, "Y1", H$t1)$residual
  
  Q1 <- function(beta, W = matrix(1, 1, 1)) {
    U1 <- rY1_H1 - beta[1] * rA1_H1_1
    m <- mean(rA1_H1_1 * U1)
    drop(t(m) %*% W %*% m)
  }
  
  beta0_1 <- c(0)
  step1_1 <- stats::optim(par = beta0_1, fn = function(b) Q1(b, W = matrix(1,1,1)), method = "BFGS")
  beta_step1_1 <- step1_1$par
  
  g_builder_1 <- function(beta) {
    U1 <- rY1_H1 - beta[1] * rA1_H1_1
    cbind(rA1_H1_1 * U1)
  }
  ow1 <- optimal_weight(g_builder_1(beta_step1_1))
  
  step2_1 <- stats::optim(
    par = beta_step1_1,
    fn  = function(b) Q1(b, W = matrix(ow1$Wopt, 1, 1)),
    method = "BFGS"
  )
  beta1 <- step2_1$par
  
  J_builder_1 <- function() {
    D11 <- mean(rA1_H1_1 * rA1_H1_1)
    matrix(-D11, 1, 1)
  }
  
  sw1 <- sandwich_eff_gest(beta1, g_builder = g_builder_1, J_builder = J_builder_1, type = type)
  se1 <- sw1$se
  
  # -----------------------------
  # Stack outputs
  # -----------------------------
  est <- c(
    "beta_1(1)" = beta1[1],
    "beta_2(1)" = beta2[1],
    "beta_2(2)" = beta2[2],
    "beta_3(1)" = beta3[1],
    "beta_3(2)" = beta3[2],
    "beta_3(3)" = beta3[3]
  )
  
  se <- c(
    "beta_1(1)" = se1[1],
    "beta_2(1)" = se2[1],
    "beta_2(2)" = se2[2],
    "beta_3(1)" = se3[1],
    "beta_3(2)" = se3[2],
    "beta_3(3)" = se3[3]
  )
  
  cc <- ci_and_cover(est, se, ci_level = ci_level, truth = truth_total)
  
  list(
    est = est,
    se = se,
    ci = cc$ci,
    cover = cc$cover,
    details = list(
      time1 = list(step1 = step1_1, step2 = step2_1, Omega = ow1$Omega, Wopt = ow1$Wopt, sandwich = sw1),
      time2 = list(step1 = step1_2, step2 = step2_2, Omega = ow2$Omega, Wopt = ow2$Wopt, sandwich = sw2),
      time3 = list(step1 = step1_3, step2 = step2_3, Omega = ow3$Omega, Wopt = ow3$Wopt, sandwich = sw3)
    )
  )
}


# -----------------------------
# IV estimation: two-parameter decay model (beta, alpha), T=3
# -----------------------------

# Helper: build lower-triangular indicator D and exponent matrix P for T
.decay_design_mats <- function(T) {
  D <- matrix(0, T, T)
  P <- matrix(0, T, T)
  for (t in 1:T) {
    for (j in 1:t) {
      D[t, j] <- 1
      P[t, j] <- (t - j)
    }
  }
  list(D = D, P = P)
}

# Objective 
# a = c(beta, alpha)
.iv_decay_objective_t3 <- function(a, AC, YC, G, cov_mat) {
  beta  <- a[1]
  alpha <- a[2]
  
  T <- ncol(AC)
  mats <- .decay_design_mats(T)
  D <- mats$D
  P <- mats$P
  
  # T x T lower-triangular coefficient matrix
  Gmat <- D * beta * (alpha ^ P)
  
  # residual_matrix[i,t] = YC[i,t] - sum_{j<=t} Gmat[t,j] * AC[i,j]
  n <- nrow(AC)
  residual_matrix <- matrix(0, nrow = n, ncol = T)
  for (t in 1:T) {
    residual_matrix[, t] <- YC[, t] - as.numeric(AC[, 1:t, drop = FALSE] %*% Gmat[t, 1:t])
  }
  
  mR <- G - mean(G)
  
  s <- solve(cov_mat, crossprod(residual_matrix, mR))
  sum(as.numeric(s)^2)
}

# Sandwich variance
.iv_decay_sandwich_t3 <- function(bb, AC, YC, G, type = c("HC0","HC1")) {
  type <- match.arg(type)
  beta  <- bb[1]
  alpha <- bb[2]
  
  n <- nrow(AC)
  T <- ncol(AC)
  mats <- .decay_design_mats(T)
  D <- mats$D
  P <- mats$P
  
  mR <- G - mean(G)
  
  # Coefficient matrix
  Gmat <- D * beta * (alpha ^ P)
  
  # Residuals Y0
  Y0 <- matrix(0, nrow = n, ncol = T)
  for (t in 1:T) {
    Y0[, t] <- YC[, t] - as.numeric(AC[, 1:t, drop = FALSE] %*% Gmat[t, 1:t])
  }
  
  # v_i,t = mR_i * Y0_i,t  (n x T)
  v <- mR * Y0
  
  # Omega_hat for moments across t
  v_ctr <- scale(v, center = TRUE, scale = FALSE)
  Omega <- crossprod(v_ctr) / n
  if (type == "HC1") {
    q <- T
    Omega <- (n / (n - q)) * Omega
  }
  
  # Derivatives of Gmat wrt beta and alpha
  dG_dbeta  <- D * (alpha ^ P)
  dG_dalpha <- (P * D) * beta * (alpha ^ pmax(P - 1, 0))
  
  # Bread M1 is T x 2
  M1 <- matrix(0, nrow = T, ncol = 2)
  for (t in 1:T) {
    dY_dbeta  <- -as.numeric(AC[, 1:t, drop = FALSE] %*% dG_dbeta[t, 1:t])
    dY_dalpha <- -as.numeric(AC[, 1:t, drop = FALSE] %*% dG_dalpha[t, 1:t])
    
    M1[t, 1] <- mean(mR * dY_dbeta)
    M1[t, 2] <- mean(mR * dY_dalpha)
  }
  
  # (T x 2) generalized inverse
  sM1 <- MASS::ginv(M1)
  
  V <- sM1 %*% Omega %*% t(sM1) / n
  SE <- sqrt(diag(V))
  
  list(Varcov = V, SE = SE, bread = M1, meat = Omega)
}


#' IV estimation under decay model (beta, alpha), T=3
#'
#' @param dat data.frame with columns: A1,A2,A3,Y1,Y2,Y3,G and baseline covariates
#' @param baseline_covars character vector of baseline covariate names (default c("F1","F2"))
#' @param instrument name of instrument column (default "G")
#' @param init initial values c(beta, alpha)
#' @param cov_mat 3x3 matrix used in objective; default identity
#' @param type "HC0" or "HC1" sandwich
#' @param ci_level CI level
#' @param truth optional named vector/list with elements beta and alpha for coverage
#' @return list(est, se, ci, cover, details)
run_iv_decay_t3 <- function(dat,
                            baseline_covars = c("F1","F2"),
                            instrument = "G",
                            init = c(-0.1, 1),
                            cov_mat = diag(1, 3),
                            type = c("HC0","HC1"),
                            ci_level = 0.95,
                            truth = NULL) {
  type <- match.arg(type)
  
  needed <- c("A1","A2","A3","Y1","Y2","Y3", instrument, baseline_covars)
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) stop("dat is missing columns: ", paste(missing, collapse = ", "))
  
  rhs <- paste(baseline_covars, collapse = " + ")
  
  AC <- cbind(
    AC_1 = stats::residuals(stats::lm(stats::as.formula(paste("A1 ~", rhs)), data = dat)),
    AC_2 = stats::residuals(stats::lm(stats::as.formula(paste("A2 ~", rhs)), data = dat)),
    AC_3 = stats::residuals(stats::lm(stats::as.formula(paste("A3 ~", rhs)), data = dat))
  )
  
  YC <- cbind(
    YC_1 = stats::residuals(stats::lm(stats::as.formula(paste("Y1 ~", rhs)), data = dat)),
    YC_2 = stats::residuals(stats::lm(stats::as.formula(paste("Y2 ~", rhs)), data = dat)),
    YC_3 = stats::residuals(stats::lm(stats::as.formula(paste("Y3 ~", rhs)), data = dat))
  )
  
  G <- dat[[instrument]]
  
  opt <- stats::optim(
    par = init,
    fn  = .iv_decay_objective_t3,
    method = "BFGS",
    hessian = TRUE,
    AC = AC, YC = YC, G = G, cov_mat = cov_mat
  )
  
  est <- opt$par
  names(est) <- c("beta", "alpha")
  
  sw <- .iv_decay_sandwich_t3(est, AC = AC, YC = YC, G = G, type = type)
  se <- sw$SE
  names(se) <- c("beta", "alpha")
  
  z <- stats::qnorm(1 - (1 - ci_level)/2)
  ci <- rbind(
    beta  = c(est["beta"]  - z * se["beta"],  est["beta"]  + z * se["beta"]),
    alpha = c(est["alpha"] - z * se["alpha"], est["alpha"] + z * se["alpha"])
  )
  colnames(ci) <- c("lower","upper")  # force lowercase
  
  cover <- NULL
  if (!is.null(truth)) {
    truth_vec <- unlist(truth)
    cover <- c(
      beta  = (truth_vec["beta"]  > ci["beta","lower"])  & (truth_vec["beta"]  < ci["beta","upper"]),
      alpha = (truth_vec["alpha"] > ci["alpha","lower"]) & (truth_vec["alpha"] < ci["alpha","upper"])
    )
  }
  
  list(
    est = est,
    se = se,
    ci = ci,
    cover = cover,
    details = list(optim = opt, sandwich = sw, AC = AC, YC = YC)
  )
}
