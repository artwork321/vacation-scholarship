source("scripts/distribution_f.R")

# Define the check function
# Source code - https://github.com/bcallaway11/BMisc/blob/master/R/BMisc.R
# https://rdrr.io/cran/qte/src/R/ciqte.R
# I implement it again for studying and trimming purpose
check_function <- function(a, tau) {
  return(a * (tau - (1 * (a <= 0))))
}


# Define the quantile objective function
q_objective_func <- function(q, y, tau, weights) {
  return(sum(weights * check_function(y - q, tau)))
}


# Estimate the `tau` quantile 
estimate_quantile <- function(tau, y, weights) {
  optim(
    par = median(y),
    fn = q_objective_func,
    y = y,
    tau = tau,
    weights = weights,
    method = "Brent",
    lower = min(y),
    upper = max(y)
  )$par
}


# Estimate QTE
estimate_QTE <- function(data, tau, method, obs=TRUE) {
  
  # QTE for experimental study
  if (obs == FALSE) {
    quantile.1 <- estimate_quantile(tau, data[data$treat==1,]$re78, rep(1, nrow(data[data$treat==1,])))
    quantile.0 <- estimate_quantile(tau, data[data$treat==0,]$re78, rep(1, nrow(data[data$treat==0,])))
    
    return(list(qte=quantile.1 - quantile.0, q1=quantile.1, q0=quantile.0))
  }
  
  if (method == "firpo") {
    w.1 <- data$treat / (data$ps * nrow(data))
    w.0 <- (1 - data$treat) / ((1 - data$ps) * nrow(data))
    
    quantile.1 <- estimate_quantile(tau, data$re78, w.1)
    quantile.0 <- estimate_quantile(tau, data$re78, w.0)
  }
  else if (method == "dist") {
    quantile.1 <- cdf2quantile(tau, estimate_CDF1, data=data)
    quantile.0 <- cdf2quantile(tau, estimate_CDF0, data=data)
  }
  else {
    print("Method is invalid")
  }
  
  return(list(qte=quantile.1 - quantile.0, q1=quantile.1, q0=quantile.0))
}


# Estimate QTE
estimate_QTT <- function(data, tau) {
  
  n.t = nrow(data[data$treat == 1,])
  
  w.1 <- data$treat / n.t
  w.0 <-  ((1 - data$treat) * data$ps) / ((1 - data$ps) * n.t)
  
  quantile.1 <- estimate_quantile(tau, data$re78, w.1)
  quantile.0 <- estimate_quantile(tau, data$re78, w.0)
  
  return(list(qte=quantile.1 - quantile.0, q1=quantile.1, q0=quantile.0))
}

# Estimate multiple QTEs
estimate_QTEs <- function(data, taus, estimand="QTE", method = "firpo", obs = TRUE) {

  if (estimand == "QTE") {
    result <- sapply(taus, function(tau) {
      estimate_QTE(data, tau = tau, method, obs)$qte
    })
  }
  else if (estimand == "QTT") {
    result <- sapply(taus, function(tau) {
      estimate_QTT(data, tau = tau)$qte
    })
  }

  return(result)
}

