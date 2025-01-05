source("scripts/pscore.R")

# Perform matching treated units to control units with replacement
ps_matching <- function(formula=NULL, data, bias_adjust=FALSE, is_plot=FALSE, estimand="ATE") {
  
  out_ps <- estimate_propensity_score(data, formula, is_plot=is_plot)
  
  data <- out_ps$filter_data
  data$ps <- out_ps$pscore
  Y <- data$re78
  treat <- data$treat
  X <- matrix(data$ps, nrow(data), 1)
  
  # Find matches
  if (bias_adjust) {
    matchest <- Match(
      Y = Y, 
      Tr = treat, 
      X = X, 
      M = 1, 
      replace = TRUE,
      estimand = estimand
    )
  }
  else {
    matchest <- Match(
      Y = Y, 
      Tr = treat, 
      X = X, 
      M = 1, 
      BiasAdjust = FALSE, 
      replace = TRUE, 
      ties = FALSE, 
      estimand = estimand
    )
  }
  
  
  # Compute the estimand and SE
  if (bias_adjust) {
    se <- matchest$se
  }
  else {
    se <- matchest$se.standard
  }
  
  list(
    m.out=matchest,
    est=matchest$est,
    se=se,
    data=data,
    estimand=matchest$estimand
  )
}


# Visualize Matching
plot_match <- function(pscore, index.treated, index.control) {
  
  # Extract propensity scores for treated and control
  pscore_treat <- pscore[index.treated]
  pscore_control <- pscore[index.control]
  
  # Sort propensity scores
  treated_sorted <- sort(pscore_treat)
  control_sorted <- sort(pscore_control)
  
  # Create plot
  plot(
    treated_sorted, 
    type = "s", 
    xlab = "Treated Units (Sorted by Propensity Score)", 
    ylab = "Estimated Propensity Score", 
    main = "Propensity Score: Treated vs Matched Control", 
    xlim = c(0, length(treated_sorted)), 
    ylim = c(0, 1),
    lwd = 2
  )
  
  # Add matched control units
  lines(
    seq_along(control_sorted), 
    control_sorted, 
    type = "s", 
    lty = 2,
    lwd = 2
  )
  
  # Add legend
  legend(
    "bottomright", 
    legend = c("Treated", "Matched Control"), 
    lty = c(1, 2), 
    lwd = 2
  )
}

# Perform every step of ps matching
ps_matching_analysis <- function(dataset, formula = NULL, bias_adjust=FALSE, estimand="ATE") {
  
  # Preprocess data
  data <- preprocess_data(dataset)
  
  # Perform matching
  match_result <- ps_matching(formula, data, bias_adjust, estimand=estimand)
  
  out <- match_result$m.out
  trim_data <- match_result$data
  matched_data <- trim_data[c(out$index.treated, out$index.control), ]
  
  # Store results
  results <- list(
    est = match_result$est,
    se = match_result$se,
    trim.data = trim_data,
    m.data = matched_data,
    estimand = match_result$estimand
  )
  
  
  return(results)
}
