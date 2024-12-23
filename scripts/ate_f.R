# Calculate ATE in experimental study
ate.experimental <- function(data) {
  matched.treated <- data$re78[data$treat==1]
  matched.control <- data$re78[data$treat==0]
  
  ate <- mean(matched.treated) - mean(matched.control)
  
  return(list(est=ate))
}


# Calculate Propensity Score using Logistic Model
calculate_pscore_logistic <- function(model, data, is_plot = FALSE) {

  # Extract fitted values
  fitted_values <- model$fitted.values
  
  # Determine minimum and maximum propensity score for treated group
  min_treat_ps <- min(fitted_values[data$treat == 1])
  max_treat_ps <- max(fitted_values[data$treat == 1])
  
  # Filter data based on propensity scores
  filter_mask <- (fitted_values > min_treat_ps & fitted_values < max_treat_ps) | data$treat == 1
  data_filter <- data[filter_mask, ]
  pscore_filter <- fitted_values[filter_mask]
  
  # Reporting
  cat("Propensity Score Summary:\n")
  cat("Control Group:\n")
  print(summary(pscore_filter[data_filter$treat == 0]))
  cat("\nTreated Group:\n")
  print(summary(pscore_filter[data_filter$treat == 1]))
  
  discarded_count <- nrow(data) - nrow(data_filter)
  cat(sprintf("\nDiscarded %d control units.\n", discarded_count))
  cat(sprintf("Remaining observations: %d\n", nrow(data_filter)))
  cat("\n")
  
  # Optional plotting
  if (is_plot) {
    plot_pscore(data_filter, pscore_filter)
  }
  
  # Return filtered data and propensity scores
  list(
    filter_data = data_filter, 
    pscore = pscore_filter
  )
}

# Plot Propensity Score Distribution
plot_pscore <- function(data, pscore) {
  # Create probability bins
  bins <- cut(pscore, breaks = seq(0, 1, by = 0.1), right = FALSE)
  
  # Count treated and control units in each bin
  counts_treat_0 <- table(bins[data$treat == 0])
  counts_treat_1 <- table(bins[data$treat == 1])
  counts_matrix <- rbind(counts_treat_0, counts_treat_1)
  
  # Generate bin labels
  group_labels <- paste(seq(0.1, 1, by = 0.1))
  
  # Create bar plot
  barplot(
    counts_matrix,
    beside = TRUE,
    col = c("white", "skyblue"),
    legend.text = c("Control (treat=0)", "Treated (treat=1)"),
    main = "Histogram of Estimated Propensity Score",
    xlab = "Probability Bin",
    ylab = "Frequency",
    names.arg = group_labels,
    ylim = c(0, 200)
  )
}

# Preprocess Data
preprocess_data <- function(data) {
  
  # Add binary unemployment indicators
  data$u74 <- as.integer(data$re74 == 0)
  data$u75 <- as.integer(data$re75 == 0)
  
  return(data)
}

# Perform matching treated units to control units with replacement
ps_matching <- function(formula=NULL, data, bias_adjust=FALSE, is_plot=FALSE, estimand="ATE") {
  
  out_ps <- estimate_propensity_score(data, formula, is_plot=is_plot)
  
  data <- out_ps$filter_data
  data$ps <- out_ps$pscore
  Y <- data$re78
  treat <- data$treat
  X <- matrix(data$ps, nrow(data), 1)
  
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
  
  
  # Compute ATT and SE
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

# Estimate Propensity Score
estimate_propensity_score <- function(data, formula = NULL, method="logistic", is_plot=FALSE) {

  if (method == "logistic") {
    # Fit logistic model
    ps_model <- glm(formula, family = binomial, data = data)
    
    # Calculate propensity scores
    return(calculate_pscore_logistic(ps_model, data, is_plot))
  }
  
  if (method == "power") {
    
  }
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

ps_IPW <- function(data, estimand="ATT", norm=FALSE) {
  
  
  if(estimand=="ATT"){
    # Calculate weights
    weights <- ifelse(data$treat==1, 1, data$ps/(1-data$ps))
    
    # Calculate estimate
    n.1 <- length(data$re78[data$treat==1])
    mean.Y.1 <- mean(data$re78[data$treat==1])
    weighted.Y.0 <- sum(weights[data$treat==0] * data$re78[data$treat==0])
    
    est <- mean.Y.1 - weighted.Y.0/n.1
  }
   
  else if(estimand=="ATE") {
    # Calculate weights
    weights <- ifelse(data$treat==1, 1/(data$ps), 1/(1-data$ps))
    
    n <- nrow(data)
    sum.weight.1 <- n
    sum.weight.0 <- n
    
    if (norm) {
      sum.weight.1 <- sum(weights[data$treat==1])
      sum.weight.0 <- sum(weights[data$treat==0])
    }
      
    # Calculate estimate
    weighted.Y.1 <- sum(weights[data$treat==1] * data$re78[data$treat==1])
    weighted.Y.0 <- sum(weights[data$treat==0] * data$re78[data$treat==0])
    
    est <- weighted.Y.1/sum.weight.1 - weighted.Y.0/sum.weight.0
  }
  
  
  return(list(est=est, weights=weights, estimand=estimand))
}