
# Calculate Propensity Score using Logistic Model
calculate_pscore_logistic <- function(model, data, method, is_plot = FALSE) {
  
  # Extract fitted values
  if (method == "logistic") {
    fitted_values <- model$fitted.values
    
    # Determine minimum and maximum propensity score for treated group
    min_treat_ps <- min(fitted_values[data$treat == 1])
    max_treat_ps <- max(fitted_values[data$treat == 1])
    
    filter_mask <- (fitted_values > min_treat_ps & fitted_values < max_treat_ps) | data$treat == 1
  }
  else if (method == "SI") {
    fitted_values <- fitted(model)
    filter_mask <- (data$treat == 0) | data$treat == 1
  }
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


# Estimate Propensity Score
estimate_propensity_score <- function(data, formula = NULL, method="logistic", optim.method="CG", is_plot=FALSE) {
  
  if (method == "logistic") {
    # Fit logistic model
    ps_model <- glm(formula, family = binomial, data = data)
  }
  
  if (method == "SI") {
    bw <- npindexbw(formula=formula, method="ichimura", data=data, optim.method=optim.method)
    
    bw$bw = bw$bw * (nrow(data)^(-1/10))
    
    ps_model <- npindex(bw, data=data)
    
    summary(ps_model)
  }
  
  
  # Calculate propensity scores
  return(calculate_pscore_logistic(ps_model, data, method=method, is_plot))
  
}
