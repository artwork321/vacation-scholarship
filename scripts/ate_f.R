# Calculate ATE in experimental study
ate.experimental <- function(data) {
  matched.treated <- data$re78[data$treat==1]
  matched.control <- data$re78[data$treat==0]
  
  ate <- mean(matched.treated) - mean(matched.control)
  
  return(list(est=ate))
}


# Preprocess Data
preprocess_data <- function(data) {
  
  # Add binary unemployment indicators
  data$u74 <- as.integer(data$re74 == 0)
  data$u75 <- as.integer(data$re75 == 0)
  
  return(data)
}
