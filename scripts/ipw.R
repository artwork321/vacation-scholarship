source("scripts/pscore.R")

# Compute ATT or ATE using IPW method
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