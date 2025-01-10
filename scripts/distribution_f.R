# Define the estimate_CDF function for a single z
estimate_CDF1 <- function(p, data, estimand="whole") {
  
  if (estimand=="whole") {
    F.1 <- sum((data$treat * (data$re78 <= p)) / data$ps) / 
      sum(data$treat / data$ps)
  }
  else if (estimand=="treated") {
    p1 <- mean(data$treat)
    F.1 <- mean(data$treat * (data$re78 <= p)) / p1
  }

  return(F.1)
}

estimate_CDF0 <- function(p, data, estimand="whole") {
  
  if (estimand=="whole") {
    F.0 <- sum(((1 - data$treat) * (data$re78 <= p)) / (1 - data$ps)) / 
      sum((1 - data$treat) / (1 - data$ps))
  }
  else if (estimand=="treated") {
    p0 <- mean(data$ps * (1 - data$treat) / (1 - data$ps))
    F.0 <- mean(((1 - data$treat) * (data$re78 <= p) * data$ps ) / (1 - data$ps)) / p0
  }
  
  
  return(F.0)
}

