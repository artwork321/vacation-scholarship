create_psm_result <- function(psm_list, 
                              row_labels=c("NSW", "PSID-1", "PSID-2", "PSID-3", "CPS-1", "CPS-2", "CPS-3"), 
                              col_labels, include.se=FALSE) {
  # Number of items in the list
  n <- length(psm_list)
  
  # Initialize result matrix
  if (include.se) {
    psm_result <- matrix(NA, n + 1, 2)
  }
  else {
    psm_result <- matrix(NA, n + 1, 1)
  }
  
  
  # Iterate through the list to populate the matrix
  for (j in 1:n) {
    out <- psm_list[[j]]
    psm_result[j + 1, 1] <- round(out$est, 0)
    
    if (include.se) {
      # Check for `se` or fallback to `se.standard`
      if (is.null(out$se) & is.null(out$se.stadard)) {
        psm_result[j + 1, 2] <- NA
      } else if (is.null(out$se)) {
        psm_result[j + 1, 2] <- round(out$se.standard, 0)
      }
      else {
        psm_result[j + 1, 2] <- round(out$se, 0)
      }
    }

  }
  
  # Assign row and column names
  rownames(psm_result) <- row_labels
  colnames(psm_result) <- col_labels
  
  return(psm_result)
}


create_qte_result <- function(qte_list, row_labels, col_labels) {
  # Number of items in the list
  n <- 7
  i = 0
  
  if (length(qte_list) < 7) {
    i <- 7 - length(qte_list)
  }
  
  # Initialize result matrix
  qte_result <- matrix(NA, n, 3)
  
  
  # Iterate through the list to populate the matrix
  for (j in 1:length(qte_list)) {
    out <- qte_list[[j]]
    qte_result[j + i, 1] <- round(out[1], 0)
    qte_result[j + i, 2] <- round(out[2], 0)
    qte_result[j + i, 3] <- round(out[3], 0)
  }
  
  # Assign row and column names
  rownames(qte_result) <- row_labels
  colnames(qte_result) <- col_labels
  
  return(qte_result)
}
