extract_hilma <- function(k, hilma_out){

  effects <- matrix(NA, 3, 3)
  colnames(effects) <- c("effect","estimate","pv")

  effects[, 1] <- c("indirect","direct","total")
  effects[1, 2] <- hilma_out$beta_hat[k, 1]
  effects[1, 3] <- hilma_out$pvalue_beta_hat[k, 1]
  effects[2, 2] <- hilma_out$alpha1_hat[k]
  effects[3, 2] <- hilma_out$alpha1_hat[k] + hilma_out$beta_hat[k, 1]

  effects <- as.data.frame(effects)
  effects$lambda <- hilma_out$lambda_used

  return(as.data.frame(effects))

}
