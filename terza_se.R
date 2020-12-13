# correct se (Terza 2016) ----
terza_se <- function(stage2, stage1, stage1_coef = "pred_stage1", robust = T) {
  # first stage
  expWalpha <- fitted.values(stage1) # xb
  alpha <- coef(stage1)
  
  # alternative terza standard errors using robust first and second stage SEs
  if(robust) {
    require(sandwich)
    covalpha <- vcovHC(stage1, type="HC1")
  } else {
    covalpha <- vcov(stage1)
  }
  
  # TSRI
  expXbeta <- fitted.values(stage2)
  beta <- coef(stage2)
  if(robust) {
    require(sandwich)
    covbeta <- vcovHC(stage2, type="HC1")
  } else {
    covbeta <- vcov(stage2)
  }
  bxu <- beta[stage1_coef]
  W <- model.matrix(stage1)
  X <- model.matrix(stage2)
  Xbeta <- X %*% beta
  
  # Compute the asymptotic covariance matrix of
  # the 2SRI estimate of beta.
  paJ <- -bxu * expXbeta * expWalpha * W
  pbJ <- expXbeta * X
  Bba <- t(pbJ) %*% paJ
  Bbb <- t(pbJ) %*% pbJ
  d22 <- solve(Bbb) %*% Bba %*% covalpha %*% t(Bba) %*% solve(Bbb) + covbeta
  # Terza standard errors
  ses <- sqrt(diag(d22))
  # t-statistics
  # tstats <- beta / ses
  # # pvalues
  # pvalues <- 2 * pnorm(-1*abs(tstats))
  # # estimates with 95% CI limits
  # cbind(beta, beta - qt(.975, stage2$df.null)*ses, beta + 1.96*ses)
  return(d22)
}
