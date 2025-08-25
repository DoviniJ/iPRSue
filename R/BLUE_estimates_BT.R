#' BLUE_estimates_BT function
#'
#' Estimates individual-level polygenic risk scores (PRS) with uncertainty using a frequentist approach
#' for binary traits. This implementation applies Firth's bias-reduced logistic regression on the discovery sample,
#' computes the coefficient covariance matrix, and uses the delta method to derive PRS variance and confidence intervals.
#'
#' @param discovery_pheno Character. Path to the phenotype file for the discovery dataset. Assumes no header and that the binary trait is in the third column.
#' @param discovery_geno_mat Character. Path to the genotype matrix file for the discovery dataset. Assumes no header.
#' @param target_pheno Character. Path to the phenotype file for the target dataset. Assumes no header and individual IDs in the second column.
#' @param target_geno_mat Character. Path to the genotype matrix file for the target dataset. Assumes no header.
#' @param significance_level Numeric. Significance level for confidence intervals (e.g., 0.05 for 95% CI). Default is 0.05.
#' @param max_iterations Integer. Maximum number of iterations allowed in Firth logistic regression. Default is 100.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{IID}{Individual identifier (from the target phenotype file).}
#'   \item{PRS}{Estimated polygenic risk score for each individual.}
#'   \item{Variance}{Estimated variance of the PRS.}
#'   \item{Lower_Limit}{Lower bound of the confidence interval.}
#'   \item{Upper_Limit}{Upper bound of the confidence interval.}
#' }
#'
#' @details
#' The function fits a Firth logistic regression model using the `logistf` package to reduce small-sample bias in the discovery set.
#' It extracts SNP effect estimates and their covariance matrix, and propagates this uncertainty through to the individual-level
#' PRS in the target dataset via the delta method. Confidence intervals are derived assuming normality.
#'
#' Missing or non-estimable coefficients and variances are set to zero.
#'
#' @import data.table
#' @import bigstatsr
#' @import logistf
#'
#' @importFrom stats binomial coef glm lm qnorm quantile rnorm var vcov
#' 
#' @examples
#' \donttest{
#'   bpd <- system.file("Bpd.txt", package = "iPRSue", mustWork = TRUE)
#'   bpt <- system.file("Bpt.txt", package = "iPRSue", mustWork = TRUE)
#'   gd  <- system.file("Gd.txt",  package = "iPRSue", mustWork = TRUE)
#'   gt  <- system.file("Gt.txt",  package = "iPRSue", mustWork = TRUE)
#'
#'   results <- BLUE_estimates_BT(
#'     discovery_pheno    = bpd,
#'     discovery_geno_mat = gd,
#'     target_pheno       = bpt,
#'     target_geno_mat    = gt,
#'     significance_level = 0.05,
#'     max_iterations     = 100
#'   )
#'   head(results)
#' }
#'
#' @export
BLUE_estimates_BT <- function(discovery_pheno, 
                                     discovery_geno_mat, 
                                     target_pheno,
                                     target_geno_mat, 
                                     significance_level = 0.05,
                                     max_iterations = 100){
  phe <- fread(discovery_pheno, header=F)
  y <- phe$V3
  
  X_discovery <- as_FBM(as.matrix(fread(discovery_geno_mat, header=F)))
  df <- data.frame(y = y, X_discovery[])
  m <- logistf(y ~ ., data=df, control=logistf.control(maxit=max_iterations), pl = FALSE)
  sigma <- vcov(m)[-1,-1]
  sigma[is.na(sigma)] <- 0
  
  X <- as_FBM(as.matrix(fread(target_geno_mat, header=F)))
  M <- sigma %*% t(X[])
  BLUE_prs_var <- diag(big_prodMat(X, M))
  
  
  B <- coef(m)[-1]
  B[is.na(B)] <- 0
  BLUE_prs <- big_prodMat(X, matrix(B, ncol = 1))[,1]
  
  z_score <- qnorm(significance_level/2, lower.tail = FALSE)
  BLUE_ci_lower <- BLUE_prs - z_score * sqrt(BLUE_prs_var)
  BLUE_ci_upper <- BLUE_prs + z_score * sqrt(BLUE_prs_var)
  
  target_y <- fread(target_pheno, header = F)
  target_IID <- target_y$V2
  
  output <- data.frame("IID" = target_IID, 
                       "PRS" = BLUE_prs, 
                       "Variance" = BLUE_prs_var,
                       "Lower_Limit" = BLUE_ci_lower,
                       "Upper_Limit" = BLUE_ci_upper)
  return(output)
}