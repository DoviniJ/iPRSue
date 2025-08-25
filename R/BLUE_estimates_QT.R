#' BLUE_estimates_QT function
#'
#' Estimates individual-level polygenic risk scores (PRS) with uncertainty using a frequentist approach
#' for quantitative traits. This implementation fits a multiple linear regression model in the discovery dataset,
#' computes the coefficient covariance matrix, and applies the delta method to propagate uncertainty to the target dataset.
#'
#' @param discovery_pheno Character. Path to the phenotype file for the discovery dataset. Assumes no header and that the quantitative trait is in the third column.
#' @param discovery_geno_mat Character. Path to the genotype matrix file for the discovery dataset. Assumes no header.
#' @param target_pheno Character. Path to the phenotype file for the target dataset. Assumes no header and individual IDs in the second column.
#' @param target_geno_mat Character. Path to the genotype matrix file for the target dataset. Assumes no header.
#' @param significance_level Numeric. Significance level for confidence intervals (e.g., 0.05 for 95% CI). Default is 0.05.
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
#' The function fits a multiple linear regression model (`lm`) using the discovery data. The estimated SNP effects
#' and their covariance matrix are used to compute PRS and associated uncertainty for each individual in the target dataset.
#' Confidence intervals are constructed using the normal approximation.
#'
#' Missing or non-estimable coefficients and variances are set to zero.
#'
#' @import data.table
#' @import bigstatsr
#' 
#' @importFrom stats binomial coef glm lm qnorm quantile rnorm var vcov
#'
#' @examples
#' \donttest{
#'   qpd <- system.file("Qpd.txt", package = "iPRSue", mustWork = TRUE)
#'   qpt <- system.file("Qpt.txt", package = "iPRSue", mustWork = TRUE)
#'   gd  <- system.file("Gd.txt",  package = "iPRSue", mustWork = TRUE)
#'   gt  <- system.file("Gt.txt",  package = "iPRSue", mustWork = TRUE)
#'
#'   results <- BLUE_estimates_QT(
#'     discovery_pheno    = qpd,
#'     discovery_geno_mat = gd,
#'     target_pheno       = qpt,
#'     target_geno_mat    = gt,
#'     significance_level = 0.05
#'   )
#'   head(results)
#' }
#'
#' @export
BLUE_estimates_QT <- function(discovery_pheno, 
                                     discovery_geno_mat, 
                                     target_pheno,
                                     target_geno_mat, 
                                     significance_level = 0.05){
  phe <- fread(discovery_pheno, header=F)
  y <- scale(phe$V3)
  
  X_discovery <- as_FBM(as.matrix(fread(discovery_geno_mat, header=F)))
  df <- data.frame(y = y, X_discovery[])
  m <- lm(y ~ ., data=df)
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