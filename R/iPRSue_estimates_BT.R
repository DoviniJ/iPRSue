#' iPRSue_estimates_BT function
#'
#' Computes individual-level polygenic risk scores (PRS) with uncertainty estimates using a simulation-based approach
#' for binary traits. This implementation follows the iPRSue framework, simulating multiple PRSs by sampling from
#' the GWAS effect size distribution and deriving individual-level confidence intervals.
#'
#' @param gwas A data frame with GWAS summary statistics for binary traits. Must contain `beta` and `se` columns representing
#' estimated SNP effect sizes and their standard errors.
#' @param target_pheno Character. Path to the target phenotype file. Assumes no header and individual IDs in the second column.
#' @param target_geno_mat Character. Path to the genotype matrix of target individuals. No header is expected; columns correspond to SNPs.
#' @param no_of_PRSs Integer. Number of simulations used to construct PRS uncertainty intervals. Default is 500.
#' @param significance_level Numeric. Significance level for confidence intervals (e.g., 0.05 gives 95% CI). Default is 0.05.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL, results may vary across runs. Default is NULL.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{IID}{Individual identifier (from target phenotype file).}
#'   \item{PRS}{Mean of simulated PRSs for each individual.}
#'   \item{Variance}{Variance across simulated PRSs.}
#'   \item{Lower_Limit}{Lower bound of the confidence interval.}
#'   \item{Upper_Limit}{Upper bound of the confidence interval.}
#' }
#'
#' @details
#' For each SNP, the function simulates `no_of_PRSs` effect sizes from a normal distribution defined by its GWAS beta and SE.
#' These sampled betas are multiplied by the genotype matrix to generate PRS replicates for each individual.
#' Confidence intervals are then calculated using the specified significance level.
#'
#' This function is designed for binary traits and should be used with GWAS summary statistics obtained from logistic regression.
#'
#' @import data.table
#' @import bigstatsr
#' 
#' @importFrom stats binomial coef glm lm qnorm quantile rnorm var vcov
#'
#' @examples
#' \dontrun{
#' # Step 1: Run GWAS on binary trait
#' gwas_res <- GWAS_BT("Bpd.txt", "Gd.txt")
#'
#' # Step 2: Estimate individual PRS with uncertainty
#' prs_estimates <- iPRSue_estimates_BT(gwas = gwas_res,
#'                                      target_pheno = "Bpt.txt",
#'                                      target_geno_mat = "Gt.txt",
#'                                      no_of_PRSs = 500,
#'                                      significance_level = 0.05,
#'                                      seed = 123)
#' head(prs_estimates)
#' }
#'
#' @export
iPRSue_estimates_BT <- function(gwas, target_pheno,
                                target_geno_mat, 
                                no_of_PRSs = 500, 
                                significance_level = 0.05,
                                seed = NULL){
  gwas <- as.data.frame(gwas)
  X <- as_FBM(as.matrix(fread(target_geno_mat, header=F)))
  
  betas <- matrix(nrow = nrow(gwas), ncol = no_of_PRSs)
  replication = seed
  for (i in 1:no_of_PRSs) {
    betas[, i] <- rnorm(n = nrow(gwas), mean = gwas$beta, sd = gwas$se)
  }
  
  matrix <- big_prodMat(X, betas)
  
  proposed_prs_mean <- c()
  proposed_prs_var <- c()
  proposed_ci_lower <- c()
  proposed_ci_upper <- c()
  
  for(j in 1:nrow(matrix)){
    proposed_prs_mean[j] <- mean(matrix[j,])
    proposed_prs_var[j] <- var(matrix[j,])
    proposed_ci_lower[j] <- quantile(matrix[j,], significance_level/2)
    proposed_ci_upper[j] <- quantile(matrix[j,], 1-(significance_level/2))
  }
  
  target_y <- fread(target_pheno, header = F)
  target_IID <- target_y$V2
  
  output <- data.frame("IID" = target_IID, 
                       "PRS" = proposed_prs_mean, 
                       "Variance" = proposed_prs_var,
                       "Lower_Limit" = proposed_ci_lower,
                       "Upper_Limit" = proposed_ci_upper)
  return(output)
}

