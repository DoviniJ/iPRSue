#' GWAS_QT function
#'
#' Performs univariate linear regression for each SNP to estimate effect sizes and standard errors
#' in a genome-wide association study (GWAS) for a quantitative trait.
#'
#' @param discovery_pheno Character. Path to the phenotype file. This file should be tab- or space-delimited,
#' with no header, and the quantitative phenotype should be located in the third column.
#' @param discovery_geno_mat Character. Path to the genotype matrix file. This file should also be
#' delimited with no header, and each column corresponds to a SNP (e.g., encoded as 0, 1, 2).
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{beta}{Estimated effect size from linear regression.}
#'   \item{se}{Standard error of the effect size estimate.}
#' }
#'
#' @details
#' The function uses linear regression (`lm`) to regress the quantitative phenotype on each SNP separately.
#' The phenotype is standardized prior to analysis. No covariates are included in the model.
#' The genotype matrix and phenotype vector are assumed to be ordered consistently.
#'
#' @import data.table
#' 
#' @importFrom stats binomial coef glm lm qnorm quantile rnorm var vcov
#' 
#' @examples
#' \donttest{
#'   # Example usage:
#'   # Phenotype file: 3rd column must contain a continuous outcome
#'   # Genotype file: SNPs in columns, rows correspond to individuals
#'
#'   qpd <- system.file("Qpd.txt", package = "iPRSue", mustWork = TRUE)
#'   gd  <- system.file("Gd.txt",  package = "iPRSue", mustWork = TRUE)
#'
#'   results <- GWAS_QT(
#'     discovery_pheno    = qpd,
#'     discovery_geno_mat = gd
#'   )
#'   head(results)
#' }
#'
#' @export
GWAS_QT <- function(discovery_pheno, discovery_geno_mat){
  
  
  phe <- fread(discovery_pheno, header=F)
  y <- scale(phe$V3)
  X_discovery <- fread(discovery_geno_mat, header=F)
  df <- as.data.frame(cbind(y,X_discovery))
  
  beta <- c()
  se <- c()
  for(i in 1:ncol(X_discovery)){
    m <- lm(y~df[,i+1], data=df)
    beta[i] <- summary(m)$coefficients[2,1]
    se[i] <- summary(m)$coefficients[2,2]
  }
  
  gwas <- as.data.frame(cbind(beta, se))
  return(gwas)
}

