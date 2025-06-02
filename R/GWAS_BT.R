#' GWAS_BT function
#'
#' Performs genome-wide association analysis for a binary trait using logistic regression.
#' It reads a phenotype file and a genotype matrix, and estimates the SNP effect sizes and standard errors.
#'
#' @param discovery_pheno A character string specifying the path to the phenotype file. The file should have no header and contain individual IDs, and the third column should contain the binary trait (0/1).
#' @param discovery_geno_mat A character string specifying the path to the genotype matrix file. The file should have no header and contain numeric genotype data (e.g., 0, 1, 2) for each SNP.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{beta}{Estimated effect size (log odds) for each SNP.}
#'   \item{se}{Standard error of the estimated effect size.}
#' }
#'
#' @details
#' The function uses logistic regression (`glm` with `binomial(link="logit")`) to regress the binary phenotype on each SNP individually. The output includes only the regression coefficient and standard error for each SNP.
#'
#' @import data.table
#' 
#' @importFrom stats binomial coef glm lm qnorm quantile rnorm var vcov
#' 
#' @examples
#' \dontrun{
#' # Phenotype file: 3rd column must contain binary outcome (0/1)
#' # Genotype file: SNPs in columns, rows correspond to individuals
#' # Run GWAS on a binary trait with discovery phenotype and genotype files
#' gwas_results <- GWAS_BT(discovery_pheno = "Bpd.txt", discovery_geno_mat = "Gd.txt")
#' head(gwas_results)
#' }
#'
#' @export
#' 
GWAS_BT <- function(discovery_pheno, discovery_geno_mat){
  
  
  phe <- fread(discovery_pheno, header=F)
  y <- phe$V3
  X_discovery <- fread(discovery_geno_mat, header=F)
  df <- as.data.frame(cbind(y,X_discovery))
  
  beta <- c()
  se <- c()
  for(i in 1:ncol(X_discovery)){
    m <- glm(y~df[,i+1], data=df, family=binomial(link="logit"))
    beta[i] <- summary(m)$coefficients[2,1]
    se[i] <- summary(m)$coefficients[2,2]
  }
  
  gwas <- as.data.frame(cbind(beta, se))
  return(gwas)
}