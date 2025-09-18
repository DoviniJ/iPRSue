#' GWAS_QT function
#'
#' Performs genome-wide association study (GWAS) using plink2 linear model and outputs 
#' the GWAS summary statistics with additive SNP effects (beta) and 
#' standard errors (se)
#'
#' @param plink_path Path to the PLINK executable application
#' @param b_file Prefix of the binary files, where all .fam, .bed and .bim files have a common prefix
#' @param discovery_pheno Name (with file extension) of the phenotype file containing family ID, individual 
#' ID and phenotype of the discovery dataset as columns, without heading
#' @param discovery_cov Name (with file extension) of the covariate file containing family ID, individual 
#' ID, and covariate(s) of the discovery dataset as columns, without heading. If no covariate is 
#' used, have a constant column (e.g. vector of 1s)
#' @param thread Number of threads used (default = 20)
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{beta}{Estimated effect size from linear regression.}
#'   \item{se}{Standard error of the effect size estimate.}
#' }
#'
#' @details
#' The function uses linear regression to regress the quantitative phenotype on each SNP separately
#' using plink 2. Then the estimated effects and standard errors are adjusted for standardization.
#' The phenotype is standardized prior to analysis. It is optional to employ covariates in the model. If no
#' covariate is used, create your covariate file with a constant in the 3rd column (e.g. vector of 1s).
#'
#' @importFrom utils read.table
#' 
#' @examples
#' \dontrun{
#'   results <- GWAS_QT(
#'     plink_path = "./plink2",
#'     b_file = "./binary_file_prefix",
#'     discovery_pheno = "./discovery_phenotype_file",
#'     discovery_cov = "./discovery_covariate_file",
#'     thread = 48
#'   )
#'   head(results)
#' }
#'
#' @export
GWAS_QT <- function(plink_path, b_file, discovery_pheno, discovery_cov, thread = 20){
  os_name <- Sys.info()["sysname"]
  if (startsWith(os_name, "Win")) {
    slash <- paste0("\\")
  } else {
    slash <- paste0("/")
  }
  cov_file <- read.table(discovery_cov)
  n_confounders = ncol(cov_file) - 2
  parameters <- 1:(n_confounders+1)
  param_vec <- paste0(parameters, collapse = ", ")
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  log_file <- runPLINK(paste0(" --bfile ", b_file, 
                              " --linear interaction --pheno ", 
                              discovery_pheno, 
                              " --covar ", discovery_cov, 
                              " --parameters ", param_vec, 
                              " --allow-no-sex --threads ", 
                              thread,
                              " --out ", tempdir(), slash, "Q_gwas"))
  first_line <- readLines(paste0(tempdir(), slash, "Q_gwas.PHENO1.glm.linear"), n = 1)
  col_names <- strsplit(first_line, "\t")[[1]]
  col_names[1] <- sub("#", "", col_names[1])
  plink_output <- read.table(paste0(tempdir(), slash, "Q_gwas.PHENO1.glm.linear"), , skip = 1, col.names = col_names, sep = "\t")
  filtered_output <- plink_output[(plink_output$TEST=="ADD"),]
  Q_out.trd.sum <- filtered_output[c("CHROM", "POS", "ID", "REF", "ALT", "A1", "OBS_CT", "BETA", "SE", "T_STAT", "P")]
  colnames(Q_out.trd.sum) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "OBS_CT", "BETA", "SE", "T_STAT", "P")
  rownames(Q_out.trd.sum) <- NULL
  gwas <- as.data.frame(Q_out.trd.sum[c("BETA", "SE")])
  runPLINK <- function(PLINKoptions = "") system(paste(plink_path, PLINKoptions))
  log_file2 <- runPLINK(paste0(" --bfile ", b_file,
                              " --freq --threads ",
                              thread,
                              " --out ", tempdir(), slash, "allelefreqs"))
  first_line2 <- readLines(paste0(tempdir(), slash, "allelefreqs.afreq"), n = 1)
  col_names2 <- strsplit(first_line2, "\t")[[1]]
  col_names2[1] <- sub("#", "", col_names2[1])
  plink_output2 <- read.table(paste0(tempdir(), slash, "allelefreqs.afreq"), , skip = 1, col.names = col_names2, sep = "\t")
  filtered_output2 <- plink_output2$ALT_FREQS
  tpq = 2*filtered_output2*(1-filtered_output2)
  beta = gwas$BETA*sqrt(tpq)
  se = gwas$SE*sqrt(tpq)
  
  adj_gwas <- data.frame(beta, se)
  return(adj_gwas)
}
