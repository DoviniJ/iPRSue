# iPRSue (Page under construction - please use CRAN)
---

Version: xxx

Authors: Dovini Jayasinghe and S. Hong Lee

---

The iPRSue is a novel method for estimating polygenic risk scores (PRSs) and their uncertainty (variance) using Genome Wide Association Studies (GWAS) summary statistics which can be constructed by ```GWAS_QT()``` or ```GWAS_BT()``` functions (see the section $\color{red} {IMPORTANT}$ for details). This method can be applied to both quantitative and binary traits through the functions ```iPRSue_estimates_QT()``` and ```iPRSue_estimates_BT()```, providing unbiased and precise estimates. In addition, we have incorporated traditionally used Best linear unbiased estimator (BLUE) method for PRS and variance estimation, available for both trait types via ```BLUE_estimates_QT()``` and ```BLUE_estimates_BT()```. While the BLUE method offers slightly but significantly higher prediction accuracy compared to iPRSue, it is limited to cases where the number of SNP genotypes is smaller than the discovery sample size. In contrast, iPRSue does not have this limitation and can be applied to any data dimension.

## Package Installation:

### CRAN 

```
install.packages("iPRSue")
```
### or

### GitHub 

Install devtools:
```
install.packages("devtools")
```
Install iPRSue:
```
devtools::install_github("DoviniJ/iPRSue")
```
## Load the library:
```
library(iPRSue)
```

### Package dependencies:
The package iPRSue depends on the following R packages. Please install them and call the respective libraries.

```
install.packages(c("data.table", "bigstatsr"))
library(data.table)
library(bigstatsr)
```

## Input files

1) Qpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
* family ID (FID) 
* individual ID (IID)  
* scaled quantitative phenotype of the discovery sample

```
ID_1 ID_1 0.967866365493657
ID_2 ID_2 -0.401199467463008
ID_3 ID_3 -0.1261351992308
ID_4 ID_4 -0.395945731035134
ID_5 ID_5 -0.131166319708339
```

2) Qpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* scaled quantitative phenotype of the target sample - optional

```
ID_801 ID_801 -0.111301037986145
ID_802 ID_802 -1.43931139178665
ID_803 ID_803 0.133752876818589
ID_804 ID_804 -1.75782675330745
ID_805 ID_805 -0.799565331602807
```

3) Bpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
* FID 
* IID  
* binary phenotype (1=controls, 2=cases) of the discovery sample

```
ID_1 ID_1 1
ID_2 ID_2 2
ID_3 ID_3 1
ID_4 ID_4 1
ID_5 ID_5 1
```

4) Bpd_0_1.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
* FID 
* IID  
* binary phenotype (0=controls, 1=cases) of the discovery sample

```
ID_1 ID_1 0
ID_2 ID_2 1
ID_3 ID_3 0
ID_4 ID_4 0
ID_5 ID_5 0
```

5) Bpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* binary phenotype (1=controls, 2=cases) of the target sample - optional

```
ID_801 ID_801 1
ID_802 ID_802 2
ID_803 ID_803 1
ID_804 ID_804 1
ID_805 ID_805 2
```

6) cd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.    
* FID 
* IID 
* 16 confounders of the discovery sample (Note: These columns are optional (at least 1 column is required). Can use any number of columns as confounders to adjust the phenotype upon user requirement.)

```
ID_1 ID_1 0.787403812314451 0.620004763647331 -3.04026 45 -12.048 2.17634 -0.940322 -0.446351 -5.45685 -2.53161 -2.13435 -1.95623 -3.82792 -0.380636 0 10
ID_2 ID_2 -0.119722138532781 0.0143333904548625 -5.1054 64 -14.5169 6.01889 -3.85803 3.62625 5.10717 -3.54574 0.393994 3.64275 4.42975 -2.26704 1 19
ID_3 ID_3 0.173372721351375 0.0300581005087816 -1.91044 59 -12.7462 5.95244 0.0738914 1.80523 4.76284 0.130369 -1.05615 0.316777 0.988783 -1.76502 1 7
ID_4 ID_4 -0.699321184695051 0.48905011936329 -1.83526 68 -10.3349 4.71264 -1.84521 -0.524855 -3.80275 0.837965 0.265233 2.10903 -0.210259 0.71504 0 20
ID_5 ID_5 3.69300366651739 13.6382760809109 -3.15649 69 -8.56737 4.78248 -1.49547 -7.49413 -5.39887 1.85316 4.07476 1.05351 0.825942 -2.09669 1 20
```
Note: If you have no confounder, then leave the 3rd column constant (e.g. vector of 1s)

cd1.txt
* FID 
* IID 
* vector of 1s

```
ID_1 ID_1 1
ID_2 ID_2 1
ID_3 ID_3 1
ID_4 ID_4 1
ID_5 ID_5 1
```  

7) Gd.txt - This is a .txt file which contains scaled (column standardized) genotype matrix with the dimension discovery_sample_size x number_of_snps (e.g. 800 x 1,000), of the discovery individuals. Note that the file has neither row nor column headings.
   
8) Gt.txt - This is a .txt file which contains scaled (column standardized) genotype matrix with the dimension target_sample_size x number_of_snps (e.g. 200 x 1,000), of the target individuals. Note that the file has neither row nor column headings.

9) mydata.fam - This is a file associated with the PLINK binary format file which contains the following columns in order. The example dataset has 1,000 individuals. Note that the file has no column headings. This follows the PLINK .fam file format.
* family ID (FID) 
* individual ID (IID) 
* father's ID 
* mother's ID 
* sex 
* phenotype value

```
ID_1 ID_1 0 0 1 -9
ID_2 ID_2 0 0 2 -9
ID_3 ID_3 0 0 2 -9
ID_4 ID_4 0 0 2 -9
ID_5 ID_5 0 0 1 -9
```
  
10) mydata.bim - This is is a file associated with the PLINK binary format file which contains the following columns in order. The example dataset has 1,000 SNPs. Note that the file has no column headings. This follows the PLINK .bim file format.
* chromosome code 
* SNP ID 
* position of centimorgans 
* base-pair coordinate 
* minor allele  
* reference allele 

```
1	SNP_1	0	768448	A	G
1	SNP_2	0	853954	C	A
1	SNP_3	0	880390	A	C
1	SNP_4	0	940203	A	G
1	SNP_5	0	987670	T	G
```

11) mydata.bed - This is the PLINK binary format file which includes genotype information. This follows the PLINK .bed file format.

## Quick start (tutorial)

### Preparation

##### Download all example data from inst directory
Download all the above example data files from inst folder to your working directory. You can copy the path to inst folder ``` https://github.com/DoviniJ/iPRSue/tree/main/inst ``` and download the directory using https://download-directory.github.io/ page, for instant downloading. Follow the commands below to obtain the expected outputs. Once practiced, you may use your own data files to generate individual PRS confidence intervals by using either the novel method iPRSue or the BLUE method.

##### Download and install plink2 in your machine. The package supports both Linux and Windows version.
Link: https://www.cog-genomics.org/plink/2.0/

##### Obtain the path to the executable plink application <plink_path>

Give the path where plink executable file is located
```
plink_path <- "<plink_path>/plink2" 
```
###### Additional note:
_Note that, all these files (user's own data files) can be placed in a separate location. It is always upto the user's choice. In that case remember to give the full path to the file location since R identifies files by name, only when they are in the same directory. For example, if your data files are in a different directory;_

* instead of ```"mydata"``` use ```"<path>/mydata"```
  
### Commands and outputs
#### When the outcome variable is quantitative,
**Commands**
```
x <- GWAS_QT(plink_path = "plink2", b_file = "mydata", discovery_pheno = "Qpd.txt", discovery_cov = "cd.txt", thread = 48)
y <- iPRSue_estimates_QT(gwas = x, target_pheno = "Qpt.txt", target_geno_mat = "Gt.txt", no_of_PRSs = 500, significance_level = 0.05, seed = set.seed(1))
```
The arguments of ```GWAS_QT()``` function, namely, ```discovery_pheno``` and ```discovery_geno_mat``` should be the text file names of phenotype and scaled genotype matrix for the individuals in the discovery dataset. The function ```iPRSue_estimates()``` use the estimated values from ```GWAS_QT()``` function as an input to the argument ```gwas```. Moreover, ```target_pheno``` and ```target_geno_mat``` should be the text file names of phenotype and scaled genotype matrix for the individuals in the target dataset. The no_of_PRSs specifies the size of the PRS distribution per target individual (default is 500), significance_level sets the level of significance for constructing PRS confidence intervals (default is 0.05), and seed is used to control random number generation for reproducibility (default is NULL).

**Outputs**
```
           beta         se
1  0.0182285503 0.03568526
2 -0.0517348657 0.03601479
3  0.0067835942 0.03327651
4 -0.0659989174 0.03386520
5 -0.0004306981 0.03517784
6  0.0284129553 0.03487432
```
Returns the dataframe ```x``` with following columns:

* x$beta : Additive effects of scaled SNP genotypes
* x$se : Standard errors of the additive effects
  
```
     IID        PRS Variance Lower_Limit Upper_Limit
1 ID_801 -0.8710867 1.326542  -3.1149828   1.2669814
2 ID_802 -2.9520295 1.194713  -5.2039386  -0.9847323
3 ID_803  0.2631768 1.382384  -2.1256907   2.5201000
4 ID_804  2.7146966 1.469916   0.2658072   5.1874619
5 ID_805  2.3592619 1.364946   0.1601901   4.6957551
6 ID_806 -1.9697418 1.210845  -4.0841447   0.2151955
```
Returns the dataframe ```y``` with following columns:

* y$IID : Target individual IDs
* y$PRS : PRS estimates of each target individual, computed using iPRSue method
* y$Variance : Variance of individual PRS, computed using iPRSue method
* y$Lower_Limit : Lower limit of the individual PRS confidence interval
* y$Upper_Limit : Lower limit of the individual PRS confidence interval

TO BE CONTINUED AFTER TEA FROMM HERE .................................................................

**Command**
```
z <- BLUE_estimates_QT(discovery_pheno = "Qpd.txt", discovery_geno_mat = "Gd.txt", target_pheno = "Qpt.txt", target_geno_mat = "Gt.txt", significance_level = 0.05)
```
The function ```BLUE_estimates_QT()``` utilizes individual level data and provides PRS and uncertainty estimates using the BLUE multiple linear regression approach. 

**Output**
```
     IID           PRS  Variance Lower_Limit Upper_Limit
1 ID_801 -0.1444712318 0.1360387  -0.8673731   0.5784306
2 ID_802 -0.3937930417 0.1255455  -1.0882554   0.3006693
3 ID_803 -0.5239819895 0.1410622  -1.2601100   0.2121460
4 ID_804  0.1379861095 0.1386364  -0.5917851   0.8677573
5 ID_805 -0.0009677077 0.1527249  -0.7669223   0.7649869
6 ID_806 -0.2428045632 0.1490595  -0.9995118   0.5139027
```
Returns the dataframe ```z``` with following columns:

* y$IID : Target individual IDs
* y$PRS : PRS estimates of each target individual, computed using BLUE method
* y$Variance : Variance of individual PRS, computed using BLUE method
* y$Lower_Limit : Lower limit of the individual PRS confidence interval
* y$Upper_Limit : Lower limit of the individual PRS confidence interval
  
#### When the outcome variable is binary,
**Commands**
```
x <- GWAS_BT(discovery_pheno = "Bpd.txt", discovery_geno_mat = "Gd.txt")
y <- iPRSue_estimates_BT(gwas = x, target_pheno = "Bpt.txt", target_geno_mat = "Gt.txt", no_of_PRSs = 500, significance_level = 0.05, seed = set.seed(1))
```
The arguments of ```GWAS_BT()``` function, namely, ```discovery_pheno``` and ```discovery_geno_mat``` should be the text file names of phenotype and scaled genotype matrix for the individuals in the discovery dataset. The function ```iPRSue_estimates()``` use the estimated values from ```GWAS_BT()``` function as an input to the argument ```gwas```. Moreover, ```target_pheno``` and ```target_geno_mat``` should be the text file names of phenotype and scaled genotype matrix for the individuals in the target dataset. The no_of_PRSs specifies the size of the PRS distribution per target individual (default is 500), significance_level sets the level of significance for constructing PRS confidence intervals (default is 0.05), and seed is used to control random number generation for reproducibility (default is NULL).

**Outputs**
```
         beta        se
1  0.02351554 0.1734839
2  0.19951530 0.1609413
3  0.01422655 0.1742960
4 -0.03176976 0.1777684
5 -0.18495543 0.1921530
6 -0.10750075 0.1847582
```
Returns the dataframe ```x``` with following columns:

* x$beta : Additive effects of scaled SNP genotypes
* x$se : Standard errors of the additive effects
  
```
     IID       PRS Variance Lower_Limit Upper_Limit
1 ID_801  2.876907 3.185277  -0.4678476   6.5059090
2 ID_802  1.802006 2.769852  -1.3031663   4.9981202
3 ID_803  2.431920 2.697261  -0.4885771   5.6501721
4 ID_804 -2.524441 3.273875  -6.4328362   0.7284884
5 ID_805  2.363432 3.063884  -1.1095964   5.6116642
6 ID_806  1.497236 3.322266  -2.0080071   5.1538390
```
Returns the dataframe ```y``` with following columns:

* y$IID : Target individual IDs
* y$PRS : PRS estimates of each target individual, computed using iPRSue method
* y$Variance : Variance of individual PRS, computed using iPRSue method
* y$Lower_Limit : Lower limit of the individual PRS confidence interval
* y$Upper_Limit : Lower limit of the individual PRS confidence interval

**Command**
```
z <- BLUE_estimates_BT(discovery_pheno = "Bpd.txt", discovery_geno_mat = "Gd.txt", target_pheno = "Bpt.txt", target_geno_mat = "Gt.txt", significance_level = 0.05, max_iterations = 100)
```
The function ```BLUE_estimates_BT()``` utilizes individual level data and provides PRS and uncertainty estimates using the BLUE multiple logistic regression approach with firth's bias correction.  

**Output**
```
     IID       PRS Variance Lower_Limit Upper_Limit
1 ID_801  2.701390 2.194493  -0.2020660   5.6048453
2 ID_802 -1.192601 2.048303  -3.9976810   1.6124780
3 ID_803  3.872388 2.347984   0.8691089   6.8756664
4 ID_804 -3.218260 2.365584  -6.2327734  -0.2037462
5 ID_805  2.755695 2.486593  -0.3349594   5.8463492
6 ID_806  1.973755 2.089906  -0.8596682   4.8071788

```
Returns the dataframe ```z``` with following columns:

* y$IID : Target individual IDs
* y$PRS : PRS estimates of each target individual, computed using BLUE method
* y$Variance : Variance of individual PRS, computed using BLUE method
* y$Lower_Limit : Lower limit of the individual PRS confidence interval
* y$Upper_Limit : Lower limit of the individual PRS confidence interval


## $\color{red} {IMPORTANT}$
We recommend producing GWAS summary statistics using scaled genotypes when applying iPRSue method. If the users have access to readily-available GWAS summary statistics which are produced using unscaled genotypes, make sure to conduct necessary adjustments to those SNP effects ($beta$) and corresponding standard errors ($se$), prior to applying ```iPRSue_estimates_QT()``` or ```iPRSue_estimates_BT()```. 

<!-- 
The adjustment can be done using minor allele frequency ($p$) or effective sample size ($n$) information as follows (ref. 1, 2):
* $beta_{adjusted} = beta . \sqrt{2p(1-p)}$
* $se_{adjusted} = se . \sqrt{2p(1-p)}$ 

or
* $beta_{adjusted} = beta / (se . \sqrt{n})$
* $se_{adjusted} = 1 / \sqrt{n}$


### References
1) Yang, J., Loos, R., Powell, J. et al. FTO genotype is associated with phenotypic variability of body mass index. Nature 490, 267–272 (2012). https://doi.org/10.1038/nature11401
2) https://github.com/euffelmann/bpc


The adjustment can be done using effective sample size ($n$) information as follows (ref. 1):

* $beta_{adjusted} = beta / (se . \sqrt{n})$
* $se_{adjusted} = 1 / \sqrt{n}$

### Reference
1) https://github.com/euffelmann/bpc
-->

The adjustment can be done using minor allele frequency ($p$) information as follows:
* $beta_{adjusted} = beta . \sqrt{2p(1-p)}$
* $se_{adjusted} = se . \sqrt{2p(1-p)}$

   
## Contact
dovini.jayasinghe@mymail.unisa.edu.au
