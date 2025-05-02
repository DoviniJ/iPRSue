# iPRSue (repository in preparation)
---

Version: 1.0.0

Authors: Dovini Jayasinghe and S. Hong Lee

---

The iPRSue is a novel method for estimating polygenic risk scores (PRSs) and their uncertainty (variance) using Genome Wide Association Studies (GWAS) summary statistics which can be constructed by ```GWAS_QT()``` or ```GWAS_BT()``` functions (see the section $\color{red} {IMPORTANT}$ for details). This method can be applied to both quantitative and binary traits through the functions ```iPRSue_estimates_QT()``` and ```iPRSue_estimates_BT()```, providing unbiased and precise estimates. In addition, we have incorporated a traditional approach for PRS and variance estimation, available for both trait types via ```traditional_estimates_QT()``` and ```traditional_estimates_BT()```. While the traditional method offers slightly but significantly higher prediction accuracy compared to iPRSue, it is limited to cases where the number of SNP genotypes is smaller than the discovery sample size. In contrast, iPRSue does not have this limitation and can be applied to any data dimension.


## GitHub installation

Install devtools:
```
install.packages("devtools")
```
Install iPRSue:
```
devtools::install_github("DoviniJ/iPRSue")
```
### Load the library
```
library(iPRSue)
```

## Input files

1) Qpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
* family ID (FID) 
* individual ID (IID)  
* quantitative phenotype of the discovery sample

```
ID_1 ID_1 31.6534
ID_2 ID_2 25.5035
ID_3 ID_3 26.7391
ID_4 ID_4 25.5271
ID_5 ID_5 26.7165
```

2) Qpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* quantitative phenotype of the target sample - optional

```
ID_801 ID_801 26.5723
ID_802 ID_802 20.2632
ID_803 ID_803 27.7365
ID_804 ID_804 18.75
ID_805 ID_805 23.3025
```

3) Bpd.txt - This is a .txt file which contains the following columns in order. The discovery dataset has 800 individuals. Note that the file has no column headings.
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

4) Bpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* binary phenotype (0=controls, 1=cases) of the target sample - optional

```
ID_801 ID_801 0
ID_802 ID_802 1
ID_803 ID_803 0
ID_804 ID_804 0
ID_805 ID_805 1
```

5) Gd.txt - This is a .txt file which contains scaled (column standardized) genotype matrix with the dimension discovery_sample_size x number_of_snps (e.g. 800 x 100), of the discovery individuals. Note that the file has neither row nor column headings.
   
6) Gt.txt - This is a .txt file which contains scaled (column standardized) genotype matrix with the dimension target_sample_size x number_of_snps (e.g. 200 x 100), of the target individuals. Note that the file has neither row nor column headings.

## Quick start (tutorial)
Download all the above example data files from inst folder to your working directory. Follow the commands below to obtain the expected outputs. Once practiced, you may use your own data files to generate individual PRS confidence intervals by using either the novel method iPRSue or the traditional method.

## Output files
#### When the outcome variable is quantitative,
**Commands**
```
x <- GWAS_QT(discovery_pheno = "Qpd.txt", discovery_geno_mat = "Gd.txt")
y <- iPRSue_estimates_QT(gwas = x, target_pheno = "Qpt.txt", target_geno_mat = "Gt.txt", no_of_PRSs = 500, significance_level = 0.05, seed = set.seed(1))
```
The arguments of ```GWAS_QT()``` function, namely, ```discovery_pheno``` and ```discovery_geno_mat``` should be the text file names of phenotype and scaled genotype matrix for the individuals in the discovery dataset. The function ```iPRSue_estimates()``` use the estimated values from ```GWAS_QT()``` function as an input to the argument ```gwas```. Moreover, ```target_pheno``` and ```target_geno_mat``` should be the text file names of phenotype and scaled genotype matrix for the individuals in the target dataset. The no_of_PRSs specifies the size of the PRS distribution per target individual (default is 500), significance_level sets the level of significance for constructing PRS confidence intervals (default is 0.05), and seed is used to control random number generation for reproducibility (default is NULL).

**Outputs**
```
          beta         se
1 -0.005625531 0.03539906
2  0.011408607 0.03539731
3 -0.024074395 0.03538936
4 -0.021290422 0.03539159
5  0.032671483 0.03538072
6  0.019838674 0.03539265
```
Returns the dataframe ```x``` with following columns:

* x$beta : Additive effects of scaled SNP genotypes
* x$se : Standard errors of the additive effects
  
```
     IID         PRS  Variance Lower_Limit Upper_Limit
1 ID_801 -0.10518765 0.1231988  -0.7764618   0.5997384
2 ID_802 -0.34199622 0.1075602  -0.9481694   0.2861330
3 ID_803 -0.35311776 0.1041257  -0.9409285   0.2779039
4 ID_804  0.07428962 0.1257936  -0.6987473   0.6993487
5 ID_805 -0.03311903 0.1191919  -0.7008964   0.6027532
6 ID_806 -0.04271522 0.1294443  -0.7324805   0.6695980
```
Returns the dataframe ```y``` with following columns:

* y$IID : Target individual IDs
* y$PRS : PRS estimates of each target individual, computed using iPRSue method
* y$Variance : Variance of individual PRS, computed using iPRSue method
* y$Lower_Limit : Lower limit of the individual PRS confidence interval
* y$Upper_Limit : Lower limit of the individual PRS confidence interval

**Command**
```
z <- traditional_estimates_QT(discovery_pheno = "Qpd.txt", discovery_geno_mat = "Gd.txt", target_pheno = "Qpt.txt", target_geno_mat = "Gt.txt", significance_level = 0.05)
```
The function ```traditional_estimates_QT()``` utilizes individual level data and provides PRS and uncertainty estimates using the traditional multiple linear regression approach. 

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
* y$PRS : PRS estimates of each target individual, computed using traditional method
* y$Variance : Variance of individual PRS, computed using traditional method
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
z <- traditional_estimates_BT(discovery_pheno = "Bpd.txt", discovery_geno_mat = "Gd.txt", target_pheno = "Bpt.txt", target_geno_mat = "Gt.txt", significance_level = 0.05)
```
The function ```traditional_estimates_BT()``` utilizes individual level data and provides PRS and uncertainty estimates using the traditional multiple logistic regression approach. 

**Output**
```
     IID        PRS Variance Lower_Limit Upper_Limit
1 ID_801   8.829390 17.43574   0.6453403   17.013440
2 ID_802  -5.406598 16.33275 -13.3275551    2.514360
3 ID_803  10.227994 15.27634   2.5674829   17.888506
4 ID_804 -11.042414 20.39588 -19.8939645   -2.190863
5 ID_805   8.724172 18.34187   0.3301530   17.118191
6 ID_806   6.061409 14.39817  -1.3756594   13.498476
```
Returns the dataframe ```z``` with following columns:

* y$IID : Target individual IDs
* y$PRS : PRS estimates of each target individual, computed using traditional method
* y$Variance : Variance of individual PRS, computed using traditional method
* y$Lower_Limit : Lower limit of the individual PRS confidence interval
* y$Upper_Limit : Lower limit of the individual PRS confidence interval


## $\color{red} {IMPORTANT}$
We recommend producing GWAS summary statistics using scaled genotypes when applying iPRSue method. If the users have access to readily-available GWAS summary statistics which are produced using unscaled genotypes, make sure to conduct necessary adjustments to those SNP effects ($beta$) and corresponding standard errors ($se$), prior to applying ```iPRSue_estimates_QT()``` or ```iPRSue_estimates_BT()```. 

The adjustment can be done using minor allele frequency ($p$) or effective number of individuals ($n$) information as follows:
* $beta_{adjusted} = beta / (se . \sqrt{2p(1-p)})$
* $se_{adjusted} = 1 / \sqrt{2p(1-p)}$

or
* $beta_{adjusted} = beta / (se . \sqrt{n})$
* $se_{adjusted} = 1 / \sqrt{n}$

  where n is the effective number of individuals (ref. 1).


### References
1) https://github.com/euffelmann/bpc

   
## Contact
dovini.jayasinghe@mymail.unisa.edu.au
