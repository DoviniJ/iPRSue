# iPRSue (repository in preparation)
---

Version: 1.0.0

Authors: Dovini Jayasinghe and S. Hong Lee

---

iPRSue is a novel method for estimating polygenic risk scores (PRSs) and their uncertainty (variance) using Genome Wide Association Studies (GWAS) summary statistics. It can be applied to both quantitative and binary traits through the functions ```iPRSue_estimates_QT()``` and ```iPRSue_estimates_BT()```, providing unbiased and precise estimates. In addition, we have incorporated a traditional approach for PRS and variance estimation, available for both trait types via ```traditional_estimates_QT()``` and ```traditional_estimates_BT()```. While the traditional method offers slightly but significantly higher prediction accuracy compared to iPRSue, it is limited to cases where the number of SNP genotypes is smaller than the discovery sample size. In contrast, iPRSue does not have this limitation and can be applied to any data dimension.


## GitHub installation

Install devtools:
```
install.packages("devtools")
```
Install iPRSue:
```
devtools::install_github("DoviniJ/iPRSue")
```
## Load the library
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

## Output files
#### When the outcome variable is quantitative
**Command**
```
x <- iPRSue_estimates_QT(discovery_pheno = "Qpd.txt", target_pheno = "Qpt.txt",
                                discovery_geno_mat = "Gd.txt", target_geno_mat = "Gt.txt", 
                                no_of_PRSs = 500, significance_level = 0.05,
                                seed = set.seed(1))
```
**Output**
```
     IID        PRS Variance Lower_Limit Upper_Limit
1 ID_801 -0.4725073 2.485959   -3.487898    2.694049
2 ID_802 -1.5362611 2.170396   -4.259216    1.285321
3 ID_803 -1.5862195 2.101094   -4.226689    1.248356
4 ID_804  0.3337120 2.538318   -3.138801    3.141503
5 ID_805 -0.1487720 2.405106   -3.148455    2.707592
6 ID_806 -0.1918785 2.611984   -3.290333    3.007862
```
Returns the dataframe ```x``` with following columns:

* x$IID : Target individual IDs
* x$PRS : PRS estimates of each target individual, computed using iPRSue method
* x$Variance : Variance of individual PRS, computed using iPRSue method
* x$Lower_Limit : Lower limit of the individual PRS confidence interval
* x$Upper_Limit : Lower limit of the individual PRS confidence interval

## Contact
dovini.jayasinghe@mymail.unisa.edu.au
