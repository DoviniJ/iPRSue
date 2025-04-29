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
# Load the library
```
library(iPRSue)
```
# Data preparation

## File formats
### Input files

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
* quantitative phenotype of the target sample

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
* binary phenotype (1=controls, 2=cases) of the discovery sample

```
ID_1 ID_1 1
ID_2 ID_2 2
ID_3 ID_3 1
ID_4 ID_4 1
ID_5 ID_5 1
```

4) Bpt.txt - This is a .txt file which contains the following columns in order. The target dataset has 200 individuals who are independent from the discovery dataset. Note that the file has no column headings.   
* FID 
* IID  
* binary phenotype (0=controls, 1=cases) of the target sample

```
ID_801 ID_801 0
ID_802 ID_802 1
ID_803 ID_803 0
ID_804 ID_804 0
ID_805 ID_805 1
```

5) Gd.txt - This is a .txt file which contains scaled (column standardized) genotype matrix with the dimension discovery_sample_size x number_of_snps (e.g. 800 x 100), of the discovery individuals. Note that the file has neither row nor column headings.
   
6) Gt.txt - This is a .txt file which contains scaled (column standardized) genotype matrix with the dimension target_sample_size x number_of_snps (e.g. 200 x 100), of the target individuals. Note that the file has neither row nor column headings.





## Contact
dovini.jayasinghe@mymail.unisa.edu.au
