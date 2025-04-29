# iPRSue
---

Version: 1.0.0

Authors: Dovini Jayasinghe and S. Hong Lee

---

iPRSue is a novel method for estimating polygenic risk scores (PRSs) and their uncertainty (variance) using GWAS summary statistics. It can be applied to both quantitative and binary traits through the functions ```iPRSue_estimates_QT()``` and ```iPRSue_estimates_BT()```, providing unbiased and precise estimates. In addition, we have incorporated the traditional approach for PRS and variance estimation, available for both trait types via ```traditional_estimates_QT()``` and ```traditional_estimates_BT()```. While the traditional method offers slightly but significantly higher prediction accuracy compared to iPRSue, it is limited to cases where the number of SNP genotypes is smaller than the discovery sample size. In contrast, iPRSue does not have this limitation and can be applied to any data dimension.


## GitHub installation

Install devtools:
```
install.packages("devtools")
```
Install iPRSue:
```
devtools::install_github("DoviniJ/iPRSue")
```
Call the library iPRSue:
```
library(iPRSue)
```

## Contact
dovini.jayasinghe@mymail.unisa.edu.au
