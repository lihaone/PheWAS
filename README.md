# PheWAS
Contains the R script used in the paper to perform systems genetics analysis, e.g. PheWAS, ePheWAS, and mediation analysis. 

## Dependencies: 

Install `emma` from http://mouse.cs.ucla.edu/emma/. 

```sh
    wget http://mouse.cs.ucla.edu/emma/emma_1.1.2.tar.gz
    R CMD INSTALL emma_1.1.2.tar.gz
```

Install other dependent packages `plyr`, `qtl`, `mlmm`, `intermediate`
```R
    install.packages("devtools")
    install.packages("plyr")
    install.packages("qtl")
    
    library(devtools)
    install_github("Gregor-Mendel-Institute/mlmm")
    install_github("simecek/intermediate")

```
