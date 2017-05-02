## Welcome to the GitHub page for the KSEAapp R Package!

### Overview

The KSEAapp R package offers a panel of functions for kinase activity inference on an input phosphoproteomics dataset using the Kinase-Substrate Enrichment Analysis (KSEA) method originally described by Casado *et al.* (see references below). 

For details on its usage and features, please read the Overview.Rmd file within the **vignettes** folder.

Note: if you are interested in performing KSEA through a user-friendly online interface, please visit [this link](https://casecpb.shinyapps.io/ksea/), or get details from our [other GitHub site](https://github.com/casecpb/KSEA/). 

-----
### Page Contents

1. **data** folder: collection of all the sample datasets in .RData format
2. **man** folder: collection of documentation files
3. **R** folder: collection of the raw R functions found within the package
4. **vignettes** folder: currently contains an Overview.Rmd file that walks readers through the functionality of the package. This is a good starting point for understanding how to use the available functions.
5. DESCRIPTION file
6. LICENSE file
7. NAMESPACE file


Elements found in 2, 5, 6, and 7 are key components for building the R package but are not crucial for understanding the basis of the KSEA calculations. 

-----
### Package Installation

There are 2 options for installing the KSEAapp package within the R console:

**(1) Via CRAN**

This gives access to the newest [CRAN release version](https://CRAN.R-project.org/package=KSEAapp). 
```
install.packages("KSEAapp")
```

**(2) Via GitHub**

This gives access to the newest development version. 
```
install.packages("devtools") # this installs the devtools R package 
devtools::install_github("casecpb/KSEAapp")
```


-----
References:

1. Wiredja D.D., Koyutürk M., Chance M.R. (2017) The KSEA App: a web-based tool for kinase activity inference from quantitative phosphoproteomics. *Submitted for review*.
2. [Casado, P., et al. (2013) Kinase-substrate enrichment analysis provides insights into the heterogeneity of signaling pathway activation in leukemia cells. *Sci. Signal*. 6, rs6-rs6](http://stke.sciencemag.org/content/6/268/rs6.long)
3. [Hornbeck P.V., et al. (2015) PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. *Nucleic Acids Res*. 43:D512-20](https://academic.oup.com/nar/article/43/D1/D512/2439467/PhosphoSitePlus-2014-mutations-PTMs-and)
4. [Horn H., et al. (2014) KinomeXplorer: an integrated platform for kinome biology studies. *Nat Methods*. 11(6):603-4](http://www.nature.com/nmeth/journal/v11/n6/full/nmeth.2968.html)
