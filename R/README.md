# User-directed unsupervised identification of cell populations
	
## DAFi: User-directed unsupervised filtering and identification of cell populations from flow cytometry data

### Author: 
Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu_cs@yahoo.com, Ivan Chang, Ph.D., ichang@jcvi.org, and Bob Sinkovits, Ph.D., sinkovit@sdsc.edu
	
### Copyright: Author and J. Craig Venter Institute
	
Paper link: https://www.biorxiv.org/content/early/2017/09/26/193912

#### install into R
For the DAFi R implementation framework, R version > 3.4 is required (https://cran.r-project.org/bin/). In addition, please have installed:
1) flowCore (https://www.bioconductor.org/packages/release/bioc/html/flowCore.html)
2) flowViz (http://bioconductor.org/packages/release/bioc/html/flowViz.html)
3) ClusterR (https://cran.r-project.org/web/packages/ClusterR/index.html)
4) FlowSOM (https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)


For automated build and install, including dependent packages listed above, please install the devtools library
```r
install.packages("devtools")
```

so you can initiate the automated install of DAFi package

```r
devtools::install_github("JCVenterInstitute/DAFi-gating", build_vignettes = TRUE)
```

then checkout the built-in vignette for the DAFi library for more documentation

```r
library(DAFi)
browseVignettes(DAFi)
```
