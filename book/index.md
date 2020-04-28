
--- 
title: "ComplexHeatmap Complete Reference"
author: "Zuguang Gu"
date: "last revised on 2020-04-28"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: jokergoo/ComplexHeatmap-reference
cover-image: complexheatmap-cover.jpg
url: 'https\://jokergoo.github.io/ComplexHeatmap-reference/book'
description: "Complex heatmaps are efficient to visualize associations between different sources of data sets and reveal potential patterns. Here the ComplexHeatmap R package provides a highly flexible way to arrange multiple heatmaps and supports various annotation graphics. This book is the complete reference to ComplexHeatmap pacakge."
---

# About {-}

This is the documentation of the
[**ComplexHeatmap**](https://github.com/jokergoo/ComplexHeatmap) package. Examples in the book
are generated under version 2.3.5.

You can get a stable Bioconductor version from http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html, but the most up-to-date version is always on Github and you can install it by:


```r
library(devtools)
install_github("jokergoo/ComplexHeatmap")
```

The [development branch on Bioconductor](http://bioconductor.org/packages/devel/bioc/html/ComplexHeatmap.html) is
basically synchronized to Github repository. 

The **ComplexHeatmap** package is inspired from the [**pheatmap**](https://CRAN.R-project.org/package=pheatmap) package. You can find many arguments in **ComplexHeatmap** have the same names as in **pheatmap**. Also you
can find [this old package](https://github.com/jokergoo/pheatmap2) that I tried to develop by modifying **pheatmap**.

**Please note, this documentation is not completely compatible with older versions (< 1.99.0, before
Oct, 2018), but the major functionality keeps the same.**

If you use **ComplexHeatmap** in your publications, I am appreciated if you can cite:

Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data.
DOI: [10.1093/bioinformatics/btw313](https://doi.org/10.1093/bioinformatics/btw313)


<img src="complexheatmap-cover.jpg" style="width:500px;border:2px solid black;" />

Session info:


```r
sessionInfo()
```

```
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS:   /usr/lib64/libblas.so.3.4.2
## LAPACK: /usr/lib64/liblapack.so.3.4.2
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] dendsort_0.3.3       dendextend_1.13.2    circlize_0.4.9      
## [4] ComplexHeatmap_2.3.5
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4          pillar_1.4.3        compiler_3.6.0     
##  [4] RColorBrewer_1.1-2  viridis_0.5.1       tools_3.6.0        
##  [7] digest_0.6.25       viridisLite_0.3.0   evaluate_0.14      
## [10] lifecycle_0.2.0     tibble_3.0.1        gtable_0.3.0       
## [13] clue_0.3-57         pkgconfig_2.0.3     png_0.1-7          
## [16] rlang_0.4.5         yaml_2.2.1          parallel_3.6.0     
## [19] xfun_0.13           gridExtra_2.3       dplyr_0.8.3        
## [22] stringr_1.4.0       cluster_2.1.0       knitr_1.28         
## [25] vctrs_0.2.4         GlobalOptions_0.1.1 tidyselect_1.0.0   
## [28] glue_1.4.0          R6_2.4.1            GetoptLong_0.1.8   
## [31] rmarkdown_1.18      bookdown_0.17       purrr_0.3.4        
## [34] ggplot2_3.2.1       magrittr_1.5        scales_1.1.0       
## [37] htmltools_0.4.0     ellipsis_0.3.0      assertthat_0.2.1   
## [40] shape_1.4.4         colorspace_1.4-1    stringi_1.4.6      
## [43] lazyeval_0.2.2      munsell_0.5.0       crayon_1.3.4       
## [46] rjson_0.2.20
```
