--- 
title: "ComplexHeatmap Complete Reference"
author: "Zuguang Gu"
date: "last revised on 2021-08-02"
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
are generated under version 2.9.3.

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
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] dendsort_0.3.4       dendextend_1.15.1    circlize_0.4.13     
## [4] ComplexHeatmap_2.9.3 BiocManager_1.30.16  colorout_1.2-2      
## 
## loaded via a namespace (and not attached):
##  [1] shape_1.4.6         bslib_0.2.5.1       GetoptLong_1.0.5   
##  [4] tidyselect_1.1.1    xfun_0.24           purrr_0.3.4        
##  [7] colorspace_2.0-2    vctrs_0.3.8         generics_0.1.0     
## [10] htmltools_0.5.1.1   stats4_4.1.0        viridisLite_0.4.0  
## [13] yaml_2.2.1          utf8_1.2.1          rlang_0.4.11       
## [16] jquerylib_0.1.4     pillar_1.6.1        glue_1.4.2         
## [19] DBI_1.1.1           BiocGenerics_0.38.0 RColorBrewer_1.1-2 
## [22] matrixStats_0.59.0  foreach_1.5.1       lifecycle_1.0.0    
## [25] stringr_1.4.0       munsell_0.5.0       gtable_0.3.0       
## [28] GlobalOptions_0.1.2 codetools_0.2-18    evaluate_0.14      
## [31] knitr_1.33          IRanges_2.26.0      doParallel_1.0.16  
## [34] parallel_4.1.0      fansi_0.5.0         scales_1.1.1       
## [37] S4Vectors_0.30.0    jsonlite_1.7.2      gridExtra_2.3      
## [40] rjson_0.2.20        ggplot2_3.3.5       png_0.1-7          
## [43] digest_0.6.27       stringi_1.6.2       bookdown_0.22      
## [46] dplyr_1.0.7         clue_0.3-59         tools_4.1.0        
## [49] sass_0.4.0          magrittr_2.0.1      tibble_3.1.2       
## [52] cluster_2.1.2       crayon_1.4.1        pkgconfig_2.0.3    
## [55] ellipsis_0.3.2      assertthat_0.2.1    rmarkdown_2.9      
## [58] iterators_1.0.13    viridis_0.6.1       R6_2.5.0           
## [61] compiler_4.1.0
```

