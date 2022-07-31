--- 
title: "ComplexHeatmap Complete Reference"
author: "Zuguang Gu"
date: "last revised on 2022-07-31"
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
are generated under version 2.13.1.

You can get a stable Bioconductor version from
http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html, but
the most up-to-date version is always on Github and you can install it by:


```r
library(devtools)
install_github("jokergoo/ComplexHeatmap")
```

The [development branch on Bioconductor](http://bioconductor.org/packages/devel/bioc/html/ComplexHeatmap.html) is
basically synchronized to the Github repository. 

The **ComplexHeatmap** package is inspired by the [**pheatmap**](https://CRAN.R-project.org/package=pheatmap) package. You can find many arguments in **ComplexHeatmap** have the same names as in **pheatmap**. Also you
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
## R version 4.2.0 (2022-04-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur/Monterey 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] dendsort_0.3.4        dendextend_1.15.2     circlize_0.4.16      
## [4] ComplexHeatmap_2.13.1 colorout_1.2-2       
## 
## loaded via a namespace (and not attached):
##  [1] png_0.1-7           assertthat_0.2.1    digest_0.6.29      
##  [4] foreach_1.5.2       utf8_1.2.2          R6_2.5.1           
##  [7] stats4_4.2.0        evaluate_0.15       ggplot2_3.3.6      
## [10] pillar_1.7.0        GlobalOptions_0.1.2 rlang_1.0.2        
## [13] jquerylib_0.1.4     S4Vectors_0.34.0    GetoptLong_1.0.5   
## [16] rmarkdown_2.14      stringr_1.4.0       munsell_0.5.0      
## [19] compiler_4.2.0      xfun_0.31           pkgconfig_2.0.3    
## [22] BiocGenerics_0.42.0 shape_1.4.6         htmltools_0.5.2    
## [25] downlit_0.4.0       tidyselect_1.1.2    tibble_3.1.7       
## [28] gridExtra_2.3       bookdown_0.27       IRanges_2.30.0     
## [31] codetools_0.2-18    matrixStats_0.62.0  fansi_1.0.3        
## [34] viridisLite_0.4.0   crayon_1.5.1        dplyr_1.0.9        
## [37] jsonlite_1.8.0      gtable_0.3.0        lifecycle_1.0.1    
## [40] DBI_1.1.2           magrittr_2.0.3      scales_1.2.0       
## [43] stringi_1.7.6       cli_3.3.0           cachem_1.0.6       
## [46] viridis_0.6.2       fs_1.5.2            doParallel_1.0.17  
## [49] xml2_1.3.3          bslib_0.3.1         ellipsis_0.3.2     
## [52] generics_0.1.2      vctrs_0.4.1         rjson_0.2.21       
## [55] RColorBrewer_1.1-3  iterators_1.0.14    tools_4.2.0        
## [58] glue_1.6.2          purrr_0.3.4         parallel_4.2.0     
## [61] fastmap_1.1.0       yaml_2.3.5          clue_0.3-61        
## [64] colorspace_2.0-3    cluster_2.1.3       memoise_2.0.1      
## [67] knitr_1.39          sass_0.4.1
```

