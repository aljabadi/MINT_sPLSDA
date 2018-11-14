# Omics Data Integration with **MINT**

This vignette book explains the functionalities of MINT tool from mixOmics package.

## How to reproduce the vignette

You need the 'bookdown' package to reproduce this book. You can clone the repository, open the Rproj file and load the Rmd files and create a 'gitbook' from the 'build' pane. 

Alternatively, you can run the standalone vignettes from the *vignettes-standalone* folder.

## Prerequisites

The 'bookdown' package needs to be installed to reproduce this vignette.

```
## install only if not installed
if (!requireNamespace('bookdown', quietly = TRUE)){
  paste('Trying to install Bookdown')
  install.packages('bookdown')
}
```

mixOmics and the following Bioconductor packages must be installed prior to proceeding. You can look them up in "packages" pane to see which ones are already installed. A full installation of TeX along with 'knitr' package are also needed to build the vignette in pdf.
```
if (!requireNamespace('BiocManager', quietly = TRUE)){
paste('Trying to install BiocManager')
install.packages('BiocManager')
}

## package installations (might take a bit of time)
BiocManager::install('mixOmics', update = F)
BiocManager::install('SingleCellExperiment', update = F) ## single-cell experiment data analysis
BiocManager::install('scran', update = F) ## sc-RNAseq data analysis
BiocManager::install('scater', update = F) ## sc gene expression analysis
BiocManager::install('vennDiagram', update = F) ## Venn diagrams
BiocManager::install('tibble', update = F) ## for data tables
```
Please keep us updated of your experience on mixOmics [issues page](https://bitbucket.org/klecao/package-mixomics/issues).

## Authors

* Al J Abadi (ajalalabadi@unimelb.edu.au)
* Dr Kim-Anh Le Cao

## Acknowledgments

* Matt Ritchie and Luyi Tian, our collaborators from WEHI who provided us with the data
* CZI who funded the project

```
sessionInfo()
```

R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] grid      parallel  stats4    stats     graphics  grDevices utils
[8] datasets  methods   base

other attached packages:
[1] tibble_1.4.2                VennDiagram_1.6.20
[3] futile.logger_1.4.3         scater_1.9.24
[5] scran_1.9.39                mixOmics_6.4.6
[7] ggplot2_3.1.0               lattice_0.20-38
[9] MASS_7.3-51.1               SingleCellExperiment_1.3.12
[11] SummarizedExperiment_1.11.6 DelayedArray_0.7.49
[13] BiocParallel_1.15.15        matrixStats_0.54.0
[15] Biobase_2.41.2              GenomicRanges_1.33.14
[17] GenomeInfoDb_1.17.4         IRanges_2.15.19
[19] S4Vectors_0.19.24           BiocGenerics_0.27.1
[21] knitr_1.20

loaded via a namespace (and not attached):
[1] viridis_0.5.1             dynamicTreeCut_1.63-1
[3] edgeR_3.23.7              tidyr_0.8.2
[5] viridisLite_0.3.0         DelayedMatrixStats_1.3.11
[7] ellipse_0.4.1             assertthat_0.2.0
[9] statmod_1.4.30            highr_0.7
[11] vipor_0.4.5               GenomeInfoDbData_1.2.0
[13] yaml_2.2.0                pillar_1.3.0
[15] backports_1.1.2           glue_1.3.0
[17] limma_3.37.10             digest_0.6.18
[19] RColorBrewer_1.1-2        XVector_0.21.4
[21] colorspace_1.3-2          htmltools_0.3.6
[23] Matrix_1.2-15             plyr_1.8.4
[25] pkgconfig_2.0.2           bookdown_0.7
[27] zlibbioc_1.27.0           purrr_0.2.5
[29] corpcor_1.6.9             scales_1.0.0
[31] HDF5Array_1.9.19          RSpectra_0.13-1
[33] withr_2.1.2               lazyeval_0.2.1
[35] magrittr_1.5              crayon_1.3.4
[37] evaluate_0.12             beeswarm_0.2.3
[39] tools_3.5.0               formatR_1.5
[41] stringr_1.3.1             Rhdf5lib_1.3.3
[43] munsell_0.5.0             locfit_1.5-9.1
[45] lambda.r_1.2.3            bindrcpp_0.2.2
[47] compiler_3.5.0            rlang_0.3.0.1
[49] rhdf5_2.25.11             RCurl_1.95-4.11
[51] BiocNeighbors_0.99.22     rstudioapi_0.8
[53] igraph_1.2.2              labeling_0.3
[55] bitops_1.0-6              rmarkdown_1.10
[57] gtable_0.2.0              codetools_0.2-15
[59] rARPACK_0.11-0            reshape2_1.4.3
[61] R6_2.3.0                  gridExtra_2.3
[63] dplyr_0.7.8               bindr_0.1.1
[65] rprojroot_1.3-2           futile.options_1.0.1
[67] ggbeeswarm_0.6.0          stringi_1.2.4
[69] Rcpp_1.0.0                tidyselect_0.2.5
[71] xfun_0.4
