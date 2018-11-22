# Single Cell Data Analysis and Integration using MINT

## Rproject

Contains the R project that produced the vignette available at [this link](https://ajabadi.github.io/MINT_sPLSDA/).
*sessioninfo.md* contains the R session information when the current version of the vignette was built.

## Rscripts

Contains the R files extracted from the Rmd files.

## vignettes-standlone

Self-contained standalone Rmd files that can be run separately.

## output

Where the run data can be saved when running the R project. See *params* section of render_book.R / Rmd files.

## data

Where the input data can be saved. Alternatively, all required data for this vignette can be loaded from GitHub (default). See *params* section of render_book.R / Rmd files.

#### data/subset

A subset of data can be stored in this location for a faster initial run. See *params* section of render_book.R.

## docs

Contains the bookdown outputs of the Rproject file to create the [this vignette](https://ajabadi.github.io/MINT_sPLSDA/)
