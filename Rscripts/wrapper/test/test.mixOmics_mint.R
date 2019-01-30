### tests the mint.splsda wrapper
### works only with the given RNAmix dataset

library(testthat)
library(mixOmics) # make sure you upload the latest version of mixOmics on gitHub: library('devtools'); install_github("mixOmicsTeam/mixOmics")
library(SingleCellExperiment)

## load a combined sce with batch data 
sce = readRDS("~/Documents/Projects/MINT_Wrapper/RNAmix_combined_scran_norm.Rds")

source("../Rscripts/wrapper/mixOmics_mint.R")

# ## subset data with Highly Variable Genes only
# sce.hvg = sce[rowData(sce)$hi_var,]

#################################################################
########################## Test the Expected Success ############
#################################################################

## a function for usual sce sanity checks
sce_sanity = function(sce){
  test_that("outputs SCE object",inherits(sce, "SingleCellExperiment"))
  test_that("MINT markers are in rowData", sum(rowData(sce)$mint_marker)>0)
  test_that("MINT components are in reducedDim", !is.null(reducedDims(sce)$mint_comps_global))
}
                       
##########  run with default parameters
test_default = mixOmics_mint(sce)
sce_sanity(test_default)

##########  tuned mint.splsda
test_tune = mixOmics_mint(sce, tune.keepX = seq(10,20,10))
sce_sanity(test_tune)

########## when hvgs are given as a:

## logical rowData
test_hvgs_rowData = mixOmics_mint(sce, hvgs = "hi_var")
sce_sanity(test_hvgs_rowData)
test_that("all genes in output are hi_var",
          all(rowData(sce[rownames(test_hvgs_rowData),])$hi_var))

## string
test_hvgs_string = mixOmics_mint(sce, hvgs = rownames(sce)[rowData(sce)$hi_var])
sce_sanity(test_hvgs_string)
test_that("all genes in output are hi_var",
          all(rowData(sce[rownames(test_hvgs_string),])$hi_var))

##########  colData.batch

## it accommodates any name
colData(sce)$new_batch = colData(sce)$batch
test_batch = mixOmics_mint(sce, colData.batch = "new_batch")
sce_sanity(test_batch)

##########  colData.class

## it accommodates any name
colData(sce)$new_class = colData(sce)$mix
test_class = mixOmics_mint(sce, colData.class = "new_class")
sce_sanity(test_class)


#################################################################
########################## Test the Expected Error ############
#################################################################

## since it is not desired to output the errors in a list
## errors cannot be checked routinely using testthat

######### invalid sce creates error
mixOmics_mint(sce = "non-sce")
mixOmics_mint(sce = something)

############ batch test
## colData.batch must be a string or number pertaining to one of colData(sce)
mixOmics_mint(sce = sce, colData.batch = sth)
## colData.class does not correspond to a valid  colData
mixOmics_mint(sce, colData.batch = "nonbatch")
## there must be more than one batch in the data to perform mint.splsda
mixOmics_mint(sce[,sce %>% colData() %>% rownames() %>% grepl("CEL",.)], colData.class = "mix")

############ Y test
## colData.class must be a string or number pertaining to one of colData(sce)
mint_tester(sce = sce, colData.class = sth)
## colData.class does not correspond to a valid  colData
mixOmics_mint(sce, colData.class = "nonclass")

############  keepX and ncomp
## the length of keepX should be at least ncomp
mixOmics_mint(sce, keepX = c(2,3), ncomp = 3)
## warning if length(keepX) > ncomp
mixOmics_mint(sce, keepX = c(2,3), ncomp = 1)
