library(testthat)
library(mixOmics) # make sure you upload the latest version of mixOmics on gitHub: library('devtools'); install_github("mixOmicsTeam/mixOmics")
library(SingleCellExperiment)

## load a combined sce with batch data 
sce = readRDS("~/Documents/Projects/MINT_Wrapper/RNAmix_combined_scran_norm.Rds")

# ## subset data with Highly Variable Genes only
# sce.hvg = sce[rowData(sce)$hi_var,]

################ Basic run
system.time({ # 3 sec - RNA mix all genes
sce.mint = mixOmics_mint(sce)
})

## tests
expect_is(object = sce.mint, class = "SingleCellExperiment")
expect_output(str(sce.mint), "mint_variates")
expect_output(str(sce.mint), "mint.markers")




################ With tuning
system.time({ # 110 sec - RNA mix all genes
sce.mint.tuned = mixOmics_mint(sce, tune.hps = T)
})

expect_is(object = sce.mint.tuned, class = "SingleCellExperiment")
expect_output(str(sce.mint.tuned), "mint_variates")
expect_output(str(sce.mint.tuned), "mint.markers")

expect(isTRUE(sum(rowData(sce.mint.tuned)$mint.markers)>0) , failure_message = "No signatures were chosen")

################ Tuning the optimal ncomp 
keepX = c(50,40,30)
system.time({ # 16 sec - RNA mix all genes
sce.mint.opt.comp = mixOmics_mint(sce, optimum.ncomp = T, ncomp = 3L, keepX = keepX)
})
expect_is(object = sce.mint.opt.comp, class = "SingleCellExperiment")
expect_output(str(sce.mint.opt.comp), "mint_variates")
expect_output(str(sce.mint.opt.comp), "mint.markers")


################ component plots

mint.ggplot = function(sce.mint=sce.mint, comps = c(1,2), colData.class = 'mix'){
  variates =as.data.frame(reducedDim(sce.mint))
  variates$Batch = colData(sce.mint)$batch
  variates$Class = factor(colData(sce.mint)[[colData.class]])
  
  p = ggplot(variates) + geom_point(aes(x=variates[,comps[1]], y=variates[,comps[2]], col =Batch, shape = Class), size=2) +
    labs(x = paste0("Variate ", comps[1]), y=paste0("Variate ", comps[2])) +
    scale_shape_manual(values = c(1:length(unique(variates$Class)))) + theme_bw()
  
  return(p)
}

## plot the basic model
mint.ggplot(sce.mint, comps=c(1,2), colData.class = 'mix')
  ## other components
  mint.ggplot(sce.mint, comps=c(2,3), colData.class = 'mix')
  
## plot the tuned model
mint.ggplot(sce.mint.tuned, comps=c(1,2), colData.class = 'mix')


################ entry checks
## expect error when length(keepX) < ncomp
mixOmics_mint(sce, ncomp = 4, keepX=c(30,10))

## expect error when ncomp =1 and opt = T
mixOmics_mint(sce, ncomp = 1, optimum.ncomp = T, keepX=c(30,10))
