params <-
list(local.input = FALSE, subset.data = FALSE, output.data = FALSE, 
    Rscripts = FALSE, recalc = TRUE)

## installing the required packages for this vignette if necessary
required.pkgs = c('mixOmics',
                  'SingleCellExperiment', ## single-cell experiment data analysis
                  'scran', ## sc-RNAseq data analysis
                  'VennDiagram', ## Venn diagrams
                  'tibble') ## for data tables

## installing BiocManager to install packages
## it can install CRAN packages as well
if (!requireNamespace('BiocManager', quietly = T)){
  paste('Trying to install BiocManager')
  install.packages('BiocManager')
}
## package installer function - only if it is not installed.
package.installer = function(pkgs=required.pkgs){
  for (package in pkgs){
    if (!requireNamespace(package, quietly = T)){
  paste0('Trying to install ', package)
  BiocManager::install(package, update = F)
    }
    }
}
## run function
package.installer(required.pkgs)
## load the required libraries
library(SingleCellExperiment)
library(mixOmics)
library(scran)
library(knitr)
library(VennDiagram)
library(tibble)
check.exists <- function(object) ## function to assess existence of objects
{
  exists(as.character(substitute(object)))
}

## input/output from parameters
io = list()

## whether or where from to locally load data - FALSE: GitHub load; or a directory
io$local.input = ifelse(check.exists(params$local.input), params$local.input, F)

## whether or where to save run data - FALSE: do not save; or a directory
io$output.data = ifelse(check.exists(params$output.data), params$output.data, F)

## whether or where to save R scripts - FALSE: do not save; or a directory
io$Rscripts=ifelse(check.exists(params$Rscripts), params$Rscripts, F)
if (isFALSE(io$local.input)){
  ## load from GitHub
  DataURL='https://tinyurl.com/sincell-with-class-RData-LuyiT'
  load(url(DataURL))
} else {
  load(file.path(io$local.input, 'sincell_with_class.RData'))
}
## make a summary of QC'ed cell line data processed by each protocol
sce10xqc_smr =  summary(as.factor(sce_sc_10x_qc$cell_line))
sce4qc_smr =    summary(as.factor(sce_sc_CELseq2_qc$cell_line))
scedropqc_smr = summary(as.factor(sce_sc_Dropseq_qc$cell_line))
## combine the summaries
celline_smr = rbind(sce10xqc_smr,sce4qc_smr,scedropqc_smr)
## produce a 'total' row as well
celline_smr = cbind(celline_smr, apply(celline_smr,1,sum))
## add the genes as well
celline_smr = cbind(celline_smr,
                    c(dim(counts(sce_sc_10x_qc))[1],
                      dim(counts(sce_sc_CELseq2_qc))[1],
                      dim(counts(sce_sc_Dropseq_qc))[1]))
## label the rows
row.names(celline_smr) = c('10X', 'CEL-seq2', 'Drop-seq')
colnames(celline_smr) = c('H1975', 'H2228', 'HCC827', 
                          'Total Cells','Total Genes')
## tabulate the summaries
kable(celline_smr,
      caption = 'Summary of cell and gene data 
      per cell line for each protocol')
## show the cell line summaries
celline_smr
## create venn diagram of genes in each protocol:
venn.plot = venn.diagram(
  x = list(Chrom.10X = rownames(sce_sc_10x_qc),
           CEL.seq2 = rownames(sce_sc_CELseq2_qc),
           Drop.seq = rownames(sce_sc_Dropseq_qc)),
  filename = NULL, label=T, margin=0.05,
  height = 1400, width = 2200,
  col = 'transparent', fill = c('cornflowerblue','green', 'red'),
  alpha = 0.60, cex = 2, fontfamily = 'serif', fontface = 'bold',
  cat.col = c('darkblue', 'darkgreen', 'black'), cat.cex = 2.2)

## save the plot, change to your own directory
png(filename = 'figures/GeneVenn.png')
grid.draw(venn.plot)
dev.off()
## normalise the QC'ed count matrices
sc10x.norm =  computeSumFactors(sce_sc_10x_qc) ## deconvolute using size factors
sc10x.norm =  normalize(sc10x.norm) ## normalise expression values
## DROP-seq
scdrop.norm = computeSumFactors(sce_sc_Dropseq_qc)
scdrop.norm = normalize(scdrop.norm)
## CEL-seq2
sccel.norm =  computeSumFactors(sce_sc_CELseq2_qc)
sccel.norm =  normalize(sccel.norm)
## pca on the normalised count matrices and find 10 PCs
pca.res.10x =     pca(t(logcounts(sc10x.norm)),  ncomp = 10,
                      center=T, scale=F)
pca.res.celseq =  pca(t(logcounts(sccel.norm)),  ncomp = 10,
                      center=T, scale=F)
pca.res.dropseq = pca(t(logcounts(scdrop.norm)), ncomp = 10,
                      center=T, scale=F)
## details of the pca output
pca.res.10x

## arrange the plots in 1 row and 3 columns
par(mfrow=c(1,3))
## find the maximum explained variance in all PCs:
ymax = max (pca.res.10x$explained_variance[1],
             pca.res.celseq$explained_variance[1],
             pca.res.dropseq$explained_variance[1])

## plot the pca objects and limit the Y axes to ymax for all
plot(pca.res.10x, main= '(A) 10X', ylim=c(0,ymax))
plot(pca.res.celseq,  main= '(B) CEL-seq2', ylim=c(0,ymax))
plot(pca.res.dropseq,  main= '(C) Drop-seq', ylim=c(0,ymax))
## define custom colours/shapes
## colour by cell line
col.cell = c('H1975'='#0000ff', 'HCC827'='grey30', 'H2228' ='#ff8000')
## shape by batch
shape.batch = c('10X' = 1, 'CEL-seq2'=2, 'Drop-seq'=3 )
## pca plots for protocols
## 10x
plotIndiv(pca.res.10x, legend = T, title = 'PCA 10X', pch = shape.batch['10X'], col = col.cell,
          group = sce_sc_10x_qc$cell_line, legend.title = 'Cell Line')
## CEL-seq2
plotIndiv(pca.res.celseq, legend = T, title = 'PCA CEL-seq2', pch = shape.batch['CEL-seq2'],
          col = col.cell, group = sce_sc_CELseq2_qc$cell_line, legend.title = 'Cell Line')
## Drop-seq
plotIndiv(pca.res.dropseq, legend = T, title = 'PCA Drop-seq', pch = shape.batch['Drop-seq'],
          col = col.cell, group = sce_sc_Dropseq_qc$cell_line, legend.title = 'Cell Line')
## find the intersect of the genes for integration
list.intersect = Reduce(intersect, list(
## the rownames of the original (un-transposed) count matrix will -
## output the genes
  rownames(logcounts(sc10x.norm)),
  rownames(logcounts(sccel.norm)),
  rownames(logcounts(scdrop.norm))
))
## combine the data at their intersection
data.combined = t( ## transpose of all 3 datasets combined
  data.frame(
    ## the genes from each protocol that match list.intersect
    logcounts(sc10x.norm)[list.intersect,],
    logcounts(sccel.norm)[list.intersect,],
    logcounts(scdrop.norm)[list.intersect,] ))
## the number of cells and genes in the intersect dataset
dim(data.combined)
## create a factor variable of cell lines
## must be in the same order as the data combination
cell.line = as.factor(c(sce_sc_10x_qc$cell_line,
                         sce_sc_CELseq2_qc$cell_line,
                         sce_sc_Dropseq_qc$cell_line))
## name the factor variable with the cell ID
names(cell.line) = rownames(data.combined)

## produce a character vector of batch names
## must be in the same order as the data combination
batch = as.factor(
  c(rep('10X',      ncol(logcounts(sc10x.norm))),
    rep('CEL-seq2',  ncol(logcounts(sccel.norm))),
    rep('Drop-seq', ncol(logcounts(scdrop.norm))) ))
## name it with corresponding cell IDs
names(batch) = rownames(data.combined)
## perform PCA on concatenated data and retrieve 2 PCs
pca.combined = pca(data.combined, ncomp = 2)
## plot the combined pca coloured by batches
plotIndiv(pca.combined, title = 'PCA Combined',
          pch = batch, ## shape by cell line
          group = cell.line, ## colour by batch
          legend = T, legend.title = 'Cell Line',
          legend.title.pch = 'Study')
## plot the combined pca coloured by protocols
plotIndiv(pca.combined, title = 'PCA Combined',
          pch = cell.line, ## shape by cell line
          group = batch, ## colour by protocol
          col.per.group = c('red', 'purple', 'green'),
          legend = T, legend.title = 'Study',
          legend.title.pch = 'Cell Line')
## create variables needed for MINT
## factor variable of cell lines
Y = as.factor(cell.line[rownames(data.combined)])
## factor variable of studies
study = batch ## defined in the combined PCA section
## MINT on the combined dataset with 5 components
mint.plsda.res = mint.plsda(X = data.combined, Y = Y,
                             study = study, ncomp = 5)
## plot the mint.plsda plots for the combined dataset
plotIndiv(mint.plsda.res, group = cell.line,
          legend  = T, subtitle     = 'MINT - Coloured by Cell Line',
          ellipse = F, legend.title = 'Cell Line',
          legend.title.pch = 'protocol',
          X.label = 'PLS-DA component 1',
          Y.label = 'PLS-DA component 2')
## perform cross validation and calculate classification error rates
set.seed(12321)  # for reproducibility of the results
perf.mint.plsda.res = perf(mint.plsda.res,
          progressBar = F)
## plot the classification error rate vs number of components
plot(perf.mint.plsda.res, col = color.mixo(5:7))
perf.mint.plsda.res$global.error$BER ## further error diagnostics 
## number of variables to select in each component
perf.mint.plsda.res$choice.ncomp
## number of variables to keep in each component
list.keepX = c(50,50) 
## perform sparse PLS-DA using MINT with 2 components
mint.splsda.res = mint.splsda(X = data.combined, Y = Y,
                              study = study, ncomp = 2, keepX = list.keepX)
plotIndiv(mint.splsda.res,
          legend  = T, subtitle = 'Sparse MINT', ellipse = T,
          X.label = 'sPLS-DA component 1', 
          Y.label = 'sPLS-DA component 2',
          group = Y, ## colour by cell line
          legend.title = 'Cell Line',
          legend.title.pch = 'Protocol')
## tune MINT for component 1 and then 2 and record the run time
## we tune individual components for visualisation purposes
## one can run only the tune.mint.c2 without already.tested.X
start.time = Sys.time()
## tune using a test set of variable numbers
## component 1
tune.mint.c1 = tune(
  X = data.combined, Y = Y, study = study, ncomp = 1,
  ## assess numbers 5,10,15...50,60,70,...100:
  test.keepX = c(seq(5,35,5),seq(40,70,10), 100), method = 'mint.splsda',
  ## use all distances to estimate the classification error rate
  dist = c('max.dist',  'centroids.dist', 'mahalanobis.dist'),
  progressBar = F
)
## component 1 to 2
tune.mint.c2 = tune(
  X = data.combined, Y = Y, study = study, ncomp = 2,
  ## already tuned component 1
  already.tested.X = tune.mint.c1$choice.keepX,
  test.keepX = c(seq(5,35,5),seq(40,70,10), 100), method = 'mint.splsda',
  dist = c('max.dist',  'centroids.dist', 'mahalanobis.dist'),
  progressBar = F
)
end.time = Sys.time()
## see how long it takes to find the optimum number of variables:
run.time =  end.time - start.time
## look at the optimal selected variables for each PC
tune.mint.c2$choice.keepX
## plot the error rates for all test variable numbers
par(mfrow=c(1,3))
plot(tune.mint.c1, col = 'darkred')
plot(tune.mint.c2, col = 'darkblue')
## run sparse mint using optimum parameters:
mint.splsda.tuned.res = mint.splsda( X =data.combined, Y = Y,
                              study = study, ncomp = 2,  
                              keepX = tune.mint.c2$choice.keepX)
## plot the tuned mint.splsda plot for the combined dataset
plotIndiv(mint.splsda.tuned.res, study = 'global', legend = T,
          title = 'MINT sPLS-DA',  subtitle = 'Global', ellipse=T, legend.title = 'Cell Line')
## tuned mint.splsda plot for each protocol
plotIndiv(mint.splsda.tuned.res, study = 'all.partial',  title = 'MINT sPLS-DA', 
          subtitle = c('10X', 'CEL-seq2', 'Drop-seq'))
set.seed(12321)  # for reproducibility of the results
## perform classification with leave-one-group-out cross validation 
perf.mint.final = perf(mint.splsda.res, progressBar = F, dist = 'max.dist')
## classification error rate
perf.mint.final$global.error
## ROC curves for both components
auc.mint.splsda1 = auroc(mint.splsda.tuned.res, roc.comp = 1, roc.study='CEL-seq2')
auc.mint.splsda2 = auroc(mint.splsda.tuned.res, roc.comp = 2, roc.study='CEL-seq2')
## if specified, save the run data
if (!isFALSE(io$output.data)){
  ## save normlaised sce objects
  sce.norm = c('sc10x.norm', 'sccel.norm', 'scdrop.norm')
  save(list=sce.norm,  file = file.path(io$output.data,'sce.norm.RData'))
  ## save the final mint object
  save(mint.splsda.tuned.res, file =file.path(io$output.data,'mint.splsda.tuned.res.RData'))
  ## save the rest of the data separately
  session01 = ls()
  session01 = session01[!session01 %in% c(sce.norm, 'mint.splsda.tuned.res')]
  save(session01, file =file.path(io$output.data,'session.integration.RData'))
}
writeLines(capture.output(sessionInfo()), "sessionInfo.md")
## session information to build this vignette
sessionInfo()
