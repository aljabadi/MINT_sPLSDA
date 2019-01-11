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
## check for valid input/output setup
check.exists <- function(object) ## function to assess existence of objects
{
  exists(as.character(substitute(object)))
}

## input/output from parameters
io = list()

## whether or where from to locally load data - FALSE: GitHub load; or a directory
io$local.input = ifelse(check.exists(params$local$input), params$local.input, F)

## whether or where to save run data - FALSE: do not save; or a directory
io$output.data = ifelse(check.exists(params$output.data), params$output.data, F)

## whether or where to save R scripts - FALSE: do not save; or a directory
io$Rscripts=ifelse(check.exists(params$Rscripts), params$Rscripts, F)

## whether the normalised sce and mint objects from output.data should be re-calculated (TRUE) or directly loaded (FALSE)
io$recalc = ifelse(check.exists(params$recalc), params$recalc, F)
if (isFALSE(io$local.input)&&io$recalc){
  ## load from GitHub
  DataURL='https://tinyurl.com/sincell-with-class-RData-LuyiT'
  load(url(DataURL))
} else if (!isFALSE(io$local.input)&&io$recalc){
  load(file.path(io$local.input, 'sincell_with_class.RData'))
  }
if (io$recalc){
  ## normalise the QC'ed count matrices
  sc10x.norm =  computeSumFactors(sce_sc_10x_qc) ## deconvolute using size factors
  sc10x.norm =  normalize(sc10x.norm) ## normalise expression values
  ## DROP-seq
  scdrop.norm = computeSumFactors(sce_sc_Dropseq_qc)
  scdrop.norm = normalize(scdrop.norm)
  ## CEL-seq2
  sccel.norm =  computeSumFactors(sce_sc_CELseq2_qc)
  sccel.norm =  normalize(sccel.norm)
} else {
  ## load locally - change to your own
  load('../output/sce.norm.RData')
}
## find the intersect of the genes
list.intersect = Reduce(intersect, list(
## the rownames of the original (un-transposed) count matrix will output the genes
  rownames(logcounts(sc10x.norm)),
  rownames(logcounts(sccel.norm)),
  rownames(logcounts(scdrop.norm))
))
## PLSDA - 10x
## extract the normalised count matrix from the SCE object (transposed)
normalised.10x = t(logcounts(sc10x.norm))
## keep the common genes only for comparability
normalised.10x = normalised.10x[,list.intersect]
## form a factor variable of cell lines
Y.10x = as.factor(sc10x.norm[list.intersect,]$cell_line)
## PLS-DA on the dataset with 5 components
plsda.10x.res = plsda(X = normalised.10x, Y = Y.10x, ncomp = 5)
## perform cross validation and find the classification error rates
start = Sys.time()
perf.plsda.10x = perf(plsda.10x.res, progressBar=F )
run.time = Sys.time()-start
## optimal number of components
plot(perf.plsda.10x, col = color.mixo(5:7))
## PLSDA - CEL-seq2 and Drop-seq 
## extract the normalised count matrix from the SCE object (transposed)
## CEL-seq2
normalised.cel = t(logcounts(sccel.norm))
normalised.cel = normalised.cel[,list.intersect]
## Drop-seq
normalised.drop= t(logcounts(scdrop.norm))
normalised.drop = normalised.drop[,list.intersect]
## factor variable of cell lines
Y.cel = as.factor(sccel.norm[list.intersect,]$cell_line)
Y.drop = as.factor(scdrop.norm[list.intersect,]$cell_line)
## PLS-DA on each dataset with 2 components
plsda.10x.res = plsda(X = normalised.10x, Y = Y.10x, ncomp = 2)
plsda.cel.res = plsda(X = normalised.cel, Y = Y.cel, ncomp = 2)
plsda.drop.res = plsda(X = normalised.drop, Y = Y.drop, ncomp = 2)
## mint.plsda plot for 10X
plotIndiv(plsda.10x.res,
          legend  = T, title     = 'PLSDA 10X', 
          ellipse = T, legend.title = 'Cell Line',
          X.label = 'PLSDA component 1', 
          Y.label = 'PLSDA component 2', pch=1)
## mint.plsda plot for CEL-seq2
plotIndiv(plsda.cel.res,
          legend  = T, title     = 'PLSDA CEL-seq2', 
          ellipse = T, legend.title = 'Cell Line',
          X.label = 'PLSDA component 1', 
          Y.label = 'PLSDA component 2', pch=2)
## run sparse PLSDA on individual studies with MINT tuned parameters
keepX = c(35,10)
splsda.10x.res = splsda( X =normalised.10x, Y = Y.10x, ncomp = 2,  
                              keepX = keepX)
splsda.cel.res = splsda( X =normalised.cel, Y = Y.cel, ncomp = 2,  
                              keepX = keepX)
splsda.drop.res = splsda( X =normalised.drop, Y = Y.drop, ncomp = 2,  
                              keepX = keepX)
## splsda plots with tuned number of variables for each sPLSDA component
## 10X
plotIndiv(splsda.10x.res, group = Y.10x,
          legend  = T, title     = 'sPLSDA - 10X',
          ellipse = F,legend.title = 'Cell Line',
          pch=1,
          X.label = 'sPLSDA component 1',
          Y.label = 'sPLSDA component 2')

## CEL-seq2
plotIndiv(splsda.cel.res, group = Y.cel,
          legend  = T, title     = 'sPLSDA - CEL-seq2',
          ellipse = F,legend.title = 'Cell Line',
          pch=2,
          X.label = 'sPLSDA component 1',
          Y.label = 'sPLSDA component 2')
## The signature genes from each sPLS-DA study
Chromium.10X.vars = unique(c(selectVar(splsda.10x.res, comp=1)$name,
                             selectVar(splsda.10x.res, comp=2)$name))
Cel.seq.vars =      unique(c(selectVar(splsda.cel.res, comp=1)$name,
                             selectVar(splsda.cel.res, comp=2)$name))
Drop.seq.vars =     unique(c(selectVar(splsda.drop.res, comp=1)$name,
                             selectVar(splsda.drop.res, comp=2)$name))
## create a venn diagram from signatures
vennProtocols <- venn.diagram(
	x = list(
		Chr.10X= Chromium.10X.vars ,
		Cel.seq= Cel.seq.vars,
		Drop.seq = Drop.seq.vars),
	filename = NULL,
	cex=1.5, cat.cex=1.5,
	fill = c('green', 'darkblue',  'yellow')
	)
png(filename = 'figures/vennProtocols.png')
grid.draw(vennProtocols)
dev.off()
if (io$recalc){
  data.combined = t( ## transpose of all 3 datasets combined
    data.frame(
      ## the genes from each protocol that match list.intersect
      logcounts(sc10x.norm)[list.intersect,],
      logcounts(sccel.norm)[list.intersect,],
      logcounts(scdrop.norm)[list.intersect,] ))
  
  ## create a factor variable of cell lines
  ## must be in the same order as the data combination
  cell.line = as.factor(c(sce_sc_10x_qc$cell_line,
                           sce_sc_CELseq2_qc$cell_line,
                           sce_sc_Dropseq_qc$cell_line))
  ## name the factor variable with the cell ID
  names(cell.line) = rownames(data.combined)
  
  ## produce a character vector of batch names
  ## must be in the same order as the data combination
  study = as.factor(
    c(rep('10X',      ncol(logcounts(sc10x.norm))),
      rep('CEL-seq2',  ncol(logcounts(sccel.norm))),
      rep('Drop-seq', ncol(logcounts(scdrop.norm))) ))
  ## name it with corresponding cell IDs
  names(study) = rownames(data.combined)
  
  ## run sparse mint using optimum parameters:
  mint.splsda.tuned.res = mint.splsda( X =data.combined, Y = Y,
                                study = study, ncomp = 2,  
                                keepX = c(35,10)) ## change for your own dataset
} else { ## if already saved
  load(file.path('../output/mint.splsda.tuned.res.RData'))
}
## Loading Plots
## 10X
plotLoadings(splsda.10x.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=F, title=NULL, 
             subtitle = '10X')
## CEL-seq2 - Comp. 1
plotLoadings(splsda.cel.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=F, title=NULL,
             subtitle = 'CEL-seq2')
## Drop-seq - Comp. 1
plotLoadings(splsda.drop.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=F, title=NULL,
             subtitle = 'Drop-seq')
## MINT - Comp. 1
plotLoadings(mint.splsda.tuned.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=F, title=NULL, 
             subtitle = c('10X', 'CEL-seq2', 'Drop-seq') )
## MINT - Comp. 2
plotLoadings(mint.splsda.tuned.res, contrib='max', method = 'mean', comp=2, 
             study='all.partial', legend=F, title=NULL, 
             subtitle = c('10X', 'CEL-seq2', 'Drop-seq') )
## MINT signature
MINT.Combined.vars = unique(c(selectVar(mint.splsda.tuned.res, comp=1)$name,
                             selectVar(mint.splsda.tuned.res, comp=2)$name))
## create venn diagram
vennMINT <- venn.diagram(
	x = list(
		Chr.10X= Chromium.10X.vars ,
		Cel.seq= Cel.seq.vars,
		Drop.seq = Drop.seq.vars,
		MINT.Combined=MINT.Combined.vars),
	filename = NULL,
	cex=1.5, cat.cex=1.5,
	fill = c('green', 'darkblue',  'yellow', 'red')
	)
png(filename = 'figures/vennMINT.png')
grid.draw(vennMINT)
dev.off()
## the signature genes identified by all studies
common.sig = Reduce(intersect, list(MINT.Combined.vars, Chromium.10X.vars, Cel.seq.vars, Drop.seq.vars))
## the 10 signature genes with highest loadings
common.sig[1:10]
## correlation circle plot
plotVar(mint.splsda.tuned.res, cex = 3)
## show genes on extreme sides

## component 1 - most positively and negatively expressed between cell lines
var.c1 = selectVar(mint.splsda.tuned.res, comp=1)$value
positive.gene.c1 = rownames(var.c1)[which.max(var.c1$value.var)]
negative.gene.c1 = rownames(var.c1)[which.min(var.c1$value.var)]

## component 2 - most positively and negatively expressed between cell lines
var.c2 = selectVar(mint.splsda.tuned.res, comp=2)$value
positive.gene.c2 = rownames(var.c2)[which.max(var.c2$value.var)]
negative.gene.c2 = rownames(var.c2)[which.min(var.c2$value.var)]
## a function to create violin + box plots for this specific dataset
violinPlot = function(mint.object, gene){
  cols = c("H2228" = "orange", "H1975" = "dodgerblue3", "HCC827" = "grey")
  ggplot() + 
    geom_boxplot(aes(mint.object$Y,  mint.object$X[,gene],
                     fill= mint.object$Y), alpha=1)+
    geom_violin(aes(mint.object$Y, mint.object$X[,gene],
                     fill= mint.object$Y), alpha=0.7) +
    labs(x = "Cell Line", y="Standardised Expression Value" )+
    guides(fill=guide_legend(title="Cell Line") ) +
    scale_fill_manual(values=cols ) 
}
## violin + box plots for the most positively expressed gene on component 1
violinPlot(mint.splsda.tuned.res, positive.gene.c1)
## violin + box plots for the most negatively expressed gene on component 1
violinPlot(mint.splsda.tuned.res, negative.gene.c1)
## violin + box plots for the most positively expressed gene on component 2
violinPlot(mint.splsda.tuned.res, positive.gene.c2)
## violin + box plots for the most negatively expressed gene on component 2
violinPlot(mint.splsda.tuned.res, negative.gene.c2)
## hierarchical clustering
cim(mint.splsda.tuned.res, comp = c(1,2), margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(mint.splsda.tuned.res$Y)), row.names = F,
    title = 'MINT sPLS-DA', save='png', name.save = 'heatmap')
if (isFALSE(io$local.input)){
  ## load from GitHub
  DE_URL <- 'https://tinyurl.com/DEtable-90cells'
  load(url(DE_URL))
} else {
  load(file.path(io$local.input, 'DEtable_90cells.RData'))
  }
## keep the genes present in the sPLS-DA analysis
HCC827_DEtable = HCC827_DEtable[row.names(HCC827_DEtable) %in% list.intersect,]
H2228_DEtable = H2228_DEtable[row.names(H2228_DEtable) %in% list.intersect,]
H1975_DEtable = H1975_DEtable[row.names(H1975_DEtable) %in% list.intersect,]

## create a column of genes for ease of merging
HCC827_DE = rownames_to_column(HCC827_DEtable, 'gene')
H2228_DE = rownames_to_column(H2228_DEtable, 'gene')
H1975_DE = rownames_to_column(H1975_DEtable, 'gene')

## create a combined DEtable from the 3 cell lines
DEtable = rbind(HCC827_DE,H2228_DE, H1975_DE)
## order by gene name and FDR (increasing)
DEtable = DEtable[order(DEtable[,'gene'],DEtable[,'FDR']),]
## keep the ones with FDR<0.05
DEtable = DEtable[DEtable$FDR<0.05,]
## remove duplicate gene names
DEtable = DEtable[!duplicated(DEtable$gene),]
## overlap with MINT
DE.MINT = DEtable[(DEtable$gene %in% MINT.Combined.vars),]
## sort in increasing FDR order
DE.MINT = DE.MINT[order(DE.MINT$FDR),]
## number of MINT signature genes that are differentially expressed
dim(DE.MINT)[1]
## geometric mean of the FDR of signature
exp(mean(log(DE.MINT$FDR)))
## if specified, save the run data
if (!isFALSE(io$output.data)){
  ## set to NULL if running this vignette alone:
  already.saved = c('sc10x.norm', 'sccel.norm', 'scdrop.norm', 'mint.splsda.tuned.res')
  ## save the files not saved already
  session02 = ls()
  session02 = session02[!session02 %in% already.saved]
  save(session02, file =file.path(io$output.data,'session.signature.RData'))
}
## session information to build this vignette
sessionInfo()
