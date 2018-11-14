## if loading the required input data locally
## input/output directories:
io = list()

## SCE data dir - FALSE for GitHub load:
# io$local.sincell = '../data/sincell_with_class.RData'
io$local.sincell = F # F or a directory

## where to save the run - FALSE for not saving:
# io$save.runs = '../output'
io$save.runs = F # F or a directory

## where to save the R scripts:
io$Rscript.dir = '../R-scripts'

## have the final mint.splsda object saved in io$save.runs for signature analyses?:
io$save.mint.spls = F # T or F

## DE tables directory for signature chapter - FALSE for GitHub load:
# io$DEtables = '../data/DEtable_90cells.RData'
io$DEtables = F # F or a directory
## load the required libraries
library(SingleCellExperiment)
library(mixOmics)
library(scran)
library(scater)
library(knitr)
library(VennDiagram)
library(tibble)
## load from github
## raw
RawURL='https://tinyurl.com/sincell-with-class-RData-LuyiT'
load(url(RawURL))
## ## or load from local directory, change to your own
## load(io$local.sincell)
## normalise the QC'ed count matrices
sc10x.norm =  computeSumFactors(sce10x_qc) ## deconvolute using size factors
sc10x.norm =  normalize(sc10x.norm) ## normalise expression values
## DROP-seq
scdrop.norm = computeSumFactors(scedrop_qc_qc)
scdrop.norm = normalize(scdrop.norm)
## CEL-seq2
sccel.norm =  computeSumFactors(sce4_qc)
sccel.norm =  normalize(sccel.norm)
## find the intersect of the genes
list.intersect = Reduce(intersect, list(
## the rownames of the original (un-transposed) count matrix will output the genes
  rownames(logcounts(sc10x.norm)),
  rownames(logcounts(sccel.norm)),
  rownames(logcounts(scdrop.norm))
))
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
perf.plsda.10x = perf(plsda.10x.res, progressBar=FALSE )
run.time = Sys.time()-start
## optimal number of components
plot(perf.plsda.10x, col = color.mixo(5:7))
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
          legend  = TRUE, title     = 'PLSDA 10X', 
          ellipse = TRUE, legend.title = 'Cell Line',
          X.label = 'PLSDA component 1', 
          Y.label = 'PLSDA component 2', pch=1)
## mint.plsda plot for CEL-seq2
plotIndiv(plsda.cel.res,
          legend  = TRUE, title     = 'PLSDA CEL-seq2', 
          ellipse = TRUE, legend.title = 'Cell Line',
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
          legend  = TRUE, title     = 'sPLSDA - 10X',
          ellipse = FALSE,legend.title = 'Cell Line',
          pch=1,
          X.label = 'sPLSDA component 1',
          Y.label = 'sPLSDA component 2')

## CEL-seq2
plotIndiv(splsda.cel.res, group = Y.cel,
          legend  = TRUE, title     = 'sPLSDA - CEL-seq2',
          ellipse = FALSE,legend.title = 'Cell Line',
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
data.combined = t( ## transpose of all 3 datasets combined
  data.frame(
    ## the genes from each protocol that match list.intersect
    logcounts(sc10x.norm)[list.intersect,],
    logcounts(sccel.norm)[list.intersect,],
    logcounts(scdrop.norm)[list.intersect,] ))

## create a factor variable of cell lines
## must be in the same order as the data combination
cell.line = as.factor(c(sce10x_qc$cell_line,
                         sce4_qc$cell_line,
                         scedrop_qc_qc$cell_line))
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
                              keepX = c(35,10))
## ## change to your own
## load(file.path(io$save.runs,'mint.splsda.tuned.res.RData'))
## 10X
plotLoadings(splsda.10x.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=FALSE, title=NULL, 
             subtitle = '10X')
## CEL-seq2
plotLoadings(splsda.cel.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=FALSE, title=NULL,
             subtitle = 'CEL-seq2')
## Drop-seq
plotLoadings(splsda.drop.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=FALSE, title=NULL,
             subtitle = 'Drop-seq')
## MINT Comp. 1
plotLoadings(mint.splsda.res, contrib='max', method = 'mean', comp=1, 
             study='all.partial', legend=FALSE, title=NULL, 
             subtitle = c('10X', 'CEL-seq2', 'Drop-seq') )
## MINT Comp. 2
plotLoadings(mint.splsda.tuned.res, contrib='max', method = 'mean', comp=2, 
             study='all.partial', legend=FALSE, title=NULL, 
             subtitle = c('10X', 'CEL-seq2', 'Drop-seq') )
## MINT signature
MINT.Combined.vars = unique(c(selectVar(mint.splsda.tuned.res, comp=1)$name,
                             selectVar(mint.splsda.tuned.res, comp=2)$name))
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
cim(mint.splsda.tuned.res, comp = c(1,2), margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(mint.splsda.tuned.res$Y)), row.names = FALSE,
    title = 'MINT sPLS-DA', save='png', name.save = 'heatmap')
## load directly from github
DE_URL <- 'https://tinyurl.com/DEtable-90cells'
load(url(DE_URL))
## ## or load the DEGs from the 'data' folder
## load(io$DEtables)
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
## ## run this code if you wish to save the RData for loading
## save.image(file = file.path(io$save.runs,'02-Signature.RData'))
## session information to build this vignette
sessionInfo()
