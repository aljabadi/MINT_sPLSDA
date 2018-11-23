### function

subset.sces = function(protocols = list(sce10x_qc=list('sce'=sce10x_qc), ## list of all sce objects that include $cell_line info
                                        sce4_qc=list('sce'=sce4_qc),
                                        scedrop_qc_qc=list('sce'=scedrop_qc_qc)),
                       n.genes = 500, ## total # of genes kept
                       pct.gene.comm = 0.8, ## percent of genes common among all
                       n.cell.tot = 90 ## total # of cells kept
                       ){
######################################### INPUT ^^^^^
cell.lines = unique(protocols[[1]]$sce$cell_line) ## cell lines
n.each = floor(n.cell.tot)/length(cell.lines) ## cells in each cell line

n.uniq = round((1 - pct.gene.comm )*n.genes) ## no. of unique genes to be kept
n.genes = round(pct.gene.comm*n.genes) ## no. of common genes

## get the gene names from all sces
rownames.list = list()
for (i in 1:length(protocols)){
  protocol = protocols[[i]]
  rownames.list[[i]] = rownames(protocol$sce)
}
#################### Reduce genes
## common genes
common.genes = data.frame(row.names =Reduce(intersect, rownames.list))
## calculate HVGs with a pseudo mean to deflate BCV near zero
common.genes$cv = apply(counts(protocols[[1]][["sce"]][rownames(common.genes),]),1,function(x) sd(x)/(mean(x)+3))
## only adding to have 2 columns in df
common.genes$HVG = common.genes$cv>quantile(common.genes$cv, 0.75)
## order genes based on BCV
common.genes=common.genes[with(common.genes, order(-cv)),]


## non-common genes
for (i in 1:length(protocols)){
  genes = as.vector(rownames(protocols[[i]]$sce[!(rownames(protocols[[i]]$sce) %in% rownames(common.genes))]))
  names(genes)= genes
  protocols[[i]]$uniq.genes = genes[1:n.uniq]
}

## cell.lines
for (i in 1:length(protocols)){
  protocol = protocols[[i]]$sce
  cells=list()
  for (j in cell.lines){
    cells[[j]] = colnames(protocol[,protocol$cell_line==j])[1:n.each]
  }
  protocols[[i]]$cells = unlist(cells)
}

## common HVGs
common.hvgs = rownames(common.genes)[1:n.genes]
names(common.hvgs) = common.hvgs

## add a $subset object to protocols
protocols$subsets = list()
for (i in 1:length(protocols)){
  ## name of the object to add
  name = paste0('subset_', names(protocols[i]))
  ## common + unique genes
  genes = c(common.hvgs, protocols[[i]]$uniq.genes)
  ## add with subsetting on genes and cells
  protocols$subsets[[name]] = protocols[[i]]$sce[genes,protocols[[i]]$cells ]
}

return(unlist(protocols$subsets))

}