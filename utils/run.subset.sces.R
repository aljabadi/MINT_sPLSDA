###### subset sce objects

## load the sces
load('../data/sincell_with_class.RData')
## list of all sce objects that include $cell_line info like below - all named 'sce'
protocols = list(sce10x_qc=list('sce'=sce10x_qc),
                 sce4_qc=list('sce'=sce4_qc),
                 scedrop_qc_qc=list('sce'=scedrop_qc_qc))

n.genes = 300 ## number of genes desired to be kept
pct.gene.comm = 0.8 ## portion of genes to be common among all
n.cell.tot = 60 ## total number of cells desired to be kept equally from all cell lines from all protocols

## run subset.sces.R
source('../utils/subset.sces.R')
subsets = subset.sces(protocols=protocols,n.genes=n.genes,
                      pct.gene.comm = pct.gene.comm, n.cell.tot=n.cell.tot)

## create variable names using protocols given
## save them in the same name as the original sce for reproducibility
for (i in 1:length(protocols)){
  name = paste0(names(protocols[i]))
  assign(name, subsets[[i]])
}
## save the combined RData file in specific name for ease of recognition - you can be more specific with the name
save(list=c(names(protocols)),file=paste0('../data/subset/subset_sces_',n.genes,'_',n.cell.tot,'.RData'))
