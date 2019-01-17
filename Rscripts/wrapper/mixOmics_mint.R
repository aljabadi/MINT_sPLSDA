mixOmics_mint = function(sce = sce, ## a combined sce with sce$batch that
                         ## icludes $batch and $cell_line (or similar)
                         colData.class = "mix", ## name of the colData that includes the (pseudo)cell types
                         use.hvgs = F, ## whether to use HVGs in rowData(sce)$hi_vars (T) or all genes (F)
                         ######### below can be left to default
                         ncomp = 5, ## at least 2 recommended - choose based on how many subpopulations there might be
                         keepX = c(50,40,30,20,10), ## a numeric (if ncomp=1) or a vector with same length as ncomp, number of features on each component
                         tune.hps = F, ## whether to tune hyperparameters (T) or use keepX number of features for all components (F)
                         tune.keepX = seq(10,100,10), ##  if tune.hps=T: span of number of features to explore for tuning
                         optimum.ncomp = F, ## whether to find the optimum number of sPLSDA components (T) or just use ncomp (F)
                         dist = "max.dist", ## if tune.hps || optimum.comp =T: distance to use for classification error rate
                         measure = "BER" # if tune.hps=T: measure to use for classification error rate; one of c("overall","BER")
){
  tp = system.time({
    try_res = try({
      
      ######################### < entry checks >
      
      Y = as.factor(colData(sce)[[colData.class]]) ## cell types/classes
      batch = as.factor(sce$batch) ## batch vector
      
      ## ncomp and keepX match - otherwise all genes will be included in excess comps
      if (length(keepX)<ncomp){ ## if keepX is a vector with size less than ncomp
        message ("keepX modifed to match ncomp")
        Diff = ncomp - length(keepX)
        ## repeat the last keepX for all the rest
        keepX = c(keepX, rep(keepX[length(keepX)],Diff))
      }
      
      ## sce checks
      if(class(sce)!="SingleCellExperiment") ## make sure it is SCE object
        stop("sce must be a SingleCellExperiment object with batch metadata")
      
      ## batch checks
      if(is.null(batch)){
        stop("$batch data do not exist in rowData(sce)")
      } else if(!isTRUE(length(unique(batch))>1)){
        stop("there must be more than one batch in the data to use mint.splsda")
      }
      
      ## cell class checks
      if(is.null(Y)){
        stop("cell type data do not exist in colData(sce)")
      } else if(!isTRUE(length(unique(Y))>1)){
        stop("there must be more than one cell type in the data to use mint.splsda")
      }
      
      ## logicals
      if(!all(is.logical(use.hvgs),is.logical(tune.hps), is.logical(optimum.ncomp)))
        stop("use.hvgs, tune.hps and optimum.ncomp must be logical")
      
      ## ensure there are no duplicate cell names
      if(any(colnames(sce)!= make.unique(colnames(sce)))){
        message("changed duplicate cell names")
        colnames(sce) = make.unique(colnames(sce))
      }
      
      ## hi_vars exist
      if(use.hvgs&is.null(rowData(sce)$hi_vars)){
        message("HVGs are not specified in rowData(sce)$hi_vars - continuing with all genes")
      }
      
      ## MINT checker checks the rest
      
      ######################### < /entry checks >
      
      if (max(logcounts(sce))>100){
        logcounts(sce) = log2(logcounts(sce)+1)
      }
      
      ## reduce sce if only HVGs are needed
      if(use.hvgs){
        sce = sce[rowData(sce)$hi_vars,]
      }
      
      ######################################
      ################## < wrapper > #######
      ######################################
      
      ## get the untuned unoptimised sparse MINT
      mint.res= mint.splsda(X=t(logcounts(sce)),
                            Y = Y,
                            study = batch,
                            ncomp = ncomp,
                            keepX = keepX
      )
     
      
      if(optimum.ncomp){
        ## evaluate method's performance at different no. of comp.s
        mint.performance = perf(mint.res, progressBar = F, dist=dist )
        ## change the number of comp's to the one proding minimum error
        ## the choice.ncomp argument has been changed
        ncomp = which.min(mint.performance$global.error[[measure]])
        keepX = keepX[1:ncomp]
      }

      if(tune.hps){
        ## tune the number of markers
        mint.tune = tune.mint.splsda(
          X = t(logcounts(sce)),
          Y = Y,
          study =  batch,
          ## get the optimum for measure and dist
          ncomp=ncomp,
          ## assess tune.keepX number of features:
          test.keepX = tune.keepX,
          ## use dist to estimate the classification error rate
          dist = dist,
          progressBar = F
        )
        
        ## change keepX to the tuned vector
        keepX = mint.tune$choice.keepX
      }
      
      
      ## final model using tuned parameters
        mint.res = mint.splsda(
        X = t(logcounts(sce)),
        Y = Y,
        study = batch,
        ncomp=ncomp,
        keepX = keepX
      )
    
      ## marker genes
      comp =1
      markers = NULL
      while(comp <= ncomp){
        markers = c(markers, selectVar(mint.res, comp = comp)$name)
        comp = comp+1
      }
      ## add a logical rowData as to whether the gene is a marker
      rowData(sce)$mint.markers <- row.names(sce) %in% markers
      ## add the sPLSDA variates for visualisation
      reducedDim(sce, "mint_variates") = mint.res$variates$X
      
      ######################################
      ################## < /wrapper > #######
      ######################################
      
    })
    if (class(try_res) == "try-error") {
      cat(paste(format(Sys.time(), "%a %b %d %X %Y. ERROR: "), print(try_res),"\n"), file = log_file, append = TRUE)
    }
  })
  
  method_name = "mixOmics_mint"
  method_type = "visualisation"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}