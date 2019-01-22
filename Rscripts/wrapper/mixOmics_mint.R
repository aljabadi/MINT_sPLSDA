#########################  Create a log file
log_file =file.path("../Rscripts/wrapper/log", paste("mint.wrapper.run",format(Sys.time(), "%y_%m_%d"),"txt",sep = "."))
file.create(log_file)
cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "start preprocessing...\n"), file = log_file, append = TRUE)

#########################  Wrapper
mixOmics_mint = function(sce = sce, ## a combined sce with batch info in sce$batch 
                         colData.class = "mix", ## name of the colData that includes the biological groups
                         use.hvgs = F, ## whether to use HVGs in rowData(sce)$hi_var (T) or all genes (F)
                         ######### below can be left to default
                         ncomp = 3L, ## integer: good rule of thumb: number of biological groups - 1
                         keepX = c(50,40,30), ## number of genes to select on each of the ncomp components - length of at least ncomp
                         tune.hps = F, ## whether to tune the number of genes on each component (T) or use keepX (F)
                         tune.keepX = seq(10,100,10), ##  if tune.hps=T: span of number of genes on each component to explore for tuning
                         optimum.ncomp = F, ## whether to find the optimum number of sPLSDA components (T) or use given ncomp (F)
                         dist = "max.dist", ## if tune.hps: distance(s) to use for classification error rate.
                         ## subset of c("min.distance", "centroids.dist","mahalanobis.dist" ); or "all"

                         measure = "BER" # if tune.hps=T: measure to use for classification error rate; one of c("overall","BER")
){
  tp = system.time({
    try_res = try({
      
      #########################  the t.test.process fuction for optimisation 
      t.test.process = function(mat.error.rate, alpha = 0.01)
      {
        ## mat.error.rate has nrep rows and ncomp columns
        ## we test successively whether adding a component improves the results
        
        max = ncol(mat.error.rate) ## number max of components included
        pval = NULL
        opt = 1 ## initialise the first optimal number of components
        for(opt in 1:max)
        {
          j=opt+1
          ## t.test of "is adding X comp improves the overall results"
          temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=T)
          ## temp can be NaN when error.keepX is constant
          if(any(class(temp) == "try-error") || is.na(temp)) 
          {
            pval = 1
          } else {
            pval = temp
          }
          
          while(pval> (alpha) & j<max)
          {
            j=j+1
            ## t.test of "is adding X comp improves the overall results"
            temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=T)
            ## temp can be NaN when error.keepX is constant
            if(any(class(temp) == "try-error") || is.na(temp))
            {
              pval = 1
            } else {
              pval = temp
            }
          }
          ## if all p-values were greater than alpha, then we do not increase opt and get out of the loop
          if( (pval> (alpha)))
            break
        }
        ncomp_opt = opt
        
        return(ncomp_opt)
      }
      
      #########################  entry checks 
      
      Y = as.factor(colData(sce)[[colData.class]]) ## biological groups (e.g. cell types)
      batch = as.factor(sce$batch) ## factor for batches
      {
      ## length of keepX is at least ncomp
      if (length(keepX)<ncomp)
        stop ("The length of keepX should be at least ncomp")
      
      ## if ncomp optimisation is required, ncomp>1
      if(optimum.ncomp & !isTRUE(ncomp>1))
        stop("ncomp must be an integer greater than 1 for optimisation")
      
      ## object is a sce
      if(class(sce)!="SingleCellExperiment") ## make sure it is SCE object
        stop("sce must be a SingleCellExperiment object with batch metadata")
      
      ## batches exist and are more than one
      if(is.null(batch)){
        stop("$batch data do not exist in rowData(sce)")
      } else if(!isTRUE(length(unique(batch))>1)){
        stop("there must be more than one batch in the data to use mint.splsda")
      }
      
      ## there are biological groups to discriminate
      if(is.null(Y)){
        stop("cell type data do not exist in colData(sce)")
      } else if(!isTRUE(length(unique(Y))>1)){
        stop("there must be more than one cell type in the data to use mint.splsda")
      }
      
      ## logical entries are actually logical
      if(!all(is.logical(use.hvgs),is.logical(tune.hps), is.logical(optimum.ncomp)))
        stop("use.hvgs, tune.hps and optimum.ncomp must be logical")
      
      ## ensure there are no duplicate cell names
      if(any(duplicated(colnames(sce)))){
        stop("There are duplicate cell names - change to make them unique")
      }
        
      ## ensure there are no duplicate gene names
      if(any(duplicated(rownames(sce)))){
        stop("There are duplicate gene names - change to make them unique")
      }
      
      ## if use.hvgs =T, hi_var exist
      if(use.hvgs&is.null(rowData(sce)$hi_var)){
        message("HVGs are not specified in rowData(sce)$hi_var - continuing with all genes")
      }
      
      if(optimum.ncomp & nlevels(batch)<3)
        warning("The number of components cannot be reliably optimised
                since there are less than 3 batches.
                Regard the results with care.
                Refer to mixOmics documentation for details.")
      
      ## MINT checker checks the rest
      
      ######################### < /entry checks >
      
      if (max(logcounts(sce))>100){
        logcounts(sce) = log2(logcounts(sce)+1)
      }
      
      ## reduce sce if only HVGs are needed
      if(use.hvgs){
        sce = sce[rowData(sce)$hi_var,]
      }
      }
      ###################################### MINT sPLSDA
      
      ## check if it need to be tuned or optimised:
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
          measure=measure,
          progressBar = F)
        
        ## change keepX to the tuned vector
        keepX = mint.tune$choice.keepX
        ## choose the optimum ncomp if required
        if(optimum.ncomp){
          ncomp = t.test.process(mint.tune$error.rate)
        }
      }
      
      ## if only ncomp is to be optimised
      if(optimum.ncomp&!tune.hps){
        
        ## get the tuned unoptimised sparse MINT
        mint.res= mint.splsda(X=t(logcounts(sce)),
                              Y = Y,
                              study = batch,
                              ncomp = ncomp,
                              keepX = keepX
        )
        
        ## evaluate method's performance at different no. of components
        mint.performance = perf(mint.res, progressBar = F)
        # use a t-test to find the optimum no. of components
        ncomp = t.test.process(mint.performance$global.error$error.rate.class[[dist]])
      }
      
      
      ## final model using tuned/optimum parameters
      mint.res = mint.splsda(
        X = t(logcounts(sce)),
        Y = Y,
        study = batch,
        ncomp=ncomp,
        keepX = keepX
      )
      
      ## get a vector of marker genes for output in rowData(sce)
      comp =1
      markers = NULL
      while(comp <= ncomp){
        markers = c(markers, selectVar(mint.res, comp = comp)$name)
        comp = comp+1
      }
      ## add a logical rowData as to whether the gene is a marker
      rowData(sce)$mint.markers <- row.names(sce) %in% markers
      ## add the sPLSDA variates for visualisation to reducedDim(sce)
      reducedDim(sce, "mint_variates") = mint.res$variates$X

      
    })
    ###################################### outputs
    
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