##
## function: gsva
## purpose: main function of the package which estimates activity
##          scores for each given gene-set
#setwd("E:/single-rna/scripts/my")
#system("R CMD SHLIB ks_test.c")

########### xiong
##src_folder = '/home/gfj/scTPA_2/src/'
suppressPackageStartupMessages(library("this.path"))
dirCurrent=this.dir()
src_folder=file.path(dirCurrent,'/')
print("###### PAS.R")
#print(dirCurrent)
########### xiong end

suppressPackageStartupMessages(library(Rcpp))
sourceCpp(file.path(src_folder,'ker.cpp'))
sourceCpp(file.path(src_folder,'ks_test.cpp'))

print('source PAS success')
setGeneric("PAS", function(expr, gset.idx.list, ...) standardGeneric("PAS"))



setMethod("PAS", signature(expr="matrix", gset.idx.list="list"),
          function(expr, 
                   gset.idx.list, 
                   weights,
                   #tmp_folder,
                   method=c("gsva", "ssgsea", "zscore", "plage"),
                   kcdf=c("Gaussian", "Poisson", "none"),
                   abs.ranking=FALSE,
                   min.sz=3,
                   max.sz=Inf,
                   parallel.sz=0, 
                   parallel.type="SOCK",
                   mx.diff=TRUE,
                   tau=switch(method, gsva=1, ssgsea=0.25, NA),
                   ssgsea.norm=TRUE,
                   verbose=TRUE)
{
  method = match.arg(method)
  kcdf = match.arg(kcdf)

  ## filter out genes with constant expression values
  sdGenes = apply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr = expr[sdGenes > 0 & !is.na(sdGenes), , drop=FALSE]   # filtering genes
    }
  } 
  
  if(method == "plage"){
    
    min.sz = max(2, min.sz)
  }

  if (nrow(expr) < 2)
    stop("Less than two genes in the input expression data matrix\n")
  
  if(class(weights) == 'character'){          
    mapped.gset.idx.list = lapply(gset.idx.list,
                                   function(x ,y) na.omit(match(x, y)),
                                   rownames(expr))
    mapped.weight.idx.list = 'none'
    #print("min cells ....................")
    #print(min.sz) 
    mapped.gset.idx.list = filterGeneSets(mapped.gset.idx.list,
                                           min.sz=max(3, min.sz),
                                           max.sz=max.sz)
  }else{
    mapped.gset.idx.list = vector(mode='list',length=length(gset.idx.list))
    mapped.weight.idx.list = vector(mode='list',length = length(gset.idx.list))
    for(i in 1:length(gset.idx.list)){
      match_i = match(gset.idx.list[[i]], rownames(expr))
      mapped.gset.idx.list[[i]] = match_i[!is.na(match_i)]
      mapped.weight.idx.list[[i]] = weights[[i]][!is.na(match_i)]
    }
    names(mapped.gset.idx.list) = names(gset.idx.list)
    names(mapped.weight.idx.list) = names(gset.idx.list)
    
    mapped.gset.idx.list = filterGeneSets(mapped.gset.idx.list,
                                           min.sz=max(3, min.sz),
                                           max.sz=max.sz)
    mapped.weight.idx.list = filterGeneSets(mapped.weight.idx.list,
                                             min.sz=max(3, min.sz),
                                             max.sz=max.sz)
  }


  if (length(unlist(mapped.gset.idx.list, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user
  
  if(length(gset.idx.list) == 0){
    stop("The gene set list is empty!  Filter may be too stringent.")
  }
  
  if (!missing(kcdf)) {
    if (kcdf == "Gaussian") {
      rnaseq = FALSE
      kernel = TRUE
    } else if (kcdf == "Poisson") {
      rnaseq = TRUE
      kernel = TRUE
    } else
      kernel = FALSE
  }
  
  

  rval = .PAS(expr, 
               mapped.gset.idx.list, 
               mapped.weight.idx.list,  
               method, 
               kcdf, 
               rnaseq, 
               abs.ranking,
               parallel.sz, 
               parallel.type, 
               mx.diff, 
               tau, 
               kernel, 
               ssgsea.norm, 
               verbose)
  print('parallel.sz........................')
  print(parallel.sz)

  rval
})

.PAS = function( expr, 
                 gset.idx.list,
                 weights,
                 #tmp_folder,
                 method=c("gsva", "ssgsea", "zscore", "plage"),
                 kcdf=c("Gaussian", "Poisson", "none"),
                 rnaseq=FALSE,
                 abs.ranking=FALSE,
                 parallel.sz=0, 
                 parallel.type="SOCK",
                 mx.diff=TRUE,
                 tau=1,
                 kernel=TRUE,
                 ssgsea.norm=TRUE,
                 verbose=FALSE)
{
	
  if (method == "ssgsea") {
	  if(verbose)
		  cat("Estimating ssGSEA scores for", length(gset.idx.list),"gene sets.\n")
    
    if(class(weights) == 'character'){
      return(ssgsea_noWeight(expr, 
                             gset.idx.list, 
                             alpha=tau, 
                             parallel.sz=parallel.sz,
                             parallel.type=parallel.type, 
                             normalization=ssgsea.norm,
                             verbose=verbose))
    }else{
      return(ssgsea_Weight(expr, 
                           gset.idx.list, 
                           weights, 
                           alpha=tau, 
                           parallel.sz=parallel.sz,
                           parallel.type=parallel.type, 
                           normalization=ssgsea.norm,
                           verbose=verbose))
    }
    
  }

  if (method == "zscore") {
    if (rnaseq)
      stop("rnaseq=TRUE does not work with method='zscore'.")

	  if(verbose)
		  cat("Estimating combined z-scores for", length(gset.idx.list),"gene sets.\n")

    return(zscore(expr, 
                  gset.idx.list, 
                  weights,  
                  parallel.sz, 
                  parallel.type, 
                  verbose))
  }

  if (method == "plage") {
    if (rnaseq)
      stop("rnaseq=TRUE does not work with method='plage'.")

	  if(verbose)
		  cat("Estimating PLAGE scores for", length(gset.idx.list),"gene sets.\n")

    return(plage(expr, 
                 gset.idx.list, 
                 weights,  
                 parallel.sz, 
                 parallel.type, 
                 verbose))
  }
  
  if (method == "gsva"){
    if(verbose)
      cat("Estimating GSVA scores for", length(gset.idx.list),"gene sets.\n")
    if(class(weights) == 'character'){
      return(gsva_noWeight(expr, 
                           gset.idx.list, 
                           parallel.sz,
                           parallel.type, 
                           abs.ranking, 
                           mx.diff, 
                           tau,
                           kernel, 
                           verbose))
    }else{
      return(gsva_Weight(expr, 
                         gset.idx.list, 
                         weights,  
                         parallel.sz,
                         parallel.type, 
                         abs.ranking, 
                         mx.diff, 
                         tau,
                         kernel, 
                         verbose))
    }
  }
  
}

gsva_noWeight = function(expr, 
                         gset.idx.list,
                         #tmp_folder,
                         parallel.sz,
                         parallel.type,
                         abs.ranking,
                         mx.diff,
                         tau,
                         kernel,
                         verbose=F,
                         rnaseq=FALSE){
  n = ncol(expr)
  
  if (verbose)
    cat("Computing observed enrichment scores\n")
  
  num_genes = nrow(expr)
  num_sample = ncol(expr)
  haveParallel = .isPackageLoaded("parallel")
  #print("is have paralle??????")
  #print(haveParallel)
  haveSnow = .isPackageLoaded("snow")
  #print("is have snow??????")
  #print(haveSnow)
  if (verbose) {
    if (kernel) {
      if (rnaseq)
        cat("Estimating ECDFs with Poisson kernels\n")
      else
        cat("Estimating ECDFs with Gaussian kernels\n")
    } else
      cat("Estimating ECDFs directly\n")
  }
  
  gene.density = compute.gene.density(expr, parallel.sz, rnaseq, kernel)
  
  compute_rank_score = function(sort_idx_vec){
    tmp = rep(0, num_genes)
    tmp[sort_idx_vec] = abs(seq(from=num_genes,to=1) - num_genes/2)
    return (tmp)
  }
  
  rank.scores = rep(0, num_genes)
  
  mclapp = get('mclapply', envir=getNamespace('parallel'))
  options(mc.cores=parallel.sz)
  
  #sort.sgn.idxs = apply(gene.density, 2, order, decreasing=TRUE) # d;\
  sort.sgn.idxs = mclapply(1:ncol(gene.density), 
	               function(i){
				      order(gene.density[,i], decreasing=TRUE)
				   })
  sort.sgn.idxs = do.call("cbind", sort.sgn.idxs)
  rank.scores = mclapp(1:ncol(sort.sgn.idxs), 
	               function(i){
				      compute_rank_score(sort.sgn.idxs[,i])
				   })	    
  rank.scores = do.call("cbind", rank.scores)
  #rank.scores = apply(sort.sgn.idxs, 2, compute_rank_score)  

  print("density success")
  haveParallel = .isPackageLoaded("parallel")
  haveSnow = .isPackageLoaded("snow")
  
  if (parallel.sz > 1 || haveParallel) {
    if (!haveParallel && !haveSnow) {
      stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
    }
    
    if (haveSnow) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl = get("makeCluster", mode="function")
      parSapp = get("parSapply", mode="function")
      #clEvalQ = get("clusterEvalQ", mode="function")
      stopCl = get("stopCluster", mode="function")
      
      if (verbose)
        cat("Allocating cluster\n")
      cl = makeCl(parallel.sz, type = parallel.type) 
      #clEvalQ(cl, library(GSVA))
      if (verbose) {
        cat("Estimating enrichment scores in parallel\n")
        if(mx.diff) {
          cat("Taking diff of max KS.\n")
        } else{
          cat("Evaluting max KS.\n")
        }
      }
      
      m = t(parSapp(cl, gset.idx.list, 
                     ks_matrix_gsva,     
                     gene.density=rank.scores, 
                     sort.idxs=sort.sgn.idxs,
                     tau=tau,
                     mx_diff=as.integer(mx.diff), 
                     abs_rnk=as.integer(abs.ranking)))
      if(verbose)
        cat("Cleaning up\n")
      stopCl(cl)
      
    } else if (haveParallel) {             ## use parallel
      #print("use parallel..........................")
      mclapp = get('mclapply', envir=getNamespace('parallel'))
      detCor = get('detectCores', envir=getNamespace('parallel'))
      nCores = detCor()
      #print(nCores)
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      
      pb = NULL
      if (verbose){
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
        assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
        assign("nGeneSets", ceiling(length(gset.idx.list) / getOption("mc.cores")), envir=globalenv())
        assign("iGeneSet", 0, envir=globalenv())
      }
      
      m = mclapp(gset.idx.list, 
                 ks_matrix_gsva,
                 rank.scores,
                 sort.sgn.idxs,
                 tau=tau,
                 mx_diff=as.integer(mx.diff), 
                 abs_rnk=as.integer(abs.ranking))

      m = do.call("rbind", m)
      colnames(m) = colnames(expr)
      
      if (verbose) {
        close(get("progressBar", envir=globalenv()))
      }
    } else
      stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
    
  } else {
    if (verbose) {
      cat("Estimating enrichment scores\n")
      if (mx.diff) {
        cat("Taking diff of max KS.\n")
      } else{
        cat("Evaluting max KS.\n")
      }
    }
    pb = NULL
    if (verbose){
      assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
      assign("nGeneSets", length(gset.idx.list), envir=globalenv())
      assign("iGeneSet", 0, envir=globalenv())
    }
    #print("no parallel........................")
    m = t(sapply(gset.idx.list, 
                 ks_matrix_gsva, 
                 rank.scores, 
                 sort.sgn.idxs,
                 tau=tau,
                 mx_diff=as.integer(mx.diff), 
                 abs_rnk=as.integer(abs.ranking)))
    
    if (verbose) {
      setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
      close(get("progressBar", envir=globalenv()))
    }
  }
  if (length(gset.idx.list) == 1)
    m = matrix(m, nrow=1)
  print(dim(m))
  colnames(m) = colnames(expr)
  #rownames(m) = names(gset.idx.list)
  
  m
}

cal_AUCell = function(counts,
                      gSets,
                      n_cores
                      # number of cores used for parallel
                      ){

  ## matrix could be integer or numeric

  tryCatch({
print("######cal_AUCell: step1")
    # @must step1: rank genes
    ac_rankings = AUCell::AUCell_buildRankings(counts[,],
                                               # matrix, dgCMatrix, SummarizedExperiment, ExpressionSet
                                               nCores=n_cores,
                                               plotStats=FALSE,
                                               # plot the expression boxplots or histograms
                                               #assayName = NULL,
                                               # slot name of assay containing expression matrix
                                               verbose = F)
    # @must step2: calculate AUC scores
print("######cal_AUCell: step2")
    sc_AUC = AUCell::AUCell_calcAUC(gSets,
                                    ac_rankings,
                                    normAUC=T,
                                    # Wether to normalize the maximum possible AUC to 1
                                    aucMaxRank=ceiling(0.05 * nrow(ac_rankings)),
                                    # the number of genes (maximum ranking) that is used to computation
                                    # default: 5% of pathways; recommend: 1%-20%
                                    verbose = F
    )
print("######cal_AUCell: step3")

    score = AUCell::getAUC(sc_AUC)

print("######cal_AUCell: finshed")

    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

cal_vision = function(counts,
                      gSets_path,
                      n_cores){


  ## matrix could be integer or numeric

  ### other parameters for Vsion:
  # proteinData = NULL, unnormalizedData = NULL, meta = NULL,
  # projection_genes = c("fano"), min_signature_genes = 5,
  # sig_gene_threshold = 0.001, threshold = 0.05, perm_wPCA = FALSE,
  # projection_methods = c("tSNE30"),
  # sig_norm_method = c("znorm_columns", "none", "znorm_rows",
  #                     "znorm_rows_then_columns", "rank_norm_columns"),
  # pool = "auto",
  # cellsPerPartition = 10, name = NULL, num_neighbors = NULL,
  # latentSpace = NULL, latentSpaceName = NULL,
  # latentTrajectory = NULL, pools = list()


  ## counts recommend for scale or normalized, but not transformed
  ## The expression data should not be log-transformed prior to loading into VISION.


  tryCatch({
    #if(ncol(counts)<100){
    #  projection_method = 'tSNE10'
    #}else{
    #  projection_method = 'tSNE30'
    #}

    vis = VISION::Vision(counts,            ## Gene X Cell
                         # data.frame; sparseMatrix; dgeMatrix; ExpressionSet; SummarizedExperiment; Seurat
                         signatures = gSets_path,
                         projection_method = 'UMAP',
                         sig_gene_threshold=0)

    options(mc.cores=n_cores)
    vis = VISION::analyze(vis)

    score = t(vis@SigScores)    ## pathway X cell
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

cal_pagoda2 = function(n_counts,
                       gSets,
                       trim = 5,
                       n_cores){


  ### must be counts matrix !!!!!

  ### other parameters for knn.error.models
  # min.nonfailed = 5, min.count.threshold = 1,
  # max.model.plots = 50,
  # min.size.entries = 2000, min.fpm = 0, cor.method = "pearson",
  # verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE,
  # local.theta.fit = linear.fit, theta.fit.range = c(0.01, 100),
  # alpha.weight.power = 1/2

  ### other parameters for pagoda.varnorm
  # batch = NULL, prior = NULL,
  # fit.genes = NULL, minimize.underdispersion = FALSE,
  # n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9,
  # verbose = 0, weight.df.power = 1, smooth.df = -1,
  # theta.range = c(0.01, 100), gene.length = NULL

  nPcs = min(round(ncol(n_counts)/5),5)
  #counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
  tryCatch({
    #print("00000000000000000000000: PAS")

    p2 = Pagoda2$new(n_counts, n.cores = n_cores,log.scale=F)
    #print("11111111111111111111111: PAS")

    p2$adjustVariance(plot=F)
    #print("22222222222222222222222: PAS")

    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    #print("333333333333333333333333: PAS")

    path_names = c()
    env = new.env(parent=globalenv())
    invisible(lapply(1:length(gSets),function(i) {
      genes = intersect(gSets[[i]],rownames(n_counts))
      name = paste0(names(gSets[i]),i)
      if(length(genes)>3){
        assign(name, genes, envir = env)
        path_names = c(path_names, name)
      }
    }))

    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = T,
                                 min.pathway.size = 1)

    path_names = names(p2@.xData$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(n_counts))
    rownames(score) = path_names
    colnames(score) = colnames(n_counts)
    for(i in 1:length(p2@.xData$misc$pwpca)){
      if(!is.null(p2@.xData$misc$pwpca[[i]]$xp$score)){
        score[i,] = as.numeric(p2@.xData$misc$pwpca[[i]]$xp$scores)
      }
    }

    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}

gsva_Weight = function(expr, 
                       gset.idx.list,
                       weights,
                       parallel.sz,
                       parallel.type,
                       abs.ranking,
                       mx.diff,
                       tau,
                       kernel,
                       verbose=F,
                       rnaseq=FALSE){
  if(verbose)
    cat("Estimating GSVA scores for", length(gset.idx.list),"gene sets.\n")
  
  if (verbose)
    cat("Computing observed enrichment scores\n")
  
  num_genes = nrow(expr)
  num_sample = ncol(expr)
  haveParallel = .isPackageLoaded("parallel")
  haveSnow = .isPackageLoaded("snow")
  if (verbose) {
    if (kernel) {
      if (rnaseq)
        cat("Estimating ECDFs with Poisson kernels\n")
      else
        cat("Estimating ECDFs with Gaussian kernels\n")
    } else
      cat("Estimating ECDFs directly\n")
  }
  #print("compute egen density")
  gene.density = compute.gene.density(expr, rnaseq, kernel, parallel.sz)
  
  compute_rank_score = function(sort_idx_vec){
    tmp = rep(0, num_genes)
    tmp[sort_idx_vec] = abs(seq(from=num_genes,to=1) - num_genes/2)
    return (tmp)
  }
  
  ks_weight = function(i, gene.density, geneSets, weights){
    gene.density[geneSets[[i]],] = gene.density[geneSets[[i]],]*weights[[i]]
    rank.scores = rep(0, num_genes)
    sort.sgn.idxs = apply(gene.density, 2, order, decreasing=TRUE) 
    rank.scores = apply(sort.sgn.idxs, 2, compute_rank_score)
    
    ks_matrix_gsva(geneSets[[i]], rank.scores, sort.sgn.idxs, tau, as.integer(mx.diff), as.integer(abs.ranking))
  }
  
  haveParallel = .isPackageLoaded("parallel")
  haveSnow = .isPackageLoaded("snow")
  
  if (parallel.sz > 1 || haveParallel) {
    if (!haveParallel && !haveSnow) {
      stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
    }
    
    if (haveSnow) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl = get("makeCluster", mode="function")
      parSapp = get("parSapply", mode="function")
      #clEvalQ = get("clusterEvalQ", mode="function")
      stopCl = get("stopCluster", mode="function")
      
      if (verbose)
        cat("Allocating cluster\n")
      cl = makeCl(parallel.sz, type = parallel.type) 
      #clEvalQ(cl, library(GSVA))
      if (verbose) {
        cat("Estimating enrichment scores in parallel\n")
        if(mx.diff) {
          cat("Taking diff of max KS.\n")
        } else{
          cat("Evaluting max KS.\n")
        }
      }
      
      m = t(parSapp(cl, 1:length(gset.idx.list), 
                    ks_weight,    
                    gene.density, 
                    gset.idx.list,
                    weights))
      if(verbose)
        cat("Cleaning up\n")
      stopCl(cl)
      
    } else if (haveParallel) {             ## use parallel
      
      mclapp = get('mclapply', envir=getNamespace('parallel'))
      detCor = get('detectCores', envir=getNamespace('parallel'))
      nCores = detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      
      pb = NULL
      if (verbose){
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
        assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
        assign("nGeneSets", ceiling(length(gset.idx.list) / getOption("mc.cores")), envir=globalenv())
        assign("iGeneSet", 0, envir=globalenv())
      }
      
      m = mclapp(1:length(gset.idx.list), 
                 ks_weight,     
                 gene.density, 
                 gset.idx.list,
                 weights)
      
      m = do.call("rbind", m)
      colnames(m) = colnames(expr)
      
      if (verbose) {
        close(get("progressBar", envir=globalenv()))
      }
    } else
      stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
    
  } else {
    if (verbose) {
      cat("Estimating enrichment scores\n")
      if (mx.diff) {
        cat("Taking diff of max KS.\n")
      } else{
        cat("Evaluting max KS.\n")
      }
    }
    pb = NULL
    if (verbose){
      assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
      assign("nGeneSets", length(gset.idx.list), envir=globalenv())
      assign("iGeneSet", 0, envir=globalenv())
    }
    
    m = t(sapply(1:length(gset.idx.list), 
                 ks_weight,     
                 gene.density, 
                 gset.idx.list,
                 weights))
    
    if (verbose) {
      setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
      close(get("progressBar", envir=globalenv()))
    }
  }
  
  colnames(m) = colnames(expr)
  #rownames(m) = names(gset.idx.list)
  
  m
}


compute.gene.density = function(expr,parallel.sz, rnaseq=FALSE, kernel=TRUE){
	
  gene.density = NA
  if (kernel) {
      #print("calculating gene density with paralle........")
      if(parallel.sz>0){
      #print("calculating gene density with paralle........")
      mclapp = get('mclapply', envir=getNamespace('parallel'))
	  options(mc.cores=parallel.sz)
	  gene.density = mclapp(1:nrow(expr), 
	               function(i, expr, rnaseq){
				      row_d(expr[i,], rnaseq)
				   },
				   expr,
	               as.integer(rnaseq))
	  gene.density = do.call("rbind", gene.density)}
	  else{
	    print("calculating gene.density using cpp matrix")
	    gene.density = matrix_d(expr, as.integer(rnaseq))
	  }
	  #print(dim(gene.density))
	  #gene.density = matrix_d(expr, as.integer(rnaseq))
  } else {
    gene.density = t(apply(expr, 1, function(x) {
                                     f = ecdf(x[sample.idxs])
                                     f(x)
                                   }, sample.idxs))
    gene.density = log(gene.density / (1-gene.density))
  }
  rownames(gene.density) = rownames(expr)
  colnames(gene.density) = colnames(expr)

	return(gene.density)	
}



ks_test_m = function(gset_idxs, gene.density, sort.idxs, mx.diff=TRUE,
                      abs.ranking=FALSE, tau=1, verbose=F){
	
	n.genes = nrow(gene.density)
	n.samples = ncol(gene.density)
	n.geneset = length(gset_idxs)

	geneset.sample.es = .C("ks_matrix_R",
			as.double(gene.density),
			R = double(n.samples),
			as.integer(sort.idxs),
			n.genes,
			as.integer(gset_idxs),
			n.geneset,
			as.double(tau),
			n.samples,
			as.integer(mx.diff),
      as.integer(abs.ranking))$R

  if (verbose) {
    assign("iGeneSet", get("iGeneSet", envir=globalenv()) + 1, envir=globalenv())
    setTxtProgressBar(get("progressBar", envir=globalenv()),
                      get("iGeneSet", envir=globalenv()) / get("nGeneSets", envir=globalenv()))
  }
	
	return(geneset.sample.es)
}


## ks-test in R code - testing only
ks_test_Rcode = function(gene.density, gset_idxs, tau=1, make.plot=FALSE){
	
	n.genes = length(gene.density)
	n.gset = length(gset_idxs)
	
	sum.gset = sum(abs(gene.density[gset_idxs])^tau)
	
	dec = 1 / (n.genes - n.gset)
	
	sort.idxs = order(gene.density,decreasing=T)
	offsets = sort(match(gset_idxs, sort.idxs))
	
	last.idx = 0
	values = rep(NaN, length(gset_idxs))
	current = 0
	for(i in seq_along(offsets)){
		current = current + abs(gene.density[sort.idxs[offsets[i]]])^tau / sum.gset - dec * (offsets[i]-last.idx-1)
		
		values[i] = current
		last.idx = offsets[i]
	}
	check_zero = current - dec * (n.genes-last.idx)
	#if(check_zero > 10^-15){ 
	#	stop(paste=c("Expected zero sum for ks:", check_zero))
	#}
	if(make.plot){ plot(offsets, values,type="l") } 
	
	max.idx = order(abs(values),decreasing=T)[1]
	mx.value = values[max.idx]
	
	return (mx.value)
}






rndWalk_1 = function(i, geneSets, sample_expr, alpha, weights) {

  if(!class(weights) == 'character'){
    sample_expr[geneSets[[i]]] = sample_expr[geneSets[[i]]] * weights[[i]]
  }
  sample_R = as.integer(rank(sample_expr))
  geneRanking = order(sample_R,decreasing = T)
  indicatorFunInsideGeneSet = match(geneRanking, geneSets[[i]])
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] = 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] = 0
  index_score = (abs(sample_R[geneRanking]) * indicatorFunInsideGeneSet)^alpha  

  stepCDFinGeneSet = cumsum(index_score) / sum(index_score)

  stepCDFoutGeneSet = cumsum(!indicatorFunInsideGeneSet) /
                       sum(!indicatorFunInsideGeneSet)
  walkStat = stepCDFinGeneSet - stepCDFoutGeneSet               
  mx_neg = min(walkStat)
  mx_pos = max(walkStat)
  if(abs(mx_neg)>mx_pos){
    mx_neg
  }else{
    mx_pos
  }
  
}

#ks_sample(sample_R, geneRanking,geneSets[[i]], 0.25, 1, 0)

rndWalk_2 = function(i, geneSets, sample_expr, alpha, weights) {

  sample_R = as.integer(rank(sample_expr)) 
order
  geneRanking = order(sample_R,decreasing = T)
  
  indicatorFunInsideGeneSet = match(geneRanking, geneSets[[i]])
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] = 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] = 0
  
  if(class(weights) == 'character'){
    index_score = (abs(sample_R[geneRanking]) * indicatorFunInsideGeneSet)^alpha
  }else{
    curr_weight = rep(0,length(sample_expr))
    curr_weight[geneSets[[i]]] = weights[[i]]
    curr_weight = curr_weight[geneRanking]
    index_score = (abs(sample_R[geneRanking]) * indicatorFunInsideGeneSet)^curr_weight
  }
  
  stepCDFinGeneSet = cumsum(index_score) / sum(index_score)
  stepCDFoutGeneSet = cumsum(!indicatorFunInsideGeneSet) /
    sum(!indicatorFunInsideGeneSet)
  walkStat = stepCDFinGeneSet - stepCDFoutGeneSet             
  
  sum(walkStat) 
}



ssgsea_noWeight = function(X, 
                           gset.idx.list, 
                           alpha, 
                           parallel.sz,
                           parallel.type, 
                           normalization, 
                           verbose) {
  n = length(gset.idx.list)
  cells = colnames(X)
  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nSets", n, envir=globalenv())
    assign("iSets", 0, envir=globalenv())
  }
  
  mclapp = get('mclapply', envir=getNamespace('parallel'))
  options(mc.cores=parallel.sz)
  R = mclapp(1:ncol(X), 
	               function(i){
				      as.integer(rank(X[,i]))
				   })
  R = do.call("cbind", R)
  rm(X)
  gc()
  O = mclapp(1:ncol(R), 
	               function(i){
				      order(R[,i],decreasing=TRUE)
				   })	    
  O = do.call("cbind", O)
  rm(mclapp)
  gc()
  #R = apply(X, 2, function(x) as.integer(rank(x)))    
  #O = apply(R, 2, order, decreasing=TRUE)

	haveParallel = .isPackageLoaded("parallel")
	haveSnow = .isPackageLoaded("snow")
	
	if (parallel.sz > 1 || haveParallel) {
	  if (!haveParallel && !haveSnow) {
	    stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
	  }
	  
	  if (haveSnow) {  ## use snow
	    ## copying ShortRead's strategy, the calls to the 'get()' are
	    ## employed to quieten R CMD check, and for no other reason
	    makeCl = get("makeCluster", mode="function")
	    parSapp = get("parSapply", mode="function")
	    #clEvalQ = get("clusterEvalQ", mode="function")
	    stopCl = get("stopCluster", mode="function")
	    
	    if (verbose)
	      cat("Allocating cluster\n")
	    cl = makeCl(parallel.sz, type = parallel.type) 
	    #clEvalQ(cl, library(GSVA))
	    if (verbose) {
	      cat("Estimating enrichment scores in parallel\n")
	    }
	    #print('using snow')
	    m = t(parSapp(cl, gset.idx.list, 
	                  ks_matrix_ssgsea,    
	                  R, 
	                  O,
	                  tau=alpha))
	    if(verbose)
	      cat("Cleaning up\n")
	    stopCl(cl)
	    
	  } else if (haveParallel) {             ## use parallel
	    
	    mclapp = get('mclapply', envir=getNamespace('parallel'))
	    detCor = get('detectCores', envir=getNamespace('parallel'))
	    nCores = detCor()
	    options(mc.cores=nCores)
	    if (parallel.sz > 0 && parallel.sz < nCores)
	    {
	      options(mc.cores=parallel.sz)
	      #print(parallel.sz)
	    }
	    #print('using parallel')
	    pb = NULL
	    if (verbose){
	      cat("Using parallel with", getOption("mc.cores"), "cores\n")
	      assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
	      assign("nGeneSets", ceiling(length(gset.idx.list) / getOption("mc.cores")), envir=globalenv())
	      assign("iGeneSet", 0, envir=globalenv())
	    }
	    
	    m = mclapp(gset.idx.list, 
	               ks_matrix_ssgsea,
	               R,
	               O,
	               tau=alpha)
	    
	    m = do.call("rbind", m)
	    #colnames(m) = colnames(X)
	    rm(mclapp,detCor,nCores,R,O)
		gc()
	    if (verbose) {
	      close(get("progressBar", envir=globalenv()))
	    }
	  } else
	    stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
	  
	}else {
	  if (verbose) {
	    cat("Estimating enrichment scores\n")
	  }
	  pb = NULL
	  if (verbose){
	    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
	    assign("nGeneSets", length(gset.idx.list), envir=globalenv())
	    assign("iGeneSet", 0, envir=globalenv())
	  }
	  #print('no paralle')
	  m = t(sapply(gset.idx.list, 
	               ks_matrix_ssgsea, 
	               R, 
	               O,
	               tau=alpha))
	  
	  if (verbose) {
	    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
	    close(get("progressBar", envir=globalenv()))
	  }
	}
	
	if (normalization) {
	  ## normalize enrichment scores by using the entire data set, as indicated
	  ## by Barbie et al., 2009, online methods, pg. 2
	  print("normlization ssgsea.....................")
	  score_range = range(m)[2] - range(m)[1]
	  m = m/score_range
      print(dim(m))
	  #m <- apply(m, 2, function(x, m) x / (range(m)[2] - range(m)[1]), m)
	}
	

  if (length(gset.idx.list) == 1)
    m = matrix(m, nrow=1)

  #rownames(m) = names(gset.idx.list)
  colnames(m) = cells

  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }

  m
}

ssgsea_Weight = function(X, 
                         gset.idx.list, 
                         weights,
                         alpha, 
                         parallel.sz,
                         parallel.type, 
                         normalization, 
                         verbose) {
  
  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nSamples", n, envir=globalenv())
    assign("iSample", 0, envir=globalenv())
  }
  
  ks_weight_ssgsea = function(i, expr, geneSets, weights){
    expr[geneSets[[i]],] = expr[geneSets[[i]],]*weights[[i]]
    R = apply(expr, 2, function(x) as.integer(rank(x))) 
    O = apply(R, 2, order, decreasing=TRUE) 
    ks_matrix_ssgsea(geneSets[[i]], R, O, alpha)
  }
  
  haveParallel = .isPackageLoaded("parallel")
  haveSnow = .isPackageLoaded("snow")
  
  if (parallel.sz > 1 || haveParallel) {
    if (!haveParallel && !haveSnow) {
      stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
    }
    
    if (haveSnow) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl = get("makeCluster", mode="function")
      parSapp = get("parSapply", mode="function")
      #clEvalQ = get("clusterEvalQ", mode="function")
      stopCl = get("stopCluster", mode="function")
      
      if (verbose)
        cat("Allocating cluster\n")
      cl = makeCl(parallel.sz, type = parallel.type) 
      #clEvalQ(cl, library(GSVA))
      if (verbose) {
        cat("Estimating enrichment scores in parallel\n")
      }
      
      m = t(parSapp(cl, 1:length(gset.idx.list), 
                    ks_weight_ssgsea,    
                    X, 
                    gset.idx.list,
                    weights))
      if(verbose)
        cat("Cleaning up\n")
      stopCl(cl)
      
    } else if (haveParallel) {             ## use parallel
      
      mclapp = get('mclapply', envir=getNamespace('parallel'))
      detCor = get('detectCores', envir=getNamespace('parallel'))
      nCores = detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      
      pb = NULL
      if (verbose){
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
        assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
        assign("nGeneSets", ceiling(length(gset.idx.list) / getOption("mc.cores")), envir=globalenv())
        assign("iGeneSet", 0, envir=globalenv())
      }
      
      m = mclapp(1:length(gset.idx.list), 
                 ks_weight_ssgsea,
                 X,
                 gset.idx.list,
                 weights)
      
      m = do.call("rbind", m)
      #colnames(m) = colnames(expr)
      
  if (verbose) {
        close(get("progressBar", envir=globalenv()))
      }
    } else
      stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
    
  } else {
    if (verbose) {
      cat("Estimating enrichment scores\n")
    }
    pb = NULL
    if (verbose){
      assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
      assign("nGeneSets", length(gset.idx.list), envir=globalenv())
      assign("iGeneSet", 0, envir=globalenv())
    }
    
    m = t(sapply(1:length(gset.idx.list), 
                 ks_weight_ssgsea, 
                 X, 
                 gset.idx.list,
                 weights))
    
    if (verbose) {
      setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
      close(get("progressBar", envir=globalenv()))
    }
  }
  
  if (normalization) {
    ## normalize enrichment scores by using the entire data set, as indicated
    ## by Barbie et al., 2009, online methods, pg. 2
    m <- apply(m, 2, function(x, m) x / (range(m)[2] - range(m)[1]), m)
  }
  
  if (length(gset.idx.list) == 1)
    m = matrix(m, nrow=1)
  
  #rownames(m) = names(gset.idx.list)
  colnames(m) = colnames(X)
  
  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }
  
  m
}



combinez_noWeight = function(i, Z, geneSets){
  if(length(geneSets[[i]]) == 1){
    Z[geneSets[[i]],]
  }else{
    colSums(Z[geneSets[[i]],]) / sqrt(length(geneSets[[i]]))
  }
}

combinez_Weight = function(i, Z, geneSets, weights){
  if(length(geneSets[[i]]) == 1){
    Z[geneSets[[i]],] * weights[[i]]
  }else{
    colSums(Z[geneSets[[i]],] * weights[[i]]) / sqrt(length(geneSets[[i]]))
  }
}


zscore = function(X, geneSets, weights, parallel.sz, parallel.type, verbose=F) {

  p = nrow(X)
  n = ncol(X)

  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nSamples", n, envir=globalenv())
    assign("iSample", 0, envir=globalenv())
  }

  Z = t(apply(X, 1, function(x) (x-mean(x))/sd(x)))
  #print(dim(Z))

	haveParallel = .isPackageLoaded("parallel")
	haveSnow = .isPackageLoaded("snow")
	
  cl = makeCl = parSapp = stopCl = mclapp = detCor = nCores = NA
	if (parallel.sz > 1 || haveParallel) {
		if (!haveParallel && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
		}

    if (!haveParallel) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl = get("makeCluster", mode="function")
      parSapp = get("parSapply", mode="function")
      stopCl = get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl = makeCl(parallel.sz, type = parallel.type) 
    } else {             ## use parallel

      mclapp = get('mclapply', envir=getNamespace('parallel'))
      detCor = get('detectCores', envir=getNamespace('parallel'))
      nCores = detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      if (verbose)
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
    }
	}
  
  if (parallel.sz == 1 || (is.na(cl) && !haveParallel)){
    if(class(weights) == 'character'){
      es = sapply(1:length(geneSets), combinez_noWeight, Z, geneSets)
    }else{
      es = sapply(1:length(geneSets), combinez_Weight, Z, geneSets, weights)
    }
    #print('no paralle')
  }else {
    if (is.na(cl)){
      print('paralle no cl')
      if(class(weights) == 'character'){
        es = mclapp(1:length(geneSets), combinez_noWeight, Z, geneSets)
      }else{
        es = mclapp(1:length(geneSets), combinez_Weight, Z, geneSets, weights)
      }
      es = do.call("rbind", es)
    }
    else
      if(class(weights) == 'character'){
        es = t(parSapp(cl, 1:length(geneSets), combinez_noWeight, Z, geneSets))
      }else{
        es = t(parSapp(cl, 1:length(geneSets), combinez_Weight, Z, geneSets, weights))
      }
  }


  if (length(geneSets) == 1)
    es = matrix(es, nrow=1)

  rownames(es) = names(geneSets)
  colnames(es) = colnames(X)

  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }

  if (!is.na(cl))
    stopCl(cl)

  es
}

rightsingularsvdvectorgset = function(i, Z, geneSets, weights) {
  if(class(weights) == 'character'){
    s = svd(Z[geneSets[[i]], ])
  }else{
    s = svd(Z[geneSets[[i]], ] * weights[[i]])
  }
  #print(length(s$v[, 1]))
  return(s$v[, 1])
}

rightsingularsvdvectorgset_noWeight = function(i, Z, geneSets){
  s = svd(Z[geneSets[[i]], ])
  return(s$v[, 1])
}


plage = function(X, geneSets, weights, parallel.sz, parallel.type, verbose=F) {

  p = nrow(X)
  n = ncol(X)
  
  print("dim of X")
  print(dim(X))

  if (verbose) {
    assign("progressBar", txtProgressBar(style=3), envir=globalenv()) ## show progress if verbose=TRUE
    assign("nGeneSets", length(geneSets), envir=globalenv())
    assign("iGeneSet", 0, envir=globalenv())
  }
  
  #mclapp = get('mclapply', envir=getNamespace('parallel'))
  #options(mc.cores=parallel.sz)
  #Z = mclapp(1:row(X), 
  #	               function(i){
  #				      (X[i,]-mean(X[i,]))/sd(X[i,])
  #				   })
  #Z = t(do.call("cbind", Z))
  
  Z = t(apply(X, 1, function(x) (x-mean(x))/sd(x)))
  
  print("dim of Z")
  print(dim(Z))

  #Z = t(apply(X, 1, function(x) (x-mean(x))/sd(x)))
  #print(dim(Z))
	haveParallel = .isPackageLoaded("parallel")
	haveSnow = .isPackageLoaded("snow")
	
  ## the masterDescriptor() calls are disabled since they are not available in windows
  ## they would help to report progress by just one of the processors. now all processors
  ## will reporting progress. while this might not be the right way to report progress in
  ## parallel it should not affect a correct execution and progress should be more or less
  ## being reported to some extent.
  cl = makeCl = parSapp = stopCl = mclapp = detCor = nCores = NA ## masterDesc = NA
	if(parallel.sz > 1 || haveParallel) {
		if(!haveParallel && !haveSnow) {
			stop("In order to run calculations in parallel either the 'snow', or the 'parallel' library, should be loaded first")
		}

    if (!haveParallel) {  ## use snow
      ## copying ShortRead's strategy, the calls to the 'get()' are
      ## employed to quieten R CMD check, and for no other reason
      makeCl = get("makeCluster", mode="function")
      parSapp = get("parSapply", mode="function")
      stopCl = get("stopCluster", mode="function")

      if (verbose)
        cat("Allocating cluster\n")
		  cl = makeCl(parallel.sz, type = parallel.type) 
    } else {             ## use parallel

      mclapp = get('mclapply', envir=getNamespace('parallel'))
      detCor = get('detectCores', envir=getNamespace('parallel'))
      ## masterDesc = get('masterDescriptor', envir=getNamespace('parallel'))
      nCores = detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)
      if (verbose)
        cat("Using parallel with", getOption("mc.cores"), "cores\n")
    }
  }

  if (parallel.sz == 1 || (is.na(cl) && !haveParallel)){
    es = sapply(1:length(geneSets), rightsingularsvdvectorgset_noWeight, Z, geneSets)
    #print('no paralle')
    #print(dim(es))
  }
  else {
    if (is.na(cl)) {
      ## firstproc = mclapp(as.list(1:(options("mc.cores")$mc.cores)), function(x) masterDesc())[[1]]
      es = mclapp(1:length(geneSets),
                  rightsingularsvdvectorgset_noWeight, 
                  Z, geneSets) ##, firstproc)
      #print(length(es))
      es = do.call("rbind", es)
      #print('paralle no cl')
      #print(es)
    } else {
      if (verbose)
        message("Progress reporting for plage with a snow cluster not yet implemented")

      es = parSapp(cl, 1:length(geneSets), rightsingularsvdvectorgset_noWeight, Z, geneSets)
      es = do.call("rbind", es)
    }
  }

  if (length(geneSets) == 1)
    es = matrix(es, nrow=1)
	
  print("dim of es")
  print(dim(es))

  rownames(es) = names(geneSets)
  colnames(es) = colnames(X)

  if (verbose) {
    setTxtProgressBar(get("progressBar", envir=globalenv()), 1)
    close(get("progressBar", envir=globalenv()))
  }

  if (!is.na(cl))
    stopCl(cl)

  es
}

setGeneric("filterGeneSets", function(gSets, ...) standardGeneric("filterGeneSets"))

setMethod("filterGeneSets", signature(gSets="list"),
          function(gSets, min.sz=1, max.sz=Inf) {
	gSetsLen = sapply(gSets,length)
	return (gSets[gSetsLen >= min.sz & gSetsLen <= max.sz])	
})





setGeneric("computeGeneSetsOverlap", function(gSets, uniqGenes=unique(unlist(gSets, use.names=FALSE)), ...) standardGeneric("computeGeneSetsOverlap"))

setMethod("computeGeneSetsOverlap", signature(gSets="list", uniqGenes="character"),
          function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  totalGenes = length(uniqGenes)

  ## map to the features requested
  gSets = lapply(gSets, function(x, y) as.vector(na.omit(match(x, y))), uniqGenes)

  lenGsets = sapply(gSets, length)
  totalGsets = length(gSets)

  gSetsMembershipMatrix = matrix(0, nrow=totalGenes, ncol=totalGsets,
                                  dimnames=list(uniqGenes, names(gSets)))
  members = cbind(unlist(gSets, use.names=FALSE), rep(1:totalGsets, times=lenGsets))
  gSetsMembershipMatrix[members] = 1

  .computeGeneSetsOverlap(gSetsMembershipMatrix, min.sz, max.sz)
})




.computeGeneSetsOverlap = function(gSetsMembershipMatrix, min.sz=1, max.sz=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  lenGsets = colSums(gSetsMembershipMatrix)

  szFilterMask = lenGsets >= max(1, min.sz) & lenGsets <= max.sz
  if (!any(szFilterMask))
    stop("No gene set meets the minimum and maximum size filter\n")

  gSetsMembershipMatrix = gSetsMembershipMatrix[, szFilterMask]
  lenGsets = lenGsets[szFilterMask]

  totalGsets = ncol(gSetsMembershipMatrix)

  M = t(gSetsMembershipMatrix) %*% gSetsMembershipMatrix

  M1 = matrix(lenGsets, nrow=totalGsets, ncol=totalGsets,
               dimnames=list(colnames(gSetsMembershipMatrix), colnames(gSetsMembershipMatrix)))
  M2 = t(M1)
  M.min = matrix(0, nrow=totalGsets, ncol=totalGsets)
  M.min[M1 < M2] = M1[M1 < M2]
  M.min[M2 <= M1] = M2[M2 <= M1]
  overlapMatrix = M / M.min

  return (overlapMatrix)
}

## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded = function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
}

##
## ARE THESE FUNCTIONS STILL NECESSARY ?????
##

##a = replicate(1000, compute.null.enrichment(10000,50,make.plot=F))

compute.null.enrichment = function(n.genes, n.geneset, make.plot=FALSE){
	ranks = (n.genes/2) - rev(1:n.genes)
	#null.gset.idxs = seq(1, n.genes, by=round(n.genes / n.geneset))
	null.gset.idxs = sample(n.genes, n.geneset)
	null.es = ks_test_Rcode(ranks, null.gset.idxs,make.plot=make.plot)
	return (null.es)
}


load.gmt.data = function(gmt.file.path){
	tmp = readLines(gmt.file.path)
	gsets = list()
	for(i in 1:length(tmp)){
		t = strsplit(tmp[i],'\t')[[1]]
		gsets[[t[1]]] = t[3:length(t)]
	}
	return (gsets)
}

compute.gset.overlap.score = function(gset.idxs){
	n = length(gset.idxs)
	mx.idx = max(unlist(gset.idxs, use.names=F))
	l = c(sapply(gset.idxs, length))
	
	gset.M = matrix(0, nrow=mx.idx, ncol=n)
	for(i in 1:n){
		gset.M[gset.idxs[[i]],i] = 1
	}
	M = t(gset.M) %*% gset.M
	
	M1 = matrix(l, nrow=n, ncol=n)
	M2 = t(M1)
	M.min = matrix(0, nrow=n, ncol=n)
	M.min[M1 < M2] = M1[M1 < M2]
	M.min[M2 <= M1] = M2[M2 <= M1]
	M.score = M / M.min
	return (M.score)
}

