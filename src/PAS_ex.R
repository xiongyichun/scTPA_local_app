cal_AUCell = function(counts,
                      gSets,
                      n_cores
                      # number of cores used for parallel
                      ){

  ## matrix could be integer or numeric
  library(AUCell)
  #tryCatch({
print("######cal_AUCell_ex: step1")
    # @must step1: rank genes
    ac_rankings = AUCell::AUCell_buildRankings(counts[,],
                                               # matrix, dgCMatrix, SummarizedExperiment, ExpressionSet
                                               nCores=n_cores,
                                               plotStats=FALSE,
                                               # plot the expression boxplots or histograms
                                               #assayName = NULL,
                                               # slot name of assay containing expression matrix
                                               verbose = T)
    # @must step2: calculate AUC scores
print("######cal_AUCell_ex: step2")

    sc_AUC = AUCell::AUCell_calcAUC(gSets,
                                    ac_rankings,
                                    normAUC=T,
                                    # Wether to normalize the maximum possible AUC to 1
                                    aucMaxRank=ceiling(0.05 * nrow(ac_rankings)),
                                    # the number of genes (maximum ranking) that is used to computation
                                    # default: 5% of pathways; recommend: 1%-20%
                                    verbose = F
    )
print("######cal_AUCell_ex: step3")

    score = AUCell::getAUC(sc_AUC)
    score = na.omit(score)
print("######cal_AUCell_ex: finished")

    return(score)
  #},error = function(e){
   # print(e)
   # return("error")
  #})
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
  library(VISION)

  #tryCatch({
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
    score = na.omit(score)
    return(score)
  #},error = function(e){
  #  print(e)
  #  return("error")
  #})
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
  library(pagoda2)
  nPcs = min(round(ncol(n_counts)/5),5)
  #counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
  print(dim(n_counts))
  print(n_counts[1:5,1:5])
  n_counts = n_counts[!duplicated(rownames(n_counts)),]
  print(dim(n_counts))
  #tryCatch({
    print("#####################PAS_ex, cal_pagoda2, start")
    p2 = Pagoda2$new(n_counts, n.cores = n_cores,log.scale=F)
    #print("11111111111111111111111: PAS_ex")
    p2$adjustVariance(plot=F)
    #print("22222222222222222222222: PAS_ex")
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    #print("333333333333333333333333: PAS_ex")
    path_names = c()
    env = new.env(parent=globalenv())
    invisible(lapply(1:length(gSets),function(i) {
      genes = intersect(gSets[[i]],rownames(n_counts))
      name = names(gSets[i])
      if(length(genes)>3){
        assign(name, genes, envir = env)
        path_names = c(path_names, name)
      }
    }))
    #print("4444444444444444444444444: PAS_ex")
    print(class(env)) 
    print(env)
    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = T,
                                 min.pathway.size = 1)
    #print("5555555555555555555555555: PAS_ex")

    path_names = names(p2@.xData$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(n_counts))
    rownames(score) = path_names
    colnames(score) = colnames(n_counts)
    for(i in 1:length(p2@.xData$misc$pwpca)){
      if(!is.null(p2@.xData$misc$pwpca[[i]]$xp$score)){
        score[i,] = as.numeric(p2@.xData$misc$pwpca[[i]]$xp$scores)
      }
    }
    print(dim(score))
    score = na.omit(score)

    print("########################: PAS_ex, cal_pagoda2, finished")

    return(score)
  #},error = function(e){
  #  print(e)
  #  return("error")
  #})
}