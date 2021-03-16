start_call = Sys.time()

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scImpute"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("dplyr"))



#########  xiong
##src_dir = '/home/gfj/scTPA_2/src/'
#pythonFile_folder = '/home/gfj/scTPA_2/src/'
#pathway_dir = '/home/gfj/scTPA_2/data/pathway/'

suppressPackageStartupMessages(library("this.path"))
dirCurrent=this.dir()
src_dir = file.path(dirCurrent,'')
pythonFile_folder = file.path(dirCurrent,'')
pathway_dir = file.path(dirCurrent,'../data/pathway/')
template_dir = file.path(dirCurrent,'../data/result_template/')
######### xiong end


source(file.path(src_dir,'visualization.R'))

source(file.path(src_dir,'PAS.R'))
source(file.path(src_dir,'function.R'))
print("function success")
source(file.path(src_dir,'markers.R'))
print("marker success")
source(file.path(src_dir,'clustering.R'))
print("clustering success")
## add by zyr 2020/0727
source(file.path(src_dir,'PAS_ex.R'))
print("pas_ex success")

options(digits=5)

option_list = list(
  make_option(c("-f", "--file"),
              type = "character",
              default = NULL,
              dest = 'file',
              help = "gene expression matrix, one row represents a gene, and one column is a cell",
              metavar = "file"),
  
              make_option(c("--cellType"),
              type = "character",
              default = NULL,
              dest = 'cellType',
              help = "cell type file. First column is cell name (same as the colnames of gene expression profile), second column is cell type. No header names.[default= %default]",
              metavar = "cellType"),
  
              make_option(c("--idType"),
              type = "character",
              default = 'symbol',
              dest = 'idType',
              help = "表达谱矩阵中基因id类型[default= %default]",
              metavar = "idType"),
  
              make_option(c("--normalize"),
              type = "character",
              default = 'none',
              dest = 'normalize_method',
              help = "methods used for normalization. 'log', 'CLR', 'RC' or 'scran'[default= %default]",
              metavar = "normalize_method"),
  
              make_option(c("--min_cells"),
              type = "integer",
              default = 3,
              dest = 'min_cells',
              help = "genes must be detected within a minimum number of cells. Used for filtering genes[default= %default]",
              metavar = "min_cells"),
  
              make_option(c("--min_features"),
              type = "integer",
              default = 200,
              dest = 'min_features',
              help = "cells must have at least the minimum number of genes. Used for filtering cells[default= %default]",
              metavar = "min_features"),
  
              make_option(c("--species"),
              type = "character",
              default = "homo",
              dest = 'species',
              help = "species. 'homo' or 'mus'[default= %default]",
              metavar = "species"),

              make_option(c("--imputation"),
              type = "logical",
              default = FALSE,
              dest = 'imputation',
              help = "Imputation method. 'scImpute' or 'none'[default= %default]",
              metavar = "imputation"),
  
              make_option(c("--data_type"),
              type = "character",
              default = 'TPM',
              dest = 'data_type',
              help = "data type of gene expression profile，'TPM' or 'count'[default= %default]",
              metavar = "file"),
  
              make_option(c("--pathway_database"),
              type = "character",
              default = "kegg",
              dest = 'pathway_database',
              help = "pathway database, details see https://github.com/sulab-wmu/scTPA#details [default= %default]",
              metavar = "pathway_database"),
  
              make_option(c("--topo"),
              type = "logical",
              default = FALSE,
              dest = 'topo',
              help = "是否使用通路拓扑结构信息[default= %default]",
              metavar = "topo"),

              make_option(c("--user_pathway"),
              type = "character",
              default = NULL,
              dest = "user_pathway",
              help = "user defined pathway file, in gmt format[default = %default]",
              metavar = "user_pathway"),
  
              make_option(c("--pas_method"),
              type = "character",
              default = "ssgsea",
              dest = 'pas_method',
              help = "method for calculating PAS (pathway activation signatures). 'gsva', 'ssgsea', 'zscore' or 'plage'. [default= %default]",
              metavar = "pas_method"),
  
              make_option(c("--add_weight_type"),
              type = "integer",
              default = 1,
              dest = 'add_weight_type',
              help = "添加通路权重的方式，1：权重用于排序，2：权重用于ks的指数[default= %default]",
              metavar = "add_weight_type"),
  
              make_option(c("--para_size"),
              type = "integer",
              default = 3,
              dest = 'para_size',
              help = "number of kernels used for parallel computation. [default= %default]",
              metavar = "para_size"),
  
              make_option(c("--cluster_method"),
              type = "character",
              default = "seurat",
              dest = 'cluster_method',
              help = "clustering method. 'seurat', 'hclust', 'simlr', 'kmedoids', 'kmeans' or 'dbscan'. [default= %default]",
              metavar = "cluster_method"),
  
              make_option(c("--seurat_dims"),
              type = "integer",
              default = 8,
              dest = 'seurat_dims',
              help = "dimensions used in Seurat FindNeighbors clustering. [default= %default]",
              metavar = "seurat_dims"),
  
              make_option(c("--seurat_resolution"),
              type = "double",
              default = 0.5,
              dest = 'seurat_resolution',
              help = "resolution used for Seurat clustering. [default= %default]",
              metavar = "seurat_resolution"),
  
              make_option(c("--k_cluster"),
              type = "integer",
              default = 5,
              dest = 'k_cluster',
              help = "number of clusters, used for clustering method except Seurat and dbscan. [default= %default]",
              metavar = "k_cluster"),
  
              make_option(c("--min_pts"),
              type = "integer",
              default = 3,
              dest = 'min_pts',
              help = "parameter in DBSCAN. [default= %default]",
              metavar = "min_pts"),

              make_option(c("--dims"),
              type = "integer",
              default = 20,
              dest = 'dims',
              help = "number of PCA dimensions used for TSNE or UMAP. [default= %default]"),
  
              make_option(c("--marker_method"),
              type = "character",
              default = "wilcox",
              dest = "marker_method",
              help = "method for finding siginificant markers[default= %default]",
              metavar="find_maker_method"),
  
              make_option(c("--logFC_thre"),
              type = "double",
              default = 0.25,
              dest = "logFC_thre",
              help = "threshold of logFC[default= %default]",
              metavar="threshold_logFC"),

              make_option(c("--min_pct"),
              type = "double",
              default = 0.1,
              dest = "min_pct",
              help = "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.[default= %default]",
              metavar="min_pct"),

              make_option(c("--shown_markers"),
              type = "integer",
              default = 3,
              dest = "shown_markers",
              help = "展示的markers的个数[default= %default]",
              metavar="find_maker_method"),
  
              make_option(c("-o","--out_dir"),
              type = "character",
              default = NULL,
              dest = "out_dir",
              help = "output folder[default= %default]",
              metavar="out_dir")

              #make_option(c("--jobID"),
              #type = "character",
              #default = NULL,
              #dest = "jobID",
              #help = "job ID[default= %default]",
              #metavar="jobID"),

              #make_option(c("--email"),
              #type = "character",
              #default = 'unacquainted',
              #dest = "email",
              #help = "email [default= %default]",
              #metavar="email")
  
              #make_option(c("--jobProces"),
              #type = "logical",
              #default = TRUE,
              #dest = "jobProces",
              #help = "用户邮箱[default= %default]",
              #metavar="jobProces")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("######opt:")
#print(opt)



msigdb = c('h.all','other', 'c2.cgp', 'c4.cgn', 'c4.cm', 'c5.bp', 'c5.cc', 'c5.mf', 'c6.all', 'c7.all', 'acsn2')

pipeline_cal = function(expr_path, 
                        cell_type_path,
                        idType='symbol',
                        data_type = c('TPM','count'),
                        normalize_method = c('log', 'CLR', 'RC', 'none','scran','sctrancform'),
                        species = c('homo', 'mus'),
                        min_cells=3,
                        min_features=200,
                        imputation=FALSE,
                        pathway_database,
                        pas_method=c('ssgsea','gsva','plage','zscore'),
                        topo=TRUE,
                        user_pathway = NULL,
                        add_weight_type = 1,
                        para_size=3,
                        cluster_method,
                        shown_method=c('tsne','umap'),
                        seurat_dims = 8,
                        seurat_resolution = 0.5,
                        k_cluster = 5,
                        min_pts = 3,
                        dims=20,   #用TSNE或UNAP展示聚类结果时，使用的pca是维度
                        marker_method = c('wilcox','bimod','MAST','t','LR'),
                        logFC_thre = 0.25,
                        min.pct = 0.1,
                        shown_markers = 3,
                        jobProces = TRUE,
                        out_dir,
                        jobID,
                        email
)
{
  print("start scEPA")
  min_cells = max(2, min_cells)
  
  ### check species
  if(species == 'mus'){
    imputation = FALSE
  }
  
  start_time = Sys.time()
  #if(cell_type_path == 'none'){
  #  cell_type = 'none'
  #}else{
  #  cell_type = read.csv(cell_type_path, row.names = 1)
  #}
  
  if(jobProces){
    
    #output_dir=paste0(out_dir,"/",jobID)
    endChar_out_dir = substr(out_dir,nchar(out_dir),nchar(out_dir));
    if(endChar_out_dir!="/" && endChar_out_dir!="\\"){
      out_dir = paste0(out_dir,"/");
    }

    if(!dir.exists(out_dir)){
      #dir.create(out_dir, recursive = TRUE)
      shell_cmd_cp=paste0("mkdir -p ",  out_dir)
      print(shell_cmd_cp)
      grep_out = system(shell_cmd_cp, intern = TRUE)
    }

    shell_cmd_cp=paste0("cp -rf ",  template_dir, "/app/* ", out_dir)
    print(shell_cmd_cp)
    grep_out = system(shell_cmd_cp, intern = TRUE)

    out_dir = paste0(out_dir,"/data/content/");
    shell_cmd_cp=paste0("cp ",  template_dir, "/result_frame.html ", out_dir,"/processing.html")
    print(shell_cmd_cp)
    grep_out = system(shell_cmd_cp, intern = TRUE)
    shell_cmd_cp=paste0("cp ", template_dir,"/result_frame.html ", out_dir,"/result.html")
    print(shell_cmd_cp)
    grep_out = system(shell_cmd_cp, intern = TRUE)

    shell_cmd = paste('python3', paste0(pythonFile_folder,'jobProcess_html.py'), out_dir, 1)
    #print(shell_cmd)
    grep_out = system(shell_cmd, intern = TRUE)
  }
  print("start loading data")
  ##################################################
  #####   step 1. load expression profile   ########
  ##################################################
  
  if(!file_test("-d", out_dir)){
    dir.create(out_dir,recursive=TRUE)
  }
  #log_file_path = file.path(out_dir, 'log')
  
  
  tryCatch({
   #print(cell_type_path)
   ret = load_expr(expr_path,
                   cell_type_path,
                   imputation,
                   data_type,
                   out_dir,
                   para_size,
                   normalize_method = normalize_method,
                   min_cells = min_cells,
                   min_features = min_features,
                   pas_method = pas_method,
                   jobProces = jobProces)
   print("load_expr success")
  expr = as.matrix(ret[[1]])
  cell_type = ret[[2]]
  rm(ret, min_cells, min_features, normalize_method, data_type)
  gc()
  rownames(expr) = toupper(rownames(expr))
  expr = conveEnsembl(expr, species)
  
  },error = function(e){
    shell_cmd = paste('python3', paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 8)
    grep_out = system(shell_cmd, intern = TRUE)
    print(e)
    q()
  })
  print("expr loading success")
   
  ##################################################
  ########    step 2. select pathways     ##########
  ##################################################
  tryCatch({
  ## edit by zyr 2020 0727
  pathway_path = file.path(pathway_dir, species, paste0(pathway_database,'.gmt'))
  path_list = load_gmt_data(pathway_path)
    
  
    gc()
    #if(topo){
    #   weight_list = path_info$weights
    #}else{
    weight_list = 'none'
    #}
  }, error = function(e){
    print("loading pathway error")
    shell_cmd = paste('python3', paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 9)
    grep_out = system(shell_cmd, intern = TRUE)
    print(file.path(out_dir, 'status'))
    write.table('error', file.path(out_dir, 'status'), quote = F, row.names=F,col.names=F)
    #write.table("Some wrong catched when selecting pathway database, please contact author.",
    #                  file.path(out_dir, 'error_unexcept'), quote=F, col.names=F, row.names=F)
    print(e)
    q()
  })
  print('pathways loading success')
  
  ##################################################
  ####  step-optional. user defined pathways   #####
  ##################################################
  tryCatch({
  if(! user_pathway == 'NULL'){
      user_pathway = load_gmt_data(user_pathway)
      aa = gsub('-|,|/| |\\(|\\)','_',names(user_pathway))
      aa = gsub('_+','_',aa)
      aa = gsub('_\\b','',aa)
      aa = gsub(',','',aa)
      names(user_pathway) = paste0('user_', names(user_pathway))
      path_list = c(path_list, user_pathway)
      weight_list = 'none'
  }}, error = function(e){
    shell_cmd = paste('python3', paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 10)
    grep_out = system(shell_cmd, intern = TRUE)
    write.table('error', file.path(out_dir, 'status'), row.names=F, col.names=F, quote=F)
    #write.table("Some wrong catched when loading user defined pathway, please check the gmt file.",
    #                      file.path(out_dir, 'error_unexcept'),  quote=F, col.names=F, row.names=F)
    q()
  })
  
  print("pathway preparation successful")
  print(Sys.time()-start_time)
  
  gc()
  ###调用python写html
  if(jobProces){
    shell_cmd = paste('python3', 
                      paste0(pythonFile_folder,'jobProcess_html.py'),
                      out_dir, 2)
    grep_out = system(shell_cmd, intern = TRUE)
  }

  
  ##################################################
  ########    step 3. calculate PAS       ##########
  ##################################################
  print(expr[1:3,1:3])
  print(class(expr))
  
  tryCatch({
  cluster_mat = switch(pas_method,
                       AUCell = cal_AUCell(expr[,],
                                           path_list,
                                           para_size),
                       pagoda2 = cal_pagoda2(expr[,],
                                             path_list,
                                             n_cores=para_size),
                       Vision = cal_vision(expr[,],
                                           pathway_path,
                                           para_size),
                       gsva = PAS(expr, 
                                  path_list, 
                                  method = 'gsva', 
                                  parallel.sz = para_size,
                                  weights = weight_list,
                                  verbose=TRUE),
                       ssgsea = PAS(expr, 
                                  path_list, 
                                  method = 'ssgsea', 
                                  parallel.sz = para_size,
                                  weights = weight_list,
                                  verbose=TRUE),
                       zscore = PAS(expr, 
                                  path_list, 
                                  method = 'zscore', 
                                  parallel.sz = para_size,
                                  weights = weight_list,
                                  verbose=TRUE),
                       plage = PAS(expr, 
                                  path_list, 
                                  method = 'plage', 
                                  parallel.sz = para_size,
                                  weights = weight_list,
                                  verbose=TRUE)
)

  }, error = function(e){
    rm(expr)
    gc()
    shell_cmd = paste('python3',paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 11)
    grep_out = system(shell_cmd, intern = TRUE)
    write.table('error', file.path(out_dir, 'status'), row.names=F, col.names=F, quote=F)
    print(e)
    q()
  })
  print(class(score))
  print(score)
  if(class(score) == 'character' && score == 'error'){
    rm(expr)
    gc()
    shell_cmd = paste('python3',paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 11)
    grep_out = system(shell_cmd, intern = TRUE)
    write.table('error', file.path(out_dir, 'status'), row.names=F, col.names=F, quote=F)
    print(e)
    q()
  }

  pas_dir = file.path(out_dir,'pas.csv')
  write.csv(cluster_mat, file = pas_dir, quote=F)
  rm(pas_dir)
  gc()
  
  print("calculating pathway activity score successful.")
  print(Sys.time()-start_time)
  
  
  
  ###调用python写html
  if(jobProces){
    shell_cmd = paste('python3',
                      paste0(pythonFile_folder,'jobProcess_html.py'),
                      out_dir, 3)
    grep_out = system(shell_cmd, intern = TRUE)
  }
  
  
  ##################################################
  ########    step 4. cluster idents      ##########
  ##################################################
  print('clustering.....................')
  print(class(cell_type))
  if(cell_type == 'none'){
    tryCatch(
      {
        cluster_idents = cal_cluster(cluster_mat, 
                                     method=cluster_method,
                                     seurat_dims = seurat_dims,
                                     seurat_resolution = seurat_resolution,
                                     k_cluster = k_cluster,
                                     min_pts = min_pts)
        cluster_idents = as.numeric(cluster_idents)
      },error = function(e){
        print("there are no clustering")
        
        #### write result html,can not clustering
        rm(cluster_mat,expr)
        gc()
        if(cell_type == 'none'){
          shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'no',
                    'noClust')
        
        }else{
          shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'yes',
                    'noClust')
         }
        #### write job result file
        grep_out = system(shell_cmd, intern = TRUE)
        #### write job processing file
        shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jobProcess_html.py'),
                    out_dir, 5)
        grep_out = system(shell_cmd, intern = TRUE)
        print(e)
        q()
      }
    )
       
    if(min(cluster_idents) == 0){
        cluster_idents = cluster_idents + 1
    }
    
    cluster_idents = factor(cluster_idents,levels = sort(unique(cluster_idents)))
    #cluster_idents = paste0('C',cluster_idents)
    names(cluster_idents) = colnames(cluster_mat)
    #print(cluster_idents[1:10])
    cluster_idents = sort(cluster_idents)
    #print(cluster_idents[1:10])
    cluster_mat = cluster_mat[,names(cluster_idents)]
    cluster_idents = factor(paste0('C',cluster_idents),levels = paste0('C',levels(cluster_idents)))
    names(cluster_idents) = colnames(cluster_mat)
    print(cluster_idents[1:10])
    ##### write cluster_idents
    ######write cell_typeS
    print("wriite cell_type")
    #print(data.frame(cluster_idents))
    
    #pas = pas[,names(cluster_idents)]
  }else{
    print("known cell type")
    cluster_idents = factor(cell_type,levels = sort(unique(cell_type)))
    names(cluster_idents) = colnames(cluster_mat)
    print(cluster_idents[1:10])
    cluster_idents = sort(cluster_idents)
    print(cluster_idents[1:10])
    cluster_mat = cluster_mat[,names(cluster_idents)]
  }  
  
  write.table(data.frame(cluster_idents), 
              file.path(out_dir,'cell_type.csv'),sep=',',
              quote = F, col.names = F)
  print("clustering cells successful")

  
  
  ##################################################
  ########    step 5. visulization        ##########
  ##################################################
 
  cells = c()
  if(ncol(cluster_mat) > 10000){
    ## extract some cells if the mat is too big
    print("extract subset cells")
    pst = 8000/ncol(cluster_mat)
    #print(dim(mat))
    #mat_new = matrix(nrow=nrow(mat))
    for(clust in unique(cluster_idents)){
      #print(clust)
      curr_cells = which(cluster_idents == clust)
      #print(length(curr_cells))
      curr_cells = sample(curr_cells, round(length(curr_cells)*pst))
      cells = c(cells, curr_cells)
      #curr_mat = mat[,curr_cells]
      #print(dim(curr_mat))
      #mat_new = cbind(mat_new, curr_mat)
      #print(dim(mat_new))
    }
    cat("extract subset of original expression profile, dimensions of new matrix is ",length(cells),"\n")
    #mat = mat_new[,2:ncol(mat_new)]
    #cat("extract subset of original expression profile, dimensions of new matrix is ",paste(dim(mat),collapse=','),"\n")
    #rm(mat_new)
    #gc()
  }else{
    cells = 1:ncol(cluster_mat)
  }
  print(cells[1:10])
  print(length(cells))  
  
  options(bitmapType='cairo')
  
  tryCatch(
    {
      plotPAS_heatmap(cluster_mat[,cells], cluster_idents[cells], paste0(out_dir, 'pas_heatmap.png'))
      
    },error= function(e){
      print("cannot  plotting PAS heatmap")
      rm(cluster_mat,cluster_idents,expr)
      gc()
      if(cell_type == 'none'){
          shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'no',
                    'noClust')
        }else{
          shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'yes',
                    'noClust')
         }
        #### write job result file
        grep_out = system(shell_cmd, intern = TRUE)
        #### write job processing file
        shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jobProcess_html.py'),
                    out_dir, 5)
        grep_out = system(shell_cmd, intern = TRUE)
        print(e)
        q()
    }
  )
  
  tryCatch(
    {
      obSeurat_pas = prepa_seuratOb(cluster_mat, cluster_idents, seurat_dims, out_dir, cells)
     },error= function(e){
      print("cannot  prepare Seurat object")
      rm(cluster_mat, expr, cluster_idents)
      gc()
      if(cell_type == 'none'){
          shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'no',
                    'noDime')
        
      }else{
        shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'yes',
                    'noDime')
      }
      #### write job result file
      grep_out = system(shell_cmd, intern = TRUE)
      #### write job processing file
      shell_cmd = paste('python3',
                  paste0(pythonFile_folder,'jobProcess_html.py'),
                  out_dir, 5)
      grep_out = system(shell_cmd, intern = TRUE)
      print(e)
      q()
      
    }
  )
  rm(cluster_mat)
  gc()
  


  
  print("ploting cluster picture successful")
  

  ###调用python写html
  if(jobProces){
    shell_cmd = paste('python3',
                      paste0(pythonFile_folder,'jobProcess_html.py'),
                      out_dir, 4)
    grep_out = system(shell_cmd, intern = TRUE)
  }

  ##################################################
  ########    step 6. finding markers      #########
  ##################################################
  
  
  # if(pas_method == 'gsva'){
  #     logfc_thre = 0.25 
  # }else{
  #     logfc_thre = 0.01
  # }
  
  errMarkers = FALSE
  tryCatch(
   {
     if(user_pathway == 'NULL'){
        pas_markers = parseAllMarkers(obSeurat_pas,
                                      NULL,
                                      method = marker_method,
                                      logfc_thre = logFC_thre,
                                      min.pct=min.pct,
                                      shown_markers = shown_markers,
                                      out_dir = out_dir,
                                      para_size = para_size,
                                      expr = expr,
                                      cell_type = cluster_idents,
                                      path_list = path_list,
                                      cells = cells)
    }else{
        pas_markers = parseAllMarkers(obSeurat_pas,
                                      names(user_pathway),
                                      method = marker_method,
                                      logfc_thre = logFC_thre,
                                      min.pct = min.pct,
                                      shown_markers = shown_markers,
                                      out_dir = out_dir,
                                      para_size = para_size,
                                      expr = expr,
                                      cell_type = cluster_idents,
                                      path_list = path_list,
                                      cells = cells)
  }
   },error = function(e){
       print("error in finding markers")
       print(e)
       errMarkers = TRUE
   }
  )
  rm(expr)
  gc()
  print("from scEPA.R")
  #print(dim(pas_markers))
  if(errMarkers || nrow(pas_markers) == 0){
    print("no diff genes")
    #### write plugin json file
    shell_json = paste('python3', paste0(pythonFile_folder,'make_json.py'),out_dir)
    grep_out = system(shell_json, intern = TRUE)
    #### write result html
    if(cell_type == 'none'){
    shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'no',
                    'no')
    print("there are no cell type")
    }else{
    shell_cmd = paste('python3',
                      paste0(pythonFile_folder,'jonResult_html.py'),
                      out_dir,
                      jobID,
                      email,
                      'yes',
                      'no')
     }
    #### write job result file
    grep_out = system(shell_cmd, intern = TRUE)
    #### write job processing file
    shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jobProcess_html.py'),
                    out_dir, 5)
    grep_out = system(shell_cmd, intern = TRUE)
    
    
    q()
  }
  

  
  #tryCatch({
  #print("write json......")
  shell_json = paste('python3', paste0(pythonFile_folder,'make_json.py'),out_dir)
  print(shell_json)
  grep_out = system(shell_json, intern = TRUE)
  print(grep_out)
  if(cell_type == 'none'){
    shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'no',
                    'yes')
    print("there are no cell type")
  }else{
    shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jonResult_html.py'),
                    out_dir,
                    jobID,
                    email,
                    'yes',
                    'yes')
    print("there are have cell type")
  }
  print(shell_cmd)
  grep_out = system(shell_cmd, intern = TRUE)
   
  shell_cmd = paste('python3',
                    paste0(pythonFile_folder,'jobProcess_html.py'),
                    out_dir, 5)
  print(shell_cmd)
  grep_out = system(shell_cmd, intern = TRUE)
 # },error = function(e){
 #   shell_cmd = paste('python3', paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 14)
 #   write.table('error', file.path(out_dir, 'status'), row.names=F, col.names=F, quote=F)
 # })
}
print(opt$pas_method)

start_time = Sys.time()
pipeline_cal(expr_path = opt$file,
             cell_type = opt$cellType,
             idType = opt$idType,
             normalize_method = opt$normalize_method,
             min_cells = opt$min_cells,
             min_features = opt$min_features,
             data_type = opt$data_type,
             imputation = opt$imputation,
             species = opt$species,
             pathway_database = opt$pathway_database,
             pas_method = opt$pas_method,
             add_weight_type = opt$add_weight_type,
             topo = opt$topo,
             user_pathway = opt$user_pathway,
             para_size = opt$para_size,
             cluster_method = opt$cluster_method,
             shown_method = 'tsne',
             seurat_dims = opt$seurat_dims,
             seurat_resolution = opt$seurat_resolution,
             k_cluster = opt$k_cluster,
             min_pts = opt$min_pts,
             dims = opt$dims,
             marker_method = opt$marker_method,
             shown_markers = opt$shown_markers,
             logFC_thre = opt$logFC_thre,
             min.pct = opt$min_pct,
             out_dir = opt$out_dir,
             jobProces = TRUE,
             jobID = opt$out_dir,
             email = 'unacquainted'
)
print(Sys.time()-start_call)




