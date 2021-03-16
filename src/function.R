

########### xiong
#data_dir = '/home/gfj/scTPA_2/data/'
#pythonFile_folder = '/home/gfj/scTPA_2/src/'
suppressPackageStartupMessages(library("this.path"))
dirCurrent=this.dir()
data_dir=file.path(dirCurrent,'../data')
pythonFile_folder=file.path(dirCurrent,'/')
print("######function.R")
#print(dirCurrent)
########### xiong end

geneLength_path = file.path(data_dir,'all_symbol_length.txt')
idConver_path = file.path(data_dir, 'ID_convert')


get_geneLis = function(paths,
                       path_lis,
                       expr_genes){
  res = sapply(paths,FUN = function(x){
    genes = path_lis[[x]]
    genes = intersect(genes,expr_genes)
    paste(genes,collapse = ' ')
  })
  as.character(res)
}


conveEnsembl = function(expr, species, thre_pst = 0.1){
  rr = as.character(sapply(rownames(expr), function(x) strsplit(x,split="\\.")[[1]][1]))
  expr = expr[!duplicated(rr),]
  rownames(expr) = rr[!duplicated(rr)]
  print(file.path(idConver_path, paste0(species, '.csv')))
  ids = read.csv(file.path(idConver_path, paste0(species, '.csv')), row.names=1)
  print(head(ids))
  inter_genes = intersect(rownames(expr), rownames(ids))
  if(length(inter_genes)> thre_pst*nrow(expr)){
    print(dim(ids))
    ids = ids[inter_genes,]
    print(length(ids))
    expr = expr[inter_genes,]
    print(head(ids))
    rownames(expr) = as.character(ids)
  }
  expr
}







load_expr = function(expr_file,
                     cell_type_path,
                     imputation,
                     data_type,
                     out_dir,
                     para_size,
                     normalize_method,
                     min_cells,
                     min_features,
                     pas_method,
                     jobProces)
  {
  # load gene expression matrix
  # file format:
  #            rownames: gene symbol; 
  #            colnames: sample (cell) ID; 
  #            sep: \t
  #if(imputation){
  #   print("impute yes")
  #}else{
  #   print('impute false')
  #}
  log_file_path = file.path(out_dir, 'log')
  if(is.null(expr_file)){
    write.table(paste0("Please upload expression matrix file!"), 
                log_file_path, quote=F, col.names=F, row.names=F)
    q()
  }
  #normalize_method = match.arg(normalize_method)
  expr_path_lis = tolower(strsplit(expr_file, split='[.]')[[1]])
  expr_path_type = expr_path_lis[length(expr_path_lis)]

print("##### start loading data!")
  tryCatch({
    if(expr_path_type == 'csv'){
      expr = read.csv(expr_file, stringsAsFactors = F,row.names=NULL,header=T,check.names=FALSE)
      expr = expr[!duplicated(expr[,1]),]
      rownames(expr) = expr[,1]
      expr = expr[,2:ncol(expr)]
    }else if(expr_path_type == 'txt'){
      tryCatch({
        expr = read.table(expr_file,sep='\t',stringsAsFactors = F,row.names=NULL,header=T,check.names=FALSE)
      },error = function(e){
        shell_cmd = paste('python3', paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 16)
        grep_out = system(shell_cmd, intern = TRUE)
        print(e)
        q()
      })
      expr = expr[!duplicated(expr[,1]),]
      rownames(expr) = expr[,1]
      expr = expr[,2:ncol(expr)]
    }else if(expr_path_type == 'rds'){
      expr = readRDS(expr_file)
    }else{
      rm(expr)
      gc()
      shell_cmd = paste('python3', file.path(pythonFile_folder,'jobProcess_html.py'), out_dir, 6)
      grep_out = system(shell_cmd, intern = TRUE)
      print("error: filename must be end in csv or txt! Please check it!")
      write.table('error', file.path(out_dir, 'status'), row.names=F, col.names=F, quote=F)
      #write.table(paste0("filename must be end in csv or txt! Please check it!"),
      #              file.path(out_dir, 'error_filetype'), quote=F, col.names=F, row.names=F)
      q()
    }
  },error = function(e){
    shell_cmd = paste('python3', paste0(pythonFile_folder, 'jobProcess_html.py'), out_dir, 16)
    grep_out = system(shell_cmd, intern = TRUE)
    print(e)
    q()
  })

  print(paste("reading file success, dimension of matrix is ",paste(dim(expr),collapse=',')))
  print(expr[1:3,1:3])
  if(! cell_type_path == 'NULL'){
    print("reading cell type file")
    print(cell_type_path)
    cell_type_path_lis = tolower(strsplit(cell_type_path, split='[.]')[[1]])
    #print(cell_type_path_lis)
    if(cell_type_path_lis[length(cell_type_path_lis)] == 'csv'){
      #print(cell_type_path)
      cell_type = read.csv(cell_type_path, stringsAsFactors = F,row.names=1,header=F,check.names=FALSE)
    }else if(cell_type_path_lis[length(cell_type_path_lis)] == 'txt'){
      cell_type = read.table(cell_type_path,sep='\t',stringsAsFactors = F,row.names=1,header=F,check.names=FALSE)
    }
    
  if(length(setdiff(colnames(expr), rownames(cell_type))) > 0){
      rm(expr)
      gc()
      shell_cmd = paste('python3', file.path(pythonFile_folder,'jobProcess_html.py'), out_dir, 7)
      print(file.path(out_dir, 'status'))
      #write.table('error', file.path(out_dir, 'status'), row.names=F, col.names=F, quote=F)
      grep_out = system(shell_cmd, intern = TRUE)
      print("error: some cells in expression matrix are not included in first line of cell type file! Please check it!")
      #write.table(paste0("some cells in expression matrix are not included in first line of cell type file! Please check it!"), 
      #            file.path(out_dir, 'error_celltypes'), quote=F, col.names=F, row.names=F)
      q()
    }
  }else{
    cell_type = 'none'
  }
  
  expr = as.matrix(expr)
  rownames(expr) = gsub(" +","",rownames(expr))
  print(expr[1:3,1:3])
   
  if (min_features > 0) {
    nfeatures = Matrix::colSums(x = expr > 0)
    #print(nfeatures[1:10])
    expr = expr[, which(x = nfeatures >= min_features)]
    print("after filter cells: ")
    print(dim(expr))
    rm(nfeatures)
    gc()
  }  

 
  if (min_cells > 0) {
    num_cells = Matrix::rowSums(x = expr > 0)
    expr = expr[which(x = num_cells >= min_cells), ]
    print("after filter genes: ")
    print(dim(expr))
    rm(num_cells)
    gc()
  }
  
    
  if(pas_method == 'gsva' && !imputation && class(expr[1,1]) == 'character'){
    warning("character type expression matrix is unuseful when pas_method is gsva, please check it! Casting data type to numeric...")
    genes = rownames(expr)
    expr = apply(expr,2,as.numeric)
    rownames(expr) = genes
  }
  
  #print(expr[1:3,1:3])
  #### 是否使用标准化
  genes = rownames(expr)
  cells = colnames(expr)
  if(normalize_method == 'log'){
    expr = NormalizeData(expr,normalization.method = 'LogNormalize',verbose=0)
  }else if(normalize_method == 'CLR'){
    expr = NormalizeData(expr,normalization.method = 'CLR',verbose=0)
  }else if(normalize_method == 'RC'){
    expr = NormalizeData(expr,normalization.method = 'RC',verbose=0)
  }else if(normalize_method == 'scran'){
    
    if(data_type == 'count'){
    expr = apply(expr,2,function(x) {storage.mode(x) <- 'integer'; x})
    sc = SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = expr)
    )
    clusters = scran::quickCluster(sc)
    sc = scran::computeSumFactors(sc, clusters=clusters)
    sc = scater::normalize(sc)
    }else{
    sc = SingleCellExperiment::SingleCellExperiment(
      assays = list(tpm = expr)
    )
    clusters = scran::quickCluster(sc,assay.type = "tpm")
    sc = scran::computeSumFactors(sc, clusters=clusters, assay.type = "tpm")
    sc = scater::normalize(sc,exprs_values  ='tpm')
    }
    print('sc normalize success')
    expr = sc@assays$data$logcounts

    rm(sc,clusters)
    gc()
  }else if(normalize_method == 'sctransform'){
    expr_ob = Seurat::CreateSeuratObject(counts=expr)
    expr_ob = Seurat::SCTransform(expr_ob,verbose=FALSE)
    expr = as.matrix(expr_ob@assays$SCT@data)
    rm(expr_ob)
    gc()
  }
  
  if(normalize_method != 'sctransform'){
    rownames(expr) = genes
    colnames(expr) = cells
  }

  ##imputation
  
  tryCatch({
    if(cell_type == 'none'){
      if(imputation){
        write.table("Start impute missing data without cell type.", 
                    file.path(out_dir, 'log_imput'), quote=F, col.names=F, row.names=F)
        gene_len = read.table(geneLength_path,stringsAsFactors = F,sep="\t",header=F,row.names=1,check.names=FALSE)
        tmp = intersect(rownames(gene_len),rownames(expr))
        if (length(tmp) != nrow(expr)){
          warning("check the length file")
          print(length(setdiff(rownames(expr),rownames(gene_len))))
          expr = expr[tmp,]
          gene_len = gene_len[tmp,]
          print("before writing")
          print(expr[1:3,1:3])
          write.csv(expr,file = file.path(out_dir, 'expression.csv'),quote=F)
          expr_file = file.path(out_dir, 'expression.csv')
        }else{
          print("write expr files without cell type")
          print("cell_type")
          gene_len = gene_len[rownames(expr), ]
          write.csv(expr,file = file.path(out_dir, 'expression.csv'),quote=F)
          expr_file = file.path(out_dir, 'expression.csv')
        }
        
        scimpute(count_path = expr_file, 
                 type = data_type,
                 infile = "csv", 
                 outfile = "csv",
                 genelen = gene_len,
                 out_dir = out_dir,
                 drop_thre = 0.5, 
                 labeled = F,
                 Kcluster = 8,
                 ncores = para_size)
        expr = read.csv(file.path(out_dir,'scimpute_count.csv'), stringsAsFactors = F,row.names=1,check.names=FALSE)
        rm(gene_len)
        gc()
      }
    }else{
      
      if(imputation){
        print("Start impute missing data with cell type.")
        #      write.table("Start impute missing data with cell type.", 
        #                  file.path(out_dir, 'log_imput'), quote=F, col.names=F, row.names=F)
        gene_len = read.table(geneLength_path,stringsAsFactors = F,sep="\t",header=F,row.names=1,check.names=FALSE)
        tmp = intersect(rownames(gene_len),rownames(expr))
        if (length(tmp) != nrow(expr)){
          warning("check the length file")
          print(length(setdiff(rownames(expr),rownames(gene_len))))
          expr = expr[tmp,]
          gene_len = gene_len[tmp,]
          print("before writing")
          print(expr[1:3,1:3])
          write.csv(expr,file = file.path(out_dir, 'expression.csv'),quote=F)
          expr_file = file.path(out_dir, 'expression.csv')
          
        }else{
          print("write expr files without cell type")
          write.csv(expr,file = file.path(out_dir, 'expression.csv'),quote=F)
          expr_file = file.path(out_dir, 'expression.csv')
          gene_len = gene_len[rownames(expr), ]
        }
        print("known cell_type")
        #print(expr_file)
        #print(out_dir)
        #print(dim(expr))
        #print(dim(cell_type))
        #print(head(cell_type))
        #print(length(as.vector(cell_type[,1])))
        #print(cell_type[1:10])
        #    imput_dir = file.path(out_dir,'imputed')
        scimpute(count_path = expr_file, 
                 type = data_type,
                 infile = "csv", 
                 outfile = "csv",
                 genelen = gene_len,
                 out_dir = out_dir, 
                 drop_thre = 0.5, 
                 labeled = T,
                 labels = as.vector(cell_type[,1]),
                 ncores = para_size)
        expr = read.csv(paste0(out_dir,'scimpute_count.csv'), stringsAsFactors = F,row.names=1,check.names=FALSE)
        #print(expr[1:5,1:5])
        #print(sum(sum(expr)))
        rm(gene_len)
        gc()
      }
    }
  },error = function(e){
    print(e)
    rm(expr)
    gc()
    shell_cmd = paste('python3', file.path(pythonFile_folder,'jobProcess_html.py'), out_dir, 15)
    grep_out = system(shell_cmd, intern = TRUE)
  }) 
 
  
  #rownames(expr) = genes
  colnames(expr) = cells
  #print(expr[1:5,1:5])
  #print(head(cell_type))
  if(cell_type != 'none'){
    cell_type = cell_type[colnames(expr),]
  }
  expr = expr[!duplicated(rownames(expr)),]
  gc()
  #print(cell_type[1:10])
  return (list(expr,cell_type))
}

load_gmt_data = function(gmt_file_path){
  # load pathway gene sets' gmt file
  # file format:
  #            first index: pathway's name/ID
  #            second index: pathway's url or others, it dosen't matter
  #            third to all: gene symbols in pathway
  #            sep: \t
  # return a list
  tmp = readLines(gmt_file_path)
  gsets = list()
  for(i in 1:length(tmp)){
    t = strsplit(tmp[i],'\t')[[1]]
    genes = t[3:length(t)]
    genes = genes[which(genes != "")]
    gsets[[t[1]]] = genes
  }
  return (gsets)
}




### do not use this method
calcu_PAS = function(expr, 
                     dataSets, 
                     method='ssgsea', 
                     topo=TRUE){
  # calculate pathway activity score for every gene set
  # expr: gene expression matrix throught normalization
  # dataSets: list format, genes in pathways
  # method: method used for calculate pathway activity score for every sample.
  #         By default is sssgsea and other option are zscore, gsva, plage.
  # topo: use topological information as weight. Default: TRUE.
  if(topo){
    weight = NA
  }else{
    weight = NA
  }
  score = PAS(as.matrix(expr), dataSets, method = method, weight = weight)
  score
}

cal_cluster = function(mat, 
                       method='seurat',
                       seurat_dims,
                       seurat_resolution,
                       k_cluster,
                       min_pts){
  ##return a character
  res_cluster = clustering(as.matrix(mat), 
                        method=method,
                        seurat_dims = seurat_dims,
                        seurat_resolution = seurat_resolution,
                        k_cluster,
                        min_pts)
  return(res_cluster)
}

getVarib = function(sc){
  tryCatch({
    sc = Seurat::FindVariableFeatures(sc,selection.method = 'vst',verbose=F)
    return('vst')
  },error=function(e){
    tryCatch({
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'disp',verbose=F)
      return('disp')
    },error=function(e){
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'mvp',verbose=F)
      return('mvp')
    })
  })
}

prepa_seuratOb = function(mat, 
                          cluster_idents,
                          dims,
                          out_dir,
                          cells){

  number_pc = min(50, ncol(mat)-5)
  number_pc = min(number_pc, nrow(mat))
  print("number of pcs...................")
  print(number_pc)
  dims = min(dims, number_pc-2)

  
  cell_name_frame = data.frame(colnames(mat))
  colnames(cell_name_frame) = "cellName"
  print("before seurat")
  print(dim(mat))
  obSeurat = CreateSeuratObject(counts=mat,min.cells = 0,min.features = 0)
  print("creat success")
  obSeurat = FindVariableFeatures(obSeurat, verbose = FALSE,selection.method = getVarib(obSeurat))
  print("varible success")
  obSeurat = ScaleData(obSeurat, verbose = FALSE)
  #print('lalalallallalalallal') 
  obSeurat = RunPCA(obSeurat,verbose=F,seed.use = 42,npcs = number_pc)
  ###### locate idents
  print("pca success")
  Idents(obSeurat) = cluster_idents[colnames(mat)]
  #print('lololololo...................................') 
  
  
  perplexity = min(30, (ncol(mat)-2)%/%3)
  ###### write 3D-tsne
  print(dims)
  if((ncol(mat)-2)%/%3 < 1){
    print('too small cells')
  }
  obSeurat = RunTSNE(obSeurat,dims=1:dims,seed.use = 42, dim.embed = 3, perplexity = perplexity)
  tsne_seurat = as.data.frame(Embeddings(obSeurat[['tsne']]))
  tsne_seurat[,'cell_type'] = Idents(obSeurat)
  #print(dim(cell_name_frame))
  print(dim(tsne_seurat))
  tsne_seurat = cbind(cell_name_frame, tsne_seurat)
  rownames(tsne_seurat) = tsne_seurat[,1]
  write.csv(format(tsne_seurat[cells,], digits = 3), file.path(out_dir,'tsne_3D.csv'),quote=F, row.names=F)
  print("tsne 3D success")
  ###### write 2D-tsne
  obSeurat = RunTSNE(obSeurat,dims=1:dims,seed.use = 42, perplexity = perplexity)
  tsne_seurat = as.data.frame(Embeddings(obSeurat[['tsne']]))
  tsne_seurat[,'cell_type'] = Idents(obSeurat)
  tsne_seurat = cbind(cell_name_frame, tsne_seurat)
  rownames(tsne_seurat) = tsne_seurat[,1]
  write.csv(format(tsne_seurat[cells,], digits = 3), file.path(out_dir,'tsne_2D.csv'), quote=F, row.names=F)
  rm(tsne_seurat)
  gc()
  print("tsne 2D success")
  ###### write 3D-umap
  obSeurat = RunUMAP(obSeurat,dims=1:dims,seed.use = 42, umap.method = "umap-learn",metric = "correlation", n.components=3)
  umap_seurat = as.data.frame(Embeddings(obSeurat[['umap']]))
  umap_seurat[,'cell_type'] = Idents(obSeurat)
  umap_seurat = cbind(cell_name_frame, umap_seurat)
  rownames(umap_seurat) = umap_seurat[,1]
  write.csv(format(umap_seurat[cells,], digits = 3), file.path(out_dir,'umap_3D.csv'), quote=F, row.names=F)
  print("umap 3D success")
  ###### write 2D-umap
  obSeurat = RunUMAP(obSeurat,dims=1:dims,seed.use = 42, umap.method = "umap-learn",metric = "correlation")
  umap_seurat = as.data.frame(Embeddings(obSeurat[['umap']]))
  umap_seurat[,'cell_type'] = Idents(obSeurat)
  umap_seurat = cbind(cell_name_frame, umap_seurat)
  rownames(umap_seurat) = umap_seurat[,1]
  write.csv(format(umap_seurat[cells,], digits = 3), file.path(out_dir,'umap_2D.csv'), quote=F, row.names=F)
  print("umap 2D success")
  return(obSeurat)
}

plot_cluster = function(obSeurat, 
                        method = c('tsne','umap'),
                        out_dir
                        ){
  #print(out_dir)
  #print(method)
  options(bitmapType='cairo')
  method = match.arg(method)
  if(method == 'tsne'){
    #print(out_dir)
    #print(dim(obSeurat))
    png(file=out_dir)
    DimPlot(object = obSeurat, group.by = 'ident', reduction='tsne')
    dev.off()
  }else if(method =='umap'){
    png(file=out_dir)
    DimPlot(object = obSeurat, group.by = 'ident', reduction='umap')
    dev.off()
  }
}

findMarker = function(obSeurat,
                      method
                      )
{
  markers = FindAllMarkers(obSeurat, test.use=method,
                           min.cells.feature = 3,
                           min.cells.group = 3)
  # return all markers, include adjusted p-Value > 0.05
  return(markers)           
}

parseAllMarkers = function(obSeurat,
                           user_pathways,
                           method,
                           logfc_thre,
                           min.pct,
                           shown_markers,
                           out_dir,
                           para_size,
                           expr,
                           cell_type,
                           path_list,
                           cells
){
  all_markers = FindAllMarkers_me(obSeurat, 
                                 test.use = method,
                                 logfc.threshold = logfc_thre,
                                 min.pct = min.pct,
                                 only.pos = T,
                                 min.cells.feature = 3,
                                 min.cells.group = 3,
                                 para_size = para_size)
  print("from function::::::")
  print(dim(all_markers))
  if(nrow(all_markers) == 0){
    return(data.frame())
  }
  
  ## subset pathways
  all_markers = subset(all_markers, select=-c(pct.1,pct.2))
  print("after subset:")
  print(dim(all_markers))
  
  ## write excel to web
  path_new_names = gsub('-','_',all_markers$gene)
  pathway_frame = as.data.frame(path_new_names)
  colnames(pathway_frame) = "pathways"
  all_markers_write = cbind(pathway_frame, all_markers[,1:(ncol(all_markers)-1)])
  #all_markers_write[,'cluster'] = paste0('cellType:', all_markers_write[,'cluster'])
  #all_markers['pathways'] = all_markers$gene
  #all_markers = subset(all_markers, select = -gene)
  genes = get_geneLis(paths = path_new_names,
                      path_lis = path_list,
                      expr_genes = rownames(expr))
  print(length(genes))
  all_markers_write['geneList'] = genes
  print(dim(all_markers_write))
  
  
  write.csv(format(all_markers_write, digits = 3), file.path(out_dir,'pas_markers.csv'),quote=F,row.names=F)
  print("write pas_markers success")
  #print(head(all_markers))
  n_top = 10
  tops = all_markers %>% group_by(cluster) %>% top_n(n = n_top, wt = avg_logFC)
  print("groupby success")
  
  while( nrow(tops) > 50 ){
    n_top = n_top - 1
    tops = all_markers %>% group_by(cluster) %>% top_n(n = n_top, wt = avg_logFC)
  }
  
  #print(dim(tops))
  p_heatmap = DoHeatmap_me(obSeurat, features = unique(tops$gene), size=6, angle = 90, group.bar.height = 0.01,label=FALSE)
  png(file = file.path(out_dir, 'Heatmap.png'),width = 1200, height = 800)
  print(p_heatmap)
  dev.off()

  pdf(file.path(out_dir, 'Heatmap.pdf'),width = 16,height = 10)  
  print(p_heatmap)
  dev.off()
  rm(p_heatmap)
  gc()

  #p_bubble = bubble_plot(obSeurat,as.data.frame(tops))
  #png(file = file.path(out_dir, 'Bubble.png'),width = 1200, height = 800)
  #print(p_bubble)
  #dev.off()

  #pdf(file.path(out_dir, 'Bubble.pdf'),width = 16,height = 10)
  #print(p_bubble)
  #dev.off()
  
  print("plot success")
  
  
  if(!is.null(user_pathways)){
    feature_pathways = unique(tops$gene, user_pathways)
  }else{
    feature_pathways = unique(tops$gene)
  }
  
  rm(tops)
  gc()
  
  expr = expr[,colnames(obSeurat)]
  #print(head(cell_type))
  expr = scale(expr[,cells])
  cell_type = as.data.frame(Idents(obSeurat)[cells]) 
  print("after seurat ...")
  print(head(cell_type))
  print(dim(expr))
  gc()
  print("zscore of expression success, dimensions: ")
  print(dim(expr))
  
  #rownames(cell_type) = colnames(expr)
  #n_celltype = length(feature_pathways)
  #print("number of feature pathways:")
  #print(length(feature_pathways))
  #if(n_celltype>7){
  #  x_angle = 10
  #}else{
  x_angle=0
  #}
  
  mclapp = get('mclapply', envir = getNamespace('parallel'))
  options(mc.cores = para_size)
  
  #expr_z = apply(expr[,cells], 2, FUN = function(x){(x-mean(x))/sd(x)})
  
  #expr_z = do.call("cbind", expr_z)
  #print(dim(expr_z))
  #dim_eval = dim(expr) == dim(expr_z)
  #if(dim_eval)
  #colnames(expr_z) = colnames(expr)
  #rownames(expr_z) = rownames(expr)
  #rm(expr)
  
  
  cat(length(feature_pathways))

  mclapp(1:length(feature_pathways),function(i){
    feature_name = gsub('-','_',feature_pathways[i])
    #print(feature_name)
    path_dir = file.path(out_dir, feature_name)
    dir.create(path_dir)
	#print(paste0("function.R, parseAllMarkers, path_dir: ",path_dir))
    #print(feature_name)
    pic_name_vil = paste0(feature_name,'.vil.png')
    pic_name_fet = paste0(feature_name,'.fet.png')
    pic_name_het = paste0(feature_name,'.het.png')
    mat_name_het = paste0(feature_name,'.csv')
    feature_vio(obSeurat, feature_pathways[i], file.path(path_dir, pic_name_vil), x_angle)
    feature_highlight(obSeurat, feature_pathways[i], file.path(path_dir, pic_name_fet))
    #if(ncol(expr_z) < 15000){
      #cat("Do not write gene matrix csv")
    #print(getTrueNames(feature_name))
    #print(feature_name %in% names(path_list))
    #print(path_list[[feature_name]])
    #print(length(path_list[[feature_name]]))
    #print(intersect(rownames(expr), path_list[[feature_name]]))
    genes_heatmap(expr, path_list[[feature_name]], cell_type, file.path(path_dir, pic_name_het), file.path(path_dir, mat_name_het))
    #}
    
  })
  #for(i in 1:length(feature_pathways)){
    #print(feature_pathways[i]) 
    #print(names(path_list)[1:30])
    #print(feature_pathways[i])
    #print(length(path_list[[feature_pathways[i]]]))
   # feature_name = gsub('-','_',feature_pathways[i])
    #feature_name = gsub(' ','_',feature_pathways[i])
    #path_dir = file.path(out_dir, feature_name)
    #dir.create(path_dir)
    #pic_name_vil = paste0(feature_name,'.vil.png')
    #pic_name_fet = paste0(feature_name,'.fet.png')
    #pic_name_het = paste0(feature_name,'.het.png')
    #mat_name_het = paste0(feature_name,'.csv')    
 
    
    #print(expr[1:3,1:3])
    #print(path_list[[feature_name]][1:5])
   # feature_vio(obSeurat, feature_pathways[i], file.path(path_dir, pic_name_vil))
   # feature_highlight(obSeurat, feature_pathways[i], file.path(path_dir, pic_name_fet))
   # genes_heatmap(expr, path_list[[feature_name]], cell_type, file.path(path_dir, pic_name_het), file.path(path_dir, mat_name_het))
  #}

  #feature_frame = as.data.frame(feature_pathways)
  #write.table(feature_frame, file = file.path(out_dir, 'pathway_names.csv'), col.names=F, row.names=F, quote=F)
  return(all_markers)
}

plot_feature = function(obSeurat, feature, reduction){
  FeaturePlot(obSeurat, features = feature, reduction = reduction)
}


#load_user_pathway()


plotPAS_heatmap = function(mat,
                           cluster_idents,
                           out_path,
                           max_path = 500){
  #print(cluster_idents[1:5])
  colors = hue_pal()(length(unique(cluster_idents)))
  #print(length(cluster_idents))
  #print(colors)
  names(colors) = levels(cluster_idents)
  print(colors)
  print(gsub('pas_heatmap.png','pas_color.csv', out_path))
  write.csv(as.data.frame(colors),file = gsub('pas_heatmap.png','pas_color.csv', out_path),quote=F)
  print("write color success")
  #names(cluster_idents) = colnames(mat)
  cluster_idents = as.data.frame(cluster_idents)
  print(head(cluster_idents))
  colnames(cluster_idents) = "cell types"
  mat = mat[,rownames(cluster_idents)]
  
  ####filter pathways
  if(dim(mat)[1] > max_path){
    #pathway_sd = apply(mat,1, sd)
    #mat = mat[names(sort(pathway_sd,decreasing = T)[0:max_path]),]
    # edited by zyr on 2020.04.07
    mat_obj = CreateAssayObject(mat)
    mat_obj = FindVariableFeatures(mat_obj, nfeatures=max_path)
    aa = mat_obj@var.features
    aa = gsub('-','_',aa)
    rownames(mat) = gsub('-','_',rownames(mat))
    mat = mat[aa,]
  }
 
   
  p = pheatmap(mat,
               color = colorRampPalette(c("forestgreen", "white", "orange"))(100),
               #color = "PiYG",
               show_rownames = F,
               show_colnames = F,
               fontsize = 20,
               annotation_names_col = F,
               annotation_colors = list("cell types" = colors),
               cluster_cols = F,
               annotation_col = cluster_idents)

  png(out_path,width = 1200,height = 800)
  print(p)
  dev.off()

  pdf(gsub('png','pdf', out_path),width = 16,height = 10)
  print(p)
  dev.off()
}



getTrueNames = function(c){
  m = regexec("[^0-9.*]+(?=[0-9.*])",c,perl=T)
  res = as.character(regmatches(c,m)[[1]])
  res
}



