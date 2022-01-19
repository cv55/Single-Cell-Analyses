##QC plots visualization function
seurat_qc<- function(sobject,nFeature_RNA_min,nFeature_RNA_max,percent_mito,percent_ribo,
                     nCount_RNA_min, nCount_RNA_max){
  
  
  sobject[["percent.mt"]] <- PercentageFeatureSet(sobject, pattern = "^MT-")
  sobject[["percent.rb"]] <- PercentageFeatureSet(sobject, pattern = "^RP[SL]")
  
  
  counts_df<- data.frame('Read_Counts'=sobject@meta.data['nCount_RNA'], 'Metric'='Read_Counts','Filter'='OG')
  colnames(counts_df)<- c('Read_Counts','Metric','Filter')
  genes_df<- data.frame('Gene_Counts'=sobject@meta.data['nFeature_RNA'],'Metric'='Gene_Counts','Filter'='OG')
  colnames(genes_df)<- c('Gene_Counts','Metric','Filter')
  mito<- data.frame("percent.mt"=sobject@meta.data['percent.mt'],'Metric'= 'MT%','Filter'='OG')
  
  
  sobject <- subset(sobject, subset = nFeature_RNA>nFeature_RNA_min & nFeature_RNA<nFeature_RNA_max & percent.mt<percent_mito & percent.rb< percent_ribo &
                      nCount_RNA < nCount_RNA_max & nCount_RNA > nCount_RNA_min )
  counts_df_qc<- data.frame('Read_Counts'=sobject@meta.data['nCount_RNA'], 'Metric'='Read_Counts','Filter'='QC')
  colnames(counts_df_qc)<- c('Read_Counts','Metric','Filter')
  genes_df_qc<- data.frame('Gene_Counts'=sobject@meta.data['nFeature_RNA'],'Metric'='Gene_Counts','Filter'='QC')
  colnames(genes_df_qc)<- c('Gene_Counts','Metric','Filter')
  mito_qc<- data.frame("percent.mt"=sobject@meta.data['percent.mt'],'Metric'= 'MT%','Filter'='QC')
  
  counts_df<-rbind(counts_df,counts_df_qc)
  genes_df<-rbind(genes_df,genes_df_qc)
  mito<- rbind(mito,mito_qc)
  
  p1<-ggplot(counts_df, aes(x=Metric, y=Read_Counts)) + geom_violin(fill='firebrick4') + 
    theme_bw() + facet_grid(cols  = vars(Metric), rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5)
  
  p2<-ggplot(genes_df, aes(x=Metric, y=Gene_Counts)) + geom_violin(fill='turquoise4') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 ) 
  
  p3<-ggplot(mito, aes(x=Metric, y=percent.mt)) + geom_violin(fill='forestgreen') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5) 
  
  print(p1+p2+p3)
  return(sobject)
}




####NORMALIZATION FUNCTION
seurat_norm<- function(sobject, method=NULL, regress.out=c()){
  
  if(method=="SCT"){
    print("SCT is applied")
    sobject=SCTransform(sobject, vars.to.regress = regress.out )
    p=VariableFeaturePlot(sobject, selection.method="sct")
    #p2 <- LabelPoints(plot = p, points = sobject@assays$RNA@var.features[0:10], repel = TRUE)
    print(p)
    return(sobject)
  }
  
  else if(method=="vst"){   
    print("vst is applied")
    sobject <- NormalizeData(sobject)
    sobject <- FindVariableFeatures(sobject, selection.method = "vst",nfeatures = 2000)
    sobject <- ScaleData(sobject, vars.to.regress = regress.out)
    p=VariableFeaturePlot(sobject, selection.method='vst')
    #p2 <- LabelPoints(plot = p, points = top10, repel = TRUE)
    #print(p)
    return(sobject)
  }
}


####INTERGRATION FUNCTION
integration<- function(sobject_list,method=NULL){
  
  if(method=='Seurat'){
    pca_list=list()
    #counter=1
    #merged_object<- merge(sobject_list[[1]], y=c(sobject_list[-1]))
    #merged_norm_object<- seurat_norm(merged_object, method="vst")
    #merged_norm_object <- RunPCA(merged_norm_object, npcs = 30)
    #p <- DimPlot(object = merged_norm_object, reduction = "pca", pt.size = .1, group.by = "Group")
    #print(p)
    #print("Freeze the function for 10 seconds  to view the figure!!!")
    #Sys.sleep(10)
    
    for(obj in sobject_list){
      pca_list[[paste0("Individual_", counter)]]<-RunPCA(obj, npcs=30)
      counter=counter+1
    }
    print("PCA list is made")
    print("Perform integration")
    features <- SelectIntegrationFeatures(object.list = sobject_list, nfeatures=2000)
    anchors <- FindIntegrationAnchors(object.list = sobject_list, anchor.features = features)
    combined <- IntegrateData(anchorset = anchors)
    DefaultAssay(combined) <- "integrated"
    #combined <- ScaleData(combined, verbose = FALSE)
    #combined <- RunPCA(combined, verbose = FALSE)
    #  p <- DimPlot(object = merged_norm_object, reduction = "pca", pt.size = .1, group.by = "Group")
    #  print(p)
    return(combined)
  }
  
  else if(method=='Harmony'){
    library(harmony)
    print("This part should be run on the raw data, first we merge the data, call the seurat_norm funvtion  on the merged data, run pca to view the batch effects and then we use harmony to integrate the data")
    library(SeuratData)
    merged_object<- merge(sobject_list[[1]], y=c(sobject_list[-1])) 
    merged_norm_object<- seurat_norm(merged_object, method="vst")
    merged_norm_object <- RunPCA(merged_norm_object, npcs = 30)
    #p <- DimPlot(object = merged_norm_object, reduction = "pca", pt.size = .1, group.by = "Group")
    #print(p)
    #print("Freeze the function for 10 seconds  to view the figure!!!")
    #Sys.sleep(10)
    
    print("Now let's quickly view the convergence")
    integrated_object <- merged_norm_object %>% RunHarmony("Group", plot_convergence = TRUE)
    Sys.sleep(5)
    
    p <- DimPlot(object = integrated_object, reduction = "harmony", pt.size = .1, group.by = "Group")
    print(p)
    return(integrated_object)
  }
}

#make umap for single object, normalized and qc filtered
makeUMAP<- function(sobject){
  sobject <- RunPCA(sobject, npcs = 30,features = VariableFeatures(object = sobject))
  print(ElbowPlot(sobject))
  Sys.sleep(10)
  sobject<-FindNeighbors(sobject, dims = 1:10)
  sobject <- FindClusters(sobject, verbose = TRUE,resolution = 0.5, algorithm =1)
  sobject <- RunUMAP(sobject, dims = 1:13)
  p=DimPlot(sobject, reduction = "umap",label=T)
  print(p)
  return(sobject)
}


makeUMAP_CCA<- function(combined_object){
  
  DefaultAssay(combined_object) <- "integrated"
  combined_object <- ScaleData(combined_object)
  combined_object <- RunPCA(combined_object, npcs = 30)
  print(ElbowPlot(combined_object))
  Sys.sleep(10)
  combined_object <- RunUMAP(combined_object, reduction = "pca", dims = 1:8)
  combined_object <- FindNeighbors(combined_object, reduction = "pca", dims = 1:8)
  combined_object <- FindClusters(combined_object, resolution = 0.5,graph.name = "integrated_snn")
  
  p1=DimPlot(combined_object, reduction = "umap", group.by = "Group")
  #Sys.sleep(10)
  p2=DimPlot(combined_object, reduction = "umap", label = TRUE, repel = TRUE)
  print(p1+p2)
  return(combined_object)
}








