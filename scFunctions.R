##QC plots visualization function
seurat_qc<- function(sobject,name,folder,nFeature_RNA_min,nFeature_RNA_max,percent_mito,percent_ribo,
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
  
  setEPS()
  postscript(paste0(folder,name,".eps"))
  plot(p1+p2+p3)
  dev.off()
  return(sobject)
}




####NORMALIZATION FUNCTION
seurat_norm<- function(sobject, method=NULL, regress.out=c()){
  DefaultAssay(sobject) <- "RNA"
  
  if(method=="SCT"){
    print("SCT is applied")
    sobject=SCTransform(sobject, vars.to.regress = regress.out, verbose=FALSE )
    p=VariableFeaturePlot(sobject, selection.method="sct")
    #p2 <- LabelPoints(plot = p, points = sobject@assays$RNA@var.features[0:10], repel = TRUE)
    print(p)
    return(sobject)
  }
  
  else if(method=="vst"){   
    sobject <- NormalizeData(sobject)
    sobject <- FindVariableFeatures(sobject, selection.method = "vst",nfeatures = 2000, verbose=FALSE)
    sobject <- ScaleData(sobject, vars.to.regress = regress.out)
    return(sobject)
  }
}










####INTERGRATION FUNCTION
integration<- function(sobject_list,method=NULL, normalization=NULL, resolution=NULL){
  if(length(resolution)>1){
    resolutions=paste(resolution, collapse = '_')
  }
  else{resolutions=resolution}
  
  if(method=='Seurat' && normalization=='vst'){
    
   print("Plot unintegrated data")
      
      norm_list <- lapply(X = sobject_list, FUN = function(x) {
      x<-NormalizeData(x)
      x<- FindVariableFeatures(x, selection.method = normalization, nfeatures = 2000)
    })
      
      
    features <- SelectIntegrationFeatures(object.list = norm_list, nfeatures = 2000)
    unintegrated_data<- merge(norm_list[[1]], y=c(norm_list[-1])) 
    unintegrated_data <- ScaleData(unintegrated_data, verbose = FALSE)
    unintegrated_data <- RunPCA(unintegrated_data, npcs = 30, verbose = FALSE, features=features)
    p1<-ElbowPlot(unintegrated_data)
    unintegrated_data <- RunUMAP(unintegrated_data, reduction = "pca", dims = 1:15,verbose = FALSE)
    unintegrated_data <- FindNeighbors(unintegrated_data, reduction = "pca", dims = 1:15,verbose = FALSE)
    unintegrated_data <- FindClusters(unintegrated_data, resolution = resolution,verbose = FALSE)  
    p2 <- DimPlot(unintegrated_data, reduction = "umap", group.by = "Group")
    print(p1+p2)
      
      
     print("Perform integration with Seurat and vst")
    
    
    features <- SelectIntegrationFeatures(object.list = norm_list)
    anchors <- FindIntegrationAnchors(object.list = norm_list, anchor.features = features)
    integrated <- IntegrateData(anchorset = anchors)
    DefaultAssay(integrated) <- "integrated"
    
    integrated <- ScaleData(integrated, verbose = FALSE)
    integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
      
    integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:15,verbose = FALSE)
    integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:15,verbose = FALSE)
    integrated <- FindClusters(integrated, resolution = resolution,verbose = FALSE)
    
    
    umap_values<-integrated@reductions$umap@cell.embeddings
    group_values<-integrated@meta.data$Group
    write.table(group_values, paste(method,normalization,resolutions,'group.txt',sep='_'))
    write.table(umap_values, paste(method,normalization,resolutions,'umap.txt',sep='_'))
    
    #p3 <- DimPlot(integrated, reduction = "umap", group.by = "Group")
    #p4 <- DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
    #print(p3 + p4)
    
    
    
    return(integrated)
    
    
  }
  
  else if(method=='Seurat' && normalization=='SCT'){
    print("Perform integration with Seurat and SCT")
    
    norm_list <- lapply(X = sobject_list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = norm_list, nfeatures = 3000)
    norm_list <- PrepSCTIntegration(object.list = norm_list, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = norm_list, normalization.method = "SCT",anchor.features = features)
    integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    
    integrated <- RunPCA(integrated, verbose = FALSE)
    integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30,verbose = FALSE)
    integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30,verbose = FALSE)
    integrated <- FindClusters(integrated, resolution = resolution)
    
    
    umap_values<-integrated@reductions$umap@cell.embeddings
    group_values<-integrated@meta.data$Group
    write.table(group_values, paste(method,normalization,resolutions,'group.txt',sep='_'))
    write.table(umap_values, paste(method,normalization,resolutions,'umap.txt',sep='_'))
    
    #p1 <- DimPlot(integrated, reduction = "umap", group.by = "Group")
    #p2 <- DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
    #print(p1 + p2)
    return(integrated)
    
  }
  
  
  
  
  else if(method=='Harmony'&& normalization=='SCT'){
    library(harmony)
    library(SeuratData)
    
    
    norm_list <- lapply(X = sobject_list, FUN = SCTransform,do.scale=T,do.center=T)
    features <- SelectIntegrationFeatures(object.list = norm_list, nfeatures = 3000)
    
    merged_object<- merge(norm_list[[1]], y=c(norm_list[-1])) 
    merged_object <- RunPCA(merged_object, npcs = 30, features=features, assay='SCT')
    
    integrated <- merged_object %>% RunHarmony("Group", plot_convergence = FALSE, assay.use='SCT')
    integrated <- integrated%>% RunUMAP(reduction = "harmony", dims = 1:10) %>% FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
      FindClusters(resolution = resolution) %>% 
      identity()
    
    
    umap_values<-integrated@reductions$umap@cell.embeddings
    group_values<-integrated@meta.data$Group
    write.table(group_values, paste(method,normalization,resolutions,'group.txt',sep='_'))
    write.table(umap_values, paste(method,normalization,resolutions,'umap.txt',sep='_'))
    
    
    #p1=DimPlot(integrated, reduction = "umap", pt.size = .1,group.by = "Group")
    #p2=DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
    #print(p1+p2)
    
    return(integrated)
  }
  
  
  
  
  else if(method=='Harmony'&& normalization=='vst'){
    library(harmony)
    library(SeuratData)
    
    merged_object<- merge(sobject_list[[1]], y=c(sobject_list[-1])) 
    merged_norm_object<- seurat_norm(merged_object, method="vst",)
    merged_norm_object <- RunPCA(merged_norm_object, npcs = 30)
    
    integrated <- merged_norm_object %>% RunHarmony("Group", plot_convergence = FALSE)
    integrated <- integrated%>% RunUMAP(reduction = "harmony", dims = 1:10) %>%FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
      FindClusters(resolution = resolution,graph.name=) %>% 
      identity()
    
    
    umap_values<-integrated@reductions$umap@cell.embeddings
    group_values<-integrated@meta.data$Group
    write.table(group_values, paste(method,normalization,resolutions,'group.txt',sep='_'))
    write.table(umap_values, paste(method,normalization,resolutions,'umap.txt',sep='_'))
    
    
    #p1=DimPlot(integrated, reduction = "umap", pt.size = .1,group.by = "Group")
    #p2=DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
    #print(p1+p2)
    
    return(integrated)
  }
  
  
  
}


      
      
      
###MULTIOME QC      

seurat_qc_multiome<- function(sobject,nCount_ATAC_min,nCount_ATAC_max,
                     nCount_RNA_min, nCount_RNA_max, nFeature_RNA_min,nFeature_RNA_max,percent_mito,percent_ribo, TSS, nucl_signal){
  
  DefaultAssay(sobject) <- "RNA"
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
  
  DefaultAssay(sobject) <- "ATAC"
  sobject <- TSSEnrichment(sobject)
  sobject <- NucleosomeSignal(sobject)
  
  atac_df<- data.frame("ATAC_Counts"=sobject@meta.data['nCount_ATAC'],'Metric'= 'ATAC_Counts','Filter'='OG')
  colnames(atac_df)<- c('ATAC_Counts','Metric','Filter')
  tss.enrich<- data.frame("TSS.enrichment"=sobject@meta.data["TSS.enrichment"], 'Metric'= "TSS.enrichment",'Filter'='OG' )
  nucl.sign<- data.frame("nucleosome_signal"=sobject@meta.data["nucleosome_signal"], 'Metric'= "nucleosome_signal",'Filter'='OG' )
  
  
  sobject <- subset(sobject, subset = nCount_ATAC>nCount_ATAC_min & nCount_ATAC<nCount_ATAC_max &
                                      TSS.enrichment > TSS & nucleosome_signal<nucl_signal)
  
  
  atac_df_qc<- data.frame("ATAC_Counts"=sobject@meta.data['nCount_ATAC'],'Metric'= 'ATAC_Counts','Filter'='QC')
  colnames(atac_df_qc)<- c('ATAC_Counts','Metric','Filter')
  tss.enrich_qc<- data.frame("TSS.enrichment"=sobject@meta.data["TSS.enrichment"], 'Metric'= "TSS.enrichment",'Filter'='QC' )
  nucl.sign_qc<- data.frame("nucleosome_signal"=sobject@meta.data["nucleosome_signal"], 'Metric'= "nucleosome_signal",'Filter'='QC' )
  
  
  atac_df<-rbind(atac_df,atac_df_qc)
  tss.enrich_df<- rbind(tss.enrich,tss.enrich_qc)
  nucl.sign_df<- rbind(nucl.sign,nucl.sign_qc)
  
  
  p1<-ggplot(counts_df, aes(x=Metric, y=Read_Counts)) + geom_violin(fill='firebrick4') + 
    theme_bw() + facet_grid(cols  = vars(Metric), rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5)
  
  p2<-ggplot(genes_df, aes(x=Metric, y=Gene_Counts)) + geom_violin(fill='turquoise4') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 ) 
  
  p3<-ggplot(mito, aes(x=Metric, y=percent.mt)) + geom_violin(fill='forestgreen') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5) 

  p4<-ggplot(atac_df, aes(x=Metric, y=ATAC_Counts)) + geom_violin(fill='turquoise4') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 ) 
  
  
  p5<-ggplot(tss.enrich_df, aes(x=Metric, y=TSS.enrichment)) + geom_violin(fill='firebrick4') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5) 
  
  p6<-ggplot(nucl.sign_df, aes(x=Metric, y=nucleosome_signal)) + geom_violin(fill='forestgreen') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5) 
  
  
    
  print(p1+p2+p3+p4+p5+p6)
  return(sobject)
}




    
    
    
    
    
    
 ####ATAC QC   
seurat_qc_atac<- function(sobject,nCount_ATAC_min,nCount_ATAC_max,pct.peaks, TSS, nucl_signal){
  
  sobject$pct_reads_in_peaks <- sobject$peak_region_fragments / sobject$passed_filters * 100
  sobject$blacklist_ratio <- sobject$blacklist_region_fragments / sobject$peak_region_fragments
  
  
  atac_df<- data.frame("ATAC_Counts"=sobject@meta.data['nCount_ATAC'],'Metric'= 'ATAC_Counts','Filter'='OG')
  colnames(atac_df)<- c('ATAC_Counts','Metric','Filter')
  pct_reads<- data.frame('Reads_in_Peaks'=sobject@meta.data['pct_reads_in_peaks'], 'Metric'= '%Reads_in_Peaks','Filter'='OG')
  colnames(pct_reads)<- c('Reads_in_Peaks','Metric','Filter')
  tss.enrich<- data.frame("TSS.enrichment"=sobject@meta.data["TSS.enrichment"], 'Metric'= "TSS.enrichment",'Filter'='OG' )
  nucl.sign<- data.frame("nucleosome_signal"=sobject@meta.data["nucleosome_signal"], 'Metric'= "nucleosome_signal",'Filter'='OG' )
  
  
  sobject <- subset(sobject, subset = nCount_ATAC>nCount_ATAC_min & nCount_ATAC<nCount_ATAC_max & pct_reads_in_peaks >pct.peaks &
                      TSS.enrichment > TSS & nucleosome_signal<nucl_signal)
  
  
  atac_df_qc<- data.frame("ATAC_Counts"=sobject@meta.data['nCount_ATAC'],'Metric'= 'ATAC_Counts','Filter'='QC')
  colnames(atac_df_qc)<- c('ATAC_Counts','Metric','Filter')
  pct_reads_qc<- data.frame('Reads_in_Peaks'=sobject@meta.data['pct_reads_in_peaks'], 'Metric'= '%Reads_in_Peaks','Filter'='QC')
  colnames(pct_reads_qc)<- c('Reads_in_Peaks','Metric','Filter')
  tss.enrich_qc<- data.frame("TSS.enrichment"=sobject@meta.data["TSS.enrichment"], 'Metric'= "TSS.enrichment",'Filter'='QC' )
  nucl.sign_qc<- data.frame("nucleosome_signal"=sobject@meta.data["nucleosome_signal"], 'Metric'= "nucleosome_signal",'Filter'='QC' )
  
  
  atac_df<-rbind(atac_df,atac_df_qc)
  pct_reads_df<- rbind(pct_reads,pct_reads_qc)
  tss.enrich_df<- rbind(tss.enrich,tss.enrich_qc)
  nucl.sign_df<- rbind(nucl.sign,nucl.sign_qc)
  
  
  
  p1<-ggplot(atac_df, aes(x=Metric, y=ATAC_Counts)) + geom_violin(fill='turquoise4') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 ) 
  
  p2<-ggplot(pct_reads_df, aes(x=Metric, y=Reads_in_Peaks)) + geom_violin(fill='dodgerblue3') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 ) 
  
  
  p3<-ggplot(tss.enrich_df, aes(x=Metric, y=TSS.enrichment)) + geom_violin(fill='firebrick4') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5) 
  
  p4<-ggplot(nucl.sign_df, aes(x=Metric, y=nucleosome_signal)) + geom_violin(fill='forestgreen') + 
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5) 
  
  
  
  print(p1+p2+p3+p4)
  return(sobject)
  
  
  
  
}



    
    
    
##PEAK CALLING UNSING MACS2
callpeaks<- function(sobject,fragmentpath){
  
  DefaultAssay(sobject) <- "ATAC"
  peaks <- CallPeaks(sobject,macs2.path = "/Users/u0149445/opt/anaconda3/bin/macs2")
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

  macs2_counts <- FeatureMatrix(
    fragments = Fragments(sobject),
    features = peaks,
    cells = colnames(sobject)
  )
  
  sobject[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragmentpath,
    annotation = annotation
  )
  return(sobject)
  
  
}










###READ ATAC OBJECT
make_atac_object<- function(h5_path,fragment_path,metadata_path, metadata, annotation, combined_peaks){
  atac_counts <- Read10X_h5(h5_path)
  
 # grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
 # grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
 # atac_counts <- atac_counts[as.vector(grange.use), ]
  
  #barcodes<- read.table(barcode_path, sep = ',', header = TRUE)
  #barcodes<- barcodes[barcodes$atac_barcode %in% colnames(atac_counts),]
  #barcodes <- barcodes[order(barcodes$atac_barcode),]
  
 # metadata <- read.csv(file = metadata_path,header = TRUE, row.names = 1)
 # metadata<- metadata[,c("peak_region_fragments","blacklist_region_fragments","passed_filters")]
 # metadata <-metadata[rownames(metadata) %in% colnames(atac_counts),] 
  
 # fragpth <-fragment_path
  
 # atac_assay <- CreateChromatinAssay( counts = atac_counts, sep = c(":", "-"), fragments=fragpth, 
                                      #annotation = annotation,min.cells = 1)
 # atac_object <- CreateSeuratObject(counts = atac_assay, assay = "ATAC",meta.data = metadata)
  #atac_object<- RenameCells(atac_object, new.names=barcodes$gex_barcode)
  
 # return(atac_object)
  
    atac_counts <- FeatureMatrix(
    fragments <- CreateFragmentObject(path = fragment_path, cells = colnames(atac_counts)),
    features = combined_peaks,
    cells = colnames(atac_counts)
        
  )

  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
 
  
  fragpth <-fragment_path
  atac  <- CreateChromatinAssay( counts =atac_counts, sep = c(":", "-"), fragments=fragpth, annotation = annotation)
    atac_object <- CreateSeuratObject(counts = atac, assay = "ATAC")

    
}


    

##Make combined peaks
combine_peaks<- function(bed_list){
        grList=list()
        for(i in bed_list){
        
            bed<-read.table(i, col.names = c("chr", "start", "end"))
            grList[[i]]<-makeGRangesFromDataFrame(bed)
            
           }
        print(grList)
      
    
    combined.peaks <-  GenomicRanges::reduce(x = c( grList[[1]], grList[[2]] ))
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
    chrs <-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

    combined.peaks<-combined.peaks[which(combined.peaks@seqnames %in% chrs),]
        return(combined.peaks)
 }
    
    
    
    
###READ Multiome OBJECT

make_multiome_object<- function(h5_path,fragment_path,annotation,combined_peaks){
  counts <- Read10X_h5(h5_path)
  
  rna_counts <- counts$`Gene Expression`
  object <- CreateSeuratObject(counts = rna_counts,assay = "RNA", min.cells=3)
  
  
  #atac_counts <- counts$Peaks
  #grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  #grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  #atac_counts <- atac_counts[as.vector(grange.use), ]

    
 
   atac_counts <- FeatureMatrix(
   fragments <- CreateFragmentObject(path = fragment_path, cells = colnames(object)),
    features = combined_peaks,
    cells = colnames(object)
  )
  
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
 
  
  fragpth <-fragment_path
  object[['ATAC']]  <- CreateChromatinAssay( counts =atac_counts, sep = c(":", "-"), fragments=fragpth, annotation = annotation)
  
    #new_counts <- FeatureMatrix(
   ## fragments = Fragments(object),
  #  features = combined_peaks,
  #  cells = colnames(object)
 # )
   
  return(object)
  
}

    
    
    
    
    
    
    
    
###NORMALIZE ATAC DATA
normalize_atac<- function(sobject){
  
  DefaultAssay(sobject) <- "peaks"
  
  sobject <- FindTopFeatures(sobject, min.cutoff = 5)
  sobject <- RunTFIDF(sobject)
  sobject <- RunSVD(sobject)
  return(sobject)
  
  }


    
    

    

make_UMAP<- function(sobject,name, folder, dims, resolution, res_to_plot){

  sobject<- RunPCA(sobject,verbose = FALSE)
  p<-ElbowPlot(sobject)
  sobject<-FindNeighbors(sobject, dims = 1:dims,verbose = FALSE)
  sobject <- FindClusters(sobject, verbose = FALSE, resolution = res_to_plot, algorithm =1)
  sobject <- RunUMAP(sobject, dims = 1:dims,verbose = FALSE)
  c=DiscretePalette(25,"alphabet")
  
  df1=data.frame(sobject@reductions$umap@cell.embeddings)
  df1['plot']="UMAP"
  df1['cluster']=sobject@meta.data$seurat_clusters
  
  p1=ggplot(df1, aes(x=UMAP_1, y=UMAP_2,color=cluster)) + geom_point() + scale_color_manual(values=c) + theme_bw() +
    facet_grid(cols  = vars(plot), scales = "free")
  
  df<- sobject@meta.data[,c('nCount_RNA','seurat_clusters')]
  df['plot']="Read_Counts"
  
  p2=ggplot(df, aes(x=seurat_clusters, y=nCount_RNA, fill=as.factor(seurat_clusters))) + geom_bar(stat = "identity") + scale_fill_manual(values=c) +
    facet_grid(cols  = vars(plot),  scales = "free")
  
  all<-data.frame()
  for(i in order(as.numeric(levels(df$seurat_clusters)))-1){
  cluster_df <- as.data.frame(df %>% subset(seurat_clusters == i))  #%>% .$Group %>% table() #/ table(df$seurat_clusters))
  cluster_df<- sum(cluster_df$nCount_RNA)/nrow(cluster_df)
    all<-rbind(all,cluster_df)  
      
  }
    colnames(all)<-'average_read'
    all['cluster']<-order(as.numeric(levels(df$seurat_clusters)))-1
    all['plot']="Average_reads"
   
    p3<- ggplot(all, aes(x=cluster, y=average_read, fill=as.factor(cluster))) + geom_bar(stat = "identity") + scale_fill_manual(values=c) +
    facet_grid(cols  = vars(plot),  scales = "free")
    
    
  print(p+p1+p2+p3)
  
  
  setEPS()
  postscript(paste0(folder,name,".eps"), width=12)
  plot(p+p1+p2+p3)
  dev.off()

  
  sobject <- FindClusters(sobject, verbose = FALSE,resolution = resolutions, algorithm =1)
  
  return(sobject)
  }




    
    
###MAKE LOOM FILE
make_loom<- function(name,sobject,marker_file_list){
  

  
    build_loom(                                                                                                                                                                                                                                                                   
      file.name=name,                                                                                                                                                                                                                                                        
      dgem=sobject@assays$RNA@counts,                                                                                                                                                                                                                                                              
      title="ASAP_integrated", # A name for your data                                                                                           
      genome="Human", # Just for user information, not used internally                                                                                                                                                                                                            
      default.embedding = sobject@reductions$umap@cell.embeddings,
      default.embedding.name = "umap_gex")   




    loom <- open_loom(name, mode = "r+")




    add_seurat_clustering(loom = loom,
                      seurat = sobject,
                      seurat.assay = "RNA",
                      seurat.clustering.prefix = "integrated_snn_res.",
                      default.clustering.resolution = "res.0.5",
                      seurat.markers.file.path.list = marker_file_list,
                      seurat.marker.metric.accessors = c("avg_log2FC", "p_val_adj"),
                      seurat.marker.metric.names = c("Avg.log2FC", "adjusted P-value"),
                      seurat.marker.metric.description = 
                        c("Average log fold change", "Adjusted p-value (BF)")
)

close_loom(loom= loom)
}
    

    
    
runSoupx<- function(filt, rawh5, resolution, name){
    filt.matrix=Read10X(filt)
    raw.matrix=Read10X_h5(rawh5,use.names = T)

    soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

    srat=CreateSeuratObject(counts = filt.matrix)
    srat    <- SCTransform(srat, verbose = F)
    srat    <- RunPCA(srat, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
    srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T, resolution = 0.5)


    soup.channel  <- setClusters(soup.channel, setNames(srat@meta.data$seurat_clusters, rownames(srat@meta.data)))
    soup.channel  <- setDR(soup.channel, srat@reductions$umap@cell.embeddings)
    soup.channel  <- autoEstCont(soup.channel,forceAccept=TRUE)
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

    return(adj.matrix)

#    DropletUtils:::write10xCounts(paste0("soupX_",name,"_filt"), adj.matrix)    

    }
    
    



runSoupx_h5<- function(filt, rawh5, resolution, name){
    filt.matrix=Read10X_h5(filt)
    raw.matrix=Read10X_h5(rawh5,use.names = T)

    soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

    srat=CreateSeuratObject(counts = filt.matrix)
    srat    <- SCTransform(srat, verbose = F)
    srat    <- RunPCA(srat, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
    srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
    srat    <- FindClusters(srat, verbose = T, resolution = 0.5)


    soup.channel  <- setClusters(soup.channel, setNames(srat@meta.data$seurat_clusters, rownames(srat@meta.data)))
    soup.channel  <- setDR(soup.channel, srat@reductions$umap@cell.embeddings)
    soup.channel  <- autoEstCont(soup.channel,forceAccept=TRUE)
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

    return(adj.matrix)

#    DropletUtils:::write10xCounts(paste0("soupX_",name,"_filt"), adj.matrix)

    }

