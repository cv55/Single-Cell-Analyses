
suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(ggplot2))
suppressWarnings(library(stringr))
suppressWarnings(library(ggrepel))
suppressWarnings(library(tidyverse))
#suppressWarnings(library(cowplot))
#suppressWarnings(library(harmony))
#suppressWarnings(library(DoubletFinder))
#suppressWarnings(library(SCopeLoomR))               
suppressWarnings(library(SoupX))
#suppressWarnings(library(DropletUtils))
suppressWarnings(source("scFunctions.R"))
#suppressWarnings(library(slingshot))
suppressWarnings(library("optparse"))
})



 
option_list = list(
make_option(c("--complete"), type="character",default=NA,
              help="Run the whole pipeline, douplets-soupx-QC-DR-UMAP, Set to True if you want this", metavar="character"),

make_option(c("--only-douplets"), type="character",default='False',dest='od',
              help="When usign the complete pipeline opt to run only the duplicate removal", metavar="character"),

make_option(c("--only-ambient"), type="character",default='False',dest='oa',
              help="When usign the complete pipeline opt to run only the ambient RNA removal", metavar="character"),

make_option(c("--typical"), type="character",default=NA,
              help="Run the typical pipeline, QC-DR-UMAP,Set to True if you want this", metavar="character"),
 
make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices. Available only if complete/typical has been set", metavar="character"),

make_option(c("--raw-matrix"), type="character", default=NA,dest='raw',
              help="Path with input raw matrices. Available only if complete/typical has been set", metavar="character"),

make_option(c("--QC"), type="character", default='100,10000,1000,50000,5', dest='QC',
              help="List of comma separated QC thresholds, in the order: min_feauture,max_feature,min_count, max_count, MT%", metavar="character"),


make_option(c("--norm-method"), type="character", default="SCT",dest='norm',
              help="Normalization method, either SCT or vst. Default is SCT", metavar="character"),



make_option(c("--douplets"), type="character", default=NA,
              help="Path with output clusters from supocell. Available only if complete has been set", metavar="character"),


make_option(c("--continue-from"), type="character", default=NA, dest='CF',
              help="Continue from specidic step, the output folder must exist", metavar="character"),

make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")



); 

parser=OptionParser(option_list=option_list) 
opt=parse_args(parser)

if( (!is.na(opt$complete)) & (opt$complete!='True')){opt$complete=NA}
if( (!is.na(opt$typical)) & (opt$typical!='True')){opt$typical=NA}
if(opt$od!='True'){opt$od='False'}
if(opt$oa!='True'){opt$oa='False'}

print(opt)



if(is.na(opt$complete) & is.na(opt$typical) & is.na(opt$CF)){stop("Choose one of complete/typical/continue-from options ")}		 
if(!is.na(opt$complete) & !is.na(opt$typical)){stop("Choose either complete or typical pipelines")}
if( ((!is.na(opt$complete)) | (!is.na(opt$typical))) & is.na(opt$matrix)){stop('Pipeline has been chosen but matrix path has not been given')}
if(!is.na(opt$complete) & is.na(opt$douplets) & (opt$oa!='True')){stop("Complete pipeline has been chosen but douplet folder is not provided")}
if(!is.na(opt$typical) & is.na(opt$matrix)){stop("Complete pipeline has been chosen but matrix folder is not provided")}
if(!is.na(opt$complete) & opt$od=='False' & is.na(opt$raw)){stop("Complete pipeline with ambient RNA removal has been chosen but raw matrix folder is not provided")}
if(opt$od=='True' & opt$oa=='True'){stop("You have set only-douplets and only-ambient both to true!.Delete them and use the --complete pipeline instead")}




path=(opt$matrix)
if(str_sub(path,-1,-1) =='/'){
path=substr(path,1,nchar(path)-1)
}

if(!is.na(opt$complete) & opt$oa=='False'){
dpath=(opt$douplets)
if(str_sub(dpath,-1,-1) =='/'){
dpath=substr(dpath,1,nchar(dpath)-1)
}
}


if(!is.na(opt$complete) & opt$od=='False'){
rawpath=(opt$raw)
print(rawpath)
if(str_sub(rawpath,-1,-1) =='/'){
rawpath=substr(rawpath,1,nchar(rawpath)-1)
}}



cd=getwd()
output=paste0('/', opt$out)
dir.create(paste0(cd, output), showWarnings = FALSE)
dir.create(paste0(cd, output, '/QC'), showWarnings = FALSE)
dir.create(paste0(cd, output, '/objects'), showWarnings = FALSE)
dir.create(paste0(cd, output, '/UMAPS'), showWarnings = FALSE)
gex_list=list()
gex_qc_list=list()
gex_norm_list=list()
gex_umap_list=list()



QC=strsplit(opt$QC, ',')
print(QC)



if(!is.na(opt$complete) | !is.na(opt$typical) ){

	for(file in list.files(path)) {
		print(file)
		name=paste(strsplit(file, '_')[[1]][1], strsplit(file, '_')[[1]][2], sep='_' )
		print(name)
		gex_list[[name]]<- Read10X_h5(paste0(path,'/',file))
}
}



if(!is.na(opt$typical)){
	
	gex_list= lapply(X=gex_list, FUN=function(x){
	x=CreateSeuratObject(counts=x, assay="RNA", min.cells=3) })

}




if(!is.na(opt$complete) ){

	if(opt$oa=='False'){
	    cd=getwd()
            dir.create(paste0(cd,output,'/noDouplets'), showWarnings = TRUE)

	
	for(file in list.files(dpath)){
		
	    douplets=grep('/', readLines(paste0(dpath,'/',file)), fixed = TRUE) #read.table(paste0(dpath,'/',file), header = F)
	    clusters=read.table(paste0(cd,'/',dpath,'/',file), header = F)
	    name=paste(strsplit(file, '_')[[1]][1], strsplit(file, '_')[[1]][2], sep='_' )
	    douplet_cluster=clusters[douplets,]	 
	    barcodes=douplet_cluster$V1

	    current <- gex_list[[name]][,!colnames(gex_list[[name]]) %in% barcodes]
	
	    setwd(paste0(cd,output,'/noDouplets'))	
	    DropletUtils:::write10xCounts(name, current)
	    setwd(cd)
	
 	if(opt$od=="True"){gex_list[[name]]<-CreateSeuratObject(counts=current, assay="RNA", min.cells=3)}
    }
}

	if(opt$od=='False'){
	cd=getwd()
	dir.create(paste0(cd,output,'/SoupX'), showWarnings = TRUE)

	res=c(0.5,0.5,0.5,0.5,0.5)

	n=0
	for(file in list.files(path, pattern='h5')){
	n=n+1	   
	   name=paste(strsplit(file, '_')[[1]][1], strsplit(file, '_')[[1]][2], sep='_' )
	

	  if(opt$oa=="False"){filtered=paste0(cd,output,'/noDouplets/', name)}
	  if(opt$oa=="True"){filtered=paste0(cd,'/',path,'/',name,'_filtered.h5')}

	  raw=paste0(cd,'/',rawpath,'/',name,'_raw.h5')
	  setwd(paste0(cd,output,'/SoupX'))

	  if(opt$oa=="False"){adj.matrix<-runSoupx(filtered, raw, res[n], name)}
       	  if(opt$oa=="True"){adj.matrix<-runSoupx_h5(filtered, raw, res[n], name)}

	  DropletUtils:::write10xCounts(name, adj.matrix)
	  setwd(cd)
}
}
	if((opt$oa=="True") | (opt$oa=="False" & opt$od=="False"))  { 
		gex_list=list()
		for(file in list.files(paste0(cd,output,'/SoupX'), pattern='MO' ) ){
		
			name=paste(strsplit(file, '_')[[1]][1], strsplit(file, '_')[[1]][2], sep='_' )
			counts<-Read10X(paste0(cd,output,'/SoupX/',file))
			gex_list[[name]]<-CreateSeuratObject(counts=counts, assay="RNA")}

}

		
}


print(gex_list)
QC=as.numeric(QC[[1]])
print(QC)
for(i in names(gex_list)){

  current<-seurat_qc(gex_list[[i]],name=i,folder=paste0(cd, output, '/QC/'), nFeature_RNA_min = QC[1], nFeature_RNA_max = QC[2], nCount_RNA_min = QC[3], nCount_RNA_max = QC[4], percent_ribo=100,percent_mito = QC[5])
  current@meta.data['Group']<- i
  rownames(current@meta.data)<- paste(i,rownames(current@meta.data), sep="_")
  current<- RenameCells(current, new.names=rownames(current@meta.data))
  gex_qc_list[[i]]<- current
  }

saveRDS(gex_qc_list, paste0(cd, output, '/objects/QC_list.rds.gz'), compress = "gzip")

print(gex_qc_list)


for(i in names(gex_qc_list)){
  current<- seurat_norm(sobject=gex_qc_list[[i]], method=opt$norm )
   gex_norm_list[[i]]<- current
}   
saveRDS(gex_norm_list, paste0(cd, output, '/objects/Norm_list.rds.gz'), compress = "gzip")



PCs=c(13,11,12,10,10)
res=c(0.5,0.5,0.5,0.5,0.5)
resolutions=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)
n=0
for(i in names(gex_norm_list)){
  n=n+1
  current<- make_UMAP(gex_norm_list[[i]],name=i,folder=paste0(cd, output, '/UMAPS/'), dims = PCs[n], resolution = resolutions, res_to_plot=res[n])
  gex_umap_list[[i]]<- current

}
saveRDS(gex_umap_list, paste0(cd, output, '/objects/UMAP_list.rds.gz'), compress = "gzip")
