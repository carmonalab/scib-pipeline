# Change directory to script dir
getScriptPath <- function(){
	    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
        script.dir <- dirname(regmatches(cmd.args, m))
        if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
	    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
	    return(script.dir)
}

setwd(getScriptPath())

library('optparse')
library(rlang)
require(Seurat)

option_list <- list(make_option(c("-m", "--method"), type="character", default=NA, help="integration method to use"),
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
		    make_option(c("-b", "--batch"), type="character", default=NA, help="batch variable"),
		    make_option(c("-n", "--ndims"), type="integer", default=30, help="number of dimensions to use for integration"),
		    make_option(c("-c", "--celltype"), type="character", default=NA, help="cell type variable"),
		    make_option(c("-v", "--hvg"), type="character", default=NA, help="hvg list for seurat"))



opt = parse_args(OptionParser(option_list=option_list))

print(opt)

source('integration.R')
sobj = loadSeuratObject(opt$i)

if(opt$method=='seurat'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}

	out = runSeurat(sobj, opt$b, opt$n, hvg)
}

if(opt$method=='STACAS'){
  if(!is.na(opt$hvg)) {
    hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
  }
  else {
    hvg <- rownames(sobj@assays$RNA)
  }
  
  out = runSTACAS(sobj, opt$b, opt$n, hvg)
}

if(opt$method=='ownFeatSTACAS'){
  if(!is.na(opt$hvg)) {
    hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
  }
  else {
    hvg <- rownames(sobj@assays$RNA)
  }
  
  out = runOwnFeaturesSTACAS(sobj, opt$b, opt$n, hvg)
}

if(opt$method=='semiSupSTACAS'){
  if(!is.na(opt$hvg)) {
    hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
  }
  else {
    hvg <- rownames(sobj@assays$RNA)
  }
  
  out = runSemiSupSTACAS(sobj, opt$b, opt$n, hvg,opt$c)
}

if(opt$method=='ownFeatSemiSupSTACAS'){
  if(!is.na(opt$hvg)) {
    hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
  }
  else {
    hvg <- rownames(sobj@assays$RNA)
  }
  
  out = runOwnFeaturesSemiSupSTACAS(sobj, opt$b, opt$n, hvg,opt$c)
}

if(opt$method=='seuratrpca'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}
	out = runSeuratRPCA(sobj, opt$b, opt$n, hvg)
}


if(opt$method=='conos'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
		sobj <- subset(sobj, features=hvg)
	}
	out = runConos(sobj, opt$b, opt$n)
}

if(opt$method=='harmony'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
		sobj <- subset(sobj, features=hvg)
	}

	out=runHarm(sobj, opt$b, opt$n)
}

if(opt$method=='liger'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}

	out = runLiger(sobj, opt$b, opt$n, hvg)
}

if(opt$method=='fastmnn'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
		sobj <- subset(sobj, features=hvg)
	}

	out=runFastMNN(sobj, opt$b, opt$n)
}

saveSeuratObject(out, opt$o)
