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
require(scIntegrationMetrics)
require(Seurat)
require(anndata)

option_list <- list(make_option(c("-m", "--method"), type="character", default=NA, help="integration method to use"),
                    make_option(c("-i", "--integrated"), type="character", default=NA, help="input integrated data"),
                    #make_option(c("-u", "--uncorrected"), type="character", default=NA, help="input uncorrected data"),
                    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
                    make_option(c("-b", "--batch"), type="character", default=NA, help="batch variable"),
                    make_option(c("-l", "--label"), type="character", default=NA, help="cell type label"),
                    make_option(c("-n", "--nComp"), type="numeric", default=30, help="number of components to use in case type is full"),
                    make_option(c("-t", "--type"), type="character", default=NA, help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony\n knn is not supported'))



opt = parse_args(OptionParser(option_list=option_list))

print(opt)

if(opt$t == "full") {
  reductionKey <- "X_pca"
} else {
  reductionKey <- "X_emb"
}


obj = read_h5ad(opt$i)
print(obj)
if (length(obj$obsm)>0) {
sobj <- CreateSeuratObject(counts = Matrix::t(obj$X),meta.data = obj$obs[,c(opt$b,opt$l)],assay = "integrated")
sobj[[reductionKey]] <- CreateDimReducObject(embeddings = obj$obsm[[reductionKey]],key = "PC",assay = "integrated")


print("computing metrics using scIntegrationMetrics package...")


batchPerLabel <- apply(table(sobj@meta.data[,opt$l],sobj@meta.data[,opt$b]),MARGIN = 1,FUN = function(x) {length(which(x>0))})

metricsLabels <- names(batchPerLabel[batchPerLabel>1])

# use all dimensions (computed with a fixed number in the pipeline with metrics.py for full outputs or during integration for embed outputs) by default 
results <- getIntegrationMetrics(object=sobj,
                                 metrics = c( 'CiLISI', 'CiLISI_means'), 
                                 meta.label = opt$l, 
                                 meta.batch = opt$b, 
                                 method.reduction = reductionKey,
                                 metricsLabels = metricsLabels)
print(results)
metricTable <- t(data.frame(results))

colnames(metricTable) <- paste0(opt$m,"_",opt$t)

} else {
  print("Setting metrics to NA because integration failed (no integrated embedding")
  metricTable <- t(data.frame(CiLISI = c(NA),CiLISI_means = c(NA)))
  colnames(metricTable) <- paste0(opt$m,"_",opt$t)
  
  }

write.csv(metricTable,opt$o)
