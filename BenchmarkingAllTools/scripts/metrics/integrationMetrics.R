

compute_silhouette <- function (X, meta_data, label_colnames) 
{
  N <- nrow(meta_data)
  sil_df <- data.frame(matrix(NA, N, length(label_colnames)))
  sil_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message(paste("ASW: Cannot compute silhouette on missing values.", 
                    "Skipping", label_colname))
      return(rep(NA, N))
    }
    else if (sum(table(labels) > 0) < 2) {
      message(paste("ASW: Cannot compute silhouette without at least 2 label levels.", 
                    "Skipping", label_colname))
      return(rep(NA, N))
    }
    else {
      dists <- vegan::vegdist(X, method = "euclidean")
      labels.num <- as.numeric(as.factor(labels))
      #sil <- as.data.frame(silhouette(labels.num, dists))
      sil <-cluster::silhouette(labels.num, dists)
      print(sil)
      return(sil[,"sil_width"])
    }
  }))
  sil_df <- as.data.frame(sil_df)
  colnames(sil_df) <- label_colnames
  row.names(sil_df) <- row.names(meta_data)
  return(sil_df)
}
# 
getIntegrationMetrics <- function (object, metrics = NULL, meta.label, meta.batch, lisi_perplexity = 30,
                                   method.reduction = "pca", metricsLabels = NULL)
{
  metricsAvailable <- c("iLISI", "norm_iLISI", "CiLISI", "CiLISI_means",
                        "norm_cLISI", "norm_cLISI_means", "celltype_ASW", "celltype_ASW_means")
  if (is.null(metrics)) {
    metrics <- metricsAvailable
  }
  if (!all(metrics %in% metricsAvailable)) {
    metrics_vector <- paste(c(metricsAvailable), collapse = ",")
    str <- sprintf("'metrics' is unknown. Please define one or more of %s, or 'metrics=NULL' to calculate them all",
                   metrics_vector)
    stop(str)
  }
  if (!is(object, "Seurat")) {
    stop("Error: 'object' must be a 'Seurat' object.")
  }
  integrationMetrics <- list()
  if (is.null(metricsLabels))
    metricsLabels <- unique(object@meta.data[[meta.label]])
  message(paste("Cell type labels:", paste(metricsLabels, collapse = ",")))
  metricsLabels_logic <- object@meta.data[[meta.label]] %in%
    metricsLabels
  batchNames <- unique(object@meta.data[[meta.batch]])
  message(paste("Batches:", paste(batchNames, collapse = ",")))
  if (any(c("iLISI", "norm_iLISI") %in% metrics)) {
    lisi.this <- compute_lisi(object@reductions[[method.reduction]]@cell.embeddings,
                              meta_data = object@meta.data, label_colnames = meta.batch,
                              perplexity = lisi_perplexity)[[1]]
    if ("iLISI" %in% metrics) {
      integrationMetrics[["iLISI"]] <- mean(lisi.this[metricsLabels_logic])
    }
    if ("norm_iLISI" %in% metrics) {
      lisi.this.normalized <- (lisi.this - 1)/(length(batchNames) -
                                                 1)
      integrationMetrics[["norm_iLISI"]] <- mean(lisi.this.normalized[metricsLabels_logic])
    }
  }
  if (any(c("CiLISI", "CiLISI_means") %in% metrics)) {
    lisi_splitByCelltype <- compute_lisi_splitBy(object@reductions[[method.reduction]]@cell.embeddings,
                                                 meta_data = object@meta.data, label_colnames = meta.batch,
                                                 perplexity = lisi_perplexity, split_by_colname = meta.label,
                                                 normalize = T)
    if ("CiLISI" %in% metrics) {
      integrationMetrics[["CiLISI"]] <- mean(unlist(lisi_splitByCelltype)[metricsLabels_logic])
    }
    if ("CiLISI_means" %in% metrics) {
      classMeans <- sapply(lisi_splitByCelltype, function(x) mean(x[,
                                                                    1]))[metricsLabels]
      message("CiLISI: ", paste(names(lisi_splitByCelltype),
                                round(classMeans, 2), " "))
      integrationMetrics[["CiLISI_means"]] <- mean(classMeans)
    }
  }
  if (any(c("norm_cLISI", "norm_cLISI_means") %in% metrics)) {
    lisi.this <- compute_lisi(object@reductions[[method.reduction]]@cell.embeddings,
                              meta_data = object@meta.data, label_colnames = meta.label,
                              perplexity = lisi_perplexity)[[1]]
    lisi.this.normalized <- (lisi.this - 1)/(length(metricsLabels) -
                                               1)
    if ("norm_cLISI" %in% metrics) {
      integrationMetrics[["norm_cLISI"]] <- 1 - mean(lisi.this.normalized[metricsLabels_logic])
    }
    if ("norm_cLISI_means" %in% metrics) {
      integrationMetrics[["norm_cLISI_means"]] <- 1 - mean(tapply(lisi.this.normalized,
                                                                  object@meta.data[[meta.label]], mean)[metricsLabels])
    }
  }
  if (any(c("celltype_ASW", "celltype_ASW_means") %in% metrics)) {
    sil.this <- compute_silhouette(object@reductions[[method.reduction]]@cell.embeddings,
                                   meta_data = object@meta.data, label_colnames = meta.label)[[1]]
    if ("celltype_ASW" %in% metrics) {
      integrationMetrics[["celltype_ASW"]] <- mean(sil.this[metricsLabels_logic])
    }
    if ("celltype_ASW_means" %in% metrics) {
      integrationMetrics[["celltype_ASW_means"]] <- mean(tapply(sil.this,
                                                                object@meta.data[[meta.label]], mean)[metricsLabels])
    }
  }
  return(integrationMetrics)
}