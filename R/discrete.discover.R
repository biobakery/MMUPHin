#' Internal evaluation and external validation of discrete clustering structures
#'
#' @param D a dissimilarity object.
#' @param batch name of the batch variable.
#' @param data data frame of metadata, must contain batch.
#' @param k.max max number of clusters to test for
#' @param cluster.function clustering function. Should be consistent with fpc, i.e.,
#' @param internal.metric internal evaluation metric for clustering structures. Currently only
#' "prediction.strength" and "silhouette.width" are supported.
#' @param classify.method method for classification (used for external validation and if
#' internal.metric="prediction.strength"). Currently only "centroid" and "knn" are supported.
#' @param M number of validation folds (for if internal.metric="prediction.strength").
#' @param nnk number of k nearest neighbors (for if classify.method="knn").
#' @param diagnostics should diagnostic plots be generated? Deafault to FALSE.
#' @param verbose should verbose modelling information should be printed? Default to TRUE.
#'
#' @return a list of four summary statistics from internal evaluation and external validation of
#' clustering structures in the data. These are means and standard errors of
#' internal and external metrics for each number of cluster-batch combination. So each component
#' (mean/standard error, internal/external) is a (k.max - 1)*(number of batches) numeric matrix.
#' @export
discrete.discover <- function(D,
                              batch,
                              data,
                              k.max = 10,
                              cluster.function = fpc::claraCBI,
                              internal.metric = "prediction.strength",
                              classify.method = "centroid",
                              M = 30,
                              nnk = 1,
                              diagnostics = FALSE,
                              verbose = TRUE) {
  # Ensure data formatts are as expected
  if(!inherits(D, "dist")) {
    stop("D must be a dissimilarity matrix!")
  }
  D.all <- as.matrix(D)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if(!(batch%in% names(data)))
    stop("Batch/covariate variable not found in data.")
  # Data dimensions need to agree with each other
  if(ncol(D.all) != nrow(data))
    stop("Dimensions of D matrix and metadata table do not agree!")

  # Check that sample names agree between the feature and metadata table
  # And assign row and column names if emppty
  if(is.null(rownames(data)))
    rownames(data) <- paste0("Sample",1:nrow(data))
  if(is.null(rownames(D.all)))
    colnames(D.all) <- rownames(D.all) <- rownames(data)
  if(any(colnames(D.all) != rownames(data)))
    stop("Sample names in feature.count and data don't agree!")

  # Check batch variable and identify groups
  batch <- data[, batch]
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  batch <- as.factor(batch)
  n.batch <- length(unique(batch))
  if(n.batch < 3) stop("Only discovery based on three or more studies is supported!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- levels(batch)

  # Number of max clusters to evaluate can't exceed half of smallest sample size
  if(k.max > floor(min(table(batch)) / 2) - 1) {
    message("k.max be less than half smallest sample size!")
    k.max <- floor(min(table(batch)) / 2) - 1
    message("Set k.max to ", k.max)
  }

  # Check that evaluation metric and classification method are supported
  if(!(internal.metric %in% c("prediction.strength", "silhouette.width")))
    stop("Only prediction.strength and silhouette.width metrics are supported!")
  if(!(classify.method %in% c("centroid", "knn")))
    stop("Only centroid and knn classification methods are supported!")

  # Clustering validation
  stats.internal <- replicate(k.max - 1, list(list()))
  stats.external <- replicate(k.max - 1, list(list()))
  for(k in 2:k.max) {
    if(verbose) message("Now evaluating clustering for k = ", k, "...")

    # Cluster each study first individually
    clusterings <- list()
    for(i in 1:n.batch) {
      i.batch <- lvl.batch[i]
      clusterings[[i]] <- cluster.function(
        data = as.dist(D.all[batch == i.batch, batch == i.batch]),
        k = k
      )
    }

    # Validation
    classifications <- list()
    for(i in 1:n.batch) {
      i.batch <- lvl.batch[i]
      if(verbose) message("Study is ", i.batch)
      if(verbose) message("Performing internal validation...")
      if(internal.metric == "prediction.strength") {
        i.result.internal <- fpc::prediction.strength(
          xdata = as.dist(D.all[batch == i.batch, batch == i.batch]),
          Gmin = k,
          Gmax = k,
          clustermethod = cluster.function,
          classification = classify.method,
          M = M,
          nnk = nnk,
          distances = TRUE,
          count = FALSE
        )
        stats.internal[[k - 1]][[i]] <- c("mean" = mean(i.result.internal$predcorr[[k]]),
                                          "sd" = sd(i.result.internal$predcorr[[k]]))
      }
      if(internal.metric == "silhouette.width") {
        sil.width <- cluster::silhouette(
          clusterings[[i]],
          dist = as.dist(D.all[batch == i.batch, batch == i.batch]))[, "sil_width"]
        stats.internal[[k - 1]][[i]] <- c("mean" = mean(sil.width),
                                          "sd" = sd(sil.width))
      }

      if(verbose) message("Performing external validation...")
      i.clustering <- rep(-1, nrow(data))
      i.clustering[batch == i.batch] <- clusterings[[i]]$partition
      classifications[[i]] <- fpc::classifdist(as.dist(D.all),
                                               clustering = i.clustering,
                                               method = classify.method,
                                               centroids = clusterings[[i]]$result$medoids,
                                               nnk = nnk)

      i.result.external <- matrix(NA, n.batch, k)
      for (j in setdiff(1:n.batch, i)) {
        j.batch <- lvl.batch[[j]]
        ctable <- fpc::xtable(clusterings[[j]]$partition,
                              classifications[[i]][batch == j.batch], k)
        for (kk in 1:k) {
          i.result.external[j, kk] <- sum(ctable[kk, ]^2 - ctable[kk,])
          cpjk <- clusterings[[j]]$partition == kk
          njk <- sum(cpjk)
          if (njk > 1)
            i.result.external[j, kk] <- i.result.external[j, kk]/(njk * (njk - 1))
          else i.result.external[j, kk] <- 1
        }
      }
      i.result.external <- apply(i.result.external[setdiff(1:n.batch, i), ],
                                 1, min)
      stats.external[[k - 1]][[i]] <- c(mean = mean(i.result.external),
                                        sd = sd(i.result.external))
    }
  }

  # If required, generate diagnostic plots
  if(diagnostics) {
    diagnostics.discrete.discover(stats.internal = stats.internal,
                                  stats.external = stats.external,
                                  lvl.batch = lvl.batch)
  }

  # compile results for output
  internal.mean <- t(sapply(stats.internal,
                            sapply(k.stats, function(i.stat) i.stat["mean"])))
  internal.se <- t(sapply(stats.internal,
                          sapply(k.stats, function(i.stat) i.stat["sd"])))
  external.mean <- t(sapply(stats.external,
                            sapply(k.stats, function(i.stat) i.stat["mean"])))
  external.se <- t(sapply(stats.external,
                          sapply(k.stats, function(i.stat) i.stat["sd"])))
  dimnames(internal.mean) <-
    dimnames(internal.se) <-
    dimnames(external.mean) <-
    dimnames(external.se) <- list(as.character(2:k.max),
                                  lvl.batch)


  return(list(internal.mean = internal.mean,
              internal.se = internal.se,
              external.mean = external.mean,
              external.se = external.se))
}

