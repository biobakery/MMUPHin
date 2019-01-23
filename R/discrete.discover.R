#' Main function for discrete structure discovery
#'
#' @param distance dissimilarity object
#' @param batch Name of batch variable
#' @param data data frame
#' @param k.max max k for testing
#' @param cluster.method method for clustering. consistent with fpc
#' @param classify.method method for classification. consistent with fpc
#' @param M number of folds validations during internal validation
#' @param nnk number of k nearest neighbors, if method is KNN
#' @param diagnostics logical. Should diagnostics be plotted?
#' @param verbose logical. Verbose?
#' @param ... additional parameters for clsutering. consistent with fpc
#' @return A list of MMUPHin continuous score discovery object
#' @export
discrete.discover <- function(distance,
                              batch,
                              data,
                              k.max = 10,
                              cluster.method = fpc::claraCBI, # must take distance and k as parameters
                              classify.method = "centroid",
                              M = 10,
                              nnk = 1,
                              diagnostics = TRUE,
                              verbose = TRUE,
                              ...) {
  ## Ensure data formatts are as expected
  if(!inherits(distance, "dist")) {
    stop("distance must be a dissimilarity matrix!")
  }
  dist.all <- as.matrix(distance)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if(!(batch%in% names(data)))
    stop("Batch/covariate variable not found in data.")
  ## Data dimensions need to agree with each other
  if(ncol(dist.all) != nrow(data))
    stop("Dimensions of distance matrix and metadata table do not agree!")

  ## Check that sample names agree between the feature and metadata table
  ## And assign row and column names if emppty
  if(is.null(rownames(data)))
    rownames(data) <- paste0("Sample",1:nrow(data))
  if(is.null(rownames(dist.all)))
    colnames(dist.all) <- rownames(dist.all) <- rownames(data)
  if(any(colnames(dist.all) != rownames(data)))
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

  # Clustering validation
  stats.internal <- replicate(k.max - 1, list(list()))
  stats.external <- replicate(k.max - 1, list(list()))
  for(k in 2:k.max) {
    if(verbose) message("Now evaluating clustering for k = ", k, "...")

    # Cluster each study first individually
    clusterings <- list()
    for(i in 1:n.batch) {
      i.batch <- lvl.batch[i]
      clusterings[[i]] <- cluster.method(
        data = as.dist(dist.all[batch == i.batch, batch == i.batch]),
        k = k
      )
    }

    # Validation
    classifications <- list()
    for(i in 1:n.batch) {
      i.batch <- lvl.batch[i]
      if(verbose) message("Study is ", i.batch)
      if(verbose) message("Performing internal validation...")
      i.result.internal <- fpc::prediction.strength(
        xdata = as.dist(dist.all[batch == i.batch, batch == i.batch]),
        Gmin = k,
        Gmax = k,
        M = M,
        clustermethod = cluster.method,
        classification = classify.method,
        nnk = nnk,
        distances = TRUE,
        count = FALSE
      )
      stats.internal[[k - 1]][[i]] <- i.result.internal$predcorr[[k]]

      if(verbose) message("Performing external validation...")
      i.clustering <- rep(-1, nrow(data))
      i.clustering[batch == i.batch] <- clusterings[[i]]$partition
      classifications[[i]] <- fpc::classifdist(as.dist(dist.all),
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
      stats.external[[k - 1]][[i]] <- i.result.external
    }
  }

  if(diagnostics) {
    data.internal <- Reduce("rbind",
                            lapply(1:(k.max - 1),
                                   function(i.k) {
                                     data.frame(k = i.k + 1,
                                                Reduce("rbind",
                                                       lapply(1:n.batch,
                                                              function(i)
                                                                data.frame(pred.strength = stats.internal[[i.k]][[i]],
                                                                           batch = lvl.batch[i]))))
                                   }))
    data.external <- Reduce("rbind",
                            lapply(1:(k.max - 1),
                                   function(i.k) {
                                     data.frame(k = i.k + 1,
                                                Reduce("rbind",
                                                       lapply(1:n.batch,
                                                              function(i)
                                                                data.frame(pred.strength = stats.external[[i.k]][[i]],
                                                                           batch = lvl.batch[i]))))
                                   }))
    data.plot <- rbind(data.frame(data.internal, validation = "Internal"),
                       data.frame(data.external, validation = "External"))
    print(ggplot2::ggplot(data.plot, ggplot2::aes(x = as.factor(k), y = pred.strength)) +
            ggplot2::geom_boxplot(ggplot2::aes(color = validation)) +
            ggplot2::facet_grid(.~batch) +
            ggplot2::theme_bw())
  }
  return(list(internal = stats.internal,
              external = stats.external))
}

