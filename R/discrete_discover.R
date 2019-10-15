#' Unsupervised meta-analytical discovery and validation of discrete clustering 
#' structures in microbial abundance data
#' 
#' \code{discrete_discover} takes as input sample-by-sample dissimilarity 
#' measurements (generated from microbial abundance profiles), and performs 
#' unsupervised clustering within each batch across a range of cluster numbers. 
#' It then evaluates the support for each cluster number with both internal 
#' (i.e., samples within the batch) and external (i.e., samples in other 
#' batches) data. Internal evaluation is realized with 
#' \code{\link[fpc]{prediction.strength}} and external evaluation is based on
#' a generalized version of the same method. \code{discrete_discover} generates 
#' as output the evaluation statistics for each cluster number. A cluster number
#' with good support from both internal and external evaluations provides 
#' meta-analytical evidence for discrete structures in the microbial abundance 
#' profiles.
#'
#' \code{control} should be provided as a named list of the following components
#' (can be a subset).
#' \describe{
#' 
#' \item{k_max}{
#' integer. Maximum number of clusters to evaluate. \code{discrete_discover} 
#' will evaluate clustering structures corresponding to cluster numbers ranging
#' from 2 to \code{k_max}. Default to 10.
#' }
#' \item{cluster_function}{
#' an interface function. This function will be used for unsupervised clustering
#' for discrete structure evaluation. This corresponds to the 
#' \code{clustermethod} parameter in 
#' \code{\link[fpc]{prediction.strength}}, and similarly, should also follow the
#' specifications as detailed in  \code{\link[fpc]{clusterboot}}. Default to
#' \code{\link[fpc:kmeansCBI]{claraCBI}}
#' }
#' \item{classify_method}{
#' character. Classification method used to assign observations in the method's
#' internal and external evaluation stage. Corresponds to the 
#' \code{classification} parameter in \code{\link[fpc]{prediction.strength}}, 
#' and can only be either \code{"centroid"} or \code{"knn"}. Default to 
#' "centroid".
#' }
#' \item{M}{
#' integer. Number of random iterations to partition the batch during method's 
#' internal evaluation. Corresponds to the \code{M} parameter in 
#' \code{\link[fpc]{prediction.strength}}. Default to 30.
#' }
#' \item{nnk}{
#' integer. Numbber of nearest neighbors if \code{classify_method="knn"}. 
#' Corresponds to the \code{nnk} parameter in 
#' \code{\link[fpc]{prediction.strength}}. Default to 1.
#' }
#' \item{diagnostic_plot}{
#' character. Name for the generated diagnostic figure file. Default to 
#' \code{"discrete_diagnostic.pdf"}. Can be set to \code{NULL} in which 
#' case no output will be generated.}
#' \item{verbose}{
#' logical. Indicates whether or not verbose information will be printed.
#' }
#' }
#'
#' @param D sample-by-sample dissimilarity measurements. Should be provided as 
#' a \code{\link[stats]{dist}} object.
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param data data frame of metadata, columns must include batch.
#' @param control a named list of additional control parameters. See details.
#' @return a list, with the following components:
#' \describe{
#' \item{internal_mean, internal_se}{
#' matrices of internal clustering structure evaluation measurements 
#' (prediction strengths). Columns and rows corresponds to different batches and
#' different numbers of clusters, respectively. \code{internal_mean} and
#' \code{internal_se}, as the names suggest, are the mean and standard error of
#' prediction strengths for each batch/cluster number.
#' }
#' \item{external_mean, external_se}{
#' same structure as \code{internal_mean} and \code{internal_se}, but records
#' external clustering structure evaluation measurements (generalized prediction
#' strength).
#' }
#' \item{control}{list of additional control parameters used in the function
#' call.
#' }
#' }
#' @importFrom stats as.dist
#' @export
#' @author Siyuan Ma, \email{siyuanma@@g.harvard.edu}
#' @examples
#' data("CRC_abd", "CRC_meta")
#' # Calculate Bray-Curtis dissimilarity between the samples
#' library(vegan)
#' D <- vegdist(t(CRC_abd))
#' fit_discrete <- discrete_discover(D = D,
#'                                   batch = "studyID",
#'                                   data = CRC_meta)
discrete_discover <- function(D,
                              batch,
                              data,
                              control) {
  control <- match_control(default = control_discrete_discover,
                           control = control)
  verbose <- control$verbose
  
  # Check data formats
  # Check dissimilarity
  D_all <- check_D(D)
  # Check metadata data frame
  data <- as.data.frame(data)
  samples <- check_samples_D(D = D_all,
                             data = data)
  # Check batch is included in metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  df_batch <- check_metadata(data = data,
                             variables = batch)
  
  # Check batch variable
  df_batch[[batch]] <- check_batch(df_batch[[batch]], min_n_batch = 2)
  n_batch <- nlevels(x = df_batch[[batch]])
  lvl_batch <- levels(df_batch[[batch]])
  if(verbose)
    message("Found ", n_batch, " batches")
  
  # Check some control variables
  k_max <- check_k(control$k_max, df_batch[[batch]])
  classify_method <- check_options(control$classify_method,
                                   "classify_method",
                                   c("centroid", "knn"))
  
  # Clustering validation
  stats_internal <- replicate(k_max - 1, list(list()))
  stats_external <- replicate(k_max - 1, list(list()))
  for(k in 2:k_max) {
    if(verbose) message("Now evaluating clustering for k = ", k, "...")
    
    # Cluster each study first individually
    clusterings <- list()
    for(i in seq_len(n_batch)) {
      i_batch <- lvl_batch[i]
      clusterings[[i]] <- control$cluster_function(
        data = as.dist(D_all[df_batch[[batch]] == i_batch, 
                             df_batch[[batch]]== i_batch]),
        k = k)}
    
    # Validation
    classifications <- list()
    for(i in seq_len(n_batch)) {
      i_batch <- lvl_batch[i]
      if(verbose) message("Study is ", i_batch)
      if(verbose) message("Performing internal validation...")
      i_result_internal <- fpc::prediction.strength(
        xdata = as.dist(D_all[df_batch[[batch]] == i_batch, 
                              df_batch[[batch]] == i_batch]),
        Gmin = k,
        Gmax = k,
        clustermethod = control$cluster_function,
        classification = classify_method,
        M = control$M,
        nnk = control$nnk,
        distances = TRUE,
        count = FALSE)
      stats_internal[[k - 1]][[i]] <- 
        c("mean" = mean(i_result_internal$predcorr[[k]]),
          "sd" = sd(i_result_internal$predcorr[[k]]))
      
      if(verbose) message("Performing external validation...")
      i_clustering <- rep(-1, nrow(data))
      i_clustering[df_batch[[batch]] == i_batch] <- clusterings[[i]]$partition
      classifications[[i]] <- 
        fpc::classifdist(as.dist(D_all),
                         clustering = i_clustering,
                         method = classify_method,
                         centroids = clusterings[[i]]$result$medoids,
                         nnk = control$nnk)
      
      i_result_external <- matrix(NA, n_batch, k)
      for (j in setdiff(seq_len(n_batch), i)) {
        j_batch <- lvl_batch[[j]]
        ctable <- fpc::xtable(
          clusterings[[j]]$partition,
          classifications[[i]][df_batch[[batch]] == j_batch], k)
        for (kk in seq_len(k)) {
          i_result_external[j, kk] <- sum(ctable[kk, ]^2 - ctable[kk,])
          cpjk <- clusterings[[j]]$partition == kk
          njk <- sum(cpjk)
          if (njk > 1)
            i_result_external[j, kk] <- 
            i_result_external[j, kk]/(njk * (njk - 1))
          else i_result_external[j, kk] <- 1
        }
      }
      i_result_external <- apply(i_result_external[setdiff(seq_len(n_batch), i), ,
                                                   drop = FALSE],
                                 1, min)
      stats_external[[k - 1]][[i]] <- c(mean = mean(i_result_external),
                                        sd = sd(i_result_external))
    }
  }
  
  # If required, generate diagnostic plots
  if(!is.null(control$diagnostic_plot)) {
    diagnostic_discrete_discover(stats_internal = stats_internal,
                                 stats_external = stats_external,
                                 lvl_batch = shorten_name(lvl_batch),
                                 output = control$diagnostic_plot)
  }
  
  # compile results for output
  internal_mean <- t(vapply(stats_internal,
                            function(k_stats)
                              vapply(k_stats, 
                                     function(i.stat) i.stat["mean"],
                                     0.0),
                            rep_len(0.0, n_batch)))
  internal_se <- t(vapply(stats_internal,
                          function(k_stats)
                            vapply(k_stats, 
                                   function(i.stat) i.stat["sd"],
                                   0.0),
                          rep_len(0.0, n_batch)))
  external_mean <- t(vapply(stats_external,
                            function(k_stats)
                              vapply(k_stats, 
                                     function(i.stat) i.stat["mean"],
                                     0.0),
                            rep_len(0.0, n_batch)))
  external_se <- t(vapply(stats_external,
                          function(k_stats)
                            vapply(k_stats, 
                                   function(i.stat) i.stat["sd"],
                                   0.0),
                          rep_len(0.0, n_batch)))
  dimnames(internal_mean) <-
    dimnames(internal_se) <-
    dimnames(external_mean) <-
    dimnames(external_se) <- list(as.character(2:k_max),
                                  lvl_batch)
  
  
  return(list(internal_mean = internal_mean,
              internal_se = internal_se,
              external_mean = external_mean,
              external_se = external_se,
              control = control))
}

