#' #' Wrapper function for unsupervised clustering strength evaluation
#' #'
#' #' @param physeq phyloseq object with otu table and metadata containing variables
#' #' to adjust for in normalizing
#' #' @param clust_method clustering method to be used
#' #' @param dist_measure distance metric to be used
#' #' @param kmax max of number of clusters (groups) to be used
#' #' @param plots whether or not visualization plots should be generated
#' #'
#' #' @return a data frame summarizing clustering strength evaluations
#' #' @export
#' #' @import phyloseq fpc cluster ggplot2
#' clust_strength <- function(physeq, clust_method = "pam", dist_measure = "bray", kmax, plots) {
#'   dist_clust <- phyloseq::distance(physeq, method = dist_measure, type = "samples")
#'   mat_ord <- phyloseq::ordinate(physeq, method = "MDS", distance = dist_clust)$vectors
#'   if(clust_method == "pam") {
#'     x_clust <- dist_clust
#'   } else x_clust <- mat_ord
#'
#'   if (clust_method == "hclust") {
#'     clust_fun_gap <- function(x, k) {
#'       list(cluster = fpc::hclustCBI(x, k, method = "complete", scaling = F)$partition)
#'     }
#'     clust_fun_ps <- function(data, k) fpc::hclustCBI(data, k, method = "complete")
#'   } else if (clust_method == "kmeans") {
#'     clust_fun_gap <- function(x, k) {
#'       list(cluster = fpc::kmeansCBI(x, k, scaling = F)$result$cluster)
#'     }
#'     clust_fun_ps <- fpc::kmeansCBI
#'   } else if (clust_method == "pam") {
#'     clust_fun_gap <- function(x, k) {
#'       list(cluster = fpc::claraCBI(x, k)$result$cluster)
#'     }
#'     clust_fun_ps <- fpc::claraCBI
#'   }
#'
#'   # Gap statistic
#'   result_gap <- cluster::clusGap(mat_ord, FUNcluster = clust_fun_gap, K.max = kmax, B = 30)$Tab
#'   result_gap <- data.frame(
#'     K = 2:kmax,
#'     stat = result_gap[-1, "gap"],
#'     se = result_gap[-1, "SE.sim"],
#'     measure = "Gap"
#'   )
#'
#'   # Prediction strength
#'   result_ps <- fpc::prediction.strength(x_clust,
#'                                         Gmax = kmax,
#'                                         M = 30,
#'                                         clustermethod = clust_fun_ps)
#'   result_ps <- data.frame(
#'     K = 2:kmax,
#'     stat = result_ps$mean.pred[-1],
#'     se = sapply(result_ps$predcorr[2:kmax], sd),
#'     measure = "Pred.Str"
#'   )
#'
#'   # Silhouette width
#'   result_sw <- data.frame(t(sapply( 2:kmax, function( k ) {
#'     clst_tmp <- clust_fun_gap(x=x_clust, k=k)$cluster
#'     sw_tmp <- cluster::silhouette(clst_tmp, dist=dist_clust)[, 'sil_width']
#'     return(c(mean(sw_tmp), sd(sw_tmp)))
#'   })))
#'   names(result_sw) <- c("stat", "se")
#'   result_sw <- data.frame(K = 2:kmax,
#'                           result_sw,
#'                           measure = "Sil.Width")
#'
#'   result <- rbind(result_gap,
#'                   result_ps,
#'                   result_sw)
#'   if(plots) print(
#'     ggplot(result, aes(x = K, y = stat)) +
#'       geom_point() +
#'       geom_line() +
#'       geom_errorbar(aes(ymax=stat+se, ymin=stat-se ), width=0.6) +
#'       facet_wrap(~ measure, scales = 'free_y') +
#'       scale_x_continuous(breaks = 2:kmax) +
#'       theme_bw( )
#'   )
#'   return()
#' }
