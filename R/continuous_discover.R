#' Unsupervised meta-analytical discovery and validation of continuous 
#' structures in microbial abundance data
#' 
#' \code{continuous_discover} takes as input a feature-by-sample matrix of 
#' microbial abundances. It first performs unsupervised continuous structure
#' discovery (PCA) within each batch. Loadings of top PCs from each batch are
#' then mapped against each other to identify "consensus" loadings that are 
#' reproducible across batches with a network community discovery approach with
#' \pkg{igraph}. The identified consensus loadings/scores can be viewed as 
#' continuous structures in microbial profiles that are recurrent across batches
#' and valid in a meta-analyitical sense. \code{continuous_discover} returns, 
#' among other output, the identified consensus scores for continuous 
#' structures in the provided microbial abundance profiles, as well as the 
#' consensus PC loadings which can be used to assign continuous scores to any 
#' sample with the same set of microbial features.
#'
#' \code{control} should be provided as a named list of the following components
#' (can be a subset).
#' \describe{
#' \item{normalization}{
#' character. Similar to the \code{normalization} parameter in 
#' \code{\link[Maaslin2]{Maaslin2}} but only \code{"TSS"} and \code{"NONE"} are 
#' allowed. Default to \code{"TSS"} (total sum scaling).
#' }
#' \item{transform}{
#' character. Similar to the \code{transform} parameter in
#' \code{\link[Maaslin2]{Maaslin2}} but only \code{"AST"} and \code{"LOG"} are 
#' allowed. Default to \code{"AST"} (arcsine square root transformation).
#' }
#' \item{pseudo_count}{
#' numeric. Pseudo count to add feature_abd before the transformation. Default 
#' to \code{NULL}, in which case pseudo count will be set automatically to 0 if
#' \code{transform="AST"}, and half of minimal non-zero values in 
#' \code{feature_abd} if \code{transform="LOG"}.
#' }
#' \item{var_perc_cutoff}{
#' numeric. A value between 0 and 1 that indicates the percentage variability 
#' explained to cut off at for selecting top PCs in each batch. Across batches, 
#' the top PCs that in total explain more than var_perc_cutoff of the total 
#' variability will be selected for meta-analytical continuous structure 
#' discovery. Default to 0.8 (PCs included need to explain at least 80% of the 
#' total variability).
#' }
#' \item{cos_cutoff}{
#' numeric. A value between 0 and 1 that indicates cutoff for absolute cosine
#' coefficients between PC loadings to construct the method's network with. Once
#' the top PC loadings from each batch are selected, cosine coefficients between
#' each loading pair are calculated which indicate their similarity. Loading 
#' pairs with absolute cosine coefficients surpassing cos_cutoff are then 
#' considered as associated with each other, and represented as an edge 
#' between the pair in a PC loading network. Network community discovery can 
#' then be performed on this network to identified densely connected "clusters"
#' of PC loadings, which represent meta-analytically recurrent continuous 
#' structures.
#' }
#' \item{cluster_function}{
#' function. \code{cluster_function} is used to perform community structure 
#' discovery in the constructed PC loading network. This can be any of the 
#' network cluster functions provided in \pkg{igraph}. Default to 
#' \code{\link[igraph]{cluster_optimal}}. Note that this option can be slow for
#' larger datasets, in which case \code{\link[igraph]{cluster_fast_greedy}} is
#' recommended.
#' }
#' \item{network_plot}{
#' character. Name for the generated network figure file. Default to 
#' \code{"clustered_network.pdf"}. Can be set to \code{NULL} in which 
#' case no output will be generated.
#' }
#' \item{plot_size_cutoff}{
#' integer. Clusters with sizes smaller than or equal to plot_size_cutoff will 
#' be excluded in the visualized network. Defaul to 2 - visualized clusters 
#' must have at least three nodes (PC loadings).
#' }
#' \item{diagnostic_plot}{
#' character. Name for the generated diagnostic figure file. Default to 
#' \code{"continuous_diagnostic.pdf"}. Can be set to \code{NULL} in which 
#' case no output will be generated.
#' }
#' \item{verbose}{
#' logical. Indicates whether or not verbose information will be printed.
#' }
#' }
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param data data frame of metadata, columns must include batch.
#' @param control a named list of additional control parameters. See details.
#' @return a list, with the following components:
#' \describe{
#' \item{consensus_scores}{
#' matrix of identified consensus continuous scores. Columns are the identified 
#' consensus scores and rows correspond to samples in feature_abd.
#' }
#' \item{consensus_loadings}{
#' matrix of identified consensus loadings. Columns are the identified 
#' consensus scores and rows correspond to features in feature_abd.
#' }
#' \item{mat_vali}{
#' matrix of validation cosine coefficients of the identified consensus 
#' loadings. Columns correspond to the identified consensus scores and rows 
#' correspond to batches.
#' }
#' \item{network, communities, mat_cos}{
#' components for the constructed PC loading network and community
#' discovery results. \code{network} is a \pkg{igraph} \code{graph}
#' object for
#' the constructed network of associated PC loadings. \code{communities} is a  
#' \code{\link[igraph:communities]{communities}} object for the identified 
#' consensus loading clusters in \code{network} (output from 
#' \code{control$cluster_function}). \code{mat_cos} is the matrix of cosine 
#' coefficients between all selected top PCs from all batches.
#' }
#' \item{control}{list of additional control parameters used in the function
#' call.
#' }
#' }
#' @importFrom stats prcomp
#' @export
#' @author Siyuan Ma, \email{siyuanma@@g.harvard.edu}
#' @examples
#' data("CRC_abd", "CRC_meta")
#' fit_continuous <- continuous_discover(feature_abd = CRC_abd,
#'                                       batch = "studyID",
#'                                       data = CRC_meta)
continuous_discover <- function(feature_abd,
                                batch,
                                data,
                                control) {
  # Check and construct controls
  control <- match_control(default = control_continuous_discover,
                           control = control)
  verbose <- control$verbose
  
  # Check control values
  normalization <- check_options(control$normalization, 
                                 "normalization",
                                 c("NONE", "TSS"))
  transform <- check_options(control$transform, 
                             "transform",
                             c("LOG", "AST"))
  var_perc_cutoff <- check_options_continuous(control$var_perc_cutoff,
                                              "var_perc_cutoff",
                                              c(0, 1))
  cos_cutoff <- check_options_continuous(control$cos_cutoff,
                                         "cos_cutoff",
                                         c(0, 1))
  plot_size_cutoff <- check_options_continuous(control$plot_size_cutoff,
                                               "plot_size_cutoff",
                                               c(0, Inf))
  
  # Check data formats
  # Check feature abundance table
  feature_abd <- as.matrix(feature_abd)
  type_feature_abd <- check_feature_abd(feature_abd = feature_abd)
  if(verbose)
    message("feature_abd is ", type_feature_abd)
  # Check metadata data frame
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd = feature_abd,
                           data = data)
  # Check batch and covariates are included in metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  df_batch <- check_metadata(data = data,
                             variables = batch)
  
  # Check batch variable
  var_batch <- check_batch(df_batch[[batch]], min_n_batch = 3)
  n_batch <- nlevels(x = var_batch)
  lvl_batch <- levels(var_batch)
  if(verbose)
    message("Found ", n_batch, " batches")
  
  # Normalize and transform feature_abd
  if(transform == "LOG") {
    if(is.null(control$pseudo_count)) {
      pseudo_count <- set_pseudo(features = feature_abd)
      if(verbose)
        message("Pseudo count is not specified and set to half of minimal ",
                "non-zero value: ",
                format(pseudo_count, digits = 3, scientific = TRUE))
    } else 
      pseudo_count <- check_pseudo_count(control$pseudo_count)
    feature_pca <- transform_features(
      features = normalize_features(
        features = feature_abd,
        normalization = normalization,
        pseudo_count = pseudo_count),
      transform = transform)
  } 
  if(transform == "AST") {
    if(is.null(control$pseudo_count)) 
      pseudo_count <- 0
    else 
      pseudo_count <- check_pseudo_count(control$pseudo_count)
    feature_pca <- transform_features(
      features = normalize_features(
        features = feature_abd,
        normalization = normalization,
        pseudo_count = pseudo_count),
      transform = transform
    )
  }
  
  # Calculate PC for the training sets
  if(verbose) message("Performing PCA in individual datasets...")
  pca_all <- lapply(lvl_batch, function(lvl)
  {
    pc <- feature_pca[, var_batch == lvl]
    dat_pca <- prcomp(t(pc))
    return(dat_pca)
  })
  
  # Specify the smallest number of PCs to include
  ind_var_perc <- lapply(pca_all, function(dat_pca) {
    (cumsum(dat_pca$sdev^2) / sum(dat_pca$sdev^2)) > var_perc_cutoff
  })
  n_pc_top <- min(vapply(ind_var_perc, length, 1))
  ind_var_perc <- vapply(ind_var_perc, 
                         function(x) x[seq_len(n_pc_top)], 
                         rep_len(TRUE, n_pc_top))
  ind_var_perc <- apply(ind_var_perc, 1, all)
  if(any(ind_var_perc)) {
    n_pc_top <- min(seq_len(n_pc_top)[ind_var_perc])
    if(verbose) 
      message("Smallest number of PCs passing the specified threshold (",
              var_perc_cutoff, ") is ", n_pc_top)
  } else
    warning("No number of PCs passed the threshold!\n",
            "Setting number of PCs per dataset to largest possible",
            " value (",n_pc_top,")")
  
  # First n_pc_top PC loadings for each training dataset
  data_loadings <- lapply(pca_all, function(x)
  {
    loadings <- x$rotation[,seq_len(n_pc_top)]
    return(loadings)
  })
  mat_data_loading <- Reduce("cbind", data_loadings)
  colnames(mat_data_loading) <- 
    unlist(lapply(lvl_batch, function(lvl) {
      paste0(lvl, ", PC", seq_len(n_pc_top))
    }))
  
  # Calculate correlation matrix of loadings across datasets
  if(verbose) message("Calculating correlation between PCs across datasets...")
  mat_cos <- t(mat_data_loading) %*% mat_data_loading
  mat_edge <- matrix(1, nrow = nrow(mat_cos), ncol = ncol(mat_cos))
  dimnames(mat_edge) <- dimnames(mat_cos)
  mat_edge[abs(mat_cos) < cos_cutoff] <- 0
  if(sum(mat_edge) == nrow(mat_edge)) {
    warning("All edges are filtered out in the PC network!\n",
            "Consider lowering the value of cos_cutoff.")
    return(NULL)
  }
  
  # Create igraph graph object
  graph_pc <- igraph::graph_from_adjacency_matrix(mat_edge,
                                                  mode = "undirected",
                                                  weighted = NULL,
                                                  diag = FALSE)
  # Perform graph community detection
  if(verbose) message("Performing network clustering...")
  cluster_pc <- control$cluster_function(graph_pc)
  size_communities <- igraph::sizes(cluster_pc)
  size_communities <- sort(size_communities, decreasing = TRUE)
  if(all(size_communities < 2)) {
    warning("No clusters found in the PC network!\n",
            "Consider lowering the value of cos_cutoff.")
    return(NULL)
  }
  ind_cons_loading <- size_communities > 1
  membership_loading <- igraph::membership(cluster_pc)
  
  # Generate consensus loadings
  if(verbose) message("Calculating consensus loadings...")
  mat_cons_loading <- 
    vapply(as.integer(names(size_communities)[ind_cons_loading]),
           function(i) {
             # reorder the nodes in a clsuter so that the highest degree one 
             # comes first
             i_degrees <- 
               igraph::degree(graph_pc,
                              v = igraph::V(graph_pc)[membership_loading == i])
             i_order_degrees <- order(i_degrees, decreasing = TRUE)
             i_loading <- 
               mat_data_loading[, membership_loading == i][, i_order_degrees]
             for(j in (2:ncol(i_loading))) {
               if (i_loading[, 1] %*% i_loading[, j] < 0)
                 i_loading[, j] <- -i_loading[, j]
             }
             i_cons_loading <- apply(i_loading,
                                     1,
                                     mean)
             i_cons_loading / sqrt(sum(i_cons_loading^2))
           },
           rep_len(0.0, nrow(mat_data_loading)))
  colnames(mat_cons_loading) <- 
    paste0("Cluster_", names(size_communities)[ind_cons_loading])
  
  # Internal validation
  mat_vali <- t(matrix(vapply(data_loadings, 
                              function(loadings) {
                                apply(abs(t(mat_cons_loading) %*% loadings), 
                                      1, max)
                              },
                              rep_len(0.0, ncol(mat_cons_loading))), 
                       nrow = ncol(mat_cons_loading)))
  colnames(mat_vali) <- names(size_communities)[ind_cons_loading]
  rownames(mat_vali) <- lvl_batch
  
  # If required, visualize the clustered network
  if(!is.null(control$network_plot)) {
    short_names <- 
      unlist(lapply(shorten_name(lvl_batch), function(lvl) {
        paste0(lvl, ", PC", seq_len(n_pc_top))
      }))
    names(short_names) <- colnames(mat_cos)
    visualize_continuous_discover(graph_pc = graph_pc,
                                  membership_loading = membership_loading,
                                  size_communities = size_communities,
                                  plot_size_cutoff = plot_size_cutoff,
                                  short_names = short_names,
                                  output = control$network_plot)
  }
  
  # If required, generate diagnostic plots
  if(!is.null(control$diagnostic_plot))
    diagnostic_continuous_discover(mat_vali = mat_vali,
                                   lvl_batch = shorten_name(lvl_batch),
                                   cos_cutoff = cos_cutoff,
                                   output = control$diagnostic_plot)
  
  return(list(consensus_scores = t(feature_pca) %*% mat_cons_loading,
              consensus_loadings = mat_cons_loading,
              mat_vali = mat_vali,
              network = graph_pc,
              communities = cluster_pc,
              mat_cos = mat_cos,
              control = control))
}
