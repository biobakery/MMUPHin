#' Main function for continuous structure discovery
#'
#' @param feature.abd feature*sample matrix of feature abundance (counts preferred if transform
#' is specified as LOG).
#' @param batch name of the batch variable.
#' @param data  data frame of metadata, must contain batch.
#' @param normalization normalization parameter.
#' @param transform transformation parameter.
#' @param pseudo.count pseudo count to add to the count table before log-transformation. Default
#' to 0.5. This should be changed for relative abundance (will automatically set to half of minimal
#' non-zero value if not specified).
#' @param var.perc.cutoff percentage variance explained cutoff to choose the top PCs.
#' @param cor.cutoff correlation cutoff to construct edges for the PC network.
#' @param plot.clustered.network should the clustered PC network be visualized? Deafault to TRUE.
#' @param plot.size.cutoff cluster size cutoff (for cluster to be included in the visualized
#' PC network.)
#' @param diagnostics should diagnostic plots be generated? Deafault to FALSE.
#' @param verbose should verbose modelling information should be printed? Default to TRUE.
#'
#' @return a list of continuous structure discovery results, including the following components:
#' consensus.loading: matrix of the identified consensus loadings, each column for one identified
#' PC loading cluster.
#' scores: matrix of the identified continuous scores for all samples, each column for one identified
#' PC loading cluster.
#' membership: list of PC loading members for each identified cluster.
#' mat.vali: matrix of validation correlations of the identified consensus loadings, each column for
#' one identified PC loading cluster.
#' clustered.network: list for the the constructed PC network and community discovery results
#'
#' @export
continuous.discover <- function(feature.abd,
                                batch,
                                data,
                                normalization = "TSS",
                                transform = "AST",
                                pseudo.count = 0.5,
                                var.perc.cutoff = 0.8,
                                cor.cutoff = 0.5,
                                plot.clustered.network = TRUE,
                                plot.size.cutoff = 2,
                                diagnostics = FALSE,
                                verbose = TRUE) {
  # Ensure data formatts are as expected
  feature.abd <- as.matrix(feature.abd)
  if(any(is.na(feature.abd)))
    stop("Found missing values in the feature table!")
  if(any(feature.abd < 0))
    stop("Found negative values in the feature table!")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if(!(batch%in% names(data)))
    stop("Batch/covariate variable not found in data.")

  # Data dimensions need to agree with each other
  if(ncol(feature.abd) != nrow(data))
    stop("Dimensions of feature table and metadata table do not agree!")

  # Check that sample names agree between the feature and metadata table
  # And assign row and column names if emppty
  if(is.null(colnames(feature.abd))) colnames(feature.abd) <-
    paste0("Sample",
           1:ncol(feature.abd))
  if(is.null(rownames(feature.abd))) rownames(feature.abd) <-
    paste0("Feature",
           1:nrow(feature.abd))
  if(is.null(rownames(data))) rownames(data) <-
    paste0("Sample",
           1:ncol(feature.abd))
  if(any(colnames(feature.abd) != rownames(data)))
    stop("Sample names in feature.abd and data don't agree!")

  # Check batch variable and identify groups
  batch <- data[, batch]
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  batch <- as.factor(batch)
  n.batch <- length(unique(batch))
  if(n.batch < 3) stop("Only discovery based on three or more studies is supported!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- levels(batch)

  # Normalize and transform sample counts
  if(transform == "LOG") {
    if(all(apply(feature.abd, 2, sum) <= 1)) {
      warning("Feature table appears to be on the relative abundance scale!")
      if(pseudo.count == 0.5) {
        warning("pseudo.count was set to 0.5 which is only appropriate to count data,",
                " setting pseudo.count to be min(feature.abd)/2.")
        pseudo.count <- min(setdiff(feature.abd, 0)) / 2
      }
      # In this case don't renormalize the feature
      feature.pca <- transformFeatures(feature.abd,
                                       transform = transform,
                                       pseudo.count = pseudo.count)
    }
    feature.pca <- transformFeatures(normalizeFeatures(feature.abd,
                                                       normalization = normalization,
                                                       pseudo.count = pseudo.count),
                                     transform = transform)
  } else {
    feature.pca <- normalizeFeatures(feature.abd, normalization = normalization)
    feature.pca <- transformFeatures(feature.pca, transform = transform)
  }

  # Calculate PC for the training sets
  if(verbose) message("Performing PCA in individual datasets...")
  pca.all <- lapply(lvl.batch, function(lvl)
  {
    pc <- feature.pca[, batch == lvl]
    dat.pca <- prcomp(t(pc))
    return(dat.pca)
  })

  # Specify the smallest number of PCs to include
  ind.var.perc <- lapply(pca.all, function(dat.pca) {
    (cumsum(dat.pca$sdev^2) / sum(dat.pca$sdev^2)) > var.perc.cutoff
  })
  n.pc.top <- min(sapply(ind.var.perc, length))
  ind.var.perc <- sapply(ind.var.perc, function(x) x[1:n.pc.top])
  ind.var.perc <- apply(ind.var.perc, 1, all)
  if(any(ind.var.perc)) {
    n.pc.top <- min((1:n.pc.top)[ind.var.perc])
    if(verbose) message("Smallest number of PCs passing the specified threshold (",
                        var.perc.cutoff,
                        ") is ",
                        n.pc.top,
                        ".")
  } else {
    warning("No number of PCs passed the threshold!\n",
            "Setting number of PCs per dataset to largest possible",
            " value (",
            n.pc.top,
            ").")
  }

  # First n.pc.top PC loadings for each training dataset
  data.loadings <- lapply(pca.all, function(x)
  {
    loadings <- x$rotation[,1:n.pc.top]
    return(loadings)
  })
  mat.data.loading <- Reduce("cbind", data.loadings)
  colnames(mat.data.loading) <- unlist(lapply(lvl.batch, function(lvl) {
    paste0(lvl, ", PC", 1:n.pc.top)
  }))

  # Calculate correlation matrix of loadings across datasets
  if(verbose) message("Calculating correlation between PCs across datasets...")
  cor.matrix <- t(mat.data.loading) %*% mat.data.loading
  edge.matrix <- matrix(1, nrow = nrow(cor.matrix), ncol = ncol(cor.matrix))
  # edge.matrix <- cor.matrix
  dimnames(edge.matrix) <- dimnames(cor.matrix)
  edge.matrix[abs(cor.matrix) < cor.cutoff] <- 0
  edge.matrix[abs(cor.matrix) >= cor.cutoff] <- 1
  if(sum(edge.matrix) == nrow(edge.matrix)) {
    warning("All edges are filtered out in the PC network!\n",
            "Consider lowering the value of cor.cutoff.")
    return(NULL)
  }

  # Create igraph graph object
  pc.graph <- igraph::graph_from_adjacency_matrix(abs(edge.matrix),
                                                  mode = "undirected",
                                                  weighted = NULL,
                                                  diag = FALSE)
  # Perform graph community detection
  if(verbose) message("Performing network clustering...")
  pc.cluster <- igraph::cluster_edge_betweenness(pc.graph)
  size.communities <- igraph::sizes(pc.cluster)
  size.communities <- sort(size.communities, decreasing = TRUE)
  if(all(size.communities < 2)) {
    warning("No clusters found in the PC network!\n",
            "Consider lowering the value of cor.cutoff.")
    return(NULL)
  }
  ind.consensus.loading <- size.communities > 1
  membership.loading <- igraph::membership(pc.cluster)

  # Generate consensus loadings
  if(verbose) message("Calculating consensus loadings...")
  mat.cons.loading <- sapply(as.integer(names(size.communities)[ind.consensus.loading]),
                             function(i) {
                               # reorder the nodes in a clsuter so that the highest degree one comes first
                               degrees.tmp <- igraph::degree(pc.graph,
                                                             v = igraph::V(pc.graph)[membership.loading == i])
                               order.degrees.tmp <- order(degrees.tmp, decreasing = TRUE)
                               loading.tmp <- mat.data.loading[, membership.loading == i][, order.degrees.tmp]
                               for(j in (2:ncol(loading.tmp))) {
                                 if (loading.tmp[, 1] %*% loading.tmp[, j] < 0)
                                   loading.tmp[, j] <- -loading.tmp[, j]
                               }
                               cons.loading <- apply(loading.tmp,
                                                     1,
                                                     mean)
                               cons.loading / sqrt(sum(cons.loading^2))
                             })
  colnames(mat.cons.loading) <- paste0("Cluster_", names(size.communities)[ind.consensus.loading])

  # Internal validation
  mat.vali <- t(matrix(sapply(data.loadings, function(loadings) {
    apply(abs(t(mat.cons.loading) %*% loadings), 1, max)
  }), nrow = ncol(mat.cons.loading)))
  colnames(mat.vali) <- names(size.communities)[ind.consensus.loading]
  rownames(mat.vali) <- lvl.batch

  # If required, visualize the clustered network
  if(plot.clustered.network) {
    visualize.continuous.discover(pc.graph = pc.graph,
                                  membership.loading = membership.loading,
                                  size.communities = size.communities,
                                  plot.size.cutoff = plot.size.cutoff)
  }
  # If required, generate diagnostic plots
  if(diagnostics) {
    diagnostics.continuous.discover(mat.vali = mat.vali,
                                    lvl.batch = lvl.batch)
  }

  return(
    list(consensus.loading = mat.cons.loading,
         scores = t(feature.pca) %*% mat.cons.loading,
         membership = lapply(as.integer(names(size.communities)[ind.consensus.loading]),
                             function(i) colnames(mat.data.loading)[
                               igraph::membership(pc.cluster) == i
                               ]),
         mat.vali = mat.vali,
         clustered.network = list(pc.graph = pc.graph,
                                  communities = pc.cluster)
    ))
}
