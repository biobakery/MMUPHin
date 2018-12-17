#' Main function for continuous structure discovery
#'
#' @param feature.count Feature x sample matrix of feature abundance
#' @param batch Name of batch variable
#' @param data data frame
#' @param normalization normalization, consistent with Maaslin2
#' @param transform transformation, consistent with Maaslin2
#' @param var.perc.cutoff percentage variance cutoff for top PCs
#' @param cor.cutoff correlation cutoff for constructing graph
#' @param directory output directory for figures
#' @param diagnostics logical. Should diagnostics be plotted?
#' @param verbose logical. Verbose?
#' @return A list of MMUPHin continuous score discovery object
#' @export
continuous.discover <- function(feature.count,
                                batch,
                                data,
                                normalization = "TSS",
                                transform = "AST",
                                var.perc.cutoff = 0.8,
                                cor.cutoff = 0.707,
                                directory = "./",
                                diagnostics = TRUE,
                                verbose = TRUE) {
  ## Ensure data formatts are as expected
  feature.count <- as.matrix(feature.count)
  if(any(is.na(feature.count)))
    stop("Found missing values in the feature table!")
  if(any(feature.count < 0))
    stop("Found negative values in the feature table!")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if(!(batch%in% names(data)))
    stop("Batch/covariate variable not found in data.")

  ## Data dimensions need to agree with each other
  if(ncol(feature.count) != nrow(data))
    stop("Dimensions of feature table and metadata table do not agree!")

  ## Check that sample names agree between the feature and metadata table
  ## And assign row and column names if emppty
  if(is.null(colnames(feature.count))) colnames(feature.count) <-
    paste0("Sample",
           1:ncol(feature.count))
  if(is.null(rownames(feature.count))) rownames(feature.count) <-
    paste0("Feature",
           1:nrow(feature.count))
  if(is.null(rownames(data))) rownames(data) <-
    paste0("Sample",
           1:ncol(feature.count))
  if(any(colnames(feature.count) != rownames(data)))
    stop("Sample names in feature.count and data don't agree!")

  # Check batch variable and identify groups
  batch <- data[, batch]
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  batch <- as.factor(batch)
  n.batch <- length(unique(batch))
  if(n.batch < 3) stop("Only discovery based on three or more studies is supported!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- levels(batch)

  ## Transform sample counts
  feature.pca <- normalizeFeatures(feature.count, normalization = normalization)
  feature.pca <- transformFeatures(feature.pca, transformation = transform)

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
  if(sum(edge.matrix) == nrow(edge.matrix)) {
    warning("No clusters found in the PC network (all edges are filtered out)!\n",
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
                               loading.tmp <- mat.data.loading[, membership.loading == i]
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

  # Visualize clustered network
  pc.graph.sub <-  igraph::delete.vertices(pc.graph,
                                           igraph::V(pc.graph)[membership.loading %in%
                                                                 as.numeric(names(size.communities)[
                                                                   size.communities <= 2])])
  list.membership.sub <- membership.loading[names(igraph::V(pc.graph.sub))]
  list.membership.sub <- lapply(unique(list.membership.sub),
                                function(x) names(list.membership.sub)[list.membership.sub == x])
  # Visualization
  if(length(igraph::V(pc.graph.sub)) > 0) {
    pdf(paste0(directory, "network_communities.pdf"),
        width = 10,
        height = 10)
    igraph::plot.igraph(pc.graph.sub, mark.groups = list.membership.sub,
                        vertex.size = igraph::degree(pc.graph.sub) / max(igraph::degree(pc.graph.sub)) * 15,
                        edge.width = igraph::E(pc.graph.sub)$weight * 10)
    dev.off()
  }

  # Internal validation
  mat.cor.vali <- sapply(data.loadings, function(loadings) {
    apply(abs(t(mat.cons.loading) %*% loadings), 1, max)
  })
  mat.cor.vali <- matrix(mat.cor.vali, nrow = ncol(mat.cons.loading))
  if(diagnostics) {
    df.cor.vali <- data.frame(mat.cor.vali)
    colnames(df.cor.vali) <- lvl.batch
    df.cor.vali$community <- factor(1:nrow(df.cor.vali))
    df.plot <- tidyr::gather(df.cor.vali,
                             key = dataset,
                             value = correlation,
                             -community)
    plot.diagnostic <- ggplot2::ggplot(df.plot,
                                       ggplot2::aes(x = community,
                                                    y = correlation)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitter()) +
      ggplot2::geom_text(ggplot2::aes(label = dataset), size = 3) +
      ggplot2::theme_bw()
    if(!is.null(cor.cutoff))
      plot.diagnostic <- plot.diagnostic +
      ggplot2::geom_hline(yintercept = cor.cutoff,
                          color = "red")
    ggsave(filename = paste0(directory, "diagnostic.pdf"),
           plot.diagnostic,
           width = min(20, sum(size.communities > 1)),
           height = 4)
  }

  list(
    consensus.loading = mat.cons.loading,
    membership = lapply((1:length(size.communities))[size.communities > 1],
                        function(i) colnames(mat.data.loading)[igraph::membership(pc.cluster) == i]),
    mat.cor = mat.cor.vali
  )
}
