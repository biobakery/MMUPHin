#' #' Main function for continuous structure discovery
#' #'
#' #' @param feature.count Feature x sample matrix of feature abundance
#' #' @param batch Variable indicating batch membership
#' #' @param pseudo.count Pseudo count to add to the count table before log-transformation. Default
#' #' to 0.5.
#' #' @return A list of MMUPHin continuous score discovery training object
#' #' @importFrom magrittr %>%
#' #' @export
#' continuous.discover <- function(feature.count,
#'                              batch,
#'                              pseudo.count = 0.5,
#'                              var.perc.cutoff = 0.3,
#'                              cor.cutoff = NULL,
#'                              directory = "./",
#'                              diagnostics = TRUE,
#'                              verbose = TRUE,
#'                              ...) {
#'   ## Ensure data formatts are as expected
#'   feature.count <- as.matrix(feature.count)
#'   if(is.null(rownames(feature.count)))
#'     rownames(feauture.count) <- paste0("Feature", 1:nrow(feature.count))
#'   if(any(feature.count < 0, na.rm = TRUE))
#'     stop("Found negative values in the feature table!")
#'   if(any(is.na(feature.count)))
#'     stop("Found missing values in the feature table!")
#'   batch <- as.factor(batch)
#'   if(any(is.na(batch))) stop("Found missing values in the batch variable!")
#'
#'   ## Data dimensions need to agree with each other
#'   if(ncol(feature.count) != length(batch)) {
#'     stop("Dimensions of feature table and batch variable do not agree!")
#'   }
#'
#'   lvl.batch <- levels(batch)
#'   n.batch <- length(lvl.batch)
#'   if(n.batch < 3) stop("Only discovery based on three or more studies is supported!")
#'   if(verbose) message("Found ", n.batch, " batches/studies.")
#'
#'   ## Transform sample counts
#'   lib.size <- apply(feature.count, 2, sum)
#'   if(all(lib.size <= 1)) {
#'     warning("Feature table appears to be on the relative abundance scale. ",
#'             "Setting pseudo.count to be min(feature.count)/2! ")
#'     pseudo.count <- min(setdiff(feature.count, 0)) / 2
#'   }
#'   log.data <- log(apply(feature.count + pseudo.count,
#'                         2,
#'                         function(x) x / sum(x)))
#'
#'
#'   # Calculate PC for the training sets
#'   if(verbose) message("Performing PCA in individual datasets...")
#'   pca.all <- lapply(lvl.batch, function(lvl)
#'   {
#'     pc <- log.data[, batch == lvl]
#'     dat.pca <- prcomp(t(pc))
#'     return(dat.pca)
#'   })
#'
#'   # Specify the smallest number of PCs to include
#'   ind.var.perc <- lapply(pca.all, function(dat.pca) {
#'     (cumsum(dat.pca$sdev^2) / sum(dat.pca$sdev^2)) > var.perc.cutoff
#'   })
#'   n.pc.top <- min(sapply(ind.var.perc, length))
#'   ind.var.perc <- sapply(ind.var.perc, function(x) x[1:n.pc.top])
#'   ind.var.perc <- apply(ind.var.perc, 1, all)
#'   if(any(ind.var.perc)) {
#'     n.pc.top <- min((1:n.pc.top)[ind.var.perc])
#'     if(verbose) message("Smallest number of PCs passing the specified threshold (",
#'                         var.perc.cutoff,
#'                         ") is ",
#'                         n.pc.top,
#'                         ".")
#'   } else {
#'     warning("No number of PCs passed the threshold!\n",
#'             "Setting number of PCs per dataset to largest possible",
#'             " value (",
#'             n.pc.top,
#'             ").")
#'   }
#'
#'   # First n.pc.top PC loadings for each training dataset
#'   data.loadings <- lapply(pca.all, function(x)
#'   {
#'     loadings <- x$rotation[,1:n.pc.top]
#'     return(loadings)
#'   })
#'   mat.data.loading <- Reduce("cbind", data.loadings)
#'   colnames(mat.data.loading) <- unlist(lapply(lvl.batch, function(lvl) {
#'     paste0(lvl, ", PC", 1:n.pc.top)
#'   }))
#'
#'   # Calculate correlation matrix of loadings across datasets
#'   if(verbose) message("Calculating correlation between PCs across datasets...")
#'   cor.matrix <- abs(t(mat.data.loading) %*% mat.data.loading)
#'   if(!is.null(cor.cutoff)) {
#'     cor.matrix[cor.matrix < cor.cutoff] <- 0
#'   }
#'
#'   # Create igraph graph object
#'   pc.graph <- igraph::graph_from_adjacency_matrix(cor.matrix,
#'                                                   mode = "undirected",
#'                                                   weighted = TRUE,
#'                                                   diag = FALSE)
#'   # Perform graph community detection
#'   if(verbose) message("Performing network clustering...")
#'   pc.cluster <- igraph::cluster_fast_greedy(pc.graph)
#'   size.communities <- igraph::sizes(pc.cluster)
#'   if(all(size.communities < 2)) {
#'     warning("No clusters found in the PC network!\n",
#'             "Consider lowering the value of cor.cutoff.")
#'     return(NULL)
#'   }
#'   ind.consensus.loading <- size.communities > 1
#'   membership.loading <- igraph::membership(pc.cluster)
#'
#'   # Generate consensus loadings
#'   if(verbose) message("Calculating consensus loadings...")
#'   mat.cons.loading <- sapply((1:max(membership.loading))[ind.consensus.loading],
#'                              function(i) {
#'                                loading.tmp <- mat.data.loading[, membership.loading == i]
#'                                for(j in (2:ncol(loading.tmp))) {
#'                                  if (loading.tmp[, 1] %*% loading.tmp[, j] < 0)
#'                                    loading.tmp[, j] <- -loading.tmp[, j]
#'                                }
#'                                apply(loading.tmp,
#'                                      1,
#'                                      mean)
#'                              })
#'   # Internal validation
#'   mat.cor.vali <- sapply(data.loadings, function(loadings) {
#'     apply(abs(cor(mat.cons.loading, loadings)), 1, max)
#'   })
#'
#'   pc.graph.tmp <- pc.graph %>%
#'     delete.vertices(
#'       V(pc.graph)[membership(pc.cluster) %in%
#'                     (names(sizes(pc.cluster))[sizes(pc.cluster) == 1] %>% as.numeric)]
#'     )
#'   list.membership.tmp <- membership(pc.cluster)[names(V(pc.graph.tmp))]
#'   list.membership.tmp <- unique(list.membership.tmp) %>%
#'     lapply(function(x) names(list.membership.tmp)[list.membership.tmp == x])
#'   plot.igraph(pc.graph.tmp, mark.groups = list.membership.tmp,
#'               vertex.size = degree(pc.graph.tmp) / max(degree(pc.graph.tmp)) * 15,
#'               edge.width = E(pc.graph.tmp)$weight * 10)
#'   pc.cluster.tmp <- igraph::cluster_fast_greedy(pc.graph.tmp)
#'   plot(pc.cluster.tmp, pc.graph.tmp)
#'   # Visualization
#'   pdf(paste0(directory, "network_communities.pdf"),
#'       width = 10,
#'       height = 10)
#'   plot(pc.cluster, pc.graph)
#'   dev.off()
#'   if(diagnostics) {
#'     df.cor.vali <- data.frame(mat.cor.vali)
#'     colnames(df.cor.vali) <- lvl.batch
#'     df.cor.vali$community <- factor(1:nrow(df.cor.vali))
#'     df.plot <- tidyr::gather(df.cor.vali,
#'                              key = dataset,
#'                              value = correlation,
#'                              -community)
#'     plot.diagnostic <- ggplot2::ggplot(df.plot,
#'                                        ggplot2::aes(x = community,
#'                                                     y = correlation)) +
#'       ggplot2::geom_boxplot(outlier.shape = NA) +
#'       ggplot2::geom_point(position = ggplot2::position_jitter()) +
#'       ggplot2::geom_text(ggplot2::aes(label = dataset), size = 3) +
#'       ggplot2::theme_bw()
#'     if(!is.null(cor.cutoff))
#'       plot.diagnostic <- plot.diagnostic +
#'       ggplot2::geom_hline(yintercept = cor.cutoff,
#'                           color = "red")
#'     ggsave(filename = paste0(directory, "diagnostic.pdf"),
#'            plot.diagnostic,
#'            width = min(20, sum(size.communities > 1)),
#'            height = 4)
#'   }
#'
#'   list(
#'     consensus.loading = mat.cons.loading,
#'     membership = lapply((1:length(size.communities))[size.communities > 1],
#'                         function(i) colnames(mat.data.loading)[igraph::membership(pc.cluster) == i]),
#'     mat.cor = mat.cor.vali
#'   )
#' }
