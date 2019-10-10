#' Visualization of the clustered network for the continuous.discover function
#'
#' @param graph_pc the full pc network constructed from correlated PCs
#' @param membership_loading membership of PC loadings from community discovery
#' @param size_communities ordered (largest to smallest) size of the identified 
#' communities
#' @param plot_size_cutoff cluster size cutoff (for cluster to be included in 
#' the visualized PC network)
#' @param short_names shorter names of the loadings
#' @param output output file name
#'
#' @return an invisible list of the subsetted network and memberships (to 
#' reproduce the plot)
#' @import igraph
#' @keywords internal
visualize_continuous_discover <- function(graph_pc,
                                          membership_loading,
                                          size_communities,
                                          plot_size_cutoff,
                                          short_names,
                                          output) {
  # Subsetting the network to exclude singletons and those that didn't pass the 
  # minimal cluster size cutoff.
  ind_visualize <- size_communities >= plot_size_cutoff
  cluster_sub <- as.integer(names(size_communities)[!ind_visualize])
  graph_pc_sub <-  delete.vertices(graph_pc,
                                   V(graph_pc)[membership_loading %in% 
                                                 cluster_sub])
  membership_sub <- membership_loading[names(V(graph_pc_sub))]
  list_membership_sub <- 
    lapply(as.integer(names(size_communities)[ind_visualize]),
           function(x) names(membership_sub)[membership_sub == x])
  
  pdf(output, width = 8, height = 8)
  plot.igraph(graph_pc_sub,
              mark.groups = list_membership_sub,
              vertex.size = 
                degree(graph_pc_sub) / max(degree(graph_pc_sub)) * 15,
              vertex.label = short_names[names(V(graph_pc_sub))],
              edge.width = E(graph_pc_sub)$weight * 10)
  dev.off()
  invisible(list(graph_pc = graph_pc_sub,
                 list_membership = list_membership_sub))
}

#' Diagnostic visualization for continuous.discover function
#'
#' @param mat_vali matrix of maximum correlations between the cluster-specific 
#' consensus loadings and top PC loadings from each batch
#' @param lvl_batch unique batches in the data
#' @param cos_cutoff the specified consine coefficient cutoff
#' @param output output file name
#'
#' @return the invisble ggplot2 plot object
#' @import ggplot2
#' @keywords internal
diagnostic_continuous_discover <- function(mat_vali,
                                           lvl_batch,
                                           cos_cutoff,
                                           output) {
  df_mali <- data.frame(mat_vali, check.names = FALSE)
  df_mali$batch <- factor(lvl_batch, levels = lvl_batch)
  df_plot <- tidyr::gather(df_mali,
                           key = community,
                           value = correlation,
                           -batch)
  plot <- ggplot(df_plot,
                 aes(x = community,
                     y = correlation)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text(aes(label = batch)) +
    geom_hline(yintercept = cos_cutoff,
               color = "red") +
    theme_bw() +
    xlab("Clusters") + ylab("Validation correlation")
  
  ggsave(plot = plot, filename = output,
         device = "pdf",
         width = 8, 
         height = 4 + ifelse(length(lvl_batch) > 4,
                             (length(lvl_batch) - 4) * 0.5,
                             0), 
         units = "in")
  invisible(plot)
}
