#' Diagnostic visualization for adj.batch function
#'
#' @param feature.abd feature*sample matrix of unadjusted abundance
#' @param feature.abd.adj feature*sample matrix of adjusted abundance
#' @param batch vector of batch variables
#' @param gamma.hat estimated per feature-batch gamma parameters
#' @param gamma.star shrinked per feature-batch gamma parameters
#'
#' @return the invisble ggplot2 plot object
#' @import ggplot2
#' @importFrom magrittr %>%
diagnostics.adjust.batch <- function(feature.abd,
                                     feature.abd.adj,
                                     batch,
                                     gamma.hat,
                                     gamma.star) {
  df.plot <- data.frame(gamma.hat = as.vector(gamma.hat),
                        gamma.star = as.vector(gamma.star))
  p.shrinkage <- ggplot(df.plot, aes(x = gamma.hat, y = gamma.star)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    theme_bw() +
    ggtitle("Shrinkage of batch mean parameters")

  df.ra <- as.data.frame(t(apply(feature.abd, 2, function(x) x / sum(x)))) %>%
    dplyr::mutate(Sample = colnames(feature.abd),
                  Adjustment = "Original")
  df.ra.norm <- as.data.frame(t(apply(feature.abd.adj, 2, function(x) x / sum(x)))) %>%
    dplyr::mutate(Sample = colnames(feature.abd.adj),
                  Adjustment = "Adjusted")
  df.plot <- rbind(df.ra, df.ra.norm) %>%
    tidyr::gather(key = "feature",
                  value = "relative abundance",
                  -Sample, -Adjustment) %>%
    dplyr::left_join(data.frame(batch = batch,
                                Sample = df.ra$Sample,
                                stringsAsFactors = FALSE),
                     by = "Sample") %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(mean_overall = mean(`relative abundance`[Adjustment == "Original"])) %>%
    dplyr::group_by(feature, batch, Adjustment, mean_overall) %>%
    dplyr::summarise(mean_batch = mean(`relative abundance`)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Adjustment = factor(Adjustment, levels = c("Original", "Adjusted")))

  p.mean <- ggplot(df.plot, aes(x = mean_overall,
                                y = mean_batch)) +
    geom_point(aes(group = paste0(Adjustment, feature),
                   color = Adjustment),
               position = position_dodge(width = max(df.plot$mean_overall) / 100)) +
    geom_line(aes(group = paste0(Adjustment, feature),
                  color = Adjustment),
              position = position_dodge(width = max(df.plot$mean_overall) / 100)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("Original" = "black", "Adjusted" = "red")) +
    theme_bw() +
    theme(legend.position=c(0, 1),
          legend.justification=c(0, 1),
          legend.direction="horizontal",
          legend.background = element_blank()) +
    ggtitle("Original/adjusted mean abundance") +
    xlab("Overal mean") + ylab("Batch mean")

  suppressWarnings(plot <- cowplot::plot_grid(p.shrinkage, p.mean, nrow = 1)) # Because missing values
  print(plot)
  invisible(plot)
}

#' Diagnostic visualization for discrete.discover function
#'
#' @param stats.internal list of internal evaluation summary statistics
#' @param stats.external list of external validation summary statistics
#' @param lvl.batch unique batches in the data
#'
#' @return the invisble ggplot2 plot object
#' @import ggplot2
#' @importFrom magrittr %>%
diagnostics.discrete.discover <- function(stats.internal,
                                          stats.external,
                                          lvl.batch) {
  k.max <- length(stats.internal) + 1
  df.internal <- Reduce("rbind",
                        lapply(2:k.max,
                               function(k)
                                 data.frame(k = k,
                                            batch = lvl.batch,
                                            Reduce("rbind", stats.internal[[k-1]])))
  )
  df.internal$batch <- factor(df.internal$batch, levels = lvl.batch)
  df.external <- Reduce("rbind",
                        lapply(2:k.max,
                               function(k)
                                 data.frame(k = k,
                                            batch = lvl.batch,
                                            Reduce("rbind", stats.external[[k-1]])))
  )
  df.external$batch <- factor(df.external$batch, levels = lvl.batch)

  p.internal <- ggplot2::ggplot(df.internal,
                                ggplot2::aes(x = k, y = mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd,
                                        ymax = mean + sd),
                           width = 0.5) +
    ggplot2::facet_grid(.~batch) +
    ggplot2::theme_bw() +
    ggplot2::xlab("K") + ggplot2::ylab("Summary statistic") +
    ggtitle("Internal evaluation")
  p.external <- ggplot2::ggplot(df.external,
                                ggplot2::aes(x = k, y = mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd,
                                        ymax = mean + sd),
                           width = 0.5) +
    ggplot2::facet_grid(.~batch) +
    ggplot2::theme_bw() +
    ggplot2::xlab("K") + ggplot2::ylab("Summary statistic") +
    ggtitle("External validation")

  plot <- cowplot::plot_grid(p.internal, p.external, nrow = 2)
  print(plot)
  invisible(plot)
}

#' Visualization of the clustered network for the continuous.discover function
#'
#' @param pc.graph the full pc network constructed from correlated PCs
#' @param membership.loading membership of PC loadings from community discovery
#' @param size.communities ordered (largest to smallest) size of the identified communities
#' @param plot.size.cutoff cluster size cutoff (for cluster to be included in the visualized
#' PC network.)
#'
#' @return an invisible list of the subsetted network and memberships (to reproduce the plot)
#' @import igraph
visualize.continuous.discover <- function(pc.graph,
                                          membership.loading,
                                          size.communities,
                                          plot.size.cutoff) {
  # Subsetting the network to exclude singletons and those that didn't pass the minimal cluster
  # size cutoff.
  ind.visualize <- size.communities >= plot.size.cutoff
  clusters.sub <- as.integer(names(size.communities)[!ind.visualize])
  pc.graph.sub <-  delete.vertices(pc.graph,
                                   V(pc.graph)[membership.loading %in% clusters.sub])
  membership.sub <- membership.loading[names(V(pc.graph.sub))]
  list.membership.sub <- lapply(as.integer(names(size.communities)[ind.visualize]),
                                function(x) names(membership.sub)[membership.sub == x])

  plot.igraph(pc.graph.sub,
              mark.groups = list.membership.sub,
              vertex.size = degree(pc.graph.sub) / max(degree(pc.graph.sub)) * 15,
              edge.width = E(pc.graph.sub)$weight * 10)
  invisible(list(pc.graph = pc.graph.sub,
                 list.membership = list.membership.sub))
}

#' Diagnostic visualization for continuous.discover function
#'
#' @param mat.vali matrix of maximum correlations between the cluster-specific consensus loadings and
#' top PC loadings from each batch
#' @param lvl.batch unique batches in the data
#'
#' @return the invisble ggplot2 plot object
#' @import ggplot2
diagnostics.continuous.discover <- function(mat.vali,
                                            lvl.batch) {
  df.mali <- data.frame(mat.vali, check.names = FALSE)
  df.mali$batch <- factor(lvl.batch, levels = lvl.batch)
  df.plot <- tidyr::gather(df.mali,
                           key = community,
                           value = correlation,
                           -batch)
  p <- ggplot(df.plot,
              aes(x = community,
                  y = correlation)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(seed = 1)) +
    ggrepel::geom_text_repel(aes(label = batch),
                             position = position_jitter(seed = 1)) +
    geom_hline(yintercept = cor.cutoff,
               color = "red") +
    theme_bw() +
    xlab("Clusters") + ylab("Validation correlation")

  print(p)
  invisible(p)
}
