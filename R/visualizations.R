#' Diagnostic visualization for adj_batch function
#'
#' @param feature_abd original feature-by-sample matrix of abundances
#' (proportions or counts).
#' @param feature_abd_adj ffeature-by-sample matrix of batch-adjusted feature
#' abundances, with covariate effects retained and scales consistent with
#' original abundance matrix.
#' @param batch the batch variable (should be a factor).
#' @param gamma_hat estimated per feature-batch gamma parameters.
#' @param gamma_star shrinked per feature-batch gamma parameters
#'
#' @return (invisbly) the ggplot2 plot object
#' @import ggplot2
diagnostics_adjust_batch <- function(feature_abd,
                                     feature_abd_adj,
                                     batch,
                                     gamma_hat,
                                     gamma_star) {
  feature_abd <- fill_dimnames(feature_abd, "Feature", "Sample")
  dimnames(feature_abd_adj) <- dimnames(feature_abd)
  if(!is.factor(batch))
    stop("batch should be a factor!")

  # Plot gamma (i.e. location) parameters before and after shrinkage
  df_plot <- data.frame(gamma_hat = as.vector(gamma_hat),
                        gamma_star = as.vector(gamma_star),
                        batch = factor(rep(levels(batch),
                                           each = nrow(gamma_hat)),
                                       levels = levels(batch)))
  df_plot <- subset(df_plot, !is.na(gamma_hat), !is.na(gamma_star))
  p_shrinkage <- ggplot(df_plot, aes(x = gamma_hat, y = gamma_star,
                                     color = batch)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_discrete(labels = shorten_name(levels(batch))) +
    ggtitle("Shrinkage of batch mean parameters") +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.direction = "horizontal",
          legend.background = element_blank()) +
    xlab("Gamma") + ylab("Gamma (shrinked)")

  # Plot each feature's per-batch and overall mean relative abundances,
  # before and after adjustment
  # matrix for relative abundances
  mat_ra <- normalize_features(feature_abd, "TSS")
  mat_ra_adj <- normalize_features(feature_abd_adj, "TSS")

  # Prepare data frame of per-batch means
  df_mean_batch <- as.data.frame(apply(mat_ra, 1,
                                 function(x) tapply(x, batch, mean)))
  df_mean_batch_adj <- as.data.frame(apply(mat_ra_adj, 1,
                                     function(x) tapply(x, batch, mean)))
  colnames(df_mean_batch) <-
    colnames(df_mean_batch_adj) <-
    rownames(feature_abd)
  df_mean_batch$Batch <-
    df_mean_batch_adj$Batch <-
    levels(batch)
  df_mean_batch$Adjustment <- "Original"
  df_mean_batch_adj$Adjustment <- "Adjusted"
  df_batch <- rbind(df_mean_batch, df_mean_batch_adj)
  df_batch$Adjustment <- factor(df_batch$Adjustment,
                              levels = c("Original", "Adjusted"))
  df_batch <- tidyr::gather(df_batch,
                            key = "Feature",
                            value = "mean_batch",
                            - Adjustment, - Batch)

  # Prepare data frame of overall means
  df_mean_overall <- data.frame(mean_overall =
                                  apply(mat_ra, 1, mean))
  df_mean_overall_adj <- data.frame(mean_overall =
                                      df_mean_overall$mean_overall +
                                      max(df_mean_overall$mean_overall) /
                                      100)
  df_mean_overall$Feature <-
    df_mean_overall_adj$Feature <-
    rownames(feature_abd)
  df_mean_overall$Adjustment <- "Original"
  df_mean_overall_adj$Adjustment <- "Adjusted"
  df_overall <- rbind(df_mean_overall, df_mean_overall_adj)
  df_overall$Adjustment <- factor(df_overall$Adjustment,
                                levels = c("Original", "Adjusted"))

  # Merge data and plot
  df_plot <- merge(df_batch, df_overall, by = c("Feature", "Adjustment"))
  p_mean <- ggplot(df_plot, aes(x = mean_overall,
                                y = mean_batch)) +
    geom_point(aes(color = Adjustment)) +
    geom_line(aes(color = Adjustment, group = paste0(Feature, Adjustment))) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("Original" = "black", "Adjusted" = "red")) +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.direction = "horizontal",
          legend.background = element_blank()) +
    ggtitle("Original/adjusted mean abundance") +
    xlab("Overal mean") + ylab("Batch mean")

  # Because missing values
  plot <- cowplot::plot_grid(p_shrinkage, p_mean, nrow = 1)
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
  df_plot <- tidyr::gather(df.mali,
                           key = community,
                           value = correlation,
                           -batch)
  p <- ggplot(df_plot,
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
