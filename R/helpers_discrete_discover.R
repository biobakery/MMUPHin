#' Check dissimilarity object
#'
#' Make sure that the input is a dissimilarity object
#'
#' @param D dissimilarity object.
#'
#' @return returns an error if D is not a dissimilarity. Otherwise D as a 
#' matrix.
#' @keywords internal
check_D <- function(D) {
  if(!inherits(D, "dist")) {
    stop("D must be a dissimilarity object!")
  }
  as.matrix(D)
}

#' Check that sample numbers and names match between a dissimilarity matrix and 
#' a metadata data frame
#'
#' Sample names (row/column names of the D matrix, row names of the metadata
#' data frame) must be matching exactly. Note that this dictates that they
#' cannot be NULL because by design data (a data frame) should have non-empty
#' row names.
#'
#' @param D sample-by-sample matrix of dissimilarities (proportions or
#' counts).
#' @param data data frame of metadata.
#'
#' @return matched sample names
#' @keywords internal
check_samples_D <- function(D, data) {
  # Sample numbers need to agree with each other
  if(ncol(D) != nrow(data) |
     nrow(D) != nrow(data))
    stop("Dimensions of D matrix and metadata table do not agree!")
  
  # By designed use case data should always have row names
  if(is.null(rownames(data)))
    stop("data should not have empty row names!")
  
  # Sample names need to agree with each other
  if(!identical(colnames(D), rownames(data)) |
     !identical(rownames(D), rownames(data)) )
    stop("Sample names in feature_abd and data don't agree!")
  
  return(rownames(data))
}

check_k <- function(k, batch) {
  # Number of max clusters to evaluate can't exceed half of smallest sample size
  if(k > floor(min(table(batch)) / 2) - 1) {
    message("k_max must be less than half smallest sample size!")
    k <- floor(min(table(batch)) / 2) - 1
    message("Set k_max to ", k)
  }
  return(k)
}

#' Diagnostic visualization for discrete.discover function
#'
#' @param stats_internal list of internal evaluation summary statistics
#' @param stats_external list of external validation summary statistics
#' @param lvl_batch unique batches in the data
#' @param output
#'
#' @return the invisble ggplot2 plot object
#' @import ggplot2
#' @keywords internal
diagnostic_discrete_discover <- function(stats_internal,
                                         stats_external,
                                         lvl_batch,
                                         output) {
  k_max <- length(stats_internal) + 1
  df_internal <- Reduce("rbind",
                        lapply(2:k_max,
                               function(k)
                                 data.frame(
                                   k = k,
                                   batch = lvl_batch,
                                   Reduce("rbind", stats_internal[[k-1]])
                                 )
                        )
  )
  df_internal$batch <- factor(df_internal$batch, levels = lvl_batch)
  df_external <- Reduce("rbind",
                        lapply(2:k_max,
                               function(k)
                                 data.frame(
                                   k = k,
                                   batch = lvl_batch,
                                   Reduce("rbind", stats_external[[k-1]])
                                 )
                        )
  )
  df_external$batch <- factor(df_external$batch, levels = lvl_batch)
  
  p_internal <- ggplot(df_internal,
                       aes(x = k, y = mean)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  width = 0.5) +
    facet_grid(.~batch) +
    theme_bw() +
    xlab("K") + ylab("Summary statistic") +
    ggtitle("Internal evaluation")
  p_external <- ggplot(df_external,
                       aes(x = k, y = mean)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  width = 0.5) +
    facet_grid(.~batch) +
    theme_bw() +
    xlab("K") + ylab("Summary statistic") +
    ggtitle("External validation")
  
  plot <- cowplot::plot_grid(p_internal, p_external, nrow = 2)
  
  ggsave(plot = plot, filename = output,
         device = "pdf",
         width = 8, height = 8, units = "in")
  invisible(plot)
}
